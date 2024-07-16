#include <set>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include "util.h"
#include "fasta/Fasta.h"


UserParams::UserParams(
    const int mapq_thresh,
    const int _base_q_thresh,
    const double lq_rate_thresh,
    const int match_score,
    const int mismatch_penal,
    const int gap_open_penal,
    const int gap_ext_penal,
    const int kmer_size,
    const int local_thresh,
    const int retarget_thresh
) : mapq_thresh(mapq_thresh),
    lq_rate_thresh(lq_rate_thresh), 
    match_score(match_score), 
    mismatch_penal(mismatch_penal),
    gap_open_penal(gap_open_penal), 
    gap_ext_penal(gap_ext_penal),
    kmer_size(kmer_size), 
    local_thresh(local_thresh),
    retarget_thresh(retarget_thresh)
{
    base_q_thresh = static_cast<char>(_base_q_thresh + 33);    
};


std::unordered_map<int, std::string_view> dictize_reference(
    std::string_view ref_seq, 
    int ref_start, 
    int ref_end
) 
{
  std::unordered_map<int, std::string_view> ref_dict;

  int i = 0, pos = ref_start;
  while (pos <= ref_end) 
  {
    ref_dict[pos] = ref_seq.substr(i, 1);
    ++i;
    ++pos;
  }

  return ref_dict;
}


LocalReference::LocalReference(
    const string& fastafile,
    const string& chrom,
    const int start,
    const int end
) : chrom(chrom),
    start(start),
    end(end)
{
    fasta.open(fastafile);
    
    _seq = fasta.getSubSequence(chrom, start - 1, end - start + 1);
    
    seq = static_cast<std::string_view>(_seq);
    
    dict = dictize_reference(seq, start, end);
};


inline bool contain_n(std::string_view sv1, std::string_view sv2)
{
    if (sv1.find('N') != std::string_view::npos) return true;
    if (sv2.find('N') != std::string_view::npos) return true;

    return false;
}


inline bool is_rotatable(std::string_view allele)
{
    return (allele.front() == allele.back());   
}


Variant::Variant 
(    int pos, 
     const std::string& ref, 
     const std::string& alt,
     bool is_clipped_segment
) : pos(pos), 
    ref(ref), 
    alt(alt),
    is_clipped_segment(is_clipped_segment),  
    ref_len(ref.size()), 
    alt_len(alt.size()),
    is_substitute(alt_len == ref_len),
    is_ins(alt_len > ref_len),
    is_del(alt_len < ref_len)
{
    has_n = contain_n(ref, alt);
     
    if (ref_len < alt_len)
    {
        if (alt.substr(0, ref_len) == ref) is_complex = false; 
        else is_complex = true;

        indel_len = alt_len - 1;
    }
    else if (ref_len > alt_len)
    {
        if (ref.substr(0, alt_len) == alt) is_complex = false;
        else is_complex = true;

        indel_len = ref_len - 1;
    }   
}


inline void to_left(
    int& pos, 
    int& variant_end_pos,
    std::string& longer_allele,
    std::string& shorter_allele,
    bool& is_failed,
    const std::unordered_map<int, std::string_view>& ref_dict
)
{
    --pos;
    --variant_end_pos;
    char prev_base = std::string{ref_dict.at(pos)}[0];
    if (prev_base == 'N')
    {
        ++pos;
        ++variant_end_pos;
        is_failed = true;
        return;
    }

    longer_allele.pop_back();
    longer_allele.insert(0, 1, prev_base);
    shorter_allele = prev_base;
    is_failed = false;
}


void left_align(
    int& pos,
    int& variant_end_pos, 
    std::string& ref, 
    std::string& alt, 
    const bool is_ins,
    const int ref_start,
    const std::unordered_map<int, std::string_view>& ref_dict
)
{
    std::string& longer_allele = (is_ins) ? alt : ref;
    std::string& shorter_allele = (is_ins) ? ref : alt;
    bool is_failed = false;
       
    while(is_rotatable(longer_allele) && ref_start < pos) 
    {
        if (ref_dict.find(pos - 1) == ref_dict.end()) return;
        
        to_left(
            pos, 
            variant_end_pos, 
            longer_allele, 
            shorter_allele, 
            is_failed, 
            ref_dict
        );
        
        if (is_failed) return;
    }
}


void Variant::_sb_leftmost_pos(const LocalReference& loc_ref)
{
    const std::string ptrn = (alt_len > 1) ? to_tandem_rep(alt) : alt;
    const size_t prtn_len = ptrn.size();
    
    lpos = pos;
    int curr_pos = pos - prtn_len;
    int idx = curr_pos - loc_ref.start;
    
    while (loc_ref.start < curr_pos)
    {
        if (alt != loc_ref.seq.substr(idx, prtn_len))
        {
            lpos = curr_pos + prtn_len;
            return;
        }
        idx -= prtn_len;
        curr_pos -= prtn_len;
    }   
}


void Variant::set_leftmost_pos(const LocalReference& loc_ref)
{
    if (is_complex 
        || is_clipped_segment 
        || has_n
        || pos <= loc_ref.start
        || loc_ref.end <= pos) 
    {    
        lpos = pos;   
        return;
    }

    if (is_substitute)
    {
        _sb_leftmost_pos(loc_ref);
        return;
    }
    
    //copy -> keep orig data
    int _pos = pos;
    int _variant_end_pos = variant_end_pos;
    std::string _ref = ref;
    std::string _alt = alt;

    left_align(
        _pos, 
        _variant_end_pos, 
        _ref, 
        _alt, 
        is_ins, 
        loc_ref.start, 
        loc_ref.dict
    );

    lpos = _pos;
}


inline void to_right(
    int& variant_end_pos, 
    std::string& longer_allele, 
    std::string& shorter_allele,
    bool& is_failed,
    const std::unordered_map<int, std::string_view>& ref_dict
)
{
    char next_base = std::string{ref_dict.at(variant_end_pos)}[0];
    if (next_base == 'N') 
    {
        is_failed = true;    
        return;
    }

    longer_allele.erase(0, 1);
    longer_allele += next_base;
    shorter_allele = next_base;
    ++variant_end_pos;
    is_failed = false;
    
}


void right_align(
    int& pos, 
    int& variant_end_pos, 
    std::string& ref, 
    std::string& alt, 
    const bool is_ins,
    const int ref_end,
    const std::unordered_map<int, std::string_view>& ref_dict
)
{
    std::string& longer_allele = (is_ins) ? alt : ref;
    std::string& shorter_allele = (is_ins) ? ref : alt;
    bool is_failed = false;

    do {
        if (ref_dict.find(variant_end_pos) == ref_dict.cend()) return;

        to_right(
            variant_end_pos, 
            longer_allele, 
            shorter_allele,
            is_failed,
            ref_dict
        );

        if (is_failed) return;
        
        ++pos;
    }
    while (is_rotatable(longer_allele) && pos < ref_end);
    
    // undo the last right shift
    if (ref_dict.find(pos - 1) != ref_dict.cend())
    {
        to_left(
            pos, 
            variant_end_pos, 
            longer_allele, 
            shorter_allele, 
            is_failed, ref_dict
        );
    }
}


void Variant::_sb_rightmost_pos(const LocalReference& loc_ref)
{
    const std::string ptrn = (alt_len > 1) ? to_tandem_rep(alt) : alt;
    const size_t prtn_len = ptrn.size();
    
    rpos = pos;
    int curr_pos = pos + prtn_len;
    int idx = curr_pos - loc_ref.start;
    
    while (curr_pos < loc_ref.end)
    {
        if (alt != loc_ref.seq.substr(idx, prtn_len))
        {
            rpos = curr_pos - prtn_len;
            return;
        }
        idx += prtn_len;
        curr_pos += prtn_len;
    }   
}


void Variant::set_rightmost_pos(const LocalReference& loc_ref)
{
    if (is_complex 
        || is_clipped_segment
        || loc_ref.start <= pos
        || loc_ref.end >= pos) 
    {    
        rpos = pos;
        return;
    }
    
    if (is_substitute)
    {
        _sb_rightmost_pos(loc_ref);
        return;
    }
    
    //copy -> keep orig info
    int _pos = pos;
    int _variant_end_pos = variant_end_pos;
    std::string _ref = ref;
    std::string _alt = alt;

    right_align(
        _pos, 
        _variant_end_pos, 
        _ref, 
        _alt, 
        is_ins,
        loc_ref.end, 
        loc_ref.dict
    );
    
    rpos = _pos;
}


bool Variant::is_equivalent(
    const Variant& other, 
    const LocalReference& loc_ref
) const
{
    
    int _pos = other.pos;
    int _varient_end_pos = other.variant_end_pos;
    std::string _ref = other.ref;
    std::string _alt = other.alt;
    
    //cplx other
    if (_ref.size() > 1 && _alt.size() > 1)
    {
        if (ref == _ref && pos == _pos) return true; //del seq match
        if (alt.substr(1) == _alt.substr(1))
        {
            //inserted seq (w/o padding base) may be anywhere in the deletion
            if (_pos <= pos 
                && pos <= _pos + static_cast<int>(_ref.size()) - 1) return true;
        }

        return false;
    }
    
    if (pos == _pos && ref == _ref && alt == _alt) return true;
    if (is_clipped_segment || other.is_clipped_segment) return false;
    if (ref_len != other.ref_len || alt_len != other.alt_len) return false;
    if (is_substitute) return false;
    
    left_align(
        _pos, 
        _varient_end_pos, 
        _ref, 
        _alt, 
        other.is_ins, 
        loc_ref.start, 
        loc_ref.dict
     );
            
     return (pos == _pos && ref == _ref && alt == _alt);
} 


bool operator == (const Variant& lhs, const Variant& rhs)
{
    return (lhs.pos == rhs.pos && lhs.ref == rhs.ref && lhs.alt == rhs.alt);
}


bool operator != (const Variant& lhs, const Variant& rhs)
{
    return !(lhs == rhs);
}


bool operator < (const Variant& lhs, const Variant& rhs)
{
    if (lhs.pos < rhs.pos) return true;
    else if (rhs.pos < lhs. pos) return false;
    else // lhs.pos == rhs.pos
    {
        if (lhs.ref.size() < rhs.ref.size()) return true;
        else if (rhs.ref.size() < lhs.ref.size()) return false;
        else return (lhs.alt < rhs.alt); //same ref deleted
    }
}


std::string to_tandem_rep(std::string_view seq)
{
    int seq_size = seq.size();
    if (seq_size < 2) 
    {
        return std::string{seq};
    }
    else 
    {
        int mid_len = (int)seq_size / 2;

        int step = 1;
        while (step <= mid_len) 
        {
            std::vector<std::string_view> tandems;
            for (int i = 0; step * i < seq_size; ++i) 
            {
                tandems.push_back(
                    seq.substr(step * i, std::min(step, seq_size - (step * i)))
                );    
            }
            std::sort(tandems.begin(), tandems.end());
            tandems.erase(
                std::unique(tandems.begin(), tandems.end()), tandems.end()
            );
            if (tandems.size() == 1) return std::string{tandems[0]};
            else ++step;    
        }
    }
    return std::string{seq};
}


std::string Variant::minimal_repeat_unit() const
{
    std::string seq = "";
    if (is_substitute) seq = alt;
    if (is_ins) seq = alt.substr(1);
    if (is_del) seq = ref.substr(1);
    
    return to_tandem_rep(seq);
}



void find_shared_variants(
    std::vector<Variant>& shared,
    std::vector<Variant>& var_vec1,
    const std::vector<Variant>& var_vec2 //pre sorted
)
{
    std::sort(var_vec1.begin(), var_vec1.end());
    
    std::set_intersection(
        var_vec1.begin(),var_vec1.end(),
        var_vec2.begin(),var_vec2.end(),
        back_inserter(shared)
    );
}


inline std::pair<char, int> split_cigar(const std::string& cigar) 
{
  size_t last_idx = cigar.size() - 1;
  return std::make_pair(cigar.substr(last_idx, 1)[0], std::stoi(cigar.substr(0, last_idx)));
}

std::vector<std::pair<char, int>> to_cigar_vector(std::string_view cigar_str) 
{
  std::vector<std::pair<char, int>> cigar_vec;

  size_t pos = 0;
  size_t new_pos;
  const size_t len = cigar_str.size();

  while (pos < len) 
  {
    new_pos = cigar_str.find_first_of("MIDNSHPX=", pos) + 1;
    cigar_vec.emplace_back(
        split_cigar(
            std::string{cigar_str.substr(pos, new_pos - pos)}
        )
    );
    pos = new_pos;
  }

  return cigar_vec;
}


std::pair<char, int> concat_gaps(const std::vector<std::pair<char, int>>& cigar_sub_vector)
{
    int tot_gap_len = 0;
    if (!cigar_sub_vector.empty())
    {
        for(const auto & cigar : cigar_sub_vector)
        {
            tot_gap_len += cigar.second;
        }

    }

    return {cigar_sub_vector[0].first, tot_gap_len};
}


//make insertion first for complex case with gap-merging
//{'=', 2}, {'D', 2}, {'I', 2}, {'D', 4}, {'=', 3}} 
// -> {{'=', 2}, {'I', 2}, {'D', 6},{'=', 3}}
//----------------------------------------------------------------------
void move_up_insertion(std::vector<std::pair<char, int>>& cigar_vector)
{
    bool gap_ended = true;
    std::vector<std::pair<char, int>> tmp, ins, del;

    for (const auto& cigar : cigar_vector)
    {
        if (cigar.first == 'I')
        {
            ins.push_back(cigar);
            gap_ended = false;   
        }
        else if (cigar.first == 'D')
        {
            del.push_back(cigar);
            gap_ended =false;   
        }
        else gap_ended = true;
        
        if (gap_ended)
        {
            if (!ins.empty())
            {
                tmp.push_back(concat_gaps(ins));
                ins.clear();
            }

            if (!del.empty())
            {
                tmp.push_back(concat_gaps(del));
                del.clear();
            }
            
            tmp.push_back(cigar);
        }   
    }

    std::swap(tmp, cigar_vector);
}


// experimental
// swap N and simple ins
// {<'=', 2>, <'N', 6>, <'I', '4'>, <'=', 2>} -> {<'=', 2>, <'I', '4'>,  <'N', 6>,  <'=', 2>}
//------------------------------------------------------------------------------------------
void make_skip_after_ins(CigarVec& cigar_vector)
{
    std::vector<int> idx_lst; 
    int last_idx = cigar_vector.size() - 1;
    for (int i = 0; i < last_idx - 1; ++i)
    {
        if (cigar_vector[i].first == 'N' && cigar_vector[i + 1].first == 'I')
        {
            if (cigar_vector[i + 2].first == '=') //simple case
            {
                idx_lst.push_back(i);
            }
        }
    } 
    
    for (const auto j : idx_lst)
    {
        std::iter_swap(cigar_vector.begin() + j, cigar_vector.begin() + j + 1);
    }

}


// find I/D*X
//----------------------------------------------------------------------
std::vector<int> find_mismatches_after_gap(const CigarVec& cigar_vector)
{
    int i = 0;
    char op = '\0';
    bool is_prev_gap = false;
    std::vector<int> indexes;
    for (const auto& cigar : cigar_vector)
    {
        op = cigar.first;
        switch (op)
        {   
            case 'X':
                if (is_prev_gap)
                {
                    indexes.push_back(i);
                }
                is_prev_gap = false;
                break;
            case 'I':
            case 'D':
                is_prev_gap = true;
                break;
            default:
                is_prev_gap = false;
                break;
        }
        ++i;
    }
    return indexes;
}


inline bool is_target_idx(const std::vector<int>& target_indexes, const int idx)
{
    if (
        std::find(
            target_indexes.begin(), target_indexes.end(), idx
        ) != target_indexes.end()
    ) 
    {
        return true;
    }
    else 
    {    
        return false;
    } 
}


void parse_to_cplx_gaps(
    const int aln_start, 
    std::string_view ref_seq,
    std::string_view read_seq,
    const CigarVec& cigar_vector,
    const std::vector<int>& target_indexes,
    std::vector<Variant>& variants 
)
{
    char op = '\0';
    int op_len = 0, ref_idx = 0, read_idx = 0, pos = aln_start, cigar_idx = 0;
    std::string_view ref, alt;
    for (const auto & cigar : cigar_vector)
    {
        op = cigar.first;
        op_len = cigar.second;
        switch (op)
        {
            case '=':
            case 'X':
                ref_idx += op_len;
                read_idx += op_len;
                pos += op_len;
                ++cigar_idx;
                break;
            case 'I':
                if (is_target_idx(target_indexes, cigar_idx + 1))
                { 
                    //del first
                    int del_len = cigar_vector[cigar_idx + 1].second;
                    ref = ref_seq.substr(ref_idx - 1, del_len + 1);
                    alt = ref_seq.substr(ref_idx - 1, 1);
                    variants.emplace_back(
                        pos - 1, 
                        std::string(ref),
                        std::string(alt)
                    );
                    
                    //ins next
                    ref = ref_seq.substr(ref_idx + del_len - 1, 1);
                    variants.emplace_back(
                        pos - 1 + del_len,                   
                        std::string(ref),
                        std::string(ref) 
                        + std::string(read_seq.substr(read_idx, op_len + del_len))
                    );
                } 
                ++cigar_idx;
                read_idx += op_len;
                break; 
            case 'D':
                if (is_target_idx(target_indexes, cigar_idx + 1)) 
                {
                    //del(not ins) first
                    int ins_len = cigar_vector[cigar_idx + 1].second;
                    ref = ref_seq.substr(ref_idx - 1, op_len + ins_len + 1);
                    alt = ref_seq.substr(ref_idx - 1, 1);
                    variants.emplace_back(
                        pos - 1,
                        std::string(ref),
                        std::string(alt)
                    );

                    //ins
                    ref =  ref_seq.substr(ref_idx + op_len - 1, 1);
                    variants.emplace_back(
                        pos - 1 + op_len,
                        std::string(ref),
                        std::string(ref)
                        + std::string(read_seq.substr(read_idx, ins_len))
                    );
                }
                ref_idx += op_len;
                pos += op_len;
                ++cigar_idx;
                break;
            case 'N':
                pos += op_len;
                ++cigar_idx;
                break;
            case 'S':
                read_idx += op_len;
                ++cigar_idx;
                break;
         }
     }
}


inline bool is_gap(const char op)
{
    return (op == 'I' || op == 'D');
}

// include skips (N) to cigar
// {<'=', 5>, <'I', 1>, <'=', 5> + {{100, 104}, {108, 109}, {112, 114}}
// -> {<'=', 5>, <'I', 1>, <'N', 3>, <'=', 2>, <'N', 2>, <'=', 3>}
//--------------------------------------------------------------------------
void splice_cigar(
    std::vector<std::pair<char, int>>& cigar_vector,
    const int start_offset, 
    const std::vector<int>& genomic_positions, 
    const Coord& coordinates
)
{
    
    // unspliced 
    if (coordinates.size() < 2)
    {
        return; 
    }
    
    std::vector<std::pair<char, int>> tmp; 
    
    char op;
    int op_len = 0;
    bool was_skip = false;
    int consumable_op_len = 0;
    int curr_idx = 0;
    int last_idx = coordinates.size() - 1;
    int curr_pos = coordinates[curr_idx].first + start_offset - 1;
    int curr_segment_start = coordinates[curr_idx].first;
    int curr_segment_end = coordinates[curr_idx].second;
    int curr_op_end = curr_pos;
    int genome_pos_idx = start_offset - 1; 
    
    for (const auto & cigar : cigar_vector)
    {
        op = cigar.first;
        op_len = cigar.second;
        consumable_op_len = 0;
        switch (op)
        {
            case 'M':
            case 'D':
            case '=':
            case 'X':
                consumable_op_len = cigar.second;
                break;
            default:
                break;
        }

        if (consumable_op_len)
        {
            genome_pos_idx += consumable_op_len;
            curr_op_end = genomic_positions[genome_pos_idx];
        }

        if (consumable_op_len)
        {
            while(curr_segment_end < curr_op_end)
            {
                tmp.emplace_back(op, (curr_segment_end - curr_pos));
                tmp.emplace_back('N', (coordinates[curr_idx + 1].first - (curr_segment_end + 1)));
                curr_pos = coordinates[curr_idx + 1].first - 1;
                
                ++curr_idx; 

                curr_segment_start = coordinates[curr_idx].first;
                curr_segment_end = coordinates[curr_idx].second; 

            }

            if (curr_op_end - curr_pos > 0) 
            {
                tmp.emplace_back(op, curr_op_end - curr_pos);
            }
            curr_pos = curr_op_end;
        }
        else 
        {
            tmp.emplace_back(op, op_len);   
        }

        if (is_gap(op) && was_skip)
        {   
            std::iter_swap(tmp.end() - 2 , tmp.end() - 1); 
        } 
        
        if (curr_op_end == curr_segment_end && curr_idx < last_idx)
        {
            ++curr_idx;
            curr_segment_start = coordinates[curr_idx].first;
            curr_segment_end = coordinates[curr_idx].second;
            tmp.emplace_back('N', (curr_segment_start - (curr_op_end + 1)));
            curr_pos = curr_segment_start - 1;
            curr_op_end = curr_pos;

            was_skip = true;
        }
        else
        {
            was_skip = false; 
        }
   }
   
   //make_skip_after_ins(tmp);
   std::swap(tmp, cigar_vector); 
}


void parse_variants(
    const int aln_start, 
    std::string_view ref_seq, 
    std::string_view read_seq,
    std::string_view base_qualities, 
    const std::vector<std::pair<char, int>>& cigar_vector,
    const std::unordered_map<int, std::string_view>& ref_dict,
    std::vector<Variant>& variants, 
    std::string& non_ref_quals
)
{
    if (read_seq == ref_seq) return;

    char op = '\0';
    int op_len = 0, ref_idx = 0, read_idx = 0, pos = aln_start;

    std::string_view ref, alt;
    bool is_padding_base_supported = false;
     
    for (const auto & cigar : cigar_vector)
    {
        op = cigar.first;
        op_len = cigar.second;

        switch (op)
        {
            case 'M':
            case 'X':
            case '=':
                for (int i = 0; i < op_len; ++i)
                {
                    ref = ref_seq.substr(ref_idx, 1);
                    alt = read_seq.substr(read_idx, 1);
                    if (ref != alt)
                    {
                        variants.emplace_back(
                            pos, 
                            std::string(ref), 
                            std::string(alt)
                        ); 
                        non_ref_quals += std::string(
                            base_qualities.substr(read_idx, 1)
                        );
                    }
                    ++ref_idx;
                    ++read_idx;
                    ++pos;
                }

                is_padding_base_supported = true;
                break;
            case 'I':
                if (is_padding_base_supported) 
                {
                    ref = ref_seq.substr(ref_idx - 1, 1);
                }
                else
                {   
                    auto it = ref_dict.find(pos - 1);
                    if (it != ref_dict.cend()) 
                    {
                        ref = it->second;
                    }
                    else ref = "N";        
                } 
                
                variants.emplace_back(pos - 1, 
                    std::string(ref), 
                    std::string(ref) 
                    + std::string(read_seq.substr(read_idx, op_len))
                );
                non_ref_quals += std::string(
                    base_qualities.substr(read_idx, op_len)
                ); 
                read_idx += op_len;
                break; 
            case 'D':
                if (is_padding_base_supported) 
                {
                     alt = ref_seq.substr(ref_idx - 1, 1);
                }
                else
                {
                    auto it = ref_dict.find(pos - 1);
                    if (it != ref_dict.cend()) 
                    {    
                        alt = it->second;
                    }
                    else alt = "N";
                }
                
                variants.emplace_back(pos - 1, 
                    std::string(alt) 
                    + std::string(ref_seq.substr(ref_idx, op_len)), 
                    std::string(alt)
                );
                ref_idx += op_len;
                pos += op_len;
                is_padding_base_supported = true;                 
                break;
            case 'N':
                pos += op_len;
                is_padding_base_supported = false;
                break;
            case 'S':
                non_ref_quals += std::string(
                    base_qualities.substr(read_idx, op_len)
                );
                read_idx += (op_len);
                break;
            case 'H':
            case 'P':
                break;
         }
     }
}

// merge ins and del at same pos to a cplx indel
// variants are sorted ins first (by move_up_insertion)
//--------------------------------------------------------------------
std::vector<Variant> merge_to_cplx(const std::vector<Variant>& variants)
{
    std::vector<Variant> merged;
    for (size_t i = 0; i < variants.size();)
    {
        if (i < variants.size() - 1 && variants[i].pos == variants[i + 1].pos)
        {
            if (variants[i].is_ins && variants[i + 1].is_del)
            {    
                merged.emplace_back(
                    variants[i].pos,
                    variants[i + 1].ref,
                    variants[i].alt
                );
                i +=2;
                continue;
            }
        }
        
        merged.push_back(variants[i]);  
        ++i; 
    }
    return merged;
}


// expand segment start/end: {{123, 125}, {502, 504}} -> {123, 124, 125, 502, 503, 504}
//-------------------------------------------------------------------------------------
std::vector<int> expand_coordinates(const Coord& coordinates)
{
    std::vector<int> expanded = {};
    for (const auto & coord_pair : coordinates)
    {
        for (int i = coord_pair.first; i <= coord_pair.second; ++i)
        {
            expanded.push_back(i);
        }
    }
    return expanded;    
}


std::string_view find_commonest_str(const std::vector<std::string_view>& v)
{
    std::unordered_map<std::string_view, int> freq;
    for (const auto& elem : v) 
    {
        freq[elem]++;
    }
                                
    auto commonest = std::max_element(
        freq.begin(), freq.end(), [] (const auto& x, const auto& y) 
        {
            return x.second < y.second;
        }
    );

   return commonest->first;
}


bool has_this(std::string_view query, std::string_view target_ptrn)
{
    return (query.find(target_ptrn) != std::string_view::npos);
}


bool has_gaps(std::string_view cigar_string)
{   
    return (cigar_string.find('I') != std::string_view::npos 
         || cigar_string.find('D') != std::string_view::npos);
}


int count_repeats(std::string_view ptrn, std::string_view seq)
{
    int ptrn_len = ptrn.size(), seq_len = seq.size();
    int i = 0, n = 0;
    while (ptrn_len * i < seq_len) 
    {
        if (ptrn == seq.substr(ptrn_len * i, std::min(ptrn_len, seq_len - (ptrn_len * i)))) 
        {
            ++n;
            ++i;
        }
        else return n;
    }

    return n;
}


void make_kmers(std::string_view seq, const size_t k, Kmers& kmers)
{
    size_t n = seq.size();
    
    if (n <= k) 
    {
        kmers.insert(seq);
        return;
    }
    
    for (size_t i = 0; i <= n - k; ++i) {
        kmers.insert(seq.substr(i, k));
    }
}


void diff_kmers(
    std::string_view query, 
    std::string_view subject, 
    const size_t k,
    Kmers& diff
)
{
    Kmers qry_kmers, sbj_kmers; 
    make_kmers(query, k, qry_kmers);
    make_kmers(subject, k, sbj_kmers);

    std::set_difference(qry_kmers.begin(), qry_kmers.end(), 
                        sbj_kmers.begin(), sbj_kmers.end(),
                        std::inserter(diff, diff.end()));
}


int count_kmer_overlap(std::string_view seq, const Kmers& kmer_set)
{
    int kmer_cnt = 0;
    //const char* _seq = std::string(seq).c_str();
    for (const auto& kmer : kmer_set)
    {
        //if(strstr(_seq, std::string(kmer).c_str())) ++kmer_cnt;
        if(seq.find(kmer) != std::string_view::npos) ++kmer_cnt;
        //if(strstr(std::string(seq).c_str(), std::string(kmer).c_str())) ++kmer_cnt;
    }

    return kmer_cnt; 
}


//split at cplx variant start index
int find_split_idx(
    const int read_start, 
    const int target_pos,
    const std::vector<std::pair<char, int>>& cigar_vector
)
{
    char op = '\0';
    int i = 0, op_len = 0, curr_pos = read_start - 1;
    for (const auto& c: cigar_vector)
    {
        op = c.first;
        op_len = c.second;
        switch (op)
        {
            case 'M':
            case 'X':
            case 'S':
                if (curr_pos + op_len < target_pos)
                {    
                    curr_pos += op_len;
                    i += op_len;    
                }
                else if (curr_pos < target_pos && target_pos <= curr_pos + op_len)
                {
                    i += target_pos - curr_pos;
                    return i;
                }
                break;
            case 'I':
                i += op_len;
                curr_pos += 1;
                break;
            case 'D':
            case 'N':
                curr_pos += op_len;
                break;
            default:
                break;      
        }      
    }
    
    return -1; //invalid case
}

