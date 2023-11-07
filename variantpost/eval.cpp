#include <cmath>
#include <string>
#include <vector> 
#include <climits>
#include <utility>
#include <string_view>

#include "eval.h"
#include "util.h"
#include "merge.h"
#include "contig.h"


void local_alignment(
    const uint8_t match_score,
    const uint8_t mismatch_penal,
    const uint8_t gap_open_penal,
    const uint8_t gap_ext_penal,
    const std::string& seq,
    const std::string& ref_seq,
    const Filter& filter,
    Alignment& alignment
)
{
    Aligner aligner(
        match_score, mismatch_penal, gap_open_penal, gap_ext_penal
    ); 
       
    int32_t mask_len = strlen(seq.c_str()) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;

    aligner.Align(
        seq.c_str(), ref_seq.c_str(), ref_seq.size(),
        filter, &alignment, mask_len
    ); 
}

/*
bool is_gapped_aln(std::string_view cigar_str)
{
    if (cigar_str.size()) return has_gaps(cigar_str);
    else return false;
}*/


void to_simple_variants(
    const UserParams& user_params,
    const Contig& contig,
    LocalReference& loc_ref,
    std::vector<Variant>& simples
)
{
    Filter filter;
    Alignment alignment;
    local_alignment(
        user_params.match_score, user_params.mismatch_penal, 
        user_params.gap_open_penal, user_params.gap_ext_penal,
        contig.mocked_seq, contig.mocked_ref, filter, alignment
    );

    std::vector<std::pair<char, int>> cigar_vec = to_cigar_vector(
        alignment.cigar_string
    );
    
    move_up_insertion(cigar_vec);
    
    std::vector<int> pos_vec = expand_coordinates(contig.mocked_coord); 
    
    int query_begin = alignment.query_begin;
    if (cigar_vec[0].first == 'S')
    {
        query_begin -= cigar_vec[0].second;
    }
    
    //bool include_clips = false;   
    std::string _tmp = ""; //won't need
    parse_variants(
        pos_vec[alignment.ref_begin],
        contig.mocked_ref.substr(alignment.ref_begin),
        contig.mocked_seq.substr(query_begin),
        contig.mocked_seq.substr(query_begin), //for mocking
        cigar_vec, loc_ref.dict, simples, _tmp
    );

    std::sort(simples.begin(), simples.end());
}


void postprocess_alignment(
    AlnResult& rslt, 
    const std::vector<int>& pos_vec,
    const Contig& contig,
    const LocalReference& loc_ref,
    const Alignment& alignment)
{
    rslt.cigar_str = alignment.cigar_string;
    rslt.cigar_vec = to_cigar_vector(rslt.cigar_str);
    
    splice_cigar(
        rslt.cigar_vec, alignment.ref_begin, 
        pos_vec, contig.coordinates
    );
    
    move_up_insertion(rslt.cigar_vec);
   
    int query_begin = alignment.query_begin;
    if (rslt.cigar_vec[0].first == 'S')
    {      
        query_begin -= rslt.cigar_vec[0].second;
    }
    
    rslt.genomic_start_pos = pos_vec[alignment.ref_begin];
    rslt.seq = contig.seq.substr(query_begin);
    rslt.ref_seq = contig.ref_seq.substr(alignment.ref_begin);
    rslt.quals = contig.quals.substr(query_begin);

    std::string _tmp = "";
    parse_variants(
        rslt.genomic_start_pos, rslt.ref_seq,
        rslt.seq, rslt.quals, rslt.cigar_vec, 
        loc_ref.dict, rslt.variants, _tmp
    );
}


size_t rt_most_idx(
    const size_t curr_idx, 
    const size_t ins_len,
    std::string_view seq
)
{
    //size_t rt_most_idx = curr_idx;
    for (size_t i = curr_idx; i + ins_len < seq.size(); ++i)
    {
        if (seq[i] != seq[i + ins_len]) return i + ins_len - 1;     
    }

    return curr_idx;
}


int get_rt_pos(
    const size_t curr_seq_idx, 
    const size_t target_idx,
    const size_t curr_cigar_idx,
    const int curr_pos,
    const CigarVec& cigar_vec        
)
{
    char op = '\0';
    size_t seq_idx = curr_seq_idx;
    int op_len = 0, pos = curr_pos;
    for (size_t i = curr_cigar_idx; i < cigar_vec.size(); ++i)
    {
        op = cigar_vec[i].first;
        op_len = cigar_vec[i].second;
        
        switch (op)
        {
            case '=':
            case 'X':
            {
                seq_idx += op_len;
                pos += op_len;
                break;
            }
            case 'I':
            {
                seq_idx += op_len;
                break;               
            }
            case 'D':
            case 'N':
            {
                pos += op_len;
                break;           
            }
        }
        
        if (seq_idx == target_idx)
        {
            if (i + 1 < cigar_vec.size())
            {
                if (cigar_vec[i+1].first != 'D' && cigar_vec[i+1].first != 'N')   
                {
                    return pos;
                }
            }
        }
    }
    
    return 0;   
}


void rt_aln_insertions(
    const AlnResult& rslt,
    std::vector<Variant>& varlst,
    LocalReference& loc_ref
)
{
   char op = '\0';
   size_t seq_idx = 0;
   int op_len = 0, curr_pos = rslt.genomic_start_pos;
   for (size_t i = 0; i < rslt.cigar_vec.size(); ++i)
   {   
        op = rslt.cigar_vec[i].first;
        op_len = rslt.cigar_vec[i].second;

        switch (op)
        {
            case '=':
            case 'X':
            {
                seq_idx += op_len;
                curr_pos += op_len;       
                break;
            }
            case 'I':
            {
                //int _curr_pos = curr_pos - 1;
                //int _seg_idx = seq_idx - 1;
                size_t rt_idx = rt_most_idx(seq_idx, op_len, rslt.seq);
                if (seq_idx + op_len < rt_idx)
                {
                    int _pos = get_rt_pos(
                        seq_idx, rt_idx, i, curr_pos, rslt.cigar_vec
                    );
                    
                    if (_pos)
                    {
                        std::string ins_seq = rslt.seq.substr(
                            rt_idx - static_cast<size_t>(op_len),  
                            op_len + 1
                        );
                        
                        varlst.emplace_back(_pos, ins_seq.substr(0, 1), ins_seq);
                    }    
                }
                seq_idx += op_len;
                break;   
            } 
            case 'D':
            case 'N':
            {    
                curr_pos += op_len;
                break;
            }
            case 'S':
            {
                seq_idx += op_len;
                break;
            }
        }
   }
}


void swappable_gaps(
    const std::vector<Variant>& varlst,
    std::vector<Variant>& augmented,
    LocalReference& loc_ref
)
{
    // varlst.size() >= 2 guaranteed
    for (size_t i = 0; i < varlst.size() - 1; ++i)
    {
        if (varlst[i].is_ins 
            && varlst[i + 1].is_del 
            && varlst[i].pos < varlst[i + 1].pos
        )
        {
            size_t mapped_len = varlst[i + 1].pos - varlst[i].pos;
            size_t del_len = varlst[i + 1].ref.size();
            if (mapped_len < del_len)
            { 
                std::string inserted = varlst[i].alt.substr(1) 
                    + loc_ref.fasta.getSubSequence(
                        loc_ref.chrom, varlst[i].pos, mapped_len);

                if (inserted.substr(0, mapped_len) 
                    == varlst[i + 1].ref.substr(del_len - mapped_len))
                {
                    int _ins_pos = varlst[i + 1].variant_end_pos - mapped_len;
                    augmented.emplace_back(
                        _ins_pos, inserted.substr(0, 1), inserted);
                    
                    int _del_pos = varlst[i + 1].pos - mapped_len;
                    std::string _del_seq = loc_ref.fasta.getSubSequence(
                        loc_ref.chrom, _del_pos - 1, del_len); 
                    augmented.emplace_back(
                        _del_pos, _del_seq, _del_seq.substr(0, 1));   
                } 
            }
        }
    }
}


void create_target_ins(
    AlnResult& rslt, 
    std::vector<Variant>& varlst,
    const int target_pos,
    const size_t event_len
)
{
    char op = '\0';
    int op_len = 0, seq_idx = 0, curr_pos = rslt.genomic_start_pos;
    
    for (const auto& c : rslt.cigar_vec)
    {
        op = c.first;
        op_len = c.second;
         
        switch (op)
        {
            case '=':
            case 'X':
            {
                if (curr_pos <= target_pos && target_pos < curr_pos + op_len)
                {
                    seq_idx += (target_pos - curr_pos); 
                    
                    varlst.emplace_back(
                        target_pos, rslt.seq.substr(seq_idx, 1), rslt.seq.substr(seq_idx, event_len)
                    );
                    return;     
                }
                else
                {
                    seq_idx += op_len;
                    curr_pos += op_len;
                }
                break;
             }
             case 'I':
             case 'S':
                 seq_idx += op_len;
                 break;
             case 'D':
             case 'N':
                curr_pos += op_len;
                break;       
        }
    }   
}


void gap_penal_grid(
    const int gap_open_penal, 
    const int gap_extension_penal, 
    std::vector<std::pair<int, int>>& grid
)
{
    std::pair<int, int> _default = {gap_open_penal, gap_extension_penal};

    uint8_t max_gap_open = 10, min_gap_open = 2; //configureab;e
    if (gap_open_penal > max_gap_open) max_gap_open = gap_open_penal;
    
    for (uint8_t i = min_gap_open; i <= max_gap_open; i += 2)
    {      
        grid.emplace_back(i, 1);
        grid.emplace_back(i, 0);
    }
    
    if (std::find(grid.begin(), grid.end(), _default) ==  grid.end())
    {
        grid.push_back(_default);
    }
}


void penal_grid(
    const int mismatch_penal,
    const int gap_open_penal,
    const int gap_ext_penal,
    size_t & basic_grid_sz,
    std::vector<std::vector<int>>& grid  
)
{
    std::vector<int> _default = {mismatch_penal, gap_open_penal, gap_ext_penal};

    uint8_t max_gap_open = 10, min_gap_open = 2; //configureable?
    if (gap_open_penal > max_gap_open) max_gap_open = gap_open_penal;
    
    for (uint8_t i = min_gap_open; i <= max_gap_open; i += 2)
    {
        grid.push_back({mismatch_penal, i, 1});
        grid.push_back({mismatch_penal, i, 0});
    }

    if (std::find(grid.begin(), grid.end(), _default) ==  grid.end()) 
    {
        grid.push_back(_default);
    }

    basic_grid_sz = grid.size();

    //encourage gaps over mismathces
    for (uint8_t i = min_gap_open; i <= max_gap_open; i += 2)
    {
        grid.push_back({mismatch_penal * 5, i, 1});
        grid.push_back({mismatch_penal * 5, i, 0});
    }
}


bool has_mapped_ends(
    const UserParams& user_params,
    const std::vector<std::pair<char, int>>& cigar_vec
)
{
    size_t last_idx = cigar_vec.size() - 1;
    if (cigar_vec[0].first == '='
        && cigar_vec[0].second >= user_params.local_thresh
        && cigar_vec[last_idx].first == '='
        && cigar_vec[last_idx].second >= user_params.local_thresh)
    {
        return true;
    }
    else return false;
}


void annot_end_mappping(AlnResult& rslt, const UserParams& user_params)
{
    size_t last_idx = rslt.cigar_vec.size() - 1;
    if (rslt.cigar_vec[0].first == '=')
    {
        rslt.lt_ref_mapped_len = rslt.cigar_vec[0].second;
        rslt.lt_well_ref_mapped = (
            rslt.lt_ref_mapped_len >= user_params.local_thresh
        ) ? true : false;
    }

    if (rslt.cigar_vec[last_idx].first == '=')
    {
        rslt.rt_ref_mapped_len = rslt.cigar_vec[last_idx].second; 
        rslt.rt_well_ref_mapped = (
            rslt.rt_ref_mapped_len >= user_params.local_thresh
        ) ? true : false;   
    }

    rslt.is_well_ref_mapped = (
        rslt.lt_well_ref_mapped && rslt.rt_well_ref_mapped
    ) ? true : false;
}


void annot_matched_segments(Coord& coord, AlnResult& rslt)
{
    char op = '\0';
    int curr_pos = rslt.genomic_start_pos, op_len = 0;
    for (const auto& c : rslt.cigar_vec)
    {
        op = c.first;
        op_len = c.second;
        if (op == '=')
        {
            coord.emplace_back(curr_pos, curr_pos + op_len -1);
            curr_pos += op_len;
        }
        else if (op == 'D' || op == 'N')
        {
            curr_pos += op_len;
        }
    }

}


bool has_enough_margin(Coord& coord, const int pos, const UserParams& user_param)
{
    bool is_enough_lt = false, is_enough_rt = false, is_first_rt = true;
    for (const auto& seg : coord)
    {
        
        if (seg.second <= pos)
        { 
            if (seg.second - seg.first + 1 >= user_param.local_thresh) 
            {    
                is_enough_lt = true;
            }
        }
        else if (pos < seg.first)
        {
            if (is_first_rt)
            {
                if (seg.second - seg.first + 1 >= user_param.local_thresh) 
                {    
                    is_enough_rt = true;
                }
            }
        }   
    }

    return (is_enough_lt && is_enough_rt);
}


int dist_from_target(
    const int lt_lim, 
    const int target_pos,
    const int rt_lim, 
    const LocalReference& loc_ref,
    Variant& variant
)
{ 
    int d;
    if (variant.pos < target_pos)
    {
        variant.set_rightmost_pos(loc_ref);
        int lt_d =  lt_lim - (variant.rpos + variant.ref_len - 1);   
        d = (lt_d < 0) ? 0 : lt_d;
    }
    else
    {
        variant.set_leftmost_pos(loc_ref);
        int rt_d = variant.lpos - rt_lim;
        d = (rt_d < 0) ? 0 : rt_d;
    } 

    return d;
}


void annot_nearest_variant(
    const int lt_lim,
    const int target_pos,
    const int rt_lim, 
    const LocalReference& loc_ref,
    std::vector<Variant>& variants,
    int& idx,
    int& dist
)
{
    int tmp_dist = INT_MAX;
    for (size_t i = 0; i < variants.size(); ++i)
    {
        tmp_dist = dist_from_target(
            lt_lim, target_pos, rt_lim, loc_ref, variants[i]
        );

        if (tmp_dist < dist)  
        {
            dist = tmp_dist;
            idx = i;
        } 
    }
}


inline bool pass_gap_check(AlnResult& rslt, const UserParams& user_params)
{
    if (rslt.variants.empty()) return false;
    
    annot_end_mappping(rslt, user_params);
    if (rslt.is_well_ref_mapped && !has_gaps(rslt.cigar_str))
    {
        return false;
    }
    
    return true;
}


bool is_overlapping_cluster(
    const Variant& target,
    const std::vector<Variant>& varlst,
    const UserParams& user_params
)
{
    int cnt = 0;
    bool is_overlapping = false;
    const int local_radius = user_params.local_thresh / 2;
    for (const auto& var : varlst)
    {
        if (var.pos <= target.pos)
        {
            if (target.pos < var.variant_end_pos)
            {
                ++cnt;
                is_overlapping = true;     
            }
            else
            {
                if (target.pos - var.variant_end_pos < local_radius)
                {
                    ++cnt;
                    if (target.pos == var.pos) is_overlapping = true;
                }
            }
        }
        else
        {
            if (var.pos < target.variant_end_pos)
            {
                ++cnt;
                is_overlapping = true;
            }
            else
            {    
                if (var.pos - target.pos < local_radius) ++cnt;
            }
        }  
    }

    if (cnt > 1 && is_overlapping) return true;
    else return false;
}


void eval_by_variant(
    AlnResult& rslt,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref,
    int& alternative_pos
)
{   
    if (!pass_gap_check(rslt, user_params)) return;
    rslt.is_passed = true;
    
    bool is_orig_aln = false;
    for (auto& v : rslt.variants)
    {
        if (target.is_equivalent(v, loc_ref))
        {
            rslt.has_target = true;
            is_orig_aln = true;
            break;
        }
    }
    
    //augment variant list for complex cases
    if (!rslt.has_target && rslt.variants.size() > 1)
    {
        //merge consecutive I/Ds
        std::vector<Variant> augmented = merge_to_cplx(rslt.variants);
        
        //right alignment on contig 
        rt_aln_insertions(rslt, augmented, loc_ref);
        
        // merge I/D*X to complx
        std::vector<int> target_idxes 
        = find_mismatches_after_gap(rslt.cigar_vec);
        if (!target_idxes.empty())
        {
            parse_to_cplx_gaps(
                rslt.genomic_start_pos, rslt.ref_seq, rslt.seq, 
                rslt.cigar_vec, target_idxes, augmented
            );
        }
        
        // test if the expectd ins-seq recovered
        //create_target_ins(rslt, augmented, target.pos, target.alt.size());
        
        //swap without chaning the number of mapped bases
        swappable_gaps(rslt.variants, augmented, loc_ref);
        
        for (auto& v : augmented)
        {
            if (target.is_equivalent(v, loc_ref))
            {
                 rslt.has_target = true;
                 is_orig_aln = false;
                 break;
            }
        }

        /* keep this experimental 
        if (!rslt.has_target)
        {
            //cluster match
            bool looks_good = is_overlapping_cluster(target, rslt.variants, user_params);
        
            if (looks_good)
            { 
                rslt.has_target = true;
                is_orig_aln = false;
            }
        }*/
    }    
    
    if (rslt.has_target)
    {
        if (is_orig_aln)
        {
            Coord mapped_seg;
            annot_matched_segments(mapped_seg, rslt);
            bool is_enough = has_enough_margin(
                mapped_seg, target.pos, user_params
            );
        
            if (is_enough) rslt.terminate_search = true;
        }
        else
        {
            annot_nearest_variant(
                target.lpos, target.pos, target.rpos + target.ref_len - 1,
                loc_ref, rslt.variants, rslt.idx_to_closest, rslt.dist_to_closest
            );

            alternative_pos = rslt.variants[rslt.idx_to_closest].pos;
            if (rslt.is_well_ref_mapped) rslt.terminate_search = true;
        }
    }
} 


void eval_by_variant_lst(
    AlnResult& rslt,
    const std::vector<Variant>* p_decomposed,
    const UserParams& user_params,
    LocalReference& loc_ref,
    int& alternative_pos 
)
{
    if (!pass_gap_check(rslt, user_params)) return;

    std::vector<Variant> shared;
    find_shared_variants(shared, rslt.variants, *p_decomposed);
   
    size_t total_match = shared.size(), indel_match = 0;
    
    //find indel with biggest change
    int alternative_indel_len = 0;
    for (const auto& v : shared)
    {
        if (!v.is_substitute)
        {
            ++indel_match;
            if (alternative_indel_len <= v.ref_len + v.alt_len)    
            {    
                alternative_pos = v.pos;
                alternative_indel_len = (v.ref_len + v.alt_len);
            }
         }
    }
    
    if (!total_match) return;
    
    if (total_match == p_decomposed->size())
    {
        rslt.has_target = true;
        rslt.is_passed = true;
        rslt.terminate_search = true;
    } 
    else if (!indel_match)
    {   
        rslt.is_passed = true;
    }
}


void annot_alignment(Contig& contig, const AlnResult& rslt)
{
    int genomic_pos = rslt.genomic_start_pos; 

    contig.seq = rslt.seq;
    contig.len = contig.seq.size();
    contig.ref_seq = rslt.ref_seq;
    contig.quals = rslt.quals;
     
    int op_len;
    char op = '\0';
    size_t seq_idx = 0, ref_idx = 0, v_idx = 0;    
    for (const auto& cigar : rslt.cigar_vec)
    {  
        op = cigar.first;
        op_len = cigar.second;
        
        switch (op)
        {
            case '=':
            case 'X':
                for (int i = 0; i < op_len; ++i)
                {                   
                    contig.pos_by_seq_idx.push_back(genomic_pos);
                    contig.positions.push_back(genomic_pos);
                    contig.ref_bases.push_back(contig.ref_seq.substr(ref_idx, 1));
                    contig.alt_bases.push_back(contig.seq.substr(seq_idx, 1));
                    contig.base_quals.push_back(contig.quals.substr(seq_idx, 1));
                    
                    ++genomic_pos;
                    ++ref_idx;
                    ++seq_idx;

                    if (op == 'X') 
                    {    
                        ++v_idx;
                    }
                }
                break;
            case 'I':
                contig.alt_bases.pop_back();
                contig.alt_bases.push_back(contig.seq.substr(seq_idx - 1, op_len + 1));
                contig.base_quals.pop_back();
                contig.base_quals.push_back(contig.quals.substr(seq_idx - 1, op_len + 1));
                
                std::fill_n(std::back_inserter(contig.pos_by_seq_idx), op_len, genomic_pos - 1);
                
                seq_idx += op_len;
                ++v_idx;
                break;
            case 'D':
                contig.ref_bases.pop_back();
                contig.ref_bases.push_back(rslt.variants[v_idx].ref);

                ref_idx += op_len;
                genomic_pos += op_len;
                ++v_idx;
                break;
            case 'N':
                contig.skip_starts.push_back(genomic_pos);
                genomic_pos += op_len;
                contig.skip_ends.push_back(genomic_pos - 1);

                contig.pos_by_seq_idx.pop_back();
                break;
            case 'S':
                std::fill_n(std::back_inserter(contig.pos_by_seq_idx), op_len, -1);
                seq_idx += op_len;
                break;
        }
    }
}


void update_contig_layout(Contig& contig, const int pos)
{
    auto it = std::find(
        contig.positions.begin(), 
        contig.positions.end(), 
        pos
    );
    
    size_t idx;
    //keep this block for furtehr cheking
    if (it == contig.positions.end()) 
    {
        int diff = INT_MAX;
        size_t tmp_idx = 0;
        for (size_t i = 0; i < contig.positions.size(); ++i)
        {
            if (contig.ref_bases[i].size() != contig.alt_bases[i].size()
                && diff > std::abs(contig.positions[i] - pos))
            { 
                diff = std::abs(contig.positions[i] - pos);
                tmp_idx = i;
            } 
        }
        
        idx = tmp_idx;    
    }    
    else
    { 
        idx = std::distance(contig.positions.begin(), it);
    } 

    int offset = 0;
    for (size_t i = 0; i < contig.pos_by_seq_idx.size(); ++i)
    {
        if (contig.pos_by_seq_idx[i] == -1) ++offset;
        else break;
    }
    
    contig.lt_len = offset;
    contig.mid_len =0;
    for (size_t i = 0; i <= idx; ++i)
    {
        if (i == idx) 
        {
            size_t _len = contig.alt_bases[i].size();
            if (_len > 1) contig.mid_len = (_len - 1);
            contig.lt_len += 1;
        }
        else
        {
            contig.lt_len += contig.alt_bases[i].size();
        }
    }
    
    contig.lt_end_idx = contig.lt_len - 1;
    contig.rt_start_idx = contig.lt_len + contig.mid_len;
    contig.rt_len = contig.len - contig.rt_start_idx; 
}


bool is_valid_aln(std::string_view cigar_str)
{
    //alignment failed
    if (!cigar_str.size()) return false;

    //some alignments only contains 'S'
    if (!has_this(cigar_str, "=")) return false; 

    return true;
}


int local_matched_len(
    const int target_pos,
    const int aln_start,
    const CigarVec& cigar_vec
)
{
    char op = '\0';
    int op_len = 0, curr_pos = aln_start;
    for (const auto& c : cigar_vec)
    {
        op = c.first;
        op_len = c.second;
        switch (op)
        {
            case '=':
                if (curr_pos <= target_pos && target_pos < curr_pos + op_len)
                {
                    if (target_pos - curr_pos 
                        > curr_pos + op_len - (target_pos + 1))
                    {
                        return curr_pos + op_len - (target_pos + 1); 
                    }
                    else
                    {
                        return target_pos - curr_pos;
                    }
                }
                curr_pos += op_len;
                break;
            case 'N':
            case 'X':
                curr_pos += op_len;
                break;
            default:
                break;
        }
    }

    return -1;
}


size_t closest_mismatch(
    const int target_pos, 
    const int aln_start,
    const CigarVec& cigar_vec
)
{
    size_t i = 0, idx = 0;
    
    char op = '\0';
    int op_len = 0, curr_pos = aln_start, dist = INT_MAX;
    for (const auto& c : cigar_vec)
    {
        op = c.first;
        op_len = c.second;
        switch (op)
        {
            case '=':
            case 'N':
                curr_pos += op_len;
                ++i;
                break;
            case 'X':
                if (std::abs(target_pos - curr_pos) < dist) 
                {
                    idx = i;
                }
                curr_pos += op_len;
                ++i;
                break;
            case 'S':  
                ++i;
                break;
            default:
                ++i;
                break;           
        }
    }
    return idx;   
}


bool is_stable_mismatch(
    const int target_pos,
    const int aln_start,
    const CigarVec& cigar_vec,
    const int local_thresh    
)
{
    bool lt_stable = false, rt_stable = false;
    size_t idx = closest_mismatch(target_pos, aln_start, cigar_vec);
    if (idx > 0 && idx < cigar_vec.size() - 1)
    {
        if (cigar_vec[idx - 1].first == '=' 
            && cigar_vec[idx - 1].second > local_thresh) lt_stable = true;
        if (cigar_vec[idx + 1].first == '='
        && cigar_vec[idx + 1].second > local_thresh) rt_stable = true;
    }
    
    return (lt_stable && rt_stable);
}


char eval_by_aln(
    Contig& contig,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref,
    const std::vector<Variant>* p_decomposed
)
{
    size_t basic_sz = 0;
    std::vector<std::vector<int>> penals;
    penal_grid(
        user_params.mismatch_penal, user_params.gap_open_penal, 
        user_params.gap_ext_penal, basic_sz, penals
    );

    std::vector<int> pos_vec = expand_coordinates(contig.coordinates);

    Filter filter;
    Alignment aln;
    std::vector<AlnResult> rslts;
    for (const auto& penal : penals)
    { 
        AlnResult rslt;
        
        local_alignment(
            user_params.match_score, penal[0], penal[1], penal[2], 
            contig.seq, contig.ref_seq, filter, aln
        );
        
        if (!is_valid_aln(aln.cigar_string)) continue;
        
        if (contig.by_kmer_suggestion && !has_gaps(aln.cigar_string))
        {
            CigarVec cigar_vec = to_cigar_vector(aln.cigar_string);
            int loc_len = local_matched_len(
                target.pos, pos_vec[aln.ref_begin], cigar_vec
            );
             
            if (loc_len > user_params.local_thresh) return 'C';
             
            if (has_this(aln.cigar_string, "X"))
            {
                if (!has_this(aln.cigar_string, "S")) return 'C';
                
                if (
                    is_stable_mismatch(
                        target.pos, pos_vec[aln.ref_begin], 
                        cigar_vec, user_params.local_thresh
                    )
                ) return 'C';
             }   
            
            continue; 
        }

        postprocess_alignment(rslt, pos_vec, contig, loc_ref, aln);
        
        int target_pos = target.pos;
        if (p_decomposed != nullptr)
        {
            eval_by_variant_lst(rslt, p_decomposed, user_params, loc_ref, target_pos);
        }
        else
        { 
            eval_by_variant(rslt, target, user_params, loc_ref, target_pos);
        }

        if (rslt.terminate_search)
        {    
            annot_alignment(contig, rslt);
            update_contig_layout(contig, target_pos);
            return 'A';   
        }
        else if (rslt.is_passed)
        {
            rslts.push_back(rslt);
        }
    }
    
    if (rslts.empty()) return 'C';

    std::vector<std::string_view> aln_ptrns;
    if (contig.by_kmer_suggestion)
    {
        for (const auto& rslt : rslts)
        {
            if (rslt.has_target)
            {
                aln_ptrns.push_back(rslt.cigar_str);
            }
        }

        if (!aln_ptrns.empty())
        {
            std::string_view common_ptrn = find_commonest_str(aln_ptrns);
            auto it = std::find(aln_ptrns.begin(), aln_ptrns.end(), common_ptrn);
            auto rslt = rslts[std::distance(aln_ptrns.begin(), it)];
            annot_alignment(contig, rslt);
            update_contig_layout(contig, target.pos);
            return 'A';
        }
        else
        {
            return 'B';
        }   
    }
    else
    { 
        //most stable aln
        for (const auto& rslt : rslts)
        {
            aln_ptrns.push_back(rslt.cigar_str);
        }    
     
        std::string_view common_ptrn = find_commonest_str(aln_ptrns);
        auto it = std::find(aln_ptrns.begin(), aln_ptrns.end(), common_ptrn);
        auto rslt = rslts[std::distance(aln_ptrns.begin(), it)];

        //spliced -> no extension
        for (const auto& c : rslt.cigar_vec)
        {
            if (c.first == 'N') 
            {    
                annot_alignment(contig, rslt);
                return 'A';
            }
        }
    
        //do extension
        if (!rslt.lt_well_ref_mapped && rslt.rt_well_ref_mapped)
        {
            return 'L';
        }
        else if (rslt.lt_well_ref_mapped && !rslt.rt_well_ref_mapped)
        {
            return 'R';
        }
        else
        {    
            return 'E';
        }
    }
}


void aln_extended_contig( 
    Contig& contig,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
   std::vector<int> pos_vec = expand_coordinates(contig.coordinates);
    
   Filter filter;
   Alignment aln;
   AlnResult rslt; 
    
   local_alignment(
       user_params.match_score, user_params.mismatch_penal,
       user_params.gap_open_penal, user_params.gap_ext_penal,
       contig.seq, contig.ref_seq, filter, aln
   );

   postprocess_alignment(rslt, pos_vec, contig, loc_ref, aln);
   annot_alignment(contig, rslt);
   update_contig_layout(contig, target.pos);
}


void switch_to_mock(
    Contig& contig,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref,
    bool& is_mocked,
    const std::vector<Variant>* p_decomposed
)
{
    contig.seq = contig.mocked_seq;
    contig.quals = contig.m_quals;
    contig.ref_seq = contig.mocked_ref;
    contig.coordinates = contig.mocked_coord;
    contig.len = contig.mock_len;
    contig.lt_end_idx = contig.mock_lt_end_idx;
    contig.lt_len = contig.mock_lt_len;
    contig.mid_len = contig.mock_mid_len;
    contig.rt_len = contig.mock_rt_len;
    contig.rt_start_idx = contig.mock_rt_start_idx;
    
    char res = eval_by_aln(contig, target, user_params, loc_ref, p_decomposed); 
    
    if (res == 'A')
    {
        is_mocked = true;
        contig.is_mocked = true;
    }
}
