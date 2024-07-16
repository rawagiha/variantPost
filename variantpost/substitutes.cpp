#include "eval.h"
#include "util.h"
#include "reads.h"
#include "contig.h"


void substitute_patterns(
    int& a_cnt,
    int& b_cnt,
    Read& read, 
    const Variant& target,
    const UserParams& user_params
)
{
    if (read.covering_ptrn != 'A') 
    {   
        read.sb_ptrn = 'C';        
        return;
    }

    char op = '\0';
    int op_len = 0, curr_pos = read.read_start;
    size_t ref_idx = 0, read_idx = 0;
    for (const auto& c : read.cigar_vector)
    {
        op = c.first;
        op_len = c.second;

        switch (op)
        {
            case 'M':
            case 'X':
            case '=':
                if (curr_pos <= target.pos && target.pos < curr_pos + op_len)
                {
                    std::string_view query = read.seq.substr(
                        read_idx + (target.pos - curr_pos), target.alt.size());
                    
                    if (query == target.ref)
                    {
                        read.sb_ptrn = 'C';
                        return;
                    }

                    std::string quals = read.base_quals.substr(
                        read_idx + (target.pos - curr_pos), target.alt.size());
                    
                    size_t lq_base_n = 0;
                    for (size_t i = 0; i < quals.size(); ++i)
                    {
                         if (quals[i] < user_params.base_q_thresh)
                         {
                            ++lq_base_n;
                         }
                    }

                    if (lq_base_n == quals.size())
                    {
                        read.sb_ptrn = 'U';
                        return; 
                    }
                    
                    //std::string_view query = read.seq.substr(
                    //    read_idx + (target.pos - curr_pos), target.alt.size()); 
                      
                    if (query == target.alt)
                    {
                        read.sb_ptrn = 'A';
                        read_idx += (target.pos - curr_pos);
                        read.sb_read_idx = read_idx;
                        ++a_cnt;
                        return;
                    }
                    else read.sb_ptrn = 'C';
                        
                    return; 
                    /*{    
                        std::string quals = read.base_quals.substr(
                            read_idx + (target.pos - curr_pos), target.alt.size());
                        
                        bool is_hq = true;
                        for (size_t i = 0; i < quals.size(); ++i)
                        { 
                            if (quals[i] < user_params.base_q_thresh)
                            {
                                is_hq = false;
                                break;
                            }
                        }
                        
                        if (is_hq)
                        {
                            read.sb_ptrn = 'C';
                        }
                        else
                        {
                            read.sb_ptrn = 'U';
                        }
                        return;
                    }*/
                }
                
                ref_idx += op_len;
                read_idx += op_len;
                curr_pos += op_len;
                break;
            case 'I':
                read_idx += op_len;
                break;     
            case 'D':
                if (curr_pos <= target.pos && target.pos < curr_pos + op_len)
                {
                    read.sb_ptrn = 'B';
                    ++b_cnt;
                    return;
                }
                ref_idx += op_len;
                curr_pos += op_len;
                break;     
            case 'S':
                if (curr_pos <= target.pos && target.pos < curr_pos + op_len) 
                {    
                    read.sb_ptrn = 'B';
                    ++b_cnt;
                    return;
                }
                read_idx += op_len;
                curr_pos += op_len; 
                break;
            case 'N':
                ref_idx += op_len;
                curr_pos += op_len;
                break;           
            default :
                break;
        }
    }    
    
    if (read.nonref_lq_rate < user_params.lq_rate_thresh
        && read.seq.find('N') == std::string_view::npos)
    {
        read.sb_ptrn = 'C';
    }
    else
    {  
        read.sb_ptrn = 'D';
    }
}   


void collect_gaps(
    const int target_pos,
    const int retarget_thresh,
    const LocalReference& loc_ref,
    Reads& reads, 
    std::unordered_map<Variant, int>& cis_gaps
)
{
    for (auto& read : reads)
    {
        if (read.variants.empty()) continue;
        
        if (read.sb_ptrn == 'A')
        {
            for (auto& v : read.variants)
            {
                if (v.ref.size() == v.alt.size()) continue;
            
                v.set_leftmost_pos(loc_ref);
                v.set_rightmost_pos(loc_ref);

                if (v.rpos + v.ref_len <= target_pos)
                {
                    if (target_pos - (v.rpos + v.ref_len) < retarget_thresh)
                    {
                         ++cis_gaps[v];
                    }
                }
                else if(v.lpos >= target_pos)
                {
                    if (v.lpos - target_pos < retarget_thresh)
                    {
                         ++cis_gaps[v];
                    }
                } 
            }        
        }
    } 
}


void commonest_cis_gap(
    Variant& _gap,
    const std::unordered_map<Variant, int>& cis_gaps
)
{
    int cnt = 0;
    for (const auto& gap : cis_gaps)
    {
        if (gap.second > cnt) 
        {
            _gap = gap.first;
            cnt = gap.second;
        }
        
        //TODO tie case
    }
}


void is_exact_cis(
    const Reads& reads,
    const std::unordered_map<Variant, int>& cis_gaps,
    const bool has_no_cis_gaps,
    Variant& _gap,
    bool& is_excl, 
    bool& is_loc_uniq
)
{
    if (!has_no_cis_gaps)
    {
        commonest_cis_gap(_gap, cis_gaps);
    }

    int cnt = 0;
    bool is_not_excl = false;
    for (const auto& read : reads)
    {
        if (read.sb_ptrn == 'B') continue;
        
        if (read.sb_ptrn == 'C')
        {    
            if (is_not_excl || has_no_cis_gaps) continue;

            if (
                std::find(read.variants.begin(), read.variants.end(), _gap) 
                != read.variants.end()
            ) is_not_excl = true;
        }
        else if (read.sb_ptrn == 'A')
        {
         // for 'A' sb_prtn, count reads covering the retarget candidate locus   
         // generally cnt >= cis_gaps.at(_gap) due to sequencing error in inserted seq
         // in candidate
            
            if (!is_loc_uniq) is_loc_uniq = (read.local_uniqueness != -1);
            
            if (has_no_cis_gaps) continue;
           
            if (read.aln_start <= _gap.pos 
                && _gap.variant_end_pos <= read.aln_end) ++cnt;
        }
    }
    
    if (!is_not_excl && !has_no_cis_gaps) 
    {    
        is_excl = (cnt < 2 * cis_gaps.at(_gap));
    } 
}


void annot_loc_variant_siganature(
    Read& read, 
    const int target_pos, 
    const int retarget_thresh
)
{
    std::string sig = "local_vars:";
    for (const auto& variant : read.variants)
    {
        if (std::abs(variant.pos - target_pos) < retarget_thresh)
        {    
            sig += (
                std::to_string(variant.pos) 
                + "_"
                + variant.ref
                + "_"
                + variant.alt
                + ";"
            );
        }
    }
    read.sb_loc_signature = sig;
}


void collect_non_ref_ptrns(
    std::vector<std::string_view>& non_ref_sigs, 
    Reads& reads,
    const int target_pos,
    const double lq_rate_thresh,
    const int retarget_thersh
)
{
    for (auto& read : reads)
    {   
        if (read.sb_ptrn == 'A' && read.nonref_lq_rate < lq_rate_thresh)
        {
            annot_loc_variant_siganature(read, target_pos, retarget_thersh);
            non_ref_sigs.push_back(read.sb_loc_signature);       
        }
    }
}


Read find_representative_read(
    Reads& reads, 
    const int target_pos,
    const UserParams& user_params
)
{      
    std::vector<std::string_view> non_ref_sigs;
    collect_non_ref_ptrns(
        non_ref_sigs, reads, target_pos, 
        user_params.lq_rate_thresh, user_params.retarget_thresh);

    if (non_ref_sigs.empty())
    {
        collect_non_ref_ptrns(
            non_ref_sigs, reads, target_pos, 1.1, user_params.retarget_thresh);
    }
    
    std::string_view common_sig = find_commonest_str(non_ref_sigs);
    
    Reads _tmp;
    for (auto& read : reads)
    {
        if (read.sb_ptrn == 'A' && read.sb_loc_signature == common_sig)
        {
            _tmp.push_back(read);
        }
    }

    std::sort(
        _tmp.begin(), _tmp.end(), [](const Read& a, const Read& b)
        {
            return a.central_score > b.central_score ? true : false;
        }
    );

    return _tmp[0];
}

 
void make_core_kmers(
    std::string_view seq, 
    const size_t k, 
    const size_t core_idx,
    Kmers& kmers
)
{
    size_t n = seq.size();
    if (n <= k)
    {
        kmers.insert(seq);
        return;
    }
    
    for (size_t i = 0; i <= n - k; ++i)
    {
        if (i <= core_idx && core_idx < i + k)
        {
            kmers.insert(seq.substr(i, k));
        }
    }    
}


void nearest_gap(
    const std::vector<Variant>& variants,
    const int target_pos,
    Variant& _gap
)
{
    int dist = INT_MAX;    
    for (const auto& v : variants)
    {
        if (v.is_substitute) continue;
        
        if (v.pos <= target_pos)
        {
            if (target_pos - v.variant_end_pos < dist) 
            {    
                dist = target_pos - v.variant_end_pos;   
                _gap = v;
            }    
        }
        else
        {
            if (v.pos - target_pos < dist)
            {
                dist = v.pos - target_pos;
                _gap = v;
            }   
        }
    }
}


int pos_to_idx(
    const int aln_start,
    const int gap_pos,
    CigarVec& cigar_vec
)
{
    char op = '\0';
    int op_len = 0, curr_pos = aln_start, read_idx = 0;
    for (const auto& c : cigar_vec)
    {
        op = c.first;
        op_len = c.second;
        
        switch (op)
        {
            case 'M':
            case 'X':
            case '=':
                curr_pos += op_len;
                read_idx += op_len;
                if (curr_pos == gap_pos + 1) 
                {    
                    --read_idx;
                    return read_idx;
                }
                break;
            case 'S':
            case 'I':
                read_idx += op_len;
                break;
            case 'D':
            case 'N':
                curr_pos += op_len;
                break;
            default:
                break;
        }         
    }

    return -1;
}


void has_compatible_trans_gap(
    Reads& reads,
    const Kmers& core_kmers,
    Variant& target,
    bool& is_retargeted
)
{
    for (auto& read : reads)
    {
        if (read.sb_ptrn == 'C')
        {
            
            read.sb_kmer_score = count_kmer_overlap(read.seq, core_kmers);
        } 
    }

    std::sort(
        reads.begin(), reads.end(), [](const Read& a, const Read& b)
        {
            return a.sb_kmer_score > b.sb_kmer_score ? true : false;
        }
    );
    
    if (reads[0].sb_kmer_score > 0 && has_gaps(reads[0].cigar_str))
    {
        Variant _gap(-1, "N", "N");       
        nearest_gap(reads[0].variants, target.pos, _gap);
        
        int _gap_idx = pos_to_idx(reads[0].aln_start, _gap.pos, reads[0].cigar_vector);
        if (_gap_idx >= 0)
        {
            size_t i = 0;
            size_t target_idx_start = static_cast<size_t>(_gap_idx);
            size_t target_idx_end = target_idx_start + _gap.alt.size() - 1;
            int cnt = 0;
            for (const auto& kmer : core_kmers)
            {
                i = reads[0].seq.find(kmer);
                
                if (i == std::string_view::npos)
                {
                    continue;;
                }
                else if (target_idx_end < i) 
                {
                    continue;
                }
                else if (i + kmer.size() - 1 < target_idx_start)
                {
                    continue;
                }
                
                ++cnt;
            }
            
            if (cnt)
            {
                target = _gap;
                is_retargeted = true;               
            }
        }   
    }
}


void mock_target_substitute(
    std::string& mock,
    std::string& mock_ref,
    int& mock_start,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    int _start = (target.pos <= int(1 + user_params.kmer_size))
                ? 0 : target.pos - 1 - user_params.kmer_size;
    int _len = (target.pos - _start < int(user_params.kmer_size))  
                ? target.pos - _start : user_params.kmer_size; 
    
    std::string lt_frag = loc_ref.fasta.getSubSequence(
        loc_ref.chrom, _start, _len
    );
    
    std::string mid_frag = target.alt;

    std::string rt_frag = loc_ref.fasta.getSubSequence(
        loc_ref.chrom, target.variant_end_pos - 1, user_params.kmer_size
    );

    mock = lt_frag + mid_frag + rt_frag;
    mock_start = target.pos - user_params.kmer_size;

    mock_ref = loc_ref.fasta.getSubSequence(
         loc_ref.chrom, _start, mock.size()
    );
}


bool has_mapped_bases(std::string_view cigar_str)
{
    if (cigar_str.size())
    { 
        if (cigar_str.find('=') != std::string_view::npos) return true;
    }
    
    return false;
}    


bool has_many_un_matches(const std::string& cigar_str)
{
    //heuristics to disqualify low confident alignment   
    std::string::difference_type n_x = std::count(
        cigar_str.cbegin(), cigar_str.cend(), 'X');
    
    if (n_x > 2) return true;
    
    std::string::difference_type n_d = std::count(
        cigar_str.cbegin(), cigar_str.cend(), 'D');

    if (n_d > 1) return true;

    std::string::difference_type n_i = std::count(
        cigar_str.cbegin(), cigar_str.cend(), 'I');

    if (n_i > 1) return true;

    return false;
}


bool is_target_sb_compatible(
    const Alignment& aln,
    const int target_idx_start,
    const int target_idx_end
) 
{
    if (!has_mapped_bases(aln.cigar_string)) return false;
    if (has_many_un_matches(aln.cigar_string)) return false;
    
    CigarVec cigar_vec = to_cigar_vector(aln.cigar_string);
        
    char op = '\0';
    int op_len = 0;
    int curr_idx = aln.ref_begin;
    for (const auto& c : cigar_vec)
    {
        if (target_idx_start < curr_idx) return false;
        
        op = c.first;
        op_len = c.second;

        switch (op)
        {
            case '=':
                if (curr_idx <= target_idx_start 
                    && target_idx_end <= curr_idx + op_len)
                {
                    return true;
                }
                else
                {
                    curr_idx += op_len;
                }
                break;
            case 'X':
            case 'D':
                curr_idx += op_len;
                break;
            default:
                break;
        }
    }
    
    return false;
}


void parse_for_nearest_gap(
    const int seq_start,
    const std::string& mock_ref,
    const std::string& seq, 
    const std::string& base_quals,
    Alignment& aln,
    LocalReference& loc_ref,
    const int target_pos,
    std::unordered_map<Variant, int>& cis_gaps   
)
{
    std::string tmp = "";
    std::vector<Variant> variants;
    parse_variants(
        seq_start + aln.ref_begin, 
        mock_ref.substr(aln.ref_begin), seq, 
        base_quals, to_cigar_vector(aln.cigar_string), 
        loc_ref.dict, variants, tmp
    ); 

    Variant _gap(-1, "N", "N");
    nearest_gap(variants, target_pos, _gap);
    ++cis_gaps[_gap];
}


void match_to_target(
    const std::string& mock_seq,
    const UserParams& user_params,
    const int target_idx_start,
    const int target_idx_end,
    Reads& reads,
    int& a_cnt
)
{
    Filter filter;
    Alignment aln;
    std::string seq = "";
    for (auto& read : reads)
    {
        if (read.sb_ptrn == 'B'
            && read.nonref_lq_rate < user_params.lq_rate_thresh)
        {                    
            seq = static_cast<std::string>(read.seq);
            
            //against mock
            local_alignment(
                user_params.match_score, user_params.mismatch_penal,
                user_params.gap_open_penal, user_params.gap_ext_penal,
                seq, mock_seq, filter, aln
            );
            
            //rough mapped rate
            double mapped_bases_rate = aln.sw_score / static_cast<double>(user_params.match_score * mock_seq.size());
            if (mapped_bases_rate < 0.75) continue;  
                        
            if (
                is_target_sb_compatible(aln, target_idx_start, target_idx_end)
            )
            {
                read.aln_score = aln.sw_score;
                read.sb_ptrn = 'A'; //tentative may change in "search_retargetable"
                ++a_cnt;               
            }
        }
    }
}


void search_retargetable(
    const std::string& mock_ref,
    const UserParams& user_params,
    const int target_pos,
    const int mock_start,
    int& a_cnt,
    Reads& reads,
    LocalReference& loc_ref,
    std::unordered_map<Variant, int>& retargetables
)
{
    Filter filter;
    Alignment aln;
    std::string seq = "";
    for (auto& read : reads)
    {
        if (read.sb_ptrn == 'A')
        {
            seq = static_cast<std::string>(read.seq);
            local_alignment(
                user_params.match_score, user_params.mismatch_penal,
                user_params.gap_open_penal, user_params.gap_ext_penal,
                seq, mock_ref, filter, aln
            );
            
            if (read.aln_score <= aln.sw_score)
            {
                read.sb_ptrn = 'C';
                --a_cnt;
                continue;
            }
            
            if (has_gaps(aln.cigar_string))
            {
                parse_for_nearest_gap(
                    mock_start, mock_ref, seq, read.base_quals,
                    aln, loc_ref, target_pos, retargetables
                );   
            }
        }
    }     
}


void fill_contig(
    const std::string& seq, 
    const std::string& base_quals,
    const std::string& ref_seq, 
    const int aln_start,
    const CigarVec& cigar_vector,
    Contig& contig
)
{
    char op = '\0';
    int op_len = 0, curr_pos = aln_start;
    size_t ref_idx = 0, seq_idx = 0;
    for (const auto& c : cigar_vector)
    {
        op = c.first;
        op_len = c.second;
       
        switch (op)
        {
            case '=':
            case 'X':
            case 'M':
                for (int i = 0; i < op_len; ++i)
                {
                    contig.positions.push_back(curr_pos);
                    contig.ref_bases.push_back(ref_seq.substr(ref_idx, 1));
                    contig.alt_bases.push_back(seq.substr(seq_idx, 1));
                    contig.base_quals.push_back(base_quals.substr(seq_idx, 1));
                    
                    ++curr_pos;
                    ++ref_idx;
                    ++seq_idx; 
                }
                break; 
            case 'I':
                if (!contig.alt_bases.empty() && !contig.base_quals.empty())
                {
                    contig.alt_bases.pop_back();
                    contig.alt_bases.push_back(seq.substr(seq_idx - 1, op_len + 1));
                    contig.base_quals.pop_back();
                    contig.base_quals.push_back(base_quals.substr(seq_idx - 1, op_len + 1));
                }
                seq_idx += op_len;
                break;
            case 'D':
                if (!contig.ref_bases.empty())
                {
                    contig.ref_bases.pop_back();
                    contig.ref_bases.push_back(ref_seq.substr(ref_idx - 1, op_len + 1));
                }
                 
                ref_idx += op_len;
                curr_pos += op_len;
                break;
            case 'N':
                contig.skip_starts.push_back(curr_pos);
                curr_pos += op_len;
                contig.skip_ends.push_back(curr_pos - 1);
                break;
            case 'S':
                seq_idx += op_len;
                break;
            default:
                break;
        }
    }    
}


void fill_contig_by_mock(
    Reads& reads, 
    const int mock_start,
    const int target_pos,
    const std::string& mock_ref,
    const UserParams& user_params,
    Contig& contig
)
{
    Read _rep = find_representative_read(reads, target_pos, user_params);
    
    Filter filter;
    Alignment aln;
    std::string seq = static_cast<std::string>(_rep.seq);
    local_alignment(
        user_params.match_score, user_params.mismatch_penal,
        user_params.gap_open_penal, user_params.gap_ext_penal,
        seq, mock_ref, filter, aln
    );

    fill_contig(
        seq, _rep.base_quals, mock_ref.substr(aln.ref_begin), 
        mock_start + aln.ref_begin, to_cigar_vector(aln.cigar_string), contig
    );
}


void retarget_to_indel(
    Reads& reads,
    Variant& target,
    Contig& contig,
    const UserParams& user_params, 
    LocalReference& loc_ref,
    bool& is_retargeted,
    bool& is_non_supporting,
    bool& is_mocked
)
{
    int a_cnt = 0, b_cnt = 0;
    for (auto& read : reads)
    {
        annot_ref_seq(read, loc_ref);
        annot_splice_pattern(read);
        annot_covering_ptrn(read, target, loc_ref, is_retargeted);
        annot_clip_pattern(read, target);
        eval_read_quality(read, user_params);
        read.local_uniqueness = eval_loc_uniq(read, user_params, loc_ref);
        substitute_patterns(a_cnt, b_cnt, read, target, user_params);
    } 
    
    if (a_cnt)
    {
        std::unordered_map<Variant, int> cis_gaps;
        collect_gaps(target.pos, user_params.retarget_thresh, loc_ref, reads, cis_gaps);

        Variant _gap(-1, "N", "N");
        bool is_exact = false, is_loc_uniq = false;
        is_exact_cis(
        reads, cis_gaps, cis_gaps.empty(), _gap, is_exact, is_loc_uniq);
        
        // linked with a gap in cis (typically complex indel)
        if (is_exact)
        {
            is_retargeted = true;
            target = _gap;
            return; 
        }

        // all target substitute reads are not uniquely mapped
        // -> gaps may be mapped as substitutes
        if (!is_loc_uniq)
        {       
            Read _rep = find_representative_read(reads, target.pos, user_params);

            Kmers core_kmers;
            make_core_kmers(
                _rep.seq, user_params.kmer_size, _rep.sb_read_idx, core_kmers
            );
           
            has_compatible_trans_gap(reads, core_kmers, target, is_retargeted);
            if (is_retargeted) return;
        }
        
        //no retargetting
        return;
    }
    
    if (b_cnt) 
    {
        //obscured cases
        int mock_start = -1;
        std::string mock_seq, mock_ref;        
        mock_target_substitute(
            mock_seq, mock_ref, mock_start, target, user_params, loc_ref
        );
        
        const int idx_start 
            = user_params.kmer_size;
        const int idx_end 
            = idx_start + static_cast<int>(target.ref.size());
        
        match_to_target(
            mock_seq, user_params, idx_start, idx_end, reads, a_cnt
        );
        
        if (a_cnt)
        {
            std::unordered_map<Variant, int> retargetables;
            search_retargetable(
                mock_ref, user_params, target.pos, 
                mock_start, a_cnt, reads, loc_ref, retargetables
            );

            
            // mock to emptu
            retargetables.clear();
            
            if (!retargetables.empty())
            {
                Variant _gap(-1, "N", "N");
                commonest_cis_gap(_gap, retargetables);
                
                int _cnt = retargetables.at(_gap);

                if (2 * _cnt  > a_cnt)
                {
                    target = _gap;
                    is_retargeted = true;
                    return;
                }
            }
            else if (a_cnt)
            {
                fill_contig_by_mock(
                    reads, mock_start, target.pos, mock_ref, user_params, contig);
                is_mocked = true;
                return;
            }
            else
            {       
                is_non_supporting = true;
            }   
        }
    }
    
    is_non_supporting = true;
} 


void from_target_substitute_reads(
    Contig& contig,
    Reads& reads,
    Reads& targets,
    Reads& non_targets,
    const int target_pos, 
    const UserParams& user_params,
    const bool is_mocked
)
{
    size_t max_size = reads.size();
    targets.reserve(max_size);
    non_targets.reserve(max_size);
    for (size_t i = 0; i < max_size; ++i)
    {    
        if (reads[i].sb_ptrn == 'A')
        {
             transfer_elem(targets, reads, i);
        }
        else if (reads[i].sb_ptrn == 'C' && reads[i].is_tight_covering)
        {
             transfer_elem(non_targets, reads, i);
        }
    }
    
    reads.clear();
    targets.shrink_to_fit();
    non_targets.shrink_to_fit();
      
    if (!is_mocked)
    { 
        Read _rep = find_representative_read(targets, target_pos, user_params);
        
        fill_contig(
            static_cast<std::string>(_rep.seq), 
            _rep.base_quals,
            static_cast<std::string>(_rep.ref_seq),    
            _rep.aln_start,
            _rep.cigar_vector,
            contig
        );
    }
} 


void from_no_substitute_reads(
    const Variant& target,
    Contig& contig,
    Reads& reads,
    Reads& non_targets
)
{
    size_t max_size = reads.size();
    non_targets.reserve(max_size);
    for (size_t i = 0; i < max_size; ++i)
    {
        if (reads[i].is_tight_covering)
        {
             transfer_elem(non_targets, reads, i);
        }
    }
    reads.clear();
    non_targets.shrink_to_fit();
    
    contig.positions.push_back(target.pos);
    contig.ref_bases.push_back(target.ref);
    contig.alt_bases.push_back(target.ref); //ref   
    contig.base_quals.push_back("F"); //pseudo
}    
