#include <set>
#include <string>
#include <vector>
#include <string.h>
#include <algorithm>
#include <string_view>

#include "util.h"
#include "merge.h"
#include "match.h"
#include "contig.h"
#include "ssw/ssw_cpp.h"


inline void to_left(std::string& allele, std::string& lt_seq)
{
    allele.pop_back();
    lt_seq.pop_back();
    if (lt_seq.size())
    {
        allele = lt_seq.substr(lt_seq.size() - 1, 1) + allele;
    }    
}


inline void to_right(std::string& allele, std::string& rt_seq)
{
    allele.erase(0, 1);
    if (rt_seq.size())
    {
        allele += rt_seq.substr(0, 1);
    }
    rt_seq.erase(0, 1);
}


std::vector<int> pos_by_read_idx(const int aln_start, const CigarVec& cigar_vec)
{
    std::vector<int> pos_by_idx;   
    char op = '\0';
    int op_len = 0;
    int genomic_pos = aln_start;
    for (const auto& c : cigar_vec)
    { 
        op = c.first;
        op_len = c.second;
        switch (op)
        {
            case '=':
            case 'M':
            case 'X':
                for (int i = 0; i < op_len; ++i)
                {                   
                    pos_by_idx.push_back(genomic_pos);
                    
                    ++genomic_pos;
                }
                break;
            case 'I':
                std::fill_n(std::back_inserter(pos_by_idx), op_len, genomic_pos - 1);
                break;
            case 'D':
                pos_by_idx.pop_back();
                genomic_pos += op_len;
                break;
            case 'N':
                pos_by_idx.pop_back();
                genomic_pos += op_len;
                break;
            case 'S':
                std::fill_n(std::back_inserter(pos_by_idx), op_len, -1);
                break;
        }
    }
    
    // to allow match in case where read starts/ends
    // in the middle of contig's shiftable region
    if (pos_by_idx[0] != -1) pos_by_idx[0] = -1;
    size_t last = pos_by_idx.size() - 1;
    if (pos_by_idx[last] != -1) pos_by_idx[last] = -1;
       
    return pos_by_idx;
}


inline bool is_shiftable(std::string_view allele)
{
    return (allele.front() == allele.back());
}


int annot_rep(std::string_view rep_unit, std::string_view ss_seq)
{
    ss_seq.remove_prefix(1);
    ss_seq.remove_suffix(1);
    
    int n_reps = 0;
    size_t seq_len = ss_seq.size(), unit_len = rep_unit.size();
    if (seq_len >= unit_len)
    {
        for (size_t k = 0; k <= (seq_len - unit_len); k += unit_len)
        {
            if (ss_seq.substr(k, unit_len) == rep_unit) ++n_reps;
            else return n_reps;
        }
    }
    
    return n_reps; 
}


void annot_shiftable_segment(
    ShiftableSegment& ss,
    const Variant& target, 
    const Contig& contig
)
{
    std::string lt_seq = contig.seq.substr(0, contig.lt_len);
    std::string rt_seq = contig.seq.substr(contig.rt_start_idx);
    
    std::string lt_allele = "", rt_allele = "";
    if (target.is_ins)
    {
        lt_allele = target.alt;
        rt_allele = target.alt;
    }
    else 
    {
        lt_allele = target.ref;
        rt_allele = target.ref;
    }
    
    size_t i = 0;
    while(i < contig.lt_len && is_shiftable(lt_allele)) 
    {
        to_left(lt_allele, lt_seq);
        ++i;
    }

    size_t j = 0;
    do
    {
        to_right(rt_allele, rt_seq);
        ++j;
    } while(j < contig.rt_len && is_shiftable(rt_allele));
    --j; //undo once
    
    ss.start = contig.lt_len < i + 1 ? 0 : contig.lt_len - i - 1;
    ss.start = static_cast<int>(ss.start);
    ss.end = contig.rt_start_idx + j;
    ss.end = static_cast<int>(ss.end);

    //shiftable seq (with boundary bases (1 each side))
    const size_t region_len = ss.end - ss.start + 1;
    ss.seq = contig.seq.substr(ss.start, region_len);
    
    ss.rep_unit = target.minimal_repeat_unit();
    ss.unit_len = ss.rep_unit.size();

    ss.n_reps = annot_rep(ss.rep_unit, ss.seq);
    
    //offset for self in simple insertions
    if (target.is_ins && ss.n_reps) ss.n_reps -=1; 
}


inline bool has_mismatches(string_view cigar_str)
{
    return (cigar_str.find('X') != std::string_view::npos);
}


inline bool has_N(std::string_view seq)
{
    return (seq.find('N') != std::string_view::npos); 
}


char ss_covered_patterns(
    const Alignment& aln, 
    const ShiftableSegment& ss, 
    const size_t contig_len,
    const size_t query_len
)
{
    //bool is_complete_cover = (aln.ref_begin <= ss.start && ss.end <= aln.ref_end);
   /* bool is_end_mapped = (
        aln.query_begin == 0 || aln.query_end == int(query_len - 1)
    );*/
    
    if (aln.cigar_string.find('X') != std::string::npos) return 'N';
    
    const int max_query_end = static_cast<int>(query_len) - 1;
    const int max_ref_end = static_cast<int>(contig_len) - 1;

    if (!aln.query_begin && aln.ref_end == max_ref_end) return 'R';        
    if (aln.query_end == max_query_end && !aln.ref_begin) return 'L';
    
    return 'M';
}


double fw_sim_score(
    std::string_view query, 
    std::string_view quals,
    std::string_view subject,
    const UserParams& user_params
)
{ 
    double miss = 0.0;
    for (size_t i = 0; i < query.size(); ++i)
    {
        //boundary must be exact
        if (i == 0 || i == query.size() - 1)
        {
            if (query[i] != subject[i]) return 0.0;
        }
        else
        {
            if (query[i] != subject[i])
            {
                //high quality mismathces -> higher penalty
                if (quals[i] >= user_params.base_q_thresh) ++miss;
                else miss += 0.5;
            }
        }
    }
    return 1.0 - miss/query.size();
}


double rv_sim_score(
    std::string_view query, 
    std::string_view quals,
    std::string_view subject,
    const UserParams& user_params
) 
{
    double miss = 0.0;
    size_t q_last = query.size() - 1, s_last = subject.size() - 1;
    for (size_t i = 0; i <= q_last; ++i) 
    {
        if (i == 0)
        {
            if (query[q_last - i] != subject[s_last - i]) return 0.0;
        }
        else
        {        
            if (query[q_last - i] != subject[s_last - i])
            {
                if (quals[q_last - i] >= user_params.base_q_thresh) ++miss;
                else miss += 0.5;
            
            }
        }
    }
    return 1.0 - miss/query.size();
}


char indel_match_pattern(
    const std::string& query, 
    std::string_view base_quals,
    const int lt_mapped_cnt,
    const int rt_mapped_cnt,
    const std::vector<int>& pos_by_idx,
    const int local_uniqueness,
    const Contig& contig,
    const ShiftableSegment& ss,
    const UserParams& user_params,
    const Filter& filter,
    const Aligner& aligner, 
    Alignment& aln
)          
{   
    int32_t mask_len = strlen(query.c_str()) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;
            
    aligner.Align(
        query.c_str(), contig.seq.c_str(), 
        contig.len, filter, &aln, mask_len
    );
    
    if (contig.is_mocked)
    {
        if (aln.ref_begin <= ss.start - 5 && ss.end + 5 <= aln.ref_end)
        {
            const char covered_ptrn 
            = ss_covered_patterns(aln, ss, contig.len, query.size());
            
            //exact match
            if (covered_ptrn != 'N') return covered_ptrn;
            
            /*
            if (is_exact_match(aln, ss, query.size()))
            {
                if (!aln.ref_begin) return 'L';
                if (aln.ref_end == static_cast<int>(contig.len) - 1) return 'R';   
                return 'M';  
            }*/       
        }
        return 'F';
    }
    
    // aln ends before ss
    if (aln.ref_end <= ss.start)
    {
        // due to rt-clip
        if (aln.query_end != static_cast<int>(query.size()) - 1) return 'F';
        // too short -> undetermined
        else return 'U';
    }    
    
    // aln starts after ss
    if (ss.end <= aln.ref_begin) 
    {    
        // due to lt-clip
        if (aln.query_begin) return 'F'; //lt-clipped
        // too short
        else return 'U'; 
    } 
    
    // aln starts between ss_start/ss_end
    if (ss.start < aln.ref_begin 
        &&  aln.ref_begin < ss.end)
    {
        if (aln.query_begin) return 'F';
        //allow if read starts after ss-start 
    }
    // aln ends between ss_start/ss_end
    else if (ss.start < aln.ref_end 
             && aln.ref_end < ss.end)
    {
        if (aln.query_end != static_cast<int>(query.size()) - 1) return 'F';  
        //allow if read ends before ss.end 
    }  
    // aln covers shiftable but less than the orig mapping
    else if (aln.ref_begin <= ss.start
             && ss.end <= aln.ref_end)
    {
        const char covered_ptrn 
            = ss_covered_patterns(aln, ss, contig.len, query.size());
        
        //exact matches
        if (covered_ptrn != 'N') return covered_ptrn;
        
        //short margin
        if (ss.start - aln.ref_begin < 10)
        {
            if (lt_mapped_cnt > 2 * (ss.start - aln.ref_begin)) return 'F';
        }

        if (aln.ref_end - ss.end < 10)
        {
            if (rt_mapped_cnt > 2 * (aln.ref_end - ss.end)) return 'F';
        }
    }
    
    const double thresh = 0.95;
    //const int  margin = 2;
    size_t crit_start = 0;
    std::string crit_seq;
    std::string_view quals;
    if (aln.ref_begin <= ss.start)
    {
        //if (int(ss.start) - aln.ref_begin + aln.query_begin >= margin)
        if (ss.start >= aln.ref_begin - aln.query_begin)
        {   
            crit_start = ss.start - aln.ref_begin + aln.query_begin;
        }

        crit_seq = query.substr(
            //crit_start, ss.seq.size() + 2 * (margin + 1)
            crit_start, ss.seq.size()
        );
        
        quals = base_quals.substr(
            //crit_start, ss.seq.size() + 2 * (margin + 1)  
            crit_start, ss.seq.size()  
        );   

        if (fw_sim_score(crit_seq, quals, ss.seq, user_params) < thresh) return 'F';
    }
    else
    { 
        if (ss.n_reps) 
        {    
            if (!aln.query_begin) return 'U';
            else return 'F';
        }
        else
        {
            crit_seq = query.substr(aln.query_begin, ss.end - aln.ref_begin + 1);
            
            quals = base_quals.substr(aln.query_begin, ss.end - aln.ref_begin);
            
            if (rv_sim_score(crit_seq, quals, ss.seq, user_params) < thresh) return 'F'; 
        }
    } 
 
    int observed_reps = 0;
    if (ss.n_reps && ss.unit_len < 4) 
    { 
        if (aln.ref_end <= ss.end) return 'U';
        
        const int query_lt_len = ss.start - aln.ref_begin + aln.query_begin + 1;
        
        //mid + rt
        std::string non_lt_fgmt = query.substr(query_lt_len);
        observed_reps = count_repeats(ss.rep_unit, non_lt_fgmt);

        if (contig.mid_len)
        {    
            int mid_rep = 0;
            std::string mid_fgmt = query.substr(query_lt_len, contig.mid_len);
            mid_rep = count_repeats(ss.rep_unit, mid_fgmt);
            
            if (mid_rep) observed_reps -= 1; // offset self
        }
        
        if (ss.n_reps !=  observed_reps) return 'F';  
    }
    
    int query_ss_start = aln.query_begin + ss.start - aln.ref_begin;
    int query_ss_end = query_ss_start + ss.end - ss.start;
    
    if (contig.pos_by_seq_idx.empty()) //extension case
    {
        if (has_mismatches(aln.cigar_string)) return 'U';       
    }
    else if (aln.ref_begin < ss.start && ss.end < aln.ref_end) 
    {   
        if (pos_by_idx[query_ss_start] == -1 
            && pos_by_idx[query_ss_end] == -1) return 'U';
        
        if (pos_by_idx[query_ss_start] == -1 
            &&  contig.pos_by_seq_idx[ss.end] == pos_by_idx[query_ss_end]) return 'M';
        if (contig.pos_by_seq_idx[ss.start] == pos_by_idx[query_ss_start] 
            &&  pos_by_idx[query_ss_end] == -1) return 'M';   
        
        if (local_uniqueness == 1)
        {
            if (contig.pos_by_seq_idx[ss.start] != pos_by_idx[query_ss_start]
                || contig.pos_by_seq_idx[ss.end] != pos_by_idx[query_ss_end]) return 'F';
        }
        else
        {
            if (contig.pos_by_seq_idx[ss.start] == pos_by_idx[query_ss_start] 
                && pos_by_idx[query_ss_end] <= contig.pos_by_seq_idx[ss.end]) return 'M';
            if (contig.pos_by_seq_idx[ss.start] <= pos_by_idx[query_ss_start] 
                && pos_by_idx[query_ss_end] == contig.pos_by_seq_idx[ss.end]) return 'M';
        }  
    } 
    
    return 'M';
}


void classify_cand_indel_reads(
    Reads& candidates,
    Reads& non_targets,

    Reads& lt_matches,
    Reads& mid_matches,
    Reads& rt_matches,

    Reads& undetermined,
   
    const Contig& contig,
    const ShiftableSegment& ss,
    const UserParams& user_params
)
{
    Filter filter;
    Alignment aln;
    Aligner aligner(
        user_params.match_score, 
        user_params.mismatch_penal,
        255, //gap_open_penalty
        255 // gap_extention_penalty
    );

    for (size_t i = 0; i < candidates.size(); ++i)
    {         
        /*if (ss.n_reps && candidates[i].incomplete_shift)
        {
            transfer_elem(non_targets, candidates, i);
            continue;   
        }*/
        
        std::vector<int> pos_by_idx = pos_by_read_idx(
            candidates[i].aln_start, 
            candidates[i].cigar_vector
        );

        char match_rslt = indel_match_pattern(
            static_cast<std::string>(candidates[i].seq), candidates[i].base_quals,
            candidates[i].lt_end_matches, candidates[i].rt_end_matches,
            pos_by_idx, candidates[i].local_uniqueness, contig, ss, user_params, filter, aligner, aln
        );
        
        switch (match_rslt)
        {
            case 'L':
                transfer_elem(lt_matches, candidates, i);
                break;
            case 'M':
                transfer_elem(mid_matches, candidates, i);
                break;
            case 'R':
                transfer_elem(rt_matches, candidates, i);
                break;
            case 'F':
                transfer_elem(non_targets, candidates, i);
                break;
            case 'U':
                transfer_elem(undetermined, candidates, i);
                break;
        }                
    }
    candidates.clear();
    candidates.shrink_to_fit();
}


void classify_simplified(
    Reads& candidates,
    Reads& targets,
    Reads& non_targets,
    Reads& undetermined,
    const Contig& contig,
    const ShiftableSegment& ss,
    const UserParams& user_params
)
{
    Filter filter;
    Alignment aln;
    Aligner aligner(
        user_params.match_score, 
        user_params.mismatch_penal,
        255, //gap_open_penalty
        255 // gap_extention_penalty
    );

    for (size_t i = 0; i < candidates.size(); ++i)
    {         
        std::vector<int> pos_by_idx = pos_by_read_idx(
            candidates[i].aln_start, candidates[i].cigar_vector
        );
        
        char match_rslt = indel_match_pattern(
            static_cast<std::string>(candidates[i].seq), candidates[i].base_quals,
            candidates[i].lt_end_matches, candidates[i].rt_end_matches,
            pos_by_idx, candidates[i].local_uniqueness, contig, ss, user_params, filter, aligner, aln
        );

        switch (match_rslt)
        {
            case 'L':
            case 'M':
            case 'R':
                transfer_elem(targets, candidates, i);
                break;
            case 'F':
                transfer_elem(non_targets, candidates, i);
                break;
            case 'U':
                transfer_elem(undetermined, candidates, i);
                break;
        }                
    } 
    
    candidates.clear();
    candidates.shrink_to_fit();
}


void classify_cand_indel_read_2(
    Reads& targets,
    Reads& candidates,
    Reads& non_targets,
    Reads& undetermined,
    const std::vector<Variant>* p_decomposed,
    const Variant& target,
    const Contig& contig,
    const UserParams& user_params
)
{
    if (p_decomposed != NULL)
    {
        Reads tmp;
        std::vector<Variant> shared;
        for (size_t i = 0; i < candidates.size(); ++i)
        {
            find_shared_variants(shared, candidates[i].variants, *p_decomposed);
            
            if (shared.size() == p_decomposed->size())
            {
                transfer_elem(targets, candidates, i);
            }
            else
            {
                transfer_elem(tmp, candidates, i);
            }
            shared.clear();
        }
        
        std::swap(candidates, tmp);
    }
    
    if (candidates.empty()) return;
    
    ShiftableSegment ss;
    annot_shiftable_segment(ss, target, contig);
    
    classify_simplified(
        candidates, targets, non_targets, undetermined,
        contig, ss, user_params
    );   
}
