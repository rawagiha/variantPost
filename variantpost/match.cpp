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


inline bool is_shiftable(std::string_view allele)
{
    return (allele.front() == allele.back());
}

void annot_shiftable_segment(
    ShiftableSegment& ss,
    const Variant& target, 
    const Contig& contig
)
{
    //std::cout << "contig prep" << std::endl;
    //std::cout << contig.seq << " " << contig.seq.size() << " " << contig.rt_start_idx << std::endl;
    const size_t lt_len = contig.lt_len;
    const size_t rt_len = contig.rt_len;
    std::string lt_seq = contig.seq.substr(0, lt_len);
    std::string rt_seq = contig.seq.substr(contig.rt_start_idx);
    
    std::string allele_for_lt = "";
    std::string allele_for_rt = "";
    if (target.is_ins)
    {
        allele_for_lt = target.alt;
        allele_for_rt = target.alt;
    }
    else // here, target must not be substitute
    {
        allele_for_lt = target.ref;
        allele_for_rt = target.ref;
    }
    
    //std::cout << "do lt" << std::endl;
    size_t i = 0;
    while(i < lt_len && is_shiftable(allele_for_lt)) 
    {
        to_left(allele_for_lt, lt_seq);
        ++i;
    }
    
    //std::cout << "do rt" << std::endl;
    size_t j = 0;
    do
    {
        to_right(allele_for_rt, rt_seq);
        ++j;
    } while(j < rt_len && is_shiftable(allele_for_rt));
    --j; //undo once for rt aln
    
    ss.boundary_start = lt_len - (i + 1);
    ss.boundary_end = contig.rt_start_idx + j - 1;
    
    //std::cout << ss.boundary_end << " " << ss.boundary_start  << std::endl;
    const size_t shiftable_len = ss.boundary_end - ss.boundary_start - 1;
    std::string shiftable = contig.seq.substr(
        ss.boundary_start, shiftable_len
    );
    
    ss.fw_repeat_unit = target.minimal_repeat_unit();
    
    ss.unit_len = ss.fw_repeat_unit.size();
    
    std::string rv_repeat_unit = ss.fw_repeat_unit;
    std::reverse(rv_repeat_unit.begin(), rv_repeat_unit.end());   
    ss.rv_repeat_unit = rv_repeat_unit;   
    
    //std::cout << " repeat "  << std::endl;
    int n_tandem_repeats = 0;
    bool is_complete_tandem_repeat = false;
    if (shiftable_len >= ss.unit_len)
    {
        for (size_t k = 0; k <= (shiftable_len - ss.unit_len); k += ss.unit_len)
        {
            if (shiftable.find(ss.fw_repeat_unit, k) != std::string::npos)
            {
                ++n_tandem_repeats;
                is_complete_tandem_repeat = true;
            }
            else
            {
                is_complete_tandem_repeat = false;
                break;
            }               
        }
        
        size_t remainder = shiftable_len % ss.unit_len;
        if (remainder) is_complete_tandem_repeat = false;
    }
    
    ss.n_tandem_repeats = n_tandem_repeats;
    ss.is_complete_tandem_repeat = is_complete_tandem_repeat;
}


inline bool has_mismatches(string_view cigar_str)
{
    return (cigar_str.find('X') != std::string_view::npos);
}


bool has_critical_mismatches(
    string_view cigar_str,
    const int ref_start,
    const int lt_end_idx, 
    const int rt_start_idx,
    const int margin = 5
)
{   
    //critical region
    const int lt_bound = lt_end_idx - margin;   
    const int rt_bound = rt_start_idx + margin;
    std::vector<std::pair<char, int>> cigar_vec{to_cigar_vector(cigar_str)};
    char op;
    int op_len = 0, i = ref_start;  
    for (const auto& c : cigar_vec) 
    {
         op = c.first;
         op_len = c.second;    
         switch (op) 
         {
            case '=':
                i += op_len;
                break;
            case 'X':
                for (int j = 0; j < op_len; ++j)
                {
                    if (lt_bound <= i && i <= rt_bound)
                    {
                        return true;
                    }
                    ++i;
                }
                break;
            case 'S':
                if (lt_bound <= i && i <= rt_bound)
                {
                    return true;
                }
                break;
            default:
                break;
         }        
    }
    return false;
}    


char indel_match_pattern
(
    const std::string& query, 
    const Contig& contig,
    const ShiftableSegment& ss,
    const Filter& filter,
    const Aligner& aligner, 
    Alignment& alignment
)          
{
    //contig layout
    const int lt_end_idx = contig.lt_end_idx;
    const int rt_start_idx = contig.rt_start_idx;
    const int contig_len = contig.len;
    const int mid_len = contig.mid_len;
    const bool is_del = (mid_len == 0);
    
    int32_t mask_len = strlen(query.c_str()) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;
            
    aligner.Align(query.c_str(), contig.seq.c_str(), contig_len, filter, &alignment, mask_len);
     
    const int ref_start = alignment.ref_begin;
    const int ref_end = alignment.ref_end;
    
    // still has gap
    // 'F': failed to match 
    if (has_gaps(alignment.cigar_string))  return 'F';
    
    // doesn't overlap the target region (might've contained target if longer)
    // 'U': undertermined 
    if (ref_end <= lt_end_idx || rt_start_idx <= ref_start) return 'U';
    
    // del needs to be covered
    if (is_del)
    {
        if ((ref_start <= lt_end_idx) && (rt_start_idx <= ref_end)) { /*pass*/ }
        else return 'F'; 
    }
    
    // checks for complete tandem repeat  
    int lt_rep = 0, mid_rep = 0, rt_rep = 0;
    if (ss.is_complete_tandem_repeat && ss.unit_len < 4) 
    { 
        //const int boundary_start = ss.boundary_start;
        const int read_start = alignment.query_begin;

        // must be bounded
        if (ss.boundary_start < ref_start 
            || ref_end <= ss.boundary_end) return 'U';
        
        std::string lt_fragment = query.substr(
            0, (ss.boundary_start + 1 - ref_start) + read_start
        );
        
        std::reverse(lt_fragment.begin(), lt_fragment.end());
        lt_rep = count_repeats(ss.rv_repeat_unit, lt_fragment);

        if (!is_del)
        {
            std::string mid_fragment = query.substr(
                (ss.boundary_start + 1 - ref_start) + read_start, mid_len
            );
            
            // -1 offset for own sequence
            mid_rep = count_repeats(ss.fw_repeat_unit, mid_fragment) - 1;   
        }

        std::string rt_fragment = query.substr(
            (ss.boundary_start + 1 - ref_start) + mid_len + read_start
        );
        rt_rep = count_repeats(ss.fw_repeat_unit, rt_fragment);

        // repeat as expected
        if (ss.n_tandem_repeats != (lt_rep + mid_rep + rt_rep)) return 'F';  
    } 
    
    if (contig.is_mocked)
    {   
        bool disqualified = has_critical_mismatches(
            alignment.cigar_string,
            ref_start,
            lt_end_idx, 
            rt_start_idx
        );
        
        if (disqualified) return 'F';
    }   
    
    // annotate match pattern 
    if (ref_start < lt_end_idx 
        && rt_start_idx < ref_end 
        && !has_mismatches(alignment.cigar_string)) 
    {
        if (ref_start == 0) 
        {
            //left extender (may be used for contig extension)
            return 'L';
        }
       
        if (ref_end == (contig_len - 1))
        {
            //right extender
            return 'R';
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
    bool is_complete_tandem_repeat = ss.is_complete_tandem_repeat;
    
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
        if (is_complete_tandem_repeat && candidates[i].incomplete_shift)
        {
            transfer_elem(non_targets, candidates, i);
            continue;   
        }

        char match_rlst = indel_match_pattern(
            std::string(candidates[i].seq),
            contig,
            ss,
            filter,
            aligner,
            aln
        );

        switch (match_rlst)
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


void classify_cand_indel_reads_2(
    Reads& candidates,
    
    Reads& targets,
    Reads& non_targets,
    Reads& undetermined,
   
    const Contig& contig,
    const ShiftableSegment& ss,
    const UserParams& user_params
)
{
    bool is_complete_tandem_repeat = ss.is_complete_tandem_repeat;
    
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
        if (is_complete_tandem_repeat && candidates[i].incomplete_shift)
        {
            transfer_elem(non_targets, candidates, i);
            continue;   
        }

        char match_rlst = indel_match_pattern(
            std::string(candidates[i].seq),
            contig,
            ss,
            filter,
            aligner,
            aln
        );

        switch (match_rlst)
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
    const Variant& target,
    const Contig& contig,
    const UserParams& user_params
)
{
    //std::cout << "shift annot" << std::endl;
    ShiftableSegment ss;
    if (!target.is_complex)
    {
        annot_shiftable_segment(ss, target, contig);
    }
       
    //std::cout << "do clsfy " << std::endl;
    classify_cand_indel_reads_2(
        candidates,
        
        targets,
        non_targets, 
        undetermined,
        
        contig,
        ss,
        user_params
    );   
}
