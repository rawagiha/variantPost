#include <vector>
#include <string>
#include <chrono>
#include <algorithm>

#include "eval.h"
#include "util.h"
#include "match.h"
#include "reads.h"
#include "search.h"
#include "contig.h"
#include "substitutes.h"

SearchResult::SearchResult() {}


void SearchResult::fill_read_info(const Reads& reads, const int target_status)
{
    for (const auto& read : reads)
    {
        read_names.push_back(std::string(read.name));
        are_reverse.push_back(read.is_reverse);
        target_statuses.push_back(target_status);
        are_from_first_bam.push_back(read.is_from_first_bam);          
    }
}


void rescue_sb_target_reads(Reads& targets, Reads& tmp, Reads& src)
{
    for (size_t i = 0; i < src.size(); ++i)
    {
        if (src[i].sb_ptrn == 'A')
        {
            transfer_elem(targets, src, i);
        }
        else
        {
            transfer_elem(tmp, src, i);
        }
    }
}


void SearchResult::report(
    const Contig& contig,
    Reads& targets,
    Reads& non_targets,
    Reads& undetermined,
    const bool _is_retargeted //default to false in header
)
{
    size_t buff_size = (
        targets.size() + non_targets.size() + undetermined.size()
    );

    read_names.reserve(buff_size);
    are_reverse.reserve(buff_size);
    target_statuses.reserve(buff_size);
    are_from_first_bam.reserve(buff_size);
    
    if (_is_retargeted)
    {
        Reads tmp_n, tmp_u;
        rescue_sb_target_reads(targets, tmp_n, non_targets);       
        rescue_sb_target_reads(targets, tmp_u, undetermined);
        fill_read_info(targets, 1);
        fill_read_info(tmp_n, 0);
        fill_read_info(tmp_u, -1);
    }
    else
    {
        fill_read_info(targets, 1);
        fill_read_info(non_targets, 0);
        fill_read_info(undetermined, -1);
    }
     
    positions = contig.positions;
    ref_bases = contig.ref_bases;
    alt_bases = contig.alt_bases;
    base_quals = contig.base_quals;
    skip_starts = contig.skip_starts;
    skip_ends = contig.skip_ends;
    is_retargeted = _is_retargeted;
}


Variant prep_target(
    const int pos, 
    const std::string& ref, 
    const std::string& alt,
    LocalReference& loc_ref
)
{   
    Variant target(pos, ref, alt);
    target.set_leftmost_pos(loc_ref);
    target.set_rightmost_pos(loc_ref); 
    target.is_shiftable = target.lpos != target.rpos ? true : false;
    
    return target;    
}


void prep_reads(
    Reads& reads,
    const std::vector<std::string>& read_names,
    const std::vector<bool>& are_reverse,
    const std::vector<std::string>& cigar_strs,
    const std::vector<int>& aln_starts,
    const std::vector<int>& aln_ends,
    const std::vector<std::string>& read_seqs,
    const std::vector<std::vector<int>>& quals,
    const std::vector<int>& mapqs,
    const std::vector<bool>& are_from_first_bam
)
{
    size_t n_reads = read_names.size();
    reads.reserve(n_reads);
    for (size_t i = 0; i < n_reads; ++i)
    {
        reads.emplace_back(
            read_names[i], are_reverse[i], cigar_strs[i],
            aln_starts[i], aln_ends[i], read_seqs[i],
            quals[i], mapqs[i], are_from_first_bam[i]
        );
    }
}


void from_target_reads(
    Contig& contig,
    const Variant& target,
    Reads& targets,
    Reads& candidates,
    Reads& non_targets,
    Reads& undetermined,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    make_contig(
        contig,
        target, 
        targets, 
        user_params,
        loc_ref
    );
    
    char _eval = eval_by_aln(contig, target, user_params, loc_ref);
    
    if (_eval == 'C')
    {
        //NO -> rather make contig with centered most.
        transfer_vector(non_targets, targets);
        transfer_vector(non_targets, candidates); 
        return;      
    }

    ShiftableSegment ss;
    annot_shiftable_segment(ss, target, contig);   

    Reads lt_matches, mid_matches, rt_matches;    
    classify_cand_indel_reads(
        candidates, non_targets, 
        lt_matches, mid_matches, rt_matches,
        undetermined, contig, ss, user_params
    );
    
    if (_eval != 'A')
    {
        extend_contig(_eval, contig, lt_matches, rt_matches, loc_ref);
        aln_extended_contig(contig, target, user_params, loc_ref);
    }
     
    transfer_vector(targets, lt_matches);
    transfer_vector(targets, mid_matches);
    transfer_vector(targets, rt_matches);
}


void from_candidate_reads(
    Contig& contig,
    const Variant& target,
    Reads& targets,
    Reads& candidates,
    Reads& non_targets,
    Reads& undetermined,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    std::vector<Variant>* p_decomposed = NULL;
    std::vector<Variant> decomposed;
    if (target.is_complex)
    {
        bool has_pass_candidates = false;
        prefilter_cplx_candidates(
            contig, candidates, non_targets, has_pass_candidates, 
            target, user_params, loc_ref, decomposed
        );

        if (!has_pass_candidates) return;
    }
    else
    {
        prefilter_candidates(
            contig, candidates, non_targets, target, 
            user_params, loc_ref
        );
    }
    
    if (candidates.empty()) return;    
    
    // for complex inputs
    if (!decomposed.empty()) p_decomposed = &decomposed;
    
    suggest_contig(target, contig, candidates, user_params, loc_ref);       
    
    //may happen if all candidates are of low-qual bases
    if (contig.seq.empty()) return;
    
    char _eval = eval_by_aln(contig, target, user_params, loc_ref, p_decomposed);
    
    if (_eval == 'B')
    {
        bool is_mocked = false;
        switch_to_mock(
            contig, target, user_params, loc_ref, is_mocked, p_decomposed
        );
        if (!is_mocked)
        {
            transfer_vector(non_targets, candidates);
            return;           
        }
    }
    else if (_eval == 'C')
    {
        transfer_vector(non_targets, candidates);
        return;
    }

    classify_cand_indel_read_2(
        targets, candidates, non_targets, undetermined,
        p_decomposed, target, contig, user_params
    );
}   


void _search_target(
    
    /* data passed from Python interface*/
    
    SearchResult& rslt,
    
    //reference fastafile
    const std::string& fastafile,
    
    //target variant info
    const std::string& chrom,
    const int pos, 
    const std::string& ref,
    const std::string& alt,
    
    //user defined 
    const int mapq_thresh,
    const int base_q_thresh,
    const double lq_base_rate_thresh,
    const int match_score,
    const int mismatch_penal,
    const int gap_open_penal,
    const int gap_ext_penal,
    const int kmer_size,
    const int local_thresh,
    const int retarget_thresh,
    const int ref_start, //local ref start/end
    const int ref_end,   //defined by user's window choice
    
    //reads mapped to the region of interest 
    const std::vector<std::string>& read_names,
    const std::vector<bool>& are_reverse,
    const std::vector<std::string>& cigar_strs,
    const std::vector<int>& aln_starts,
    const std::vector<int>& aln_ends,
    const std::vector<std::string>& read_seqs,
    const std::vector<std::vector<int>>& quals,
    const std::vector<int>& mapqs,
    const std::vector<bool>& are_from_first_bam)
{  
    
    //auto t1 = std::chrono::high_resolution_clock::now();
        
    // do input validation at python ends
    // 1) no fetched reads, 
    // 2) undefined variants...
    // 3) ref_start/ref_end must be non-N region
    
    // packing data from Python 
    UserParams user_params(
        mapq_thresh, base_q_thresh, lq_base_rate_thresh,
        match_score, mismatch_penal, gap_open_penal, gap_ext_penal, 
        kmer_size, local_thresh, retarget_thresh 
    );

    LocalReference loc_ref(fastafile, chrom, ref_start, ref_end);   

    Variant target = prep_target(pos, ref, alt, loc_ref);
   
    Reads reads, targets, candidates, non_targets, undetermined;
    
    // read parsing
    prep_reads(
        reads, read_names, are_reverse, cigar_strs,
        aln_starts, aln_ends, read_seqs, quals, mapqs, are_from_first_bam
    ); 
       
    Contig contig;
    bool is_retargeted = false, is_non_supporting = false, is_mocked = false;
    if (target.is_substitute)
    {
        retarget_to_indel(
            reads, target, contig, user_params, loc_ref, 
            is_retargeted, is_non_supporting, is_mocked
        );
                
        if (is_retargeted) // process as indel
        {    
            target = prep_target(target.pos, target.ref, target.alt, loc_ref);   
        }
        else 
        {    
            if (is_non_supporting)
            {
                from_no_substitute_reads(target, contig, reads, non_targets);    
            }
            else
            {
                from_target_substitute_reads(
                    contig, reads, targets, non_targets, 
                    target.pos, user_params, is_mocked
                );
            }      
            
            rslt.report(contig, targets, non_targets, undetermined);
            
            return;
        }
    }
    
    // indel read annotation
    annotate_reads(reads, target, user_params, loc_ref, is_retargeted);  
    
    classify_reads(reads, targets, candidates, non_targets, user_params);
    
    // indel contig processing
    if (!targets.empty())
    {
        from_target_reads(
            contig, target,
            targets, candidates, non_targets, 
            undetermined, user_params, loc_ref
        );
    }
    else if (!candidates.empty()) 
    {   
        from_candidate_reads(
            contig, target,
            targets, candidates, non_targets,
            undetermined, user_params, loc_ref
        );
    }
    
    rslt.report(contig, targets, non_targets, undetermined, is_retargeted);
    if (is_retargeted) rslt.retarget_pos = target.pos;
}   
