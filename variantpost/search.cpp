#include <vector>
#include <string>
#include <chrono>

#include "eval.h"
#include "util.h"
#include "match.h"
#include "reads.h"
#include "search.h"
#include "contig.h"




void from_target_reads(
    Contig& contig,
    const Variant& target,
    Reads& targets,
    Reads& candidates,
    Reads& non_targets,
    const UserParams& user_param,
    LocalReference& loc_ref
);


void from_candidate_reads(
    Contig& contig,
    const Variant& target,
    Reads& targets,
    Reads& candidates,
    Reads& non_targets,
    const UserParams& user_params,
    LocalReference& loc_ref
);

/*
ProcessedPileup prepare_processed_rslt(const Contig & contig,
                                       int target_pos,
                                       std::string & target_ref,
                                       std::string & target_alt,
                                       const std::vector<Read> & targets,
                                       const std::vector<Read> & candidates,
                                       const std::vector<Read> & non_targets);
*/

SearchResult::SearchResult() {}

/*
ProcessedPileup::ProcessedPileup
(
    const std::vector<int> & positions,
    const std::vector<std::string> & ref_bases,
    const std::vector<std::string> & alt_bases,
    const std::vector<std::string> & base_quals,
    const std::vector<int> & skip_starts,
    const std::vector<int> & skip_ends,
    const int target_pos,
    const std::string ref,
    const std::string alt,
    std::vector<std::string> & read_names,
    std::vector<bool> & are_reverse,
    std::vector<int> & target_statuses,
    std::vector<bool> & are_from_first_bam
) : positions(positions), ref_bases(ref_bases), alt_bases(alt_bases), base_quals(base_quals),
    skip_starts(skip_starts), skip_ends(skip_ends),
    target_pos(target_pos), ref(ref), alt(alt), 
    read_names(read_names), are_reverse(are_reverse),
    target_statuses(target_statuses), 
    are_from_first_bam(are_from_first_bam)                
{}
*/

void pack_user_params(
    UserParams& user_params,

    const int mapq_thresh,
    const int base_q_thresh,
    const double low_q_base_rate_thresh,
    const int match_score,
    const int mismatch_penal,
    const int gap_open_penal,
    const int gap_ext_penal,
    const int kmer_size,
    const int local_thresh
)
{
    user_params.mapq_thresh = mapq_thresh;
    user_params.base_q_thresh = static_cast<char>(base_q_thresh + 33);
    user_params.lq_rate_thresh = low_q_base_rate_thresh;
    user_params.match_score = match_score;
    user_params.mismatch_penal = mismatch_penal;
    user_params.gap_open_penal = gap_open_penal;
    user_params.gap_ext_penal = gap_ext_penal;
    user_params.kmer_size = kmer_size;
    user_params.local_thresh = local_thresh;
}
      

void pack_reference_info(
    LocalReference& loc_ref,

    const std::string& fastafile,
    const std::string& chrom,
    const int ref_start,
    const int ref_end
)
{
    loc_ref.chrom = chrom;
    loc_ref.start = ref_start;
    loc_ref.end = ref_end;
    loc_ref.set_up(fastafile);  
}


Variant pack_target_info(
    const int pos, 
    const std::string& ref, 
    const std::string& alt,
    LocalReference& loc_ref
)
{   
    Variant target = Variant(pos, ref, alt);
    target.set_leftmost_pos(loc_ref);
    target.set_rightmost_pos(loc_ref); 
    target.is_shiftable = (target.lpos != target.rpos) ? true : false;
    
    return target;    
}


void pack_read_info(
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
            read_names[i], 
            are_reverse[i], 
            cigar_strs[i], 
            aln_starts[i], 
            aln_ends[i], 
            read_seqs[i], 
            quals[i], 
            mapqs[i], 
            are_from_first_bam[i]
        ); 
    }             
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
    const double low_q_base_rate_thresh,
    const int match_score,
    const int mismatch_penal,
    const int gap_open_penal,
    const int gap_ext_penal,
    const int kmer_size,
    const int local_thresh,
    const int ref_start, //local ref start/end
    const int ref_end,   //defined by user's window choice
    
    //reads mapped to the region of interest 
    const std::vector<std::string>& read_names,
    const std::vector<bool>& are_reverse,
    const std::vector<std::string> & cigar_strs,
    const std::vector<int>& aln_starts,
    const std::vector<int>& aln_ends,
    const std::vector<std::string>& read_seqs,
    const std::vector<std::vector<int>>& quals,
    const std::vector<int>& mapqs,
    const std::vector<bool>& are_from_first_bam)
{  
    
    // do input validation at python ends
    // 1) no fetched reads, 
    // 2) undefined variants...
    // 3) ref_start/ref_end must be non-N region
    
    // packing data from Python 
    UserParams user_params;
    pack_user_params(
        user_params, 
        mapq_thresh, 
        base_q_thresh, 
        low_q_base_rate_thresh,
        match_score, 
        mismatch_penal, 
        gap_open_penal, 
        gap_ext_penal, 
        kmer_size, 
        local_thresh 
    );

    LocalReference loc_ref;   
    pack_reference_info(
        loc_ref,
        fastafile, 
        chrom, 
        ref_start, 
        ref_end
    );

    Variant target = pack_target_info(pos, ref, alt, loc_ref);
    
    Reads reads = {};
    pack_read_info(
        reads, 
        read_names, 
        are_reverse, 
        cigar_strs, 
        aln_starts, 
        aln_ends, 
        read_seqs, 
        quals, 
        mapqs, 
        are_from_first_bam
    );
       
    
    // read processing
    annotate_reads(
        reads, 
        target, 
        user_params, 
        loc_ref
    );  
    
    Reads targets, candidates, non_targets;
    classify_reads(reads, targets, candidates, non_targets, user_params);
    
    // contig processing
    Contig contig;
    if (!targets.empty())
    {
        from_target_reads(
            contig,
            target,
            targets,
            candidates,
            non_targets,
            user_params,
            loc_ref
        );
        
        std::cout << targets.size() << " " << candidates.size() << " " << non_targets.size() << std::endl;
        return;    
    }
    else if (!candidates.empty()) 
    {   
        from_candidate_reads(
            contig,
            target,
            targets,
            candidates,
            non_targets,
            user_params,
            loc_ref
        );

        return;  
    }
        
    //no resutl 
     
    return;
}   



void from_target_reads(
    Contig& contig,
    const Variant& target,
    Reads& targets,
    Reads& candidates,
    Reads& non_targets,
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
        transfer_vector(non_targets, targets);
        transfer_vector(non_targets, candidates); 
        return;      
    }

    ShiftableSegment ss;
    annot_shiftable_segment(ss, target, contig);   
     
    Reads lt_matches, mid_matches, rt_matches, undetermined;    
    classify_cand_indel_reads(
        candidates,
        non_targets, 
            
        lt_matches,
        mid_matches,
        rt_matches,

        undetermined,
        contig,
        ss,
        user_params
    );
        
    if (_eval == 'B')
    {
        //ext
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
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    prefilter_candidates(
        contig,
        candidates,
        non_targets,
        target,
        user_params,
        loc_ref
    );
        
    if (candidates.empty()) return;
        
    if (target.is_complex)
    {
        //cplx input -> reduce to from_target_reads
    }
    
    suggest_contig(contig, candidates, user_params, loc_ref); 
       
    //may happen if all candidates are of low-qual bases
    if (contig.seq.empty()) return;
        
    char _eval = eval_by_aln(contig, target, user_params, loc_ref);
        
    if (_eval == 'A')
    {
        //-> do count again?
    }
    else
    {
        //exact match only? 
    }
}   
