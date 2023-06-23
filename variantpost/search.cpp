#include <vector>
#include <string>
#include <chrono>

#include "eval.h"
#include "util.h"
#include "match.h"
#include "reads.h"
#include "search.h"
#include "contig.h"


SearchResult::SearchResult() {}



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

void SearchResult::report(
    const Contig& contig,
    const Reads& targets,
    const Reads& non_targets
    //const Reads& undetermined
)
{
    //size_t buff_size = targets.size() + non_targets.size() + undetermined.size();
    size_t buff_size = targets.size() + non_targets.size();
    
    read_names.reserve(buff_size);
    are_reverse.reserve(buff_size);
    target_statuses.reserve(buff_size);
    are_from_first_bam.reserve(buff_size);
    
    fill_read_info(targets, 1);
    fill_read_info(non_targets, 0);
    //fill_read_info(undetermined, -1);
    
    positions = contig.positions;
    ref_bases = contig.ref_bases;
    alt_bases = contig.alt_bases;
    skip_starts = contig.skip_starts;
    skip_ends = contig.skip_ends;
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
    const double lq_base_rate_thresh,
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
    const std::vector<std::string>& cigar_strs,
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
    
    UserParams user_params(
        mapq_thresh, 
        base_q_thresh, 
        lq_base_rate_thresh,
        match_score, 
        mismatch_penal, 
        gap_open_penal, 
        gap_ext_penal, 
        kmer_size, 
        local_thresh 
    );

    LocalReference loc_ref(fastafile, chrom, ref_start, ref_end);   

    Variant target = prep_target(pos, ref, alt, loc_ref);
    
    Reads reads;
    prep_reads(
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
    }
    
    std::cout << targets.size() << " " << non_targets.size() << std::endl;    
    
    rslt.report(contig, targets, non_targets); 
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
    std::cout << _eval << " " << std::endl;
     
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
    
    // 'A'
    transfer_vector(targets, lt_matches);
    transfer_vector(targets, mid_matches);
    transfer_vector(targets, rt_matches);
}



void swith_to_mock_layout(Contig& contig)
{
    contig.seq = contig.mocked_seq;
    contig.ref_seq = contig.mocked_ref;
    contig.coordinates = contig.mocked_coord;
    contig.len = contig.mock_len;
    contig.lt_end_idx = contig.mock_lt_end_idx;
    contig.lt_len = contig.mock_lt_len;
    contig.mid_len = contig.mock_mid_len;
    contig.rt_len = contig.mock_rt_len;
    contig.rt_start_idx = contig.mock_rt_start_idx;
       
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
        
    Reads undetermined;
    if (_eval != 'A')
    {
        swith_to_mock_layout(contig);
    }    
        
    classify_cand_indel_read_2(
        targets,
        candidates,
        non_targets,
        undetermined,
        target,
        contig,
        user_params
    );
}   
