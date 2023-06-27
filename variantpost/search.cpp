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
    const Reads& non_targets,
    const Reads& undetermined
)
{
    size_t buff_size = 
        targets.size() + non_targets.size() + undetermined.size();
    //size_t buff_size = targets.size() + non_targets.size();
    
    read_names.reserve(buff_size);
    are_reverse.reserve(buff_size);
    target_statuses.reserve(buff_size);
    are_from_first_bam.reserve(buff_size);
    
    fill_read_info(targets, 1);
    fill_read_info(non_targets, 0);
    fill_read_info(undetermined, -1);
    
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
    //std::cout << "making" << std::endl;
    make_contig(
        contig,
        target, 
        targets, 
        user_params,
        loc_ref
    );
    
    //std::cout << contig.seq << " " << contig.seq.size() << std::endl;
    //std::cout << "eval" << std::endl;             
    char _eval = eval_by_aln(contig, target, user_params, loc_ref);
    //std::cout << "eval reult " << _eval  << std::endl;

    if (_eval == 'C')
    {
        transfer_vector(non_targets, targets);
        transfer_vector(non_targets, candidates); 
        return;      
    }

    //std::cout << "repeat" << std::endl;
    ShiftableSegment ss;
    annot_shiftable_segment(ss, target, contig);   
    
    //std::cout << "cand classify" << std::endl; 
    Reads lt_matches, mid_matches, rt_matches;    
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
    
    if (_eval != 'A')
    {
        extend_contig(_eval, contig, lt_matches, rt_matches, loc_ref);
        aln_extended_contig(contig, target, user_params, loc_ref);
    }
     
    transfer_vector(targets, lt_matches);
    transfer_vector(targets, mid_matches);
    transfer_vector(targets, rt_matches);
}


inline void switch_to_mock_layout(Contig& contig)
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

    contig.is_mocked = true;
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
    //std::cout << "prefilter" << std::endl;
    
    std::vector<Variant>* p_decomposed = NULL;
    std::vector<Variant> decomposed;
    if (target.is_complex)
    {
        prefilter_cplx_candidates(
            contig,
            candidates,
            non_targets,
            target,
            user_params,
            loc_ref,
            decomposed
        );
        
        //for (const auto& h: *p_decomposed) std::cout << h.pos << " " << h.ref << " " << h.alt << std::endl; 
    }
    else
    {
        prefilter_candidates(
            contig,
            candidates,
            non_targets,
            target,
            user_params,
            loc_ref
        );
    }
        
    if (!decomposed.empty())
    {
        p_decomposed = &decomposed;
    }
    
    if (candidates.empty()) return;
        
    
    //std::cout << "suggest" << std::endl;
    suggest_contig(contig, candidates, user_params, loc_ref);       
    
    //may happen if all candidates are of low-qual bases
    if (contig.seq.empty()) return;
    
    //std::cout << contig.seq << " " << contig.ref_seq << std::endl;
    
    //std::cout << "eval contig" << std::endl;
    char _eval = eval_by_aln(contig, target, user_params, loc_ref, p_decomposed);

    if (_eval != 'A')
    {
        switch_to_mock_layout(contig);
    }    
    
    //std::cout << _eval << std::endl;
    //std::cout << contig.seq << " " << contig.ref_seq << std::endl;

    //std::cout << "classify" << std::endl;
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
    
    //auto t1 = std::chrono::high_resolution_clock::now();
        
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
    
    //std::cout << "read prep " <<std::endl;
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
    
    //std::cout << "read anot " <<std::endl;
    // read processing
    annotate_reads(
        reads, 
        target, 
        user_params, 
        loc_ref
    );  
    
    Reads targets, candidates, non_targets, undetermined;
    classify_reads(reads, targets, candidates, non_targets, user_params);
    
    // contig processing
    //Contig contig(target);
    Contig contig;
    if (!targets.empty())
    {
        //std::cout << "conting from target " <<std::endl;
        from_target_reads(
            contig,
            target,
            targets,
            candidates,
            non_targets,
            undetermined,
            user_params,
            loc_ref
        );
    }
    else if (!candidates.empty()) 
    {   
        //std::cout << "conting from kmer " <<std::endl;
        from_candidate_reads(
            contig,
            target,
            targets,
            candidates,
            non_targets,
            undetermined,
            user_params,
            loc_ref
        );
    }
    
    rslt.report(contig, targets, non_targets, undetermined); 
    
    /*
    for (auto& read : non_targets)
    {
        std::cout << read.name << " " << read.cigar_str << std::endl;   
    }*/
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    //std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    //std::cout << ms_int.count() << "ms\n";
    //std::cout << ms_double.count() << "ms\n";
    
    /*
    size_t r = contig.positions.size();
    for (size_t i = 0; i < r; ++i)
    {
        std::cout << contig.positions[i] << " " << contig.ref_bases[i] << " " << contig.alt_bases[i]<< std::endl; 
    }*/
    
    return;
}   
