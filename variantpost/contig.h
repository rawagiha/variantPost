#ifndef CONTIG_H
#define CONTIG_H

#include "util.h"
#include "merge.h"
#include "reads.h"

struct Contig
{
    std::string seq;
    std::string quals;
    std::string ref_seq; 
    Coord coordinates;
    
    bool by_kmer_suggestion = false;
    size_t n_seeds = 0;
    
    //mocked from input
    bool is_mocked = false;
    std::string mocked_seq;
    std::string m_quals;
    std::string mocked_ref;
    Coord mocked_coord;
    
    size_t mock_len;
    size_t mock_lt_end_idx;
    size_t mock_lt_len;
    size_t mock_mid_len;
    size_t mock_rt_start_idx;
    size_t mock_rt_len;
     
    
    //contig lyaout
    size_t len;

    size_t lt_end_idx;
    size_t lt_len;
   
    size_t mid_len;

    size_t rt_start_idx;
    size_t rt_len;
    
    //alignments
    std::vector<int> pos_by_seq_idx{};
    std::vector<int> positions{};
    std::vector<int> skip_starts{};
    std::vector<int> skip_ends{};
    std::vector<std::string> ref_bases{};
    std::vector<std::string> alt_bases{};
    std::vector<std::string> base_quals{}; 
    
    Contig();

    Contig(const Variant& target);

    void furnish(
        const Seq& merged_target_reads,
        const Variant& target
    );
};


void make_contig(
    Contig& contig,
    const Variant& target, 
    Reads& targets, 
    const UserParams& user_param,
    LocalReference& loc_ref
);


void prefilter_candidates(
    Contig& contig,
    Reads& candidates,
    Reads& non_targets,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
);


void prefilter_cplx_candidates(
    Contig& contig,
    Reads& candidates,
    Reads& non_targets,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref,
    std::vector<Variant>& decomposed
);


void suggest_contig(
    const Variant& target,
    Contig& contig,
    Reads& candidates,
    const UserParams& user_params,
    LocalReference& loc_ref
);


void extend_contig(
    const char eval,
    Contig& contig,
    Reads& lt_matches,
    Reads& rt_matches,
    LocalReference& loc_ref
);


#endif
