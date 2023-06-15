#ifndef CONTIG_H
#define CONTIG_H

#include "util.h"
#include "merge.h"
#include "reads.h"

struct UnalignedContig
{
    std::string seq;
    std::string quals;
    std::string ref_seq; 
    Coord coordinates;
     
    size_t len;

    size_t lt_end_idx;
    size_t lt_len;
   
    size_t mid_len;

    size_t rt_start_idx;
    size_t rt_len;
    
    void furnish(
        const Seq& merged_target_reads,
        const Variant& target
    );
};


void make_unaln_contig(
    UnalignedContig& u_contig,
    const Variant& target, 
    Reads& targets, 
    const UserParams& user_param,
    LocalReference& loc_ref
);


void prefilter_candidates(
    Reads& candidates,
    Reads& non_targets,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
);


void suggest_unaln_contig(
    UnalignedContig& u_contig,
    Reads& candidates,
    const UserParams& user_params,
    LocalReference& loc_ref
);


#endif
