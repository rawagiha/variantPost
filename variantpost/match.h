#ifndef MATCH_H
#define MATCH_H

#include <string>

#include "util.h"
#include "contig.h"

struct ShiftableSegment
{
    bool is_complete_tandem_repeat = false; 
    int n_tandem_repeats = 0; 
    std::string fw_repeat_unit;
    std::string rv_repeat_unit;
    size_t unit_len = 0;
    int boundary_start = -1;
    int boundary_end = -1;
};

void match_to_unaln_contig();

void annot_shiftable_segment(
    ShiftableSegment& ss,
    const Variant& target, 
    const UnalignedContig& u_contig
);

void classify_cand_indel_reads(
    Reads& candidates,
    Reads& non_targets,

    Reads& lt_matches,
    Reads& mid_matches,
    Reads& rt_matches,

    Reads& undetermined,
   
    const UnalignedContig& u_contig,
    const ShiftableSegment& ss,
    const UserParams& user_params
);

#endif
