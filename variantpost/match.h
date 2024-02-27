#ifndef MATCH_H
#define MATCH_H

#include <string>

#include "util.h"
#include "contig.h"

struct ShiftableSegment
{
    bool is_complete_tandem_repeat = false; 
    int n_tandem_repeats = 0; 
    std::string seq;
    std::string rep_unit = "";
    std::string fw_repeat_unit;
    std::string rv_repeat_unit;
    size_t unit_len = 0;
    int start = -1;
    int end = -1;
    int n_reps = 0;
};

void annot_shiftable_segment(
    ShiftableSegment& ss,
    const Variant& target, 
    const Contig& contig
);

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
);

void classify_cand_indel_read_2(
    Reads& targets,
    Reads& candidates,
    Reads& non_targets,
    Reads& undetermined,
    const std::vector<Variant>* p_decomposed,
    const Variant& target,
    const Contig& contig,
    const UserParams& user_params
);

#endif
