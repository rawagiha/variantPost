#ifndef SUBSTITUTES_H
#define SUBSTITUTES_H

#include "reads.h"
#include "util.h"

void retarget_to_indel(
    Reads& reads,
    Variant& target,
    Contig& contig,
    const UserParams& user_params,
    LocalReference& loc_ref,
    bool& is_retargeted,
    bool& is_non_supporing,
    bool& is_mocked
);

void from_no_substitute_reads(
    const Variant& target,
    Contig& contig,
    Reads& reads,
    Reads& non_targets
);

void from_target_substitute_reads(
    Contig& contig,
    Reads& reads,
    Reads& targets,
    Reads& non_targets,
    const bool is_mocked
);


#endif 
