#ifndef EVAL_H
#define EVAL_H

#include "util.h"
#include "contig.h"



void eval_by_aln(
    const UnalignedContig& u_contig,
    const UserParams& user_params,
    LocalReference& loc_ref
);

#endif
