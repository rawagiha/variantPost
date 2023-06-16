#ifndef EVAL_H
#define EVAL_H

#include <climits>

#include "util.h"
#include "contig.h"

struct AlnResult
{
    //basic alignment info
    int genomic_start_pos;
    std::string seq;
    std::string ref_seq;
    std::string quals;
    
    std::string cigar_str;
    std::vector<std::pair<char, int>> cigar_vec;
    std::vector<Variant> variants;   
    
    //evaluations 
    int lt_ref_mapped_len = 0;
    int rt_ref_mapped_len = 0;
    bool lt_well_ref_mapped = false;
    bool rt_well_ref_mapped = false;
    bool is_well_ref_mapped = false;

    bool has_target = false;
    int confidence_level = 0;
    
    int idx_to_closest = 0;
    int dist_to_closest = INT_MAX;

    bool is_simple = true;
    
    bool is_passed = false;
    
    bool terminate_search = false;
};   



void eval_by_aln(
    const Contig& contig,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
);

#endif
