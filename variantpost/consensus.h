#ifndef CONSENSUS_H
#define CONSENSUS_H

#include "util.h"
/*
struct LocusAlignments {
    std::vector<std::string_view> refs; 
    std::vector<std::string_view> alts;   
};*/

struct Consensus {
    int start = INT_MAX; int end = INT_MIN;
    std::string ref_dot = ".";
    std::map<int, Strs> aln;
    
    // atcga
    // ..A..  pos = [1, 2, 3, 4, 5]
    //  .A... cov = [1, 2, 2, 2, 1]
    //        ref = [A, T, C, G, A]
    //        alt = [A, T, A, G, A]
    Ints pos; // 1-based genomic coordinates
    Ints cov; // how many reads covered the locus to make consensus 
    Strs ref, alt; // consensus ref/alt. ref and alt may be identical
     
    // make consensus data from variant list
    void _from_variants(const int start_, const int end_, const Vars& vars, 
                        UserParams& params, LocalReference& loc_ref);
};


#endif
