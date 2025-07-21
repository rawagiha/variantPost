#ifndef MATCH_H
#define MATCH_H

#include <bitset>
#include "util.h"
#include "pileup.h"

//------------------------------------------------------------------------------
struct NonMatch {
    NonMatch(int i_, char op_, int len_) : i(i_), op(op_), len(len_) { }; 
    int i = -1; char op = '\0'; int len = -1; 
};

//------------------------------------------------------------------------------
struct Consensus {
    int i = -1; int pos = -1; // index and genomic pos
    int cnt = 0; // times mapped
    CigarVec non_matches;
};


//------------------------------------------------------------------------------
// fs/fe: flanking start/end defined by 2-mer complexity
// es/ee: event start/end defined by left/right realignment and variant length
// 
// ............fs....es...****...ee.....fe............
// 
// es-ee is a subregion in fs-fe
//
// COVERAGE FLAGS (CF)
// flag 0: fs is covered 
// flag 1: es is covered 
// flag 2: ee is covered
// flag 3: fe is covered
// 
// MATCH FLAGS (MF)
// flag 0: non-matches (X/I/D/S) occured in fs-es
// flag 1: non-matches occurred in es-ee
// flag 2: non-matches occurred in ee-fe
typedef std::bitset<4> CF;
typedef std::bitset<3> MF;
void match_flag(Alignment& aln, const int fl_start, const int event_start,
                const int event_end, const int fl_end, CF& cf, MF& mf);


void match2haplotypes(Pileup& pileup, const std::vector<std::string>& read_seqs,
                      //const std::string& seq_t, const std::string& seq_nt0,
                      //const std::string& seq_nt1, const std::string& seq_nt2,
                      const UserParams& params);

#endif
