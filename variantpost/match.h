#ifndef MATCH_H
#define MATCH_H

#include <bitset>
#include "util.h"
//#include "pileup.h"

//------------------------------------------------------------------------------
struct NonMatch {
    NonMatch(const int i_, const string& ref, const string& alt);
    
    int i = -1; 
    std::string ref, alt; 

    bool is_snv = false, is_mnv = false, is_ins = false, is_del = false;
};

typedef std::vector<NonMatch> NMS;

void gap_grid(const UserParams& params, std::vector<Ints>& grid);

bool search_over_grid(const int start, LocalReference& loc_ref,
                      const std::string& refseq, const std::string& query,
                      const std::vector<Ints>& grid, const Variant& target);

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
/*
typedef std::bitset<4> CF;
typedef std::bitset<3> MF;
void match_flag(Alignment& aln, const int fl_start, const int event_start,
                const int event_end, const int fl_end, CF& cf, MF& mf);

void rerank_by_realn(Pileup& pileup, const Strs& read_seqs, 
                     const UserParams& params, Variant& target);

void match2haplotypes(Pileup& pileup, const Strs& read_seqs,
                      //const std::string& seq_t, const std::string& seq_nt0,
                      //const std::string& seq_nt1, const std::string& seq_nt2,
                      const UserParams& params);
*/
#endif
