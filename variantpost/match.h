#ifndef MATCH_H
#define MATCH_H

#include <bitset>
#include "util.h"
#include "pileup.h"

struct SearchResult;

void gap_grid(const UserParams& params, std::vector<Ints>& grid);

bool search_over_grid(const int start, LocalReference& loc_ref,
                      const UserParams& params, const size_t n_vars,
                      const std::string& refseq, const std::string& query,
                      const std::vector<Ints>& grid, const Variant& target);

void check_match_pattern(Alignment& aln, std::bitset<3>& check_points,
                         const int fss, const int fse, 
                         const int ts, const int te,
                         const int fes, const int fee);

void match2haplotypes(Pileup& pileup, const Strs& read_seqs,
                      const UserParams& params);

void personalize(const Pileup& pileup, 
                 LocalReference& loc_ref, 
                 const UserParams& params, 
                 const Variant& target, SearchResult& rslt);

#endif
