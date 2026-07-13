#ifndef CONSENSUS_H
#define CONSENSUS_H

#include "util.h"
#include "reads.h"
#include "variant_types.h"

struct SearchResult;

void from_consensus_variant_list(SearchResult& rslt, LocalReference& loc_ref, 
                                 const int start, const int end, const Vars& consensus);

void variant_consensus(const std::vector<Read>& reads, const Ints& idx, Vars& consensus);

#endif
