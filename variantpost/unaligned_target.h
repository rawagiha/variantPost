#ifndef UNALIGNED_TARGET
#define UNALIGNED_TARGET

#include <string>
#include <vector>
#include "util.h"
#include "read_classifier.h"
#include "fasta/Fasta.h"


void process_unaligned_target(Variant & target, FastaReference & fr, 
                              const int base_quality_threshold,
                              const double low_quality_base_rate_threshold, 
                              const int match_score,
                              const int mismatch_penalty,
                              const int gap_open_penalty,
                              const int gap_extention_penalty,
                              const size_t kmer_size,
                              const int unspl_loc_ref_start,
                              const std::unordered_map<int, char> & indexed_local_reference,
                              Contig & contig,
                              Reads & targets, 
                              Reads & candidates, 
                              Reads & non_targets);
#endif 
