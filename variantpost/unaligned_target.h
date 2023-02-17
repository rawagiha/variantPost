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
                              std::vector<Read> & targets, 
                              std::vector<Read> & candidates, 
                              std::vector<Read> & non_targets,
                              const size_t kmer_size);
#endif 
