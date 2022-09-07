#ifndef ALIGNED_VARIANT
#define ALIGNED_VARIANT

#include <string>
#include <vector>

#include "pileup_parser.h"
#include "fasta/Fasta.h"

typedef std::vector<ParsedRead> Reads;
typedef std::pair<std::string, std::string> SeqAndQual;

void process_aligned_target(const std::string & chrom,
                            FastaReference & fr,
                            const int base_quality_threshold,
                            const double low_quality_base_rate_threshold,
                            const int kmer_size,
                            std::string & contig,
                            int & target_pos,
                            std::string & target_ref,
                            std::string & target_alt,
                            std::string & repeat_unit, 
                            Reads & targets,
                            Reads & candidates,
                            Reads & non_targets);

#endif
