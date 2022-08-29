#ifndef ALIGNED_VARIANT
#define ALIGNED_VARIANT

#include <string>
#include <vector>

#include "pileup_parser.h"
#include "fasta/Fasta.h"

void process_aligned_target(const std::string & chrom,
                            FastaReference & fr,
                            std::string & contig,
                            int & target_pos,
                            std::string & target_ref,
                            std::string & target_alt,
                            std::string & repeat_unit, 
                            std::vector<ParsedRead> & targets,
                            std::vector<ParsedRead> & candidates,
                            std::vector<ParsedRead> & non_targets);

#endif
