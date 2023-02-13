#ifndef ALIGNED_TARGET
#define ALIGNED_TARGET

#include <string>
#include <vector>

#include "read_classifier.h"
#include "fasta/Fasta.h"

typedef std::vector<Read> Reads;
typedef std::pair<std::string, std::string> SeqAndQual;
/*
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
*/
void process_aligned_target(Variant & target,
                            FastaReference & fr, 
                            const int base_quality_threshold,
                            const double low_quality_base_rate_threshold, 
                            const int kmer_size,
                            std::string & _contig,
                            Reads & targets, Reads & candidates, Reads & non_targets);


#endif
