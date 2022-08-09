#include <vector>
#include <string>

#include "processor.h"
#include "pileup_parser.h"

std::string  pp::process_pileup(
    const std::string & fastafile,
    const std::string & chrom,
    int pos, 
    const std::string & ref,
    const std::string & alt,
    int base_quality_threshold,
    int unspliced_local_reference_start,
    int unspliced_local_reference_end,
    const std::vector<std::string> & read_names,
    const std::vector<bool> & are_reverse,
    const std::vector<std::string> & cigar_strings,
    const std::vector<int> & aln_starts,
    const std::vector<int> & aln_ends,
    const std::vector<std::string> & read_seqs,
    const std::vector<std::vector<int>> & quals,
    const std::vector<int> & mapqs,
    const std::vector<bool> & is_from_first_bam)
{
    std::vector<ParsedRead> targets;
    std::vector<ParsedRead> candidates;
    std::vector<ParsedRead> non_targets; 

    parse_pileup(targets, candidates, non_targets,
                         fastafile,
                         chrom, pos, ref, alt,
                         base_quality_threshold,
                         unspliced_local_reference_start,
                         unspliced_local_reference_end,
                         read_names,
                         are_reverse,
                         cigar_strings,
                         aln_starts,
                         aln_ends,
                         read_seqs,
                         quals,
                         mapqs,
                         is_from_first_bam);
    
    //[[read_names], [orientations], [are_countable], [are_targets], [are_from_bam1], [tar_pos], [tar_alt], [tar_ref], [contig]] 
    
    return "done";
}

