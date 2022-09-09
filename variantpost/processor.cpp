#include <vector>
#include <string>

#include "util.h"
#include "processor.h"
#include "pileup_parser.h"
#include "aligned_variant.h"

pp::ProcessedPileup prepare_processed_rslt(std::string & contig,
                                       int target_pos,
                                       std::string & target_ref,
                                       std::string & target_alt,
                                       const std::vector<ParsedRead> & targets,
                                       const std::vector<ParsedRead> & non_targets);


pp::ProcessedPileup::ProcessedPileup() {}

pp::ProcessedPileup::ProcessedPileup
(
    const std::string & contig,
    const int target_pos,
    const std::string ref,
    const std::string alt,
    std::vector<std::string> & read_names,
    std::vector<bool> & are_reverse,
    std::vector<bool> & are_target,
    std::vector<bool> & are_from_first_bam
) : contig(contig), target_pos(target_pos), ref(ref), alt(alt), 
    read_names(read_names), are_reverse(are_reverse),
    are_target(are_target), 
    are_from_first_bam(are_from_first_bam)                
{}

pp::ProcessedPileup  pp::process_pileup(
    const std::string & fastafile,
    const std::string & chrom,
    const int pos, 
    const std::string & ref,
    const std::string & alt,
    const int mapping_quality_threshold,
    const int base_quality_threshold,
    const double low_quality_base_rate_threshold,
    const int kmer_size,
    const int unspliced_local_reference_start,
    const int unspliced_local_reference_end,
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

    FastaReference fr;
    fr.open(fastafile);

    parse_pileup(targets, candidates, non_targets,
                 fr,
                 chrom, pos, ref, alt,
                 mapping_quality_threshold,
                 base_quality_threshold,
                 unspliced_local_reference_start,
                 unspliced_local_reference_end,
                 is_from_first_bam,
                 read_names,
                 are_reverse,
                 cigar_strings,
                 aln_starts,
                 aln_ends,
                 read_seqs,
                 quals,
                 mapqs);
    
    std::string contig = "";
    
    int target_pos = pos;
    std::string target_ref = ref;
    std::string target_alt = alt;
    Variant target(pos, ref, alt);
    
    std::string minimal_repeat = target.minimal_repeat_unit();
     
    if (targets.size() > 0) {
         process_aligned_target(chrom,
                                fr,
                                base_quality_threshold,
                                low_quality_base_rate_threshold,
                                kmer_size,
                                contig,
                                target_pos,
                                target_ref,
                                target_alt,
                                minimal_repeat,
                                targets,
                                candidates,
                                non_targets);
    }
    else if (candidates.size() > 0) {
        //candidate processor
    }
    
    ProcessedPileup prp = prepare_processed_rslt(contig,
                                                 target_pos,
                                                 target_ref, 
                                                 target_alt, 
                                                 targets, 
                                                 non_targets);
    //output preparer
    
    
    /*std::string _contig = "AATTCCGG";
    int _pos = 123;
    std::string _ref = "AT";
    std::string _alt = "A";
    std::vector<std::string> _read_names {"read1", "read2", "read3"};
    std::vector<bool> _are_rev {true, false, false};
    std::vector<bool> _are_target {false, false, true};
    std::vector<bool> _are_from {false, true, false};
    pp::ProcessedPileup prp {_contig, _pos, _ref, _alt, _read_names, _are_rev, _are_target, _are_from};
    */
    return prp;
}

pp::ProcessedPileup prepare_processed_rslt(std::string & contig,
                                       int target_pos,
                                       std::string & target_ref,
                                       std::string & target_alt,
                                       const std::vector<ParsedRead> & targets,
                                       //const std::vector<ParsedRead> & candidates, //this shouldn't exit at this stage
                                       const std::vector<ParsedRead> & non_targets)
{
    std::vector<std::string> read_names;
    std::vector<bool> are_reverse;
    std::vector<bool> are_target;
    std::vector<bool> are_from_first_bam;   
       
    for (const auto & read : targets) {
        read_names.push_back(read.read_name);
        are_reverse.push_back(read.is_reverse);
        are_target.push_back(true);
        are_from_first_bam.push_back(read.is_from_first);
    }
    
    for (const auto & read : non_targets) {
        read_names.push_back(read.read_name);
        are_reverse.push_back(read.is_reverse);
        are_target.push_back(false);
        are_from_first_bam.push_back(read.is_from_first);   
    }

    pp::ProcessedPileup prp {contig, target_pos, target_ref, target_alt, read_names, are_reverse, are_target, are_from_first_bam};
    
    return prp;
}        
