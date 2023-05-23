#include <vector>
#include <string>
#include <unordered_map>

#include "util.h"
#include "pileup_processor.h"
#include "read_classifier.h"
#include "aligned_target.h"
#include "unaligned_target.h"

ProcessedPileup prepare_processed_rslt(const Contig & contig,
                                       int target_pos,
                                       std::string & target_ref,
                                       std::string & target_alt,
                                       const std::vector<Read> & targets,
                                       const std::vector<Read> & non_targets);


ProcessedPileup::ProcessedPileup() {}

ProcessedPileup::ProcessedPileup
(
    const std::vector<int> & positions,
    const std::vector<std::string> & ref_bases,
    const std::vector<std::string> & alt_bases,
    const std::vector<std::string> & base_quals,
    const int target_pos,
    const std::string ref,
    const std::string alt,
    std::vector<std::string> & read_names,
    std::vector<bool> & are_reverse,
    std::vector<bool> & are_target,
    std::vector<bool> & are_from_first_bam
) : positions(positions), ref_bases(ref_bases), alt_bases(alt_bases),
    target_pos(target_pos), ref(ref), alt(alt), 
    read_names(read_names), are_reverse(are_reverse),
    are_target(are_target), 
    are_from_first_bam(are_from_first_bam)                
{}

ProcessedPileup  cpp_process_pileup(
    const std::string & fastafile,
    const std::string & chrom,
    const int pos, 
    const std::string & ref,
    const std::string & alt,
    const int mapping_quality_threshold,
    const int base_quality_threshold,
    const double low_quality_base_rate_threshold,
    const int match_score,
    const int mismatch_penalty,
    const int gap_open_penalty,
    const int gap_extention_penalty,
    const int kmer_size,
    const int unspl_loc_ref_start,
    const int unspl_loc_ref_end,
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
    std::vector<Read> targets = {};
    std::vector<Read> candidates = {};
    std::vector<Read> non_targets = {}; 

    //prep preference
    FastaReference fr;
    fr.open(fastafile); //closed in destructor
    std::string unspl_loc_ref = fr.getSubSequence(
                                    chrom,
                                    unspl_loc_ref_start - 1, 
                                    unspl_loc_ref_end - unspl_loc_ref_start + 1                            
                                );
    std::unordered_map<int, char> ref_dict = reference_by_position(    
                                                unspl_loc_ref,
                                                unspl_loc_ref_start, unspl_loc_ref_end
                                             );
    
    Variant target = Variant(pos, ref, alt, chrom);
    
    classify_reads(
        targets, candidates, non_targets,
        fr,
        target,
        mapping_quality_threshold,
        base_quality_threshold,
        unspl_loc_ref_start,
        unspl_loc_ref_end,
        unspl_loc_ref,
        ref_dict,
        is_from_first_bam,
        read_names,
        are_reverse,
        cigar_strings,
        aln_starts,
        aln_ends,
        read_seqs,
        quals,
        mapqs
    );
    
    
    if (target.is_substitute)
    {
        cout << "This SNV is found: " << targets.size() << " alt: " << candidates.size() << " cand: " << non_targets.size() << " non-tar" << std::endl;
        
        ProcessedPileup _prp_snv;
        return _prp_snv;
    }

    // no targets/candidates
    if (targets.empty() && candidates.empty())
    {
        std::cout <<"zero yannke" << std::endl;
        // return result here
    }
    
    
    //std::string contig = "";
    
    Contig contig;
    
    //std::string minimal_repeat = target.minimal_repeat_unit();
     
    
    if (targets.size() > 0) 
    {
         /*
         process_aligned_target(chromi,
                                fr,
                                base_quality_threshold,
                                low_quality_base_rate_threshold,
                                kmer_size,
                                contig,
                                target.pos,
                                target.ref,
                                target.alt,
                                minimal_repeat,
                                targets,
                                candidates,
                                non_targets);
       */
       //pass variables needed for equ check, 
        process_aligned_target(target, fr, 
                               base_quality_threshold, low_quality_base_rate_threshold, 
                               match_score, mismatch_penalty, gap_open_penalty, gap_extention_penalty,
                               kmer_size, unspl_loc_ref_start, ref_dict, 
                               contig, targets, candidates, non_targets);
    }
    else if (candidates.size() > 0) 
    {
        process_unaligned_target(target, fr, 
                                 base_quality_threshold, low_quality_base_rate_threshold, 
                                 kmer_size, unspl_loc_ref_start, ref_dict,
                                 targets, candidates, non_targets);
    }
    
    
    ProcessedPileup prp = prepare_processed_rslt(contig,
                                                 target.pos,
                                                 target.ref, 
                                                 target.alt, 
                                                 targets, 
                                                 non_targets);
    
    
    std::cout << "target N: " << targets.size() << std::endl;
    std::cout << "candidate N: " << candidates.size() << std::endl; 
    std::cout << "non_target N: " << non_targets.size() << std::endl;  
    /*
    for (auto & c : non_targets)
    {
        std::cout << c.read_name << " " << c.cigar_string << std::endl;
    }    
    */
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

ProcessedPileup prepare_processed_rslt(const Contig & contig,
                                       int target_pos,
                                       std::string & target_ref,
                                       std::string & target_alt,
                                       const std::vector<Read> & targets,
                                       //const std::vector<ParsedRead> & candidates, //this shouldn't exit at this stage
                                       const std::vector<Read> & non_targets)
{
    std::vector<int> positions;
    std::vector<std::string> ref_bases;
    std::vector<std::string> alt_bases;
    std::vector<std::string> base_quals;
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

    for (const auto & base_level_aln : contig.alignment)
    {
        positions.push_back(base_level_aln.genomic_pos);
        ref_bases.push_back(base_level_aln.ref_base);
        alt_bases.push_back(base_level_aln.alt_base);
        base_quals.push_back(base_level_aln.base_qual);
    }
    
    ProcessedPileup prp {positions, ref_bases, alt_bases, base_quals, target_pos, target_ref, target_alt, read_names, are_reverse, are_target, are_from_first_bam};
    
    return prp;
}        
