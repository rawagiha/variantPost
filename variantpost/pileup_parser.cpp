#include <map>
#include <cmath>
#include <random>
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <climits>
#include <algorithm>
#include <unordered_map>

#include "util.h"
#include "pileup_parser.h"
#include "fasta/Fasta.h"


//@function
//    classify read overlap with target locus
//@return
//    'A' : covering
//    'D' : non covering (but may contain target (e.g. del))
//    'X' : disqulified (within spliced-out interval)
//
//    applicable shiftable target only
//    'B' : covering with non-ref bases
//       -> may support the target
//    'C' : covering with ref-bases only
//       -> undetermined if this read supports the target
//-------------------------------------------------------------------------------
char classify_covering(
    int & covering_start, int & covering_end,
    const int lpos, const int pos, const int rpos,
    const bool is_shiftable,
    const int start_offset, const int end_offset,
    const std::vector<std::pair<int, int>> & aligned_segments,
    const std::vector<std::pair<int, int>> & skipped_segments,
    const std::vector<Variant> & variants
)
{
    
    // skipping reads
    if (!skipped_segments.empty()) 
    {
        for (const auto & segment : skipped_segments) 
        {
            if (segment.first < lpos && rpos < segment.second) 
            {
                return 'X';
            }
        }
    }
    
    for (const auto & segment : aligned_segments) 
    {
        covering_start = segment.first;
        covering_end = segment.second;
 
        if (!is_shiftable) 
        {
            if (is_ascending(covering_start, pos, covering_end)) 
            {
                return 'A';
            }
        } 
        else 
        {
            // both bounded
            if (covering_start <= lpos && rpos + 1 <= covering_end) 
            {
                return 'A';
            }
            // left unbounded
            else if (is_ascending(lpos + 1, covering_start, rpos)) 
            {
                // first & lt-clipped
                if ((&segment == &aligned_segments.front()) && start_offset) 
                {
                    return 'B';
                }
                // variants exist between start and rpos
                for (const auto & variant : variants) 
                {
                    if (is_ascending(covering_start, variant.pos, rpos)) 
                    {
                        return 'B'; // repeats broken
                    }
                }
                return 'C';
            }
            // right unbounded
            else if (is_ascending(lpos, covering_end, rpos)) 
            {
                // last & rt-clipped
                if ((&segment == &aligned_segments.back()) && end_offset) 
                {
                    return 'B';
                }
                // variants exist between lpos and stop
                for (const auto & variant : variants) 
                {
                    if (is_ascending(lpos, variant.pos, covering_end)) 
                    {
                        return 'B'; // repeats broken
                    }
                }
                return 'C';
            }
        }
    }
    return 'D';
}


//@function 
//   measure distance to nearest variant that is not the-already-aligned target.
//   if the read has the-already-aligned target, read is flagged so.
//@return
//   dist to the nearest non-target variant (may be zero->multiallelelic)
//   INT_MAX if no variants in the read 
//----------------------------------------------------------------------------------------
int dist_to_non_target_variant(
    bool & has_target,
    bool & has_gteq_five_indel,
    bool & has_non_target_in_critical_region,
    int & target_aligned_pos,
    std::string & target_aligned_ref,
    std::string & target_aligned_alt,
    const int lpos, const int pos, const int rpos,
    const Variant & target,
    const std::vector<Variant> & variants,
    const int unspliced_local_reference_start,
    const int unspliced_local_reference_end,
    const std::unordered_map<int, char> & indexed_local_reference
)
{    
    if (variants.empty()) return INT_MAX;
       
    int target_start = target.pos;
    int target_end = target_start + (target.ref_len - 1);
    
    std::vector<int> dist;
    for (const auto & variant : variants) 
    {
        bool is_target = target.is_equivalent(variant, 
                                              unspliced_local_reference_start, 
                                              indexed_local_reference);
        
        if (!has_gteq_five_indel) 
        {
            has_gteq_five_indel = (std::max({variant.ref_len, variant.alt_len}) >= 5);
        }

        if (is_target) 
        {
            has_target = true;
            
            //as aligned (may differ from input)
            target_aligned_pos = variant.pos;
            target_aligned_ref = variant.ref;
            target_aligned_alt = variant.alt;
        }
        else 
        {
            if (lpos <= variant.pos && variant.pos <= rpos) 
            {
                has_non_target_in_critical_region = true;
            }

            int event_region_start = variant.get_leftmost_pos(
                                        unspliced_local_reference_start, 
                                        indexed_local_reference
                                    );
            int event_region_end = variant.get_rightmost_pos(
                                        unspliced_local_reference_end,
                                        indexed_local_reference
                                    ) + (variant.ref_len - 1);

            if (event_region_end <= target_start) 
            {
                dist.push_back(target_start - event_region_end);
            }
            else if (target_end <= event_region_start) 
            {
                dist.push_back(event_region_start - target_end);
            }
            else 
            {
                int ld = std::abs(pos - event_region_start);
                int rd = std::abs(event_region_end - pos);
            
                if (ld < rd) dist.push_back(ld);
                else dist.push_back(rd);
            }
        }
    }
   
    if (dist.empty()) return INT_MAX;
    else return *std::min_element(dist.begin(), dist.end());
    
}

//@function
//   classify softclip location relative to target pos
//   clipping in non-target exon will be classified as unclipped
//@return 
//   'U': unclipped (non-target exon may be clipped)
//   'L': target pos is closer to the left clipping
//   'R': target pos is closer to the right clippping
//------------------------------------------------------------------------   
char classify_clip_pattern (
    int & dist_to_clip,
    const int pos, 
    const int aln_start,
    const int aln_end,
    const int start_offset,
    const int end_offset,
    const int covering_start,
    const int covering_end
)
{
    dist_to_clip = INT_MAX;

    //upclipped
    if (!start_offset && !end_offset) return 'U';
    
    // spliced side check
    bool is_left_spliced = (aln_start < covering_start) ? true : false;
    bool is_right_spliced = (covering_end < aln_end) ? true : false;
    
    int lt_dist = std::abs(pos - aln_start);
    int rt_dist = std::abs(aln_end - pos);
    bool is_lefty = (lt_dist < rt_dist) ? true : false;

    // both side clipped
    if (start_offset && end_offset)
    {
        // target is middle exon
        if (is_left_spliced && is_right_spliced) 
        {
            return 'U';
        }
        else if (is_lefty) 
        {
            dist_to_clip = lt_dist;
            return 'L';            
        }
        else 
        {
            dist_to_clip = rt_dist;
            return 'R';
        }
    }
    //left-only clipped
    else if (start_offset) 
    {
        if (is_left_spliced) 
        {
            return 'U';
        }
        else 
        {
            dist_to_clip = lt_dist;
            return 'L';
        }
     }
    //right-only clipped
    else  
    {
        if (is_right_spliced) 
        {
            return 'U';
        }
        else 
        {
            dist_to_clip = rt_dist;
            return 'R';
        }
    }
}


//@function
//  classify if the read is worth further analysis
//  also check for possible involvement of complex event
//@return
//  'A': has target variant (already aligned)
//  'B': may have target (worth analyzing)
//  'N': not worth further check
//----------------------------------------------------------------------------
char classify_local_pattern(
    bool & may_be_complex,
    bool & has_non_target_in_critical_region,
    char & clip_ptrn,
    int & target_aligned_pos,
    std::string & target_aligned_ref,
    std::string & target_aligned_alt,
    const char covering_ptrn,
    const std::string & ref_seq,
    const std::string & read_seq,
    const int lpos,
    const int pos,
    const int rpos,
    const int aln_start,
    const int aln_end,
    const int start_offset, 
    const int end_offset,
    const int covering_start, 
    const int covering_end,
    const int ref_allele_len,
    const Variant & target,
    const std::vector<Variant> & variants,
    const int unspliced_local_reference_start,
    const int unspliced_local_reference_end,
    const std::unordered_map<int, char> & indexed_local_reference
)
{
    //trivial cases
    if (ref_seq == read_seq) return 'N';
    if (covering_ptrn == 'X') return 'N'; 
    if (covering_ptrn == 'B') return 'B';
       
    // covering but repeats are not bounded -> can't characterize target indel
    if (covering_ptrn == 'C') return 'N';

    // no overlapping with target lt-most and rt-most positions  
    if (covering_ptrn == 'D') 
    {
        if ((rpos + ref_allele_len) < covering_start || covering_end < lpos) 
        {
            return 'N';
        }
    }
    
    bool has_target = false;
    bool has_gteq_five_indel = false;
    int d2var = dist_to_non_target_variant(
                    has_target,
                    has_gteq_five_indel,
                    has_non_target_in_critical_region,
                    target_aligned_pos, target_aligned_ref, target_aligned_alt,
                    lpos, pos, rpos,  
                    target, variants,
                    unspliced_local_reference_start, 
                    unspliced_local_reference_end,
                    indexed_local_reference
                );
    
    int d2clip;
    clip_ptrn = classify_clip_pattern(
                    d2clip,
                    pos, 
                    aln_start, aln_end,
                    start_offset, end_offset,
                    covering_start, covering_end
                );
    
    int read_len = read_seq.length();
    double cplx_thresh = 10;
    double variant_free_dist = read_len / 5; //20%
    double clip_free_dist = read_len / 2; //50% 
    
    if (has_gteq_five_indel || d2var <= cplx_thresh || d2clip <= 2 * cplx_thresh) may_be_complex = true;
                   
    if (has_target) return 'A';
    if (!d2var) return 'B';  // partial match to target for cplx or multiallelic
    
    if (covering_ptrn == 'D') 
    {
        if (covering_end < pos) 
        {
            if (clip_ptrn == 'R' || d2var <= variant_free_dist) 
            {
                return 'B';
            }
        }
        else if (pos < covering_start) 
        {
            if (clip_ptrn == 'L' || d2var <= variant_free_dist) 
            {
                return 'B';
            }
        }
       
        return 'N';
    }
    else 
    {
        if (clip_ptrn == 'U') 
        {
            if (variant_free_dist < d2var) return 'N';
            else return 'B';
        }
        else 
        {
            if (variant_free_dist < d2var && clip_free_dist < d2clip) return 'N';
            else return 'B';
        }
    }
}


double low_qual_base_rate(
    const char base_qual_thresh, 
    const std::string & non_ref_quals, 
    const int read_len
)
{
    double dirty_base_cnt = 0.0;
    for (const char & q : non_ref_quals) 
    {
        if (q <= base_qual_thresh) dirty_base_cnt += 1.0;
    }
    return dirty_base_cnt / read_len; 
} 


std::string non_ref_pattern_str(
    const std::vector<Variant> & variants, 
    const std::vector<std::pair<int, int>> & skipped_segments, 
    const int aln_start, int start_offset, 
    const int aln_end, int end_offset
)
{
    std::string str;
    for (auto & variant : variants) 
    {
        str += (std::to_string(variant.pos) + "_" + variant.ref + "_" + variant.alt + ";");
    } 
    
    for (auto & segment : skipped_segments) 
    {
        str += ("spl=" + std::to_string(segment.first) 
                 + "-" + std::to_string(segment.second) + ";");
    }
    
    if (start_offset) str += ("lt_clip_end=" + std::to_string(aln_start - 1) + ";");
    if (end_offset) str += ("rt_clip_start=" + std::to_string(aln_end + 1));

    return str;
}
                                                                            

// ParsedRead constuctor
ParsedRead::ParsedRead() {};

ParsedRead::ParsedRead
( 
    const int unspliced_local_reference_start,
    const int unspliced_local_reference_end,
    const std::string & unspliced_local_reference,
    const char base_qual_thresh,
    const bool is_from_first,
    const std::string & read_name,
    const bool is_reverse,
    const std::string & cigar_string,
    const int aln_start,
    const int aln_end,
    const std::string & read_seq,
    const std::vector<int> & q,
    const int mapq,
    const Variant & target,
    const int ref_allele_len,
    const int lpos,
    const int pos,
    const int rpos,
    const bool is_shiftable,
    const std::unordered_map<int, char> & indexed_local_reference,
    const std::string & chrom,
    FastaReference & fr
) : is_from_first(is_from_first), read_name(read_name), is_reverse(is_reverse), 
    cigar_string(cigar_string), aln_start(aln_start),  aln_end(aln_end), 
    read_seq(read_seq), mapq(mapq)
{
    cigar_vector = to_cigar_vector(cigar_string);

    std::pair<char, int> first_cigar = cigar_vector[0];
    std::pair<char, int> last_cigar = cigar_vector.back();
    
    start_offset = (first_cigar.first == 'S') ? first_cigar.second : 0;
    read_start = aln_start - start_offset;
    end_offset = (last_cigar.first == 'S') ? last_cigar.second : 0;
    read_end = aln_end + end_offset;

    //set reference seq
    if (cigar_string.find('N') != std::string::npos) 
    {
        ref_seq = get_spliced_ref_seq(chrom, aln_start, cigar_vector, fr);
    
    } else 
    {
        ref_seq = get_unspliced_ref_seq(aln_start, aln_end,
                                        unspliced_local_reference_start, 
                                        unspliced_local_reference);
    }

    parse_splice_pattern(aligned_segments, skipped_segments, 
                         cigar_vector, read_start, read_end);

    
    if (ref_seq.empty()) 
    {
        covering_ptrn = 'X'; // X: disqualified flag 
    }
    else 
    {
        base_qualities = to_fastq_qual(q);
        
        //variants already aligned
        
        /*
        variants = find_mapped_variants(
                        aln_start, aln_end, 
                        ref_seq, read_seq, base_qualities,
                        cigar_vector, non_ref_quals
                    );
        */

        variants = find_variants(aln_start, ref_seq, read_seq, base_qualities, cigar_vector, non_ref_quals);
        
        
        //check covering pattern (indel position shift considered)
        covering_ptrn = classify_covering(
                            covering_start, covering_end,
                            lpos, pos, rpos, 
                            is_shiftable,
                            start_offset, end_offset,
                            aligned_segments, skipped_segments, variants
                        );
    }
    
    local_ptrn = classify_local_pattern(
                        may_be_complex,
                        has_non_target_in_critical_region,
                        clip_ptrn,
                        target_aligned_pos, target_aligned_ref, target_aligned_alt,
                        covering_ptrn, 
                        ref_seq, read_seq,
                        lpos, pos, rpos, 
                        aln_start, aln_end,
                        start_offset, end_offset, 
                        covering_start, covering_end,
                        ref_allele_len, target, variants,
                        unspliced_local_reference_start,
                        unspliced_local_reference_end,
                        indexed_local_reference
                    );          
    
    if (local_ptrn != 'N') 
    {
        //string to summarize local non-reference base & splice pattern 
        non_ref_ptrn_str = non_ref_pattern_str(
                                variants,
                                skipped_segments,
                                aln_start, start_offset, aln_end, end_offset
                            );  
        
        dirty_base_rate = low_qual_base_rate(
                                base_qual_thresh,
                                non_ref_quals, 
                                read_seq.size()
                            );  
    
    }
}


void parse_pileup
(
    std::vector<ParsedRead> & targets,
    std::vector<ParsedRead> & candidates,
    std::vector<ParsedRead> & non_targets,         
    FastaReference & fr,
    const std::string & chrom,
    const int pos, 
    const std::string & ref,
    const std::string & alt,
    const int mapping_quality_threhold,
    const int base_quality_threshold,
    const int unspl_loc_ref_start,
    const int unspl_loc_ref_end,
    const std::vector<bool> & is_from_first_bam,
    const std::vector<std::string> & read_names,
    const std::vector<bool> & are_reverse,
    const std::vector<std::string> & cigar_strings,
    const std::vector<int> & aln_starts,
    const std::vector<int> & aln_ends,
    const std::vector<std::string> & read_seqs,
    const std::vector<std::vector<int>> & base_quals,
    const std::vector<int> & mapqs
)
{
    // prep reference 
    std::string unspl_loc_ref = fr.getSubSequence(
                                    chrom, 
                                    unspl_loc_ref_start - 1, 
                                    unspl_loc_ref_end - unspl_loc_ref_start + 1
                                ); 
    
    std::unordered_map<int, char> ref_dict = reference_by_position(
                                                unspl_loc_ref,
                                                unspl_loc_ref_start, unspl_loc_ref_end
                                            );
    
    // find target's shiftable range 
    Variant target = Variant(pos, ref, alt);
    target.left_aln(unspl_loc_ref_start, ref_dict);
    int ref_allele_len = target.ref_len;
    int lpos = target.get_leftmost_pos(unspl_loc_ref_start, ref_dict);
    int rpos = target.get_rightmost_pos(unspl_loc_ref_end, ref_dict);
    bool is_shiftable = (lpos != rpos) ? true : false;
    
    // parse reads
    char base_qual_thresh = static_cast<char>(base_quality_threshold + 33);
    size_t pileup_size = read_names.size(); 
    std::vector<ParsedRead> parsed_reads;
    for (size_t i = 0; i < pileup_size; ++i) 
    {
        parsed_reads.emplace_back(
                        unspl_loc_ref_start, unspl_loc_ref_end, unspl_loc_ref, base_qual_thresh,
                        is_from_first_bam[i], read_names[i], are_reverse[i], cigar_strings[i],
                        aln_starts[i], aln_ends[i], read_seqs[i], base_quals[i], mapqs[i], 
                        target, ref_allele_len, lpos, pos, rpos, is_shiftable, 
                        ref_dict, chrom, fr
                    );

    }

    // first pass read annotation
    for (const auto & read : parsed_reads) 
    {
        if (read.local_ptrn == 'A') 
        {   targets.push_back(read);
        }
        else
        {
            if (read.mapq >= mapping_quality_threhold)
            {
                if (read.local_ptrn == 'B') 
                {
                    candidates.push_back(read);
                }
                else if (read.local_ptrn == 'N' && read.covering_ptrn == 'A') 
                {
                    non_targets.push_back(read);
                }
            } 
        }
    }

}
