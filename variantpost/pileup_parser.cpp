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
#include "swlib.h"
#include "pileup_parser.h"
#include "fasta/Fasta.h"


// check covering pattern
char classify_covering(const int lpos, const int pos, const int rpos,
                       const bool is_shiftable, 
                       const int start_offset, const int end_offset,
                       int & covering_start, int & covering_end,
                       const std::vector<std::pair<int, int>> & unspl_segments,
                       const std::vector<std::pair<int, int>> & spl_segments,
                       const std::vector<Variant> & variants);



// check local non-ref base pattern
char classify_local_pattern(bool & may_be_complex,
                            char & clip_ptrn,
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
                            const std::unordered_map<int, char> & indexed_local_reference);


// check base quals
double dirty_rate(const char base_qual_thresh, const std::string & base_qualities);

// non-reference base patterns to string
std::string non_ref_pattern_str(std::vector<Variant> & variants, 
                                std::vector<std::pair<int, int>> & spliced_segments, 
                                int aln_start, int start_offset, 
                                int aln_end, int end_offset);


// select seed reads for assembly
std::vector<std::pair<std::string, std::string>> get_seed_reads
(
    const std::vector<ParsedRead> & targets, 
    size_t n = 5
);

// ParsedRead constuctor
ParsedRead::ParsedRead
( 
    const int unspliced_local_reference_start,
    const int unspliced_local_reference_end,
    const std::string & unspliced_local_reference,
    const char base_qual_thresh,
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
) : read_name(read_name), is_reverse(is_reverse), cigar_string(cigar_string),
    aln_start(aln_start),  aln_end(aln_end), read_seq(read_seq), mapq(mapq)
{
    cigar_vector = to_cigar_vector(cigar_string);

    std::pair<char, int> first_cigar = cigar_vector[0];
    std::pair<char, int> last_cigar = cigar_vector.back();
    
    start_offset = (first_cigar.first == 'S') ? first_cigar.second : 0;
    read_start = aln_start - start_offset;
    end_offset = (last_cigar.first == 'S') ? last_cigar.second : 0;
    read_end = aln_end + end_offset;

    //ref seq
    if (cigar_string.find('N') != std::string::npos) {
        ref_seq = get_spliced_ref_seq(chrom, aln_start, cigar_vector, fr);
    
    } else {
        ref_seq = get_unspliced_ref_seq(aln_start, aln_end,
                                        unspliced_local_reference_start, 
                                        unspliced_local_reference);
    }

    //splice pattern
    std::vector<std::pair<int, int>> un_spliced_segments;
    std::vector<std::pair<int, int>> spliced_segments;
    parse_splice_pattern(un_spliced_segments, spliced_segments, 
                         cigar_vector, read_start, read_end);

    covering_start = 0;
    covering_end = 0;
    if (ref_seq.empty()) {
        covering_ptrn = 'X'; // X: disqualified flag 
    }
    else {
        //variants already aligned
        variants = find_mapped_variants(aln_start, aln_end, 
                                        ref_seq, read_seq, 
                                        cigar_vector, 
                                        unspliced_local_reference_start, 
                                        unspliced_local_reference_end,
                                        indexed_local_reference);
        
        //check covering pattern (indel position shift considered)
        covering_ptrn = classify_covering(lpos, pos, rpos, is_shiftable,
                                          start_offset, end_offset,
                                          covering_start, covering_end,
                                          un_spliced_segments, spliced_segments, 
                                          variants);
    }
    
    //classify local non-reference base pattern 
    may_be_complex = false;
    local_ptrn = classify_local_pattern(may_be_complex,
                                        clip_ptrn,
                                        covering_ptrn, 
                                        ref_seq, read_seq,
                                        lpos, pos, rpos, 
                                        aln_start, aln_end,
                                        start_offset, end_offset, 
                                        covering_start, covering_end,
                                        ref_allele_len,
                                        target,
                                        variants,
                                        unspliced_local_reference_start,
                                        unspliced_local_reference_end,
                                        indexed_local_reference);          
                                                                                                                                                 
    base_qualities = to_fastq_qual(q);
    
    if (local_ptrn != 'N') {
        //string to summarize local non-reference base & splice pattern 
        non_ref_ptrn_str = non_ref_pattern_str(variants,
                                               spliced_segments,
                                               aln_start, start_offset,
                                               aln_end, end_offset);    
    }
    
    if (local_ptrn == 'B') {
        dirty_base_rate = dirty_rate(base_qual_thresh, base_qualities);
    }
    else dirty_base_rate = 0.0;
}


void parse_pileup
(
    std::vector<ParsedRead> & targets,
    std::vector<ParsedRead> & candidates,
    std::vector<ParsedRead> & non_targets,         
    const std::string & fastafile,
    const std::string & chrom,
    int pos, 
    const std::string & ref,
    const std::string & alt,
    int base_quality_threshold,
    int unspliced_local_reference_start,
    int unspliced_local_reference_end,
    //const std::string & unspliced_local_reference,
    const std::vector<std::string> & read_names,
    const std::vector<bool> & are_reverse,
    const std::vector<std::string> & cigar_strings,
    const std::vector<int> & aln_starts,
    const std::vector<int> & aln_ends,
    const std::vector<std::string> & read_seqs,
    //const std::vector<std::string> & ref_seqs,
    const std::vector<std::vector<int>> & quals,
    const std::vector<int> & mapqs,
    const std::vector<bool> & is_from_first_bam 
)
{
    // prep reference 
    FastaReference fr;
    fr.open(fastafile);
    std::string unspliced_local_reference = fr.getSubSequence(chrom, 
                                                              unspliced_local_reference_start - 1, 
                                                              unspliced_local_reference_end - unspliced_local_reference_start + 1); 
    
    std::unordered_map<int, char> aa = reference_by_position(unspliced_local_reference,
                                                             unspliced_local_reference_start,
                                                             unspliced_local_reference_end);
    
    // target variant
    Variant target = Variant(pos, ref, alt);
    target.left_aln(unspliced_local_reference_start, aa);
    int ref_allele_len = target.ref_len;
    int lpos = target.get_leftmost_pos(unspliced_local_reference_start, aa);
    int rpos = target.get_rightmost_pos(unspliced_local_reference_end, aa);
    bool is_shiftable = (lpos != rpos) ? true : false;
                
    
    // parse reads
    char base_qual_thresh = static_cast<char>(base_quality_threshold + 33);
    size_t pileup_size = read_names.size(); 
    std::vector<ParsedRead> parsed_reads;
    for ( size_t i = 0; i < pileup_size; ++i ) {
        parsed_reads.emplace_back( unspliced_local_reference_start,
                                   unspliced_local_reference_end,
                                   unspliced_local_reference,
                                   base_qual_thresh,
                                   read_names[i],
                                   are_reverse[i],
                                   cigar_strings[i],
                                   aln_starts[i],
                                   aln_ends[i],
                                   read_seqs[i],
                                   //ref_seqs[i],
                                   quals[i],
                                   mapqs[i],
                                   target,
                                   ref_allele_len,
                                   lpos,
                                   pos,
                                   rpos,
                                   is_shiftable,
                                   aa,
                                   chrom,
                                   fr
                                 );

    }

    //std::vector<pileup::ParsedRead> targets;
    //std::vector<pileup::ParsedRead> candidates;
    //std::vector<pileup::ParsedRead> non_targets; 
    
    for (auto & read : parsed_reads) {
        if (read.local_ptrn == 'A') {
            targets.push_back(read);
        }
        else if (read.local_ptrn == 'B') {
            candidates.push_back(read);
        }
        else if ((read.covering_ptrn == 'A') & (read.local_ptrn == 'N')) {
            non_targets.push_back(read);
        } 
    }

    /*
    std::cout << targets.size() << " num of target reads ---" << std::endl;
    std::cout << candidates.size() << " num of worth-cheking reads ---" << std::endl;
    std::cout << non_targets.size() << " num of nontarget reads ---" << std::endl;
    
   
    int dirty_cnt = 0;
    if (candidates.size() > 0) {
        for (auto & read : candidates) {
            if (is_dirty(base_qual_thresh, read.base_qualities)) ++dirty_cnt;
        }
        
        std::cout << "dirties: " << dirty_cnt << std::endl; 
    } 
     
      
    
    std::string _contig;
    size_t target_pileup_size = targets.size();
    
    if (target_pileup_size > 0) {
        if (target_pileup_size == 1) {
            _contig = targets[0].read_seq;
        } 
        else {
            std::vector<std::pair<std::string, std::string>> seeds = get_seed_reads(targets);
            if (seeds.size() > 1) {
                std::vector<std::pair<std::string, std::string>> _others = {seeds.begin() + 1, seeds.end()};
                _contig = sw::flatten_reads(seeds[0], _others);
            }
            else {
                _contig = seeds[0].first; 
            }
        }

        //get_ref N case
        
        // _contig to ref
        // match candidates to _contig
            

    }
    else if (candidates.size() > 0) {
        // not straightforward case  
    }
    else {
        // not found case -> return analysis result
    }


    //[[read_names], [orientations], [are_countable], [are_targets], [are_from_bam1], [tar_pos], [tar_alt], [tar_ref], [contig]] 
    
    return "done";
    */   
}

// select gapped_seed

std::vector<std::pair<std::string, std::string>> get_seed_reads
(
    const std::vector<ParsedRead> & targets, 
    size_t n
)
{
    
    std::vector<std::string> non_ref_ptrns;
    for (auto & read : targets) {
        non_ref_ptrns.push_back(read.non_ref_ptrn_str);
    }
    const std::string commonest_ptrn = find_commonest_str(non_ref_ptrns);
    
    std::vector<std::pair<std::string, std::string>> seed_candidates;
    for (auto & read : targets) {
        if (read.non_ref_ptrn_str == commonest_ptrn) {
            seed_candidates.emplace_back(read.read_seq, read.base_qualities);
        }    
    }

    if (seed_candidates.size() > n) { 
        std::shuffle(seed_candidates.begin(), 
                     seed_candidates.end(), 
                     std::default_random_engine(123));
        std::vector<std::pair<std::string, std::string>> tmp(seed_candidates.begin(), 
                                                             seed_candidates.begin() + n);
        return tmp;
    }
    else return seed_candidates;   
}        


//@function
//    check how the read covers the locus of interest
//@return
//    'A' : covering
//    'D' : non covering (but may contain target (e.g. del))
//    'X' : disqulified (completely within intron)
//
//    applicable shiftable target only
//    'B' : covering with non-ref bases
//       -> may support the target
//    'C' : covering with ref-bases only
//       -> undetermined if this read supports the target
//-------------------------------------------------------------------------------
char classify_covering(const int lpos, const int pos, const int rpos,
                       const bool is_shiftable,
                       const int start_offset, const int end_offset,
                       int & covering_start, int & covering_end,
                       const std::vector<std::pair<int, int>> & unspl_segments,
                       const std::vector<std::pair<int, int>> & spl_segments,
                       const std::vector<Variant> & variants)
{
    
    // fly through (shiftable case considered)
    if (!spl_segments.empty()) {
        for (auto & segment : spl_segments) {
            if (segment.first < lpos && rpos < segment.second) {
                return 'X';
            }
        }
    }
   
    size_t i = 0;
    size_t last = unspl_segments.size() - 1;

    for (auto & segment : unspl_segments) {

        covering_start = segment.first;
        covering_end = segment.second;

        if (!is_shiftable) {
            if (ordered(covering_start, pos, covering_end)) {
                return 'A';
            }
        } 
        else {
            // repeats bounded
            if (covering_start <= lpos && rpos + 1 <= covering_end) {
                return 'A';
            }
            // repeats left open
            else if (ordered(lpos + 1, covering_start, rpos)) {
                // 1st segment & lt-clipped
                if ( i == 0 && start_offset) {
                    return 'B';
                }
                // variants exist between start and rpos
                for (auto & variant : variants) {
                    if (ordered(covering_start, variant.pos, rpos)) {
                        return 'B'; // repeats broken
                    }
                }
                return 'C';
            }
            // repeats right open
            else if (ordered(lpos, covering_end, rpos)) {
                // last segment & rt-clipped
                if (i == last && end_offset) {
                    return 'B';
                }
                // variants exist between lpos and stop
                for (auto & variant : variants) {
                    if (ordered(lpos, variant.pos, covering_end)) {
                        return 'B'; // repeats broken
                    }
                }
                return 'C';
            }

            i += 1;
        }

    }
    return 'D';
}


//@function 
//   measure the distance to the nearest variant that is not the-already-aligned target.
//   check if the read has the target that is already aligned
//@return
//   dist to the nearest non-target variant (may be zero->multiallelelic)
//   INT_MAX if no variants in the read 
//----------------------------------------------------------------------------------------
int dist_to_non_target_variant(bool & has_target,
                               bool & has_gteq_five_indel,
                               const int pos, 
                               const Variant & target,
                               const std::vector<Variant> & variants,
                               const int unspliced_local_reference_start,
                               const int unspliced_local_reference_end,
                               const std::unordered_map<int, char> & indexed_local_reference)
{    
    if (variants.empty()) return INT_MAX;
    
    
    // adjust for del and mnv
    int offset = (target.ref_len - 1);
    int p = target.pos;
    int _p = p + offset;
    
    std::vector<int> dist;
    for (auto & variant : variants) {
        bool is_target = target.is_equivalent(variant, 
                                              unspliced_local_reference_start, 
                                              indexed_local_reference);
        
        if (!has_gteq_five_indel) {
            has_gteq_five_indel = (std::max({variant.ref_len, variant.alt_len}) >= 5);
        }

        if (is_target) {
            has_target = true;
        }
        else {
            int lp = variant.get_leftmost_pos(unspliced_local_reference_start, 
                                          indexed_local_reference);
            int rp = variant.get_rightmost_pos(unspliced_local_reference_end, 
                                           indexed_local_reference);
            int _rp = rp + (variant.ref_len - 1); 
        
            if (_rp <= p) {
                dist.push_back(p - _rp);
            }
            else if (_p <= lp) {
                dist.push_back(lp - _p);
            }
            else {
                int ld = std::abs(pos - lp);
                int rd = std::abs(rp - pos);
            
                if (ld < rd) {
                    dist.push_back(ld);
                }
                else {
                    dist.push_back(rd);
                }
            }
        }
    }
   
    if (dist.empty()) { 
        return INT_MAX;
    }
    else {
        return *std::min_element(dist.begin(), dist.end());
    }
}

//@function
//   classify softclip location relative to target pos
//   clipping in non-target exon are annotate as unclipped
//@return 
//   'U': unclipped (non-target exon may be clipped)
//   'L': target pos is closer to the left clipping
//   'R': target pos is closer to the right clippping
//------------------------------------------------------------------------   
char classify_clip_pattern (int & dist_to_clip,
                            int pos, 
                            int aln_start,
                            int aln_end,
                            int start_offset,
                            int end_offset,
                            int covering_start,
                            int covering_end)
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
    if (start_offset && end_offset){
        // target is middle exon
        if (is_left_spliced && is_right_spliced) {
            return 'U';
        }
        else if (is_lefty) {
            dist_to_clip = lt_dist;
            return 'L';            
        }
        else {
            dist_to_clip = rt_dist;
            return 'R';
        }
    }
    //left-only clipped
    else if (start_offset) {
        if (is_left_spliced) {
            return 'U';
        }
        else {
            dist_to_clip = lt_dist;
            return 'L';
        }
     }
    //right-only clipped
    else  {
        if (is_right_spliced) {
            return 'U';
        }
        else {
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
char classify_local_pattern(bool & may_be_complex,
                            char & clip_ptrn,
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
    if (covering_ptrn == 'D') {
        if (((rpos + ref_allele_len) < covering_start) 
            || (covering_end < lpos)) {
            return 'N';
        }
    }
    
    bool has_target = false;
    bool has_gteq_five_indel = false;
    int d2var = dist_to_non_target_variant(has_target,
                                           has_gteq_five_indel,
                                           pos,
                                           target, 
                                           variants,
                                           unspliced_local_reference_start, 
                                           unspliced_local_reference_end,
                                           indexed_local_reference);
    
    int d2clip;
    clip_ptrn = classify_clip_pattern (d2clip,
                                       pos, 
                                       aln_start,
                                       aln_end,
                                       start_offset,
                                       end_offset,
                                       covering_start,
                                       covering_end);
    
    int read_len = read_seq.length();
    double cplx_thresh = 10;
    double variant_free_dist = read_len / 5; //20%
    double clip_free_dist = read_len / 2; //50% 
    
    if (has_gteq_five_indel 
        || d2var <= cplx_thresh 
        || d2clip <= 2 * cplx_thresh) may_be_complex = true;
                   
    if (has_target) return 'A';
    if (!d2var) return 'B';  // partial match to target for cplx or multiallelic
    
    if (covering_ptrn == 'D') {
        if (covering_end < pos) {
            if (clip_ptrn == 'R' || (d2var <= variant_free_dist)) {
                return 'B';
            }
        }
        else if (pos < covering_start) {
            if (clip_ptrn == 'L' || (d2var <= variant_free_dist)) {
                return 'B';
            }
        }
        
        return 'N';
    }
    else {
        if (clip_ptrn == 'U') {
            if (variant_free_dist < d2var) {
                return 'N';
            }
            else {
                return 'B';
            }
        }
        else {
            if ((variant_free_dist < d2var) && (clip_free_dist < d2clip)) {
                return 'N';
            }
            else {
                return 'B';
            }
        }
    }
}


double dirty_rate(const char base_qual_thresh, const std::string & base_qualities){
    double dirty_base_cnt = 0.0;
    for (const char & c : base_qualities) {
        if (c <= base_qual_thresh) dirty_base_cnt += 1.0;
    }
    return dirty_base_cnt / base_qualities.length(); 
} 

std::string non_ref_pattern_str(std::vector<Variant> & variants, 
                                std::vector<std::pair<int, int>> & spliced_segments, 
                                int aln_start, int start_offset, 
                                int aln_end, int end_offset)
{
    std::string str;
    for (auto & variant : variants) {
        str += (std::to_string(variant.pos) + "_" + variant.ref + "_" + variant.alt + ";");
    } 
    
    for (auto & segment : spliced_segments) {
        str += ("spl=" + std::to_string(segment.first) 
                 + "-" + std::to_string(segment.second) + ";");
    }
    
    if (start_offset) str += ("lt_clip_end=" + std::to_string(aln_start - 1) + ";");
    if (end_offset) str += ("rt_clip_start=" + std::to_string(aln_end + 1));

    return str;
}
                                                                            
