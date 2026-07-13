#ifndef READ_H
#define READ_H

#include <climits>
#include <bitset>
#include "util.h"


enum class CoveringPattern {
    // Target's left-aligned (lpos) and right-aligned pos (rpos) 
    // is within mapped segments + SOFTCLIP extension (read-start/end)
    Full,

    // One of target's lpos/rpos is  within mapped segment + SOFTCLIP extension 
    // with non-reference pattern (variants or softclipping)
    Partial,

    // Otherwise (including spliced reads completely skipping the region) 
    NoCover
};

//------------------------------------------------------------------------------
// features at read level
struct Read {   
    // constructor from inputs passed by Cython wrapper
    Read(std::string_view name, bool is_reverse, const std::string& cigar_str,
         int aln_start, int aln_end, std::string_view seq,
         std::string_view quals, bool is_from_first_bam);
     
    //--------------------------------------------------------------------------
    // setup reference sequence and coordinates 
    void setReference(LocalReference& loc_ref); 

    //--------------------------------------------------------------------------
    // parse CIGAR for variants
    void setVariants(LocalReference& loc_ref);
    
    //--------------------------------------------------------------------------
    // annotate how the read covers the target locus
    void parseCoveringPattern(LocalReference& loc_ref, const Variant& target);
        
    //--------------------------------------------------------------------------
    // find target and non-target variations
    void parseLocalPattern(LocalReference& loc_ref, 
                           const Variant& target, const UserParams& params);
    
    void checkByRepeatCount(const Variant& target, bool& has_excess_ins_hap);
    
    //--------------------------------------------------------------------------
    // find freq. of low quality bases within start/end region
    void qualityCheck(const int start, const int end, 
                      const char qual_thresh, const double freq_thresh); 
   
    //--------------------------------------------------------------------------
    // find strech of low qual bases from both ends
    void trimLowQualBases(const char qual_thresh);  
    
    //--------------------------------------------------------------------------
    // annotate variants/clip/splice patterns in std::string
    void setSignatureStrings(const UserParams& params);

    //--------------------------------------------------------------------------
    // flag unambigous non-ref alignments to be used as template
    void isStableNonReferenceAlignment(LocalReference& loc_ref);
    
    //--------------------------------------------------------------------------
    // test if the target locus is middle of read mapping
    void isCenterMapped(const Variant& target);
    
    //-------------------------------------------------------------------------
    // test if read has the target by pattern matching (complex indel/MNVs only)
    void hasTargetComplexSubstitute(const Variant& target);
    void hasTargetComplexIndel(LocalReference& loc_ref, const Variant& target);
    
    //-------------------------------------------------------------------------
    bool IsRefAt(const int pos) const;
    
    // 8-byte members (Pointers, string_view, std::string, std::vector)
    std::string_view name;
    std::string_view base_quals;
    std::string_view seq;
    std::string_view ref_seq;
    std::string_view cigar_str;

    std::string spliced_ref_seq;
    std::string target_ref = "N";
    std::string target_alt = "N";
    std::string splice_sig;
    std::string non_ref_sig;

    Ints idx2pos;
    Coord aligned_segments;  // excl N
    Coord mapped_segments;  // excl N + D
    Coord skipped_segments;
    Ints var_idx;
    Ints var_pos;
    CigarVec cigar_vector;
    Vars variants;

    // 4-byte members (int)
    int aln_start = 0;
    int aln_end = 0;
    int start_offset = 0;
    int end_offset = 0;
    int read_start = 0;
    int read_end = 0;
    int qs = -1, qe = -1;
    int covering_start = 0;
    int covering_end = 0;
    int target_idx = -1;
    int target_pos = -1;
    int flnk_v_cnt = 0;
    int dist_to_non_target = INT_MAX;
    int smer = 0, nmer = 0;
   
    CoveringPattern covering_pattern = CoveringPattern::NoCover;

    // 1-byte members (bool, char)
    bool is_reverse = false;
    bool is_control = false;
    bool is_na_ref = false;
    bool covered_in_clip = false;
    bool is_ref = false;
    bool has_target = false;
    bool has_local_events = false;
    bool has_local_clip = false;
    bool has_smaller_change = false;
    bool has_anti_pattern = false;
    bool has_positional_overlap = false;
    bool qc_passed = false;
    bool fail_to_cover_flankings = false;
    bool is_stable_non_ref = false;
    bool is_central_mapped = false;
    bool is_quality_map = false;
    bool ineffective_kmer = false;
    bool high_ambiguity = false;

    char rank = '\0';
};

#endif
