#ifndef READ_H
#define READ_H

#include <climits>
#include <bitset>
#include "util.h"

//------------------------------------------------------------------------------
// features at read level
struct Read {   
    //--------------------------------------------------------------------------
    // constructor from inputs passed by Cython wrapper
    Read(std::string_view name, const bool is_reverse, const std::string&  cigar_str,
         const int aln_start, const int aln_end, std::string_view seq,
         const std::string_view quals, const bool is_from_first_bam);
     
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
    void parseLocalPattern(LocalReference& loc_ref, const Variant& target);
    
    //--------------------------------------------------------------------------
    // find freq. of low quality bases within start/end region
    void qualityCheck(const int start, const int end, const UserParams& param); 
     
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
    
    //--------------------------------------------------------------------------
    // read identifier
    std::string_view name;

    //--------------------------------------------------------------------------
    // quality scores 
    std::string_view base_quals; // base quality string 
    //std::vector<Qual> non_ref_quals; // base quality for variant and clipped bases

    //--------------------------------------------------------------------------
    // genome coordinates (1-based) 
    int aln_start = 0; int aln_end = 0; // alignment start/end 
    int start_offset = 0; int end_offset = 0; // softclip len 
    int read_start = 0; int read_end = 0; // alignmet start/end extented by softclip len
    Coord idx2pos; // maps read idx to aligned genomic position
    Coord aligned_segments; // vector<pair<int, int>> start/end of aligned segments 
    Coord skipped_segments; // vector<pair<int, int>> start/end of skipped segments
    int covering_start = 0; // start of unspliced segment covering the target locus
    int covering_end = 0; // end of above segment
                                        
    //--------------------------------------------------------------------------
    // sequences 
    std::string_view seq; // reference to read seq data (std::string) from input
    std::string_view ref_seq; // reference to refseq data
    std::string spliced_ref_seq; // store refseq data if spliced as string (not view)
                                    
    //--------------------------------------------------------------------------
    // cigar 
    std::string_view cigar_str; // reference to cigar string from input
    CigarVec cigar_vector; // CigarVec vector<pair<char, int>>
    
    //-------------------------------------------------------------------------- 
    // boolean flags
    bool is_reverse = false; // true if aligned to the complementary stran
    bool is_control = false; // true for second BAM reads. always true for single BAM
    bool is_na_ref = false; // true if no ref_seq is available for this read
    bool is_ref = false; // true if seq is identical to ref_seq
    bool has_target = false; // true if supporitng the target
    bool has_local_events = false; // true if non_ref events with event_len (Variant)
    bool qc_passed = false; // true if local freq of dirty bases < user thresh
    bool fail_to_cover_flankings = false; // if true, classify as ambigous read
    bool is_stable_non_ref = false; // true if bounded by enough 2-mer diversity
    bool is_central_mapped = false; // true if mapped in 2nd of read len tertile

    //--------------------------------------------------------------------------
    // pattern keys
    char covering_ptrn = 'C'; // 'A':complete coverage, 'B":partial, 'C':none (default)   

    //--------------------------------------------------------------------------
    // variant info
    Vars variants;
    int target_idx = -1; // target variant idx in variants (vector<Variants>)
    int target_pos = -1;  // genomic pos in the input (possibly non-normalized)
    std::string target_ref = "N"; // ref allele (possibly non-normalized)
    std::string target_alt = "N"; // alt allele (possibly non-normalized)

    //--------------------------------------------------------------------------
    // signtures to compare variations    
    //std::string flanking_sig; // for variants in flanking regions
    //std::string variant_sig; // for all variants in read
    //std::string clip_sig; // for clipping pattern
    std::string splice_sig; // for splicing pattern
    std::string non_ref_sig; // for variants and clippings
    
    //--------------------------------------------------------------------------
    // metrics
    int dist_to_non_target = INT_MAX; // distance to nearest variant or clip
    int smer = 0, nmer = 0; // supporting and non-supporting kmer count 
      
    //--------------------------------------------------------------------------
    // rank: 's' supporting target, 'n': non_supprting 'u': undetermined
    char rank = '\0'; 
};

#endif
