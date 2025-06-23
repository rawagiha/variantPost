#ifndef READ_H
#define READ_H

#include <climits>
#include <bitset>
#include "util.h"

struct Read; //forward declaration
typedef std::vector<Read> Reads;

//------------------------------------------------------------------------------
// receive inputs from Python and create Read 
void prep_reads(const std::vector<std::string>& read_names,
                const std::vector<bool>& are_reverse, 
                const std::vector<std::string>& cigar_strs,
                const std::vector<int>& aln_starts,
                const std::vector<int>& aln_ends,
                const std::vector<std::string>& read_seqs,
                const std::vector<std::vector<int>>& quals,
                const std::vector<int>& mapqs,
                const std::vector<bool>& are_from_first_bam,
                LocalReference& loc_ref, Variant& target,
                const UserParams& params, Reads& reads);

//------------------------------------------------------------------------------
void triage_reads(Reads& reads, Reads& targets, 
                  Reads& candidates, Reads& non_targets, UserParams& params);

struct Read
{   
    //--------------------------------------------------------------------------
    // constructor from inputs passed by Cython wrapper
    Read(std::string_view name, const bool is_reverse, std::string_view cigar_str,
         const int aln_start, const int aln_end, std::string_view seq,
         const std::vector<int>& quals, const int mapq, const bool is_from_first_bam);
     
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
    // annotate variants/clip/splice patterns in std::string
    void setSignatureStrings();

    //--------------------------------------------------------------------------
    // read identifier
    std::string_view name;

    //--------------------------------------------------------------------------
    // quality scores 
    int mapq = -1;  // mapping quality 
    std::string base_quals; // base quality string 
    std::vector<Qual> non_ref_quals; // base quality for variant and clipped bases

    //--------------------------------------------------------------------------
    // genome coordinates (1-based) 
    int aln_start = 0; int aln_end = 0; // alignment start/end 
    int start_offset = 0; int end_offset = 0; // softclip len 
    int read_start = 0; int read_end = 0; // alignmet start/end extented by softclip len
    Coord pos2idx; //TODO consider if this is really necessary TODO 
    Coord aligned_segments; // vector<pair<int, int>> start/end of aligned segments 
    Coord skipped_segments; // vector<pair<int, int>> start/end of skipped segments
    int covering_start = 0; // start of unspliced segment covering the target locus
    int covering_end = 0; // end of above segment
                                        
    //--------------------------------------------------------------------------
    // sequences 
    std::string_view seq; // reference to read seq data (std::string) from input
    std::string_view ref_seq; // reference to refseq data
    std::string spliced_ref_seq; // store refseq data if spliced as string (not view).  
                                    
    //--------------------------------------------------------------------------
    // cigar 
    std::string_view cigar_str; // reference to cigar string from input
    CigarVec cigar_vector; // CigarVec vector<pair<char, int>>
    
    //-------------------------------------------------------------------------- 
    // boolean flags
    bool is_reverse = false; // true if aligned to the complementary stran
    bool is_from_first_bam = false; // true if read is from the first BAM file
    bool is_na_ref = false; // true if no ref_seq is available for this read
    bool is_ref = false; // true if seq is identical to ref_seq
    bool has_target = false; // true if supporitng the target
    bool fail_to_cover_flankings = false; // if true, classify as ambigous read
    bool is_analyzable = false; // true if passed all filters
   
    //--------------------------------------------------------------------------
    // pattern keys
    char covering_ptrn = '\0'; // 'A':complete coverage, 'B":partial, 'C':none   

    //--------------------------------------------------------------------------
    // variant info
    std::vector<Variant> variants;
    int target_idx = -1; // target variant idx in variants (vector<Variants>)
    int target_pos = -1;  // genomic pos in the input (possibly non-normalized)
    std::string target_ref = "N"; // ref allele (possibly non-normalized)
    std::string target_alt = "N"; // alt allele (possibly non-normalized)

    //--------------------------------------------------------------------------
    // signtures to compare variations    
    std::string non_ref_sig; // for variants and clippings
    std::string splice_sig; // for splicing  
    
    //--------------------------------------------------------------------------
    // metrics
    int dist_to_non_target = INT_MAX; // distance to nearest variant or clip 
      
    /*
    //non reference event info
    //std::vector<Variant> variants;
    //int variants_target_idx = -1;
    //dist to closest non-target overlapping target shiftable segmen
    //int dist_to_non_target = INT_MAX; 
    int dist_to_clip = INT_MAX;
    //int target_pos = -1;  //actual pos (may not be normalized)
    //std::string target_ref = "N"; //actual alt
    //std::string target_alt = "N"; //actual ref
    //std::string non_ref_signature;
    //std::string splice_signature;
    bool incomplete_shift= false;
    bool may_be_complex = false; 
    int local_uniqueness = 0;
    int lt_end_matches = -1;
    int rt_end_matches = -1;
    
    //annotated patterns    
    //char covering_ptrn = '\0';
    char clip_ptrn = '\0';
    char local_ptrn = '\0';   
    
    // substitute-related
    char sb_ptrn = '\0';
    bool has_sb_target = false;
    size_t sb_read_idx = 0; 
    int sb_kmer_score = -1;
    int aln_score = 0;
    std::string sb_loc_signature;
    
    //metrics   
    double central_score = -1.0;
    double overall_lq_rate = 0.0;
    double nonref_lq_rate = 0.0;
    double local_lq_rate = -1.0;
    int kmer_score = -1;
    
    //other flags
    bool is_tight_covering = false;
    bool is_deprioritized = false;
    //bool is_contig_member = false;
    */
};

/*
typedef std::vector<Read> Reads;

void sort_by_start(Reads & reads);

void sort_by_kmer(Reads & reads);

//void annot_ref_seq(Read& read, LocalReference& loc_ref);

//void annot_splice_pattern(Read& read);

void annot_clip_pattern(Read& read, const Variant& target);

void annot_non_ref_signature(Read& read);

void eval_read_quality(Read& read, const UserParams& user_params);

int eval_loc_uniq(
    Read& read, 
    const UserParams& user_params, 
    const LocalReference& loc_ref
);

void annot_covering_ptrn(
    Read& read, 
    const Variant& target, 
    LocalReference& loc_ref,
    bool is_retargeted
);

void annotate_reads(
    Reads& reads, 
    const Variant& target, 
    const UserParams& user_params, 
    LocalReference& loc_ref,
    bool is_retargeted
);
*/

#endif
