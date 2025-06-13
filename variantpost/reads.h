#ifndef READ_H
#define READ_H

#include <climits>
#include "util.h"

struct Read
{   
    //-------------------------------------------------------------
    // constructor from inputs passed by Cython
    Read(std::string_view name, const bool is_reverse, std::string_view cigar_str,
        const int aln_start, const int aln_end, std::string_view seq,
        const std::vector<int>& quals, const int mapq, const bool is_from_first_bam);
    
    //-------------------------------------------------------------
    // basic attributes 
    std::string_view name;  // read identifier 
    bool is_reverse;        // true if read is aligned to the complementary strand  
    bool is_from_first_bam; // true if read is from the first BAM file 

    //-------------------------------------------------------------
    // quality scores 
    int mapq;  // mapping quality 
    std::string base_quals; // base quality string 
    std::string non_ref_quals; // base quality for non-ref bases including softclippings

    //-------------------------------------------------------------
    // genome coordinates (1-based) 
    int aln_start; int aln_end; // alignment start/end 
    int start_offset; int end_offset; // softclip len 
    int read_start; int read_end; // alignmet start/end extented by softclip len
    Coord aligned_segments = {}; // vector<pair<int, int>> start/end of aligned segments 
    Coord skipped_segments= {}; // vector<pair<int, int>> start/end of skipped segments
    int covering_start; int covering_end; // start/end of aligned segment covering the target locus
                                          // may differ from aln_start/aln_end if spliced

    //-------------------------------------------------------------
    // sequences 
    std::string_view seq; // reference to read seq data (std::string) from input
    std::string_view ref_seq; // reference to refseq data 
    std::string spliced_ref_seq; // store refseq data if spliced.  
                                 // if unspliced, stored in LocalReference in util.h
                                    
    //-------------------------------------------------------------
    // cigar 
    std::string_view cigar_str; // reference to cigar string from input
    CigarVec cigar_vector; // CigarVec vector<pair<char, int>>
    

    //non reference event info
    std::vector<Variant> variants;
    int variants_target_idx = -1;
    /*dist to closest non-target overlapping target shiftable segment*/
    int dist_to_non_target = INT_MAX; 
    int dist_to_clip = INT_MAX;
    int target_pos = -1;  //actual pos (may not be normalized)
    std::string target_ref = "N"; //actual alt
    std::string target_alt = "N"; //actual ref
    std::string non_ref_signature;
    std::string splice_signature;
    bool has_target = false;
    bool incomplete_shift= false;
    bool may_be_complex = false; 
    int local_uniqueness = 0;
    int lt_end_matches = -1;
    int rt_end_matches = -1;
    
    //annotated patterns    
    char covering_ptrn = '\0';
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
    bool is_ref = false;
    bool is_na_ref = false;
    bool is_tight_covering = false;
    bool is_deprioritized = false;
    //bool is_contig_member = false;
    
};


typedef std::vector<Read> Reads;

void sort_by_start(Reads & reads);

void sort_by_kmer(Reads & reads);

void annot_ref_seq(Read& read, LocalReference& loc_ref);

void annot_splice_pattern(Read& read);

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


void classify_reads(
    Reads& reads, 
    Reads& targets, 
    Reads& candidates, 
    Reads& non_targets, 
    const UserParams& user_params
);


#endif
