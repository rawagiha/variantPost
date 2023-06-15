#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <utility>
#include <climits>
#include <unordered_map>

#include "util.h"
#include "fasta/Fasta.h"


struct Read{
    // basic info
    std::string name;
    bool is_reverse;
    bool is_from_first_bam;
   
    //qualities
    int mapq;
    std::string base_quals;
    std::string non_ref_quals; //base qual for non reference events incl. clips
    
    //coordinates
    int aln_start;
    int aln_end;
    int read_start;
    int read_end;
    int start_offset;
    int end_offset;
    int covering_start;
    int covering_end;

    //sequences
    std::string seq;
    std::string ref_seq;
    
    //CIGARs
    std::string cigar_str;
    std::vector<std::pair<char, int>> cigar_vector;
    
    //mapsskips(splice)
    Coord aligned_segments;
    Coord skipped_segments;
    
    //non reference event info
    std::vector<Variant> variants;
    int variants_target_idx = -1;
    int dist_to_non_target = INT_MAX; //distance to closest non-target 
    int dist_to_clip = INT_MAX;
    int target_pos = -1;  //actual pos (may not be normalized)
    std::string target_ref = "N"; //actual alt
    std::string target_alt = "N"; //actual ref
    std::string non_ref_signature;
    std::string splice_signature;
    bool has_target = false;
    bool incomplete_shift= false;
    bool may_be_complex = false;

    //annotated patterns    
    char covering_ptrn;
    char clip_ptrn;
    char local_ptrn;
    
    //metrics   
    double central_score = -1.0;
    double overall_lq_rate;
    double nonref_lq_rate;
    int kmer_score = -1;

    //other flags
    bool is_ref = false;
    bool is_na_ref = false;
    bool is_tight_covering = false;
    bool is_contig_member = false;
    
    Read();

    Read(
        const std::string& name, 
        const bool is_reverse, 
        const std::string& cigar_str, 
        const int aln_start, 
        const int aln_end, 
        const std::string seq, 
        const std::vector<int>& quals,
        const int mapq, 
        const bool is_from_first_bam
    );  
    
    bool operator == (const Read& rhs) const
    {
        return (name == rhs.name && is_reverse == rhs.is_reverse);
    }
};


typedef std::vector<Read> Reads;


void sort_by_start(Reads & reads);


void sort_by_kmer(Reads & reads);


void annotate_reads(
    Reads& reads, 
    const Variant& target,
    const UserParams& user_params, 
    LocalReference& loc_ref
);


void classify_reads(
     Reads& reads, 
     Reads& targets, 
     Reads& candiates, 
     Reads& non_targets, 
     const UserParams& user_params
);


#endif
