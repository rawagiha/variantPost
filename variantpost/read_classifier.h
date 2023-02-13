#ifndef READ_CLASSIFIER_H
#define READ_CLASSIFIER_H

//#include <map>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>

#include "util.h"
#include "fasta/Fasta.h"

struct Read
{
    bool is_from_first;
    std::string read_name;
    bool is_reverse;
    std::string cigar_string;
    std::vector<std::pair<char, int>> cigar_vector;
    int aln_start;
    int aln_end;
    int read_start;
    int read_end;
    int covering_start = 0;
    int covering_end = 0;
    int start_offset;
    int end_offset;
    std::string seq;
    std::string ref_seq;
    std::string base_qualities = "";
    std::string non_ref_quals = "";
    int mapq;
    std::vector<Variant> variants;
    std::vector<std::pair<int, int>> aligned_segments;
    std::vector<std::pair<int, int>> skipped_segments;
    bool may_be_complex = false;
    bool has_non_target_in_critical_region = false;
    int target_aligned_pos = -1;
    std::string target_aligned_ref = "N";
    std::string target_aligned_alt = "N";
    char covering_ptrn;
    char local_ptrn;
    char clip_ptrn;
    double dirty_base_rate = 0.0;
    std::string non_ref_ptrn_str;
    int kmer_score = -1;

    Read();

    Read(const int unspliced_local_reference_start, 
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
         const std::vector<int> & qualities, 
         const int mapq,
         const Variant & target, // target
         const int ref_allele_len,
         const int lpos,
         const int pos,
         const int rpos,
         const bool is_shiftable,
         const std::unordered_map<int, char> & indexed_local_reference, //used for indel eq. check
         const std::string & chrom,
         FastaReference & fr
    );
};

void classify_reads(std::vector<Read> & targets,
                    std::vector<Read> & candidates,
                    std::vector<Read> & non_targets,
                    FastaReference & fr,
                    //Variant & target,
                    //const std::string & chrom,
                    //const int,
                    //const std::string &,
                    //const std::string &,
                    Variant & target,
                    
                    const int mapping_quality_threhold,
                    const int base_quality_threshold, 
                    const int unspl_loc_ref_start,
                    const int unspl_loc_ref_end,
                    const std::string & unspl_loc_ref,
                    const std::unordered_map<int, char> & ref_dict,
                    const std::vector<bool> & is_from_first_bam,
                    const std::vector<std::string> & read_names,
                    const std::vector<bool> & are_reverse, 
                    const std::vector<std::string> & cigar_strings,
                    const std::vector<int> & aln_starts, 
                    const std::vector<int> & aln_ends,
                    const std::vector<std::string> & read_seqs, 
                    const std::vector<std::vector<int>> & base_quals, 
                    const std::vector<int> & mapqs
);

#endif
