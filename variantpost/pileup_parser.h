#ifndef PILEUP_PARSER_H
#define PILEUP_PARSER_H

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>

#include "util.h"
#include "fasta/Fasta.h"

struct ParsedRead {
    bool is_from_first;
    std::string read_name;
    bool is_reverse;
    std::string cigar_string;
    std::vector<std::pair<char, int>> cigar_vector;
    int aln_start;
    int aln_end;
    int read_start;
    int read_end;
    int covering_start;
    int covering_end;
    int start_offset;
    int end_offset;
    std::string read_seq;
    std::string ref_seq;
    std::string base_qualities;
    int mapq;
    std::vector<Variant> variants;
    std::vector<std::pair<int, int>> un_spliced_segments;
    std::vector<std::pair<int, int>> spliced_segments;
    bool may_be_complex;
    int target_aligned_pos;
    std::string target_aligned_ref;
    std::string target_aligned_alt;
    char covering_ptrn;
    char local_ptrn;
    char clip_ptrn;
    double dirty_base_rate;
    std::string non_ref_ptrn_str;

    ParsedRead();

    ParsedRead(const int  unspliced_local_reference_start, 
               const int  unspliced_local_reference_end,
               const std::string & unspliced_local_reference, 
               const char base_qual_thresh,
               const bool is_from_first, 
               const std::string & read_name,
               const bool is_reverse, 
               const std::string & cigar_string,
               const int  aln_start, 
               const int  aln_end, 
               const std::string & read_seq,
               const std::vector<int> & qualities, 
               const int  mapq,
               const Variant & target, // target
               const int ref_allele_len,
               const int  lpos,
               const int  pos,
               const int  rpos,
               const bool  is_shiftable,
               const std::unordered_map<int, char> &  indexed_local_reference, //std::map<int, char> & indexed_local_reference
               const std::string & chrom,
               FastaReference & fr
    );

};

void parse_pileup( std::vector<ParsedRead> &,
                   std::vector<ParsedRead> &,
                   std::vector<ParsedRead> &,
                   FastaReference &,
                   const std::string &,
                   int,
                   const std::string &,
                   const std::string &,
                   int, 
                   int,
                   int,
                   //const std::string &, 
                   const std::vector<std::string> &,
                   const std::vector<bool> &, 
                   const std::vector<std::string> &,
                   const std::vector<int> &, 
                   const std::vector<int> &,
                   const std::vector<std::string> &, 
                   //const std::vector<std::string> &,
                   const std::vector<std::vector<int>> &, 
                   const std::vector<int> &,
                   const std::vector<bool> &
);

#endif
