#ifndef PILEUP_PARSER_H
#define PILEUP_PARSER_H

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>

#include "util.h"

namespace pileup {

struct ParsedRead {
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
    //bool is_spliced;
    std::vector<Variant> variants;
   // bool is_clipped;
    //bool is_ref_seq;
    //bool is_target;
    bool may_be_complex;
    char covering_ptrn;
    char local_ptrn;
    std::string non_ref_ptrn_str;

    ParsedRead(const int  unspliced_local_reference_start, 
               const int  unspliced_local_reference_end,
               const std::string & unspliced_local_reference, 
               const std::string & read_name,
               const bool  is_reverse, 
               const std::string & cigar_string,
               const int  aln_start, 
               const int  aln_end, 
               const std::string & read_seq,
               const std::string & ref_seq, 
               const std::vector<int> & qualities, 
               const int  mapq,
               const Variant & target, // target
               const int ref_allele_len,
               const int  lpos,
               const int  pos,
               const int  rpos,
               const bool  is_shiftable,
               const std::unordered_map<int, char> &  indexed_local_reference //std::map<int, char> & indexed_local_reference

    );

};

void parse_pileup( const std::string &,
                   int,
                   const std::string &,
                   const std::string &,
                   int, 
                   int,
                   const std::string &, 
                   const std::vector<std::string> &,
                   const std::vector<bool> &, 
                   const std::vector<std::string> &,
                   const std::vector<int> &, 
                   const std::vector<int> &,
                   const std::vector<std::string> &, 
                   const std::vector<std::string> &,
                   const std::vector<std::vector<int>> &, 
                   const std::vector<int> &,
                   const std::vector<bool> &
);

} // end of namespace "pileup"

#endif
