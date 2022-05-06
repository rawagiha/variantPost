#ifndef PILEUP_PARSER_H
#define PILEUP_PARSER_H

#include <map>
#include <string>
#include <vector>
#include <utility>
#include "util.h"

namespace pileup {

struct ParsedRead {
    std::string read_name_;
    bool is_reverse_;
    std::string cigar_string_;
    std::vector<std::pair<char, int>> cigar_vector_;
    int aln_start_;
    int aln_end_;
    std::string read_seq_;
    std::string ref_seq_;
    std::string base_qualities_;
    int mapq_;
    bool is_spliced_;
    std::vector<Variant> variants;

    ParsedRead(int, 
               int,
               const std::string &, 
               const std::string &,
               bool, 
               const std::string &,
               int, 
               int, 
               const std::string &,
               const std::string &, 
               const std::vector<int> &, 
               int,
               const std::map<int, char> & //std::map<int, char> & indexed_local_reference
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
                   const std::vector<int> & 
);

} // end of pileup

#endif
