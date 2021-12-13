#ifndef PILEUP_PARSER_H
#define PILEUP_PARSER_H

#include <string>
#include <vector>

namespace pileup {

struct ParsedRead {
    std::string read_name_;
    bool is_reverse_;
    std::string cigar_string_;
    std::vector<std::string> cigar_vector_;
    int aln_start_;
    int aln_end_;
    std::string read_seq_;
    std::string ref_seq_;
    std::string base_qualities_;
    int mapq_;
    bool is_spliced_;

    ParsedRead(int, 
               const std::string &, 
               const std::string &,
               bool, 
               const std::string &,
               int, 
               int, 
               const std::string &,
               const std::string &, 
               const std::vector<int> &, 
               int 
    );

};

void parse_pileup( int, 
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
