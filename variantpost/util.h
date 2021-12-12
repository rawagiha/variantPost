#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>
//#include "util.h"
#include <regex>

std::string to_fastq_qual(const std::vector<int> & ); 

std::vector<std::string> to_cigar_vector( const std::regex &, const std::string & );

std::string get_read_wise_ref_seq(int, int, int, const std::string &);

#endif
