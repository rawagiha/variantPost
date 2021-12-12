#include <string>
#include <vector>
#include <regex>
#include <iostream>
#include <algorithm>

#include "util.h"

inline char to_base_qual(int q) {
    return static_cast<char>(q + 33);
}

std::string to_fastq_qual(const std::vector<int> & qvec) 
{
    std::vector<char> res(qvec.size());
    std::transform(qvec.begin(), qvec.end(), res.begin(), to_base_qual);
    std::string fq(res.begin(), res.end());
    return fq;
}


std::vector<std::string> to_cigar_vector( const std::regex & cigar_pattern, const std::string & cigar_string )
{
    std::vector<std::string> cigarette;
    std::sregex_token_iterator itr(cigar_string.begin(), cigar_string.end(), cigar_pattern, -1);
    std::sregex_token_iterator end_itr;
    
    while (itr != end_itr) {
        cigarette.push_back(*itr);
        ++itr;
    }
    
    return cigarette;
}

/*
std::vector<std::string> to_cigar_vector( const std::regex & cigar_pattern, const std::string & cigar_string )
{
    std::regex_iterator<std::string::const_iterator> cigar_itr(
        cigar_string.begin(),
        cigar_string.end(), cigar_pattern );
    std::regex_iterator<std::string::const_iterator> end_of_itr;
    
    std::vector<std::string> cigarette;
    while ( cigar_itr != end_of_itr ) {
        cigarette.push_back( ( *cigar_itr ).str() );
        ++cigar_itr;
    }
    
    return cigarette;
}
*/

int get_clipped_len( const std::vector<std::string> & cigar_vector,
                     bool for_start )
{
    int clip_len = 0;

    if ( for_start ) {
        std::string first_op = cigar_vector.front();
        if ( first_op.find( 'S' ) != std::string::npos ) {
            clip_len = std::stoi( first_op.substr( 0, first_op.size() - 1 ) );
        }
    }
    else {
        std::string last_op = cigar_vector.back();
        if ( last_op.find( 'S' ) != std::string::npos ) {
            clip_len = std::stoi( last_op.substr( 0, last_op.size() - 1 ) );
        }
    }

    return clip_len;
}


std::string get_read_wise_ref_seq( int aln_start, int aln_end,
                                   int unspliced_local_reference_start,
                                   const std::string & unspliced_local_reference )
{
    int start_idx = aln_start - 1 - unspliced_local_reference_start;
    return unspliced_local_reference.substr( start_idx,
            ( aln_end - aln_start + 1 ) );
}
