#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>

std::string to_fastq_qual( const std::vector<int> & );

std::vector<std::string> to_cigar_vector( const std::string & );

std::string get_read_wise_ref_seq( int, int, int, const std::string & );


struct Variant {
    std::string chrom_;
    int pos_;
    std::string ref_;
    std::string alt_;

    int unspliced_local_reference_start_;
    int unspliced_local_reference_end_;
    std::map<int, char> indexed_local_reference_;

    const int ref_len_;
    const int alt_len_;
    const int variant_end_pos_ = pos_ + ref_len_

                                 const bool is_substitute_;
    const bool is_ins_;
    const bool is_del_;

    Variant( const std::string & chrom,
             const int pos,
             const std::string & ref,
             const std::string & alt,
             const int unspliced_local_reference_start,
             const int unspliced_local_reference_end,
             const std::map<int, char> & indexed_local_reference;
           );


    bool is_shiftable();
    // get_leftmost_pos
    // get_rightmost_pos
    //bool operator == (const Variant & rhs) const;

};

#endif
