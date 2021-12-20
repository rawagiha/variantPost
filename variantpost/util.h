#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <string>
#include <vector>

std::string to_fastq_qual( const std::vector<int> & );

std::vector<std::string> to_cigar_vector( const std::string & );

std::string get_read_wise_ref_seq( int, int, int, const std::string & );

std::map<int, char> pos_index_reference( const std::string &, int, int);

struct Variant {
    std::string chrom_;
    int pos_;
    std::string ref_;
    std::string alt_;

    int unspliced_local_reference_start_;
    int unspliced_local_reference_end_;
    std::map<int, char> indexed_local_reference_;

    int ref_len_;
    int alt_len_;
    int variant_end_pos_ = pos_ + ref_len_;

    bool is_substitute_;
    bool is_ins_;
    bool is_del_;

    Variant( const std::string & chrom,
             const int pos,
             const std::string & ref,
             const std::string & alt,
             const int unspliced_local_reference_start,
             const int unspliced_local_reference_end,
             const std::map<int, char> & indexed_local_reference
           );


    bool is_shiftable() const;
    int get_leftmost_pos() const;
    int get_rightmost_pos() const;

    bool operator == (const Variant & rhs) const;

};

#endif
