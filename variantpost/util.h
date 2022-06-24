#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <string>
#include <vector>
#include <utility>

std::string to_fastq_qual( const std::vector<int> & arr_qual );

std::vector<std::pair<char, int>> to_cigar_vector( const std::string &
                               cigar_string );

std::string get_read_wise_ref_seq( int aln_start, int aln_end,
                                   int unspliced_local_reference_start,
                                   const std::string & unspliced_local_reference );

std::map<int, char> reference_by_position( const std::string &
        unspliced_local_reference, int unspliced_local_reference_start,
        int unspliced_local_reference_end );

struct Variant {
    std::string chrom_;
    int pos_;
    std::string ref_;
    std::string alt_;
    
    int ref_len_;
    int alt_len_;
    int variant_end_pos_ = pos_ + ref_len_;
    
    bool is_substitute_;
    bool is_ins_;
    bool is_del_;
    
    Variant(const std::string & chrom, const int pos, const std::string & ref, const std::string & alt);
    
    //int ref_len_;
    //int alt_len_;
    //int variant_end_pos_ = pos_ + ref_len_;

    //bool is_substitute_;
    //bool is_ins_;
    //bool is_del_;

    /*
    Variant( const std::string & chrom,
             const int pos,
             const std::string & ref,
             const std::string & alt,
             const int unspliced_local_reference_start,
             const int unspliced_local_reference_end,
             const std::map<int, char> & indexed_local_reference
           );


    */
    bool is_shiftable(const std::map<int, char> & indexed_local_reference) const;
    int get_leftmost_pos(const int unspliced_local_reference_start, const std::map<int, char> & indexed_local_reference) const;
    int get_rightmost_pos(const int unspliced_local_reference_end, const std::map<int, char> & indexed_local_reference) const;
    bool is_equivalent(const Variant & v, const int unspliced_local_reference_start, const std::map<int, char> & indexed_local_reference) const;
    //void say_hi(const Variant & j) const;

    //bool operator == ( const Variant & rhs ) const;
    
};


std::vector<Variant> find_mapped_variants( const int aln_start,
        const int aln_end, const std::string & ref_seq, const std::string & read_seq,
        const std::vector<std::pair<char, int>> & cigar_vector,
        const std::string & chrom,
        const int & unspliced_local_reference_start,
        const int & unspliced_local_reference_end,
        const std::map<int, char> & indexed_local_reference );


#endif