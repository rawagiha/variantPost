#ifndef UTIL_H
#define UTIL_H

#include <set>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "fasta/Fasta.h"

template<typename A, typename B, typename C>
bool is_ascending(const A & a, const B & b, const C & c)
{
    return (a <= b) && (b <= c);
}

std::string find_commonest_str(const std::vector<std::string> & arr_str);

std::string to_fastq_qual( const std::vector<int> & arr_qual );

//cigar string to vec: 10M4D3M2S -> {<'M', 10>, <'D', 4>, <'M', 3>, <'S', 2>}
//----------------------------------------------------------------------------
std::vector<std::pair<char, int>> to_cigar_vector(const std::string & cigar_string);

//cigar vec to string: {<'M', 10>, <'D', 4>, <'M', 3>, <'S', 2>} -> 10M4D3M2S
//----------------------------------------------------------------------------- 
std::string to_cigar_string(const std::vector<std::pair<char, int>> & cigar_vector);

std::string get_unspliced_ref_seq(const int aln_start, const int aln_end,
                                  const int unspliced_local_reference_start,
                                  const std::string & unspliced_local_reference);

std::string get_spliced_ref_seq(const std::string & chrom, const int aln_start,
                                const std::vector<std::pair<char, int>> & cigar_vector,
                                FastaReference & fr);

void parse_splice_pattern(std::vector<std::pair<int, int>> & exons,
                          std::vector<std::pair<int, int>> & introns,
                          const std::vector<std::pair<char, int>> & cigar_vector,
                          const int start,
                          const int end);

//expand segment start/end: {{123, 125}, {502, 504}} -> {0(offset), 123, 124, 125, 502, 503, 504}
//-----------------------------------------------------------------------------------------------
std::vector<int> expand_coordinates(const std::vector<std::pair<int, int>> & coordinates, 
                                    bool with_offset = true);

std::unordered_map<int, char> reference_by_position( const std::string &
        unspliced_local_reference, int unspliced_local_reference_start,
        int unspliced_local_reference_end );

struct RefSeq {
    std::string seq;
    int start;      //1-based first base pos
    int stop;       //1-based last base pos
};

struct Variant {
    //std::string chrom;
    int pos;
    std::string ref;
    std::string alt;
    
    int ref_len;
    int alt_len;
    int variant_end_pos = pos + ref_len;
    
    bool is_substitute;
    bool is_ins;
    bool is_del;
    
    Variant(const int pos, const std::string & ref, const std::string & alt);
    
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
    void left_aln(const int unspliced_local_reference_start, const std::unordered_map<int, char> & indexed_local_reference); 
    bool is_shiftable(const std::unordered_map<int, char> & indexed_local_reference) const;
    int get_leftmost_pos(const int unspliced_local_reference_start, const std::unordered_map<int, char> & indexed_local_reference) const;
    int get_rightmost_pos(const int unspliced_local_reference_end, const std::unordered_map<int, char> & indexed_local_reference) const;
    bool is_equivalent(const Variant & v, const int unspliced_local_reference_start, const std::unordered_map<int, char> & indexed_local_reference) const;
    std::string minimal_repeat_unit() const;
    //void say_hi(const Variant & j) const;

    //bool operator == ( const Variant & rhs ) const;
    
};


std::vector<Variant> find_mapped_variants(const int aln_start, const int aln_end, 
                                          const std::string & ref_seq, 
                                          const std::string & read_seq,
                                          const std::string & base_qualities,
                                          const std::vector<std::pair<char, int>> & cigar_vector,
                                          std::string & non_ref_quals);


int count_repeats(const std::string & ptrn, const std::string & seq);

std::set<std::string> make_kmers(const std::string & seq, const size_t k);

std::unordered_map<std::string, int> generate_kmer(const std::string & seq, 
                                                      const size_t k,
                                                      std::unordered_set<std::string> & kmers);

double euclidean_dist(const std::string & query,
                      const size_t k,
                      const std::unordered_map<std::string, int> & subject_kmer_cnt,
                      const std::unordered_set<std::string> & subject_kmers);


#endif
