#ifndef UTIL_H
#define UTIL_H

#include <set>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "fasta/Fasta.h"

template<typename A, typename B, typename C>
bool is_ascending(const A & a, const B & b, const C & c)
{
    return (a <= b) && (b <= c);
}

//transfer all elems to destination vector
//source vector will be cleared
//---------------------------------------------------------------
template<typename A>
void transfer_vector(std::vector<A> & dest, std::vector<A> & src)
{
    dest.insert(
        dest.end(),
        std::make_move_iterator(src.begin()),
        std::make_move_iterator(src.end())
       );

    src.clear();
}

//transfer i-th elem to destination vector
//source vector will NOT be resized
//------------------------------------------------------------
template<typename A>
void transfer_elem(std::vector<A> & dest, std::vector<A> & src, const size_t i)
{
    dest.insert(
        dest.end(),
        std::make_move_iterator(src.begin() + i),
        std::make_move_iterator(src.begin() + i + 1)
       );
}


std::string find_commonest_str(const std::vector<std::string> & arr_str);

std::string to_fastq_qual( const std::vector<int> & arr_qual );

//cigar string to vec: 10M4D3M2S -> {<'M', 10>, <'D', 4>, <'M', 3>, <'S', 2>}
//----------------------------------------------------------------------------
std::vector<std::pair<char, int>> to_cigar_vector(const std::string & cigar_string);


//cigar vec to string: {<'M', 10>, <'D', 4>, <'M', 3>, <'S', 2>} -> 10M4D3M2S
//----------------------------------------------------------------------------- 
std::string to_cigar_string(const std::vector<std::pair<char, int>> & cigar_vector);


// include skips (N) to cigar vec
// {<'=', 5>, <'I', 1>, <'=', 5> + {{100, 104}, {108, 109}, {112, 114}}
// -> {<'=', 5>, <'I', 1>, <'N', 3>, <'=', 2>, <'N', 2>, <'=', 3>}
//--------------------------------------------------------------------------
void splice_cigar(std::vector<std::pair<char, int>> & cigar_vector,
                  const int start_offset, 
                  const std::vector<int> & genomic_positions, 
                  const std::vector<std::pair<int, int>> & extended_coordinates);


//make insertion first for complex case with gap-merging
//{'=', 2}, {'D', 2}, {'I', 2}, {'D', 4}, {'=', 3}} 
// -> {{'=', 2}, {'I', 2}, {'D', 6},{'=', 3}}
//-----------------------------------------------------------------------
void move_up_insertion(std::vector<std::pair<char, int>> & cigar_vector);

std::string get_unspliced_ref_seq(const int aln_start, const int aln_end,
                                  const int unspliced_local_reference_start,
                                  const std::string & unspliced_local_reference);

//fasta file kept open???
std::string get_spliced_ref_seq(const std::string & chrom, const int aln_start,
                                const std::vector<std::pair<char, int>> & cigar_vector,
                                FastaReference & fr);

void parse_splice_pattern(std::vector<std::pair<int, int>> & exons,
                          std::vector<std::pair<int, int>> & introns,
                          const std::vector<std::pair<char, int>> & cigar_vector,
                          const int start,
                          const int end);

//expand segment start/end: {{123, 125}, {502, 504}} -> {123, 124, 125, 502, 503, 504}
//-----------------------------------------------------------------------------------------------
std::vector<int> expand_coordinates(const std::vector<std::pair<int, int>> & coordinates);

//void make_skip_after_ins(std::vector<std::pair<char, int>> & cigar_vec);

std::unordered_map<int, char> reference_by_position( const std::string &
        unspliced_local_reference, int unspliced_local_reference_start,
        int unspliced_local_reference_end );

struct RefSeq {
    std::string seq;
    int start;      //1-based first base pos
    int stop;       //1-based last base pos
};

struct Variant {
    int pos;
    std::string ref;
    std::string alt;
    std::string chrom;
    
    int ref_len;
    int alt_len;
    int variant_end_pos = pos + ref_len;
    
    bool is_substitute;
    bool is_ins;
    bool is_del;
    bool is_complex;
    
    
    Variant(const int pos, const std::string & ref, const std::string & alt, const std::string & chrom = "N");
    
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

struct RealignedGenomicSegment
{
    int start = 0; //1-based
    int end = 0;
    int target_pos = 0;
    bool has_target = false;
    std::string ref_seq = "";
    std::string seq = "";
    std::string base_qualities = "";
    std::vector<std::pair<char, int>> cigar_vec = {};
    std::vector<Variant> variants = {};
    int seq_len = 0;
    std::pair<char, int> first_cigar;
    std::pair<char, int> last_cigar;

    RealignedGenomicSegment() {}   
    
    RealignedGenomicSegment(
        const int start,
        const int end,
        int target_pos, 
        bool has_target,
        const std::string & ref_seq,
        const std::string & seq, 
        const std::string & base_qualities,
        const std::vector<std::pair<char, int>> & cigar_vec,
        const std::vector<Variant> & variants 
    ) : start(start), end(end), target_pos(target_pos), 
        has_target(has_target),
        ref_seq(ref_seq), seq(seq), base_qualities(base_qualities),
        cigar_vec(cigar_vec), variants(variants) 
        {
            seq_len = seq.size();
            first_cigar = cigar_vec[0];
            last_cigar = cigar_vec[cigar_vec.size() - 1];   
        }
};

struct PairwiseBaseAlignmnent
{
    int genomic_pos = 0;
    std::string ref_base = "";
    std::string alt_base = "";
    std::string base_qual = "";

    PairwiseBaseAlignmnent() {}

    PairwiseBaseAlignmnent(
        const int genomic_pos,
        const std::string & ref_base,
        const std::string & alt_base,
        const std::string &  base_qual
    ) : genomic_pos(genomic_pos), 
        ref_base(ref_base), alt_base(alt_base),
        base_qual(base_qual) {}
};

struct Contig
{
    std::vector<PairwiseBaseAlignmnent> alignment;
    std::vector<std::pair<int, int>> skips;
    
    Contig(
        const std::vector<PairwiseBaseAlignmnent> & alignment,
        const std::vector<std::pair<int, int>> & skips
    ) : alignment(alignment), skips(skips) {}
};

void make_contig(const std::vector<RealignedGenomicSegment> & realns);

std::vector<Variant> find_mapped_variants(const int aln_start, const int aln_end, 
                                          const std::string & ref_seq, 
                                          const std::string & read_seq,
                                          const std::string & base_qualities,
                                          const std::vector<std::pair<char, int>> & cigar_vector,
                                          std::string & non_ref_quals);


std::vector<Variant> find_variants(const int aln_start, 
                                   const std::string & ref_seq, 
                                   const std::string & read_seq,
                                   const std::string & base_qualities,
                                   const std::vector<std::pair<char, int>> & cigar_vector,
                                   std::string & non_ref_quals);

int count_repeats(const std::string & ptrn, const std::string & seq);



std::set<std::string> diff_kmers(const std::string & query,
                                 const std::string & subject,
                                 const size_t k);

int count_kmer_overlap(const string & seq, const std::set<std::string> & kmer_set);

/*
std::set<std::string> make_kmers(const std::string & seq, const size_t k);

std::unordered_map<std::string, int> generate_kmer(const std::string & seq, 
                                                      const size_t k,
                                                      std::unordered_set<std::string> & kmers);
*/

double euclidean_dist(const std::string & query,
                      const size_t k,
                      const std::unordered_map<std::string, int> & subject_kmer_cnt,
                      const std::unordered_set<std::string> & subject_kmers);


#endif
