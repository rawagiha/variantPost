#ifndef UTIL_H
#define UTIL_H

#include <set>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iterator>
#include <algorithm>
#include <string_view>
#include <unordered_map>

#include "fasta/Fasta.h"
#include "ssw/ssw_cpp.h"

typedef std::set<std::string_view> Kmers;
typedef std::vector<std::pair<int, int>> Coord;
typedef std::vector<std::pair<char, int>> CigarVec;
typedef std::unordered_map<int, std::string_view> Dict;
typedef StripedSmithWaterman::Alignment Alignment;
typedef StripedSmithWaterman::Aligner Aligner;
typedef StripedSmithWaterman::Filter Filter;

struct Variant;

struct Qual
{
    Qual(const int idx, const int pos, const int len);
    int idx; int pos; int len;
};

struct UserParams
{ 
    int mapq_thresh;
    char base_q_thresh;
    double lq_rate_thresh;
    int match_score;
    int mismatch_penal;
    int gap_open_penal;
    int gap_ext_penal;
    int kmer_size; 
    int local_thresh;
    int retarget_thresh;
    int min_dimer_cnt;

    UserParams
    (
        const int mapq_thresh,
        const int base_q_thresh,
        const double lq_rate_thresh,
        const int match_score,
        const int mismatch_penal,
        const int gap_open_penal,
        const int gap_ext_penal,
        const int kmer_size,
        const int local_thresh,
        const int retarget_thresh
    );
}; 

//------------------------------------------------------------------------------
struct LocalReference
{ 
    LocalReference(const std::string& fastafile, 
                   const std::string& chrom, const int start, const int end);
    
    void setFlankingBoundary(const Variant& target, const size_t window);
     
    FastaReference fasta;
    std::string chrom;
    
    int start;
    int end;
    
    // flanking region defined by 2-mer diversity
    int flanking_start = -1; int flanking_end = -1;
    bool has_flankings = false;

    std::string_view seq;
    Dict dict;

private:
    std::string _seq;    
};


//------------------------------------------------------------------------------
struct Variant
{
    //--------------------------------------------------------------------------
    // base qualities may not be supplied (e.g.,deletions) 
    Variant(const int pos, 
            const std::string& ref, const std::string& alt, std::string_view qual = ""); 
    
    //--------------------------------------------------------------------------
    // perform left or right alignmemt 
    void setLeftPos(const LocalReference& loc_ref);
    void setRightPos(const LocalReference& loc_ref);
    void setEndPos(const LocalReference& loc_ref);
    
    //--------------------------------------------------------------------------
    // set left-flanking, inserted seq (middle), right flanking
    void setFlankingSequences(const LocalReference& loc_ref);
    
    //--------------------------------------------------------------------------
    // test variant identity after normalization
    bool isEquivalent(const Variant& v, const LocalReference& loc_ref) const;
    
    //--------------------------------------------------------------------------
    // inputs
    int pos; // 1-based genomic pos
    std::string ref; std::string alt;
    std::string_view qual; // base qualities, NOT including padding (aln may begin with insertion) 
    
    //--------------------------------------------------------------------------
    // flanking sequences
    std::string_view lt_seq; // left flanking
    std::string_view mid_seq; // middle (inserted seq)
    std::string_view rt_seq; // right flanking

    //--------------------------------------------------------------------------
    // numeric data
    int ref_len = 0, alt_len = 0; // allele len
    int indel_len = 0; // len of inserted or deleted sequence
    int lpos = -1, rpos = -1; // left and right aligned positions
    int end_pos = rpos + ref_len; // end postion of event after right-aligned
    
    //--------------------------------------------------------------------------
    // boolean flags
    bool has_n = false; // true if "N" base in alleles
    bool has_flankings = false; // true if flanking sequences set
    bool is_substitute = false, is_ins = false, is_del = false; // variant class
    bool is_complex = false; // true if complex event
    bool is_shiftable = false; // true if left- and right-positions are different
    
    //bool is_overlapping = false;
    
    void _sb_leftmost_pos(const LocalReference& loc_ref);
    //void set_leftmost_pos(const LocalReference& loc_ref);
    
    void _sb_rightmost_pos(const LocalReference& loc_ref);
    //void set_rightmost_pos(const LocalReference& loc_ref);

    //bool is_equivalent(const Variant& v, const LocalReference& loc_ref) const;
    
    std::string minimal_repeat_unit() const;
};


template<> struct std::hash<Variant>
{
    size_t operator () (const Variant& v) const noexcept
    {
        size_t h1 = std::hash<int>()(v.pos);
        size_t h2 = std::hash<std::string>()(v.ref);
        size_t h3 = std::hash<std::string>()(v.alt);
        return h1 ^ h2 ^ h3;
    }
};


bool operator == (const Variant& lhs, const Variant& rhs); 


bool operator != (const Variant& lhs, const Variant& rhs);


bool operator < (const Variant& lhs, const Variant& rhs);


void find_shared_variants(
    std::vector<Variant>& shared,
    std::vector<Variant>& var_vec1,
    const std::vector<Variant>& var_vec2
);


CigarVec to_cigar_vector(std::string_view cigar_string);

//overloading
CigarVec to_cigar_vector(std::vector<uint32_t>& cigar);


void move_up_insertion(std::vector<std::pair<char, int>>& cigar_vector);


std::vector<int> find_mismatches_after_gap(const CigarVec& cigar_vector);


void splice_cigar(
    CigarVec& cigar_vector,
    const int start_offset,
    const std::vector<int>& genomic_positions,
    const Coord& coordinates
);


void parse_to_cplx_gaps(
    const int aln_start, 
    std::string_view ref_seq,
    std::string_view read_seq,
    const CigarVec& cigar_vector,
    const std::vector<int>& target_indexes,
    std::vector<Variant>& variants 
);


void read2variants(
    const int aln_start, 
    std::string_view ref_seq, 
    std::string_view read_seq,
    std::string_view base_qualities, 
    const CigarVec& cigar_vector,
    const Dict& ref_dict,
    std::vector<Variant>& variants, 
    Coord& idx2pos
);


std::vector<Variant> merge_to_cplx(const std::vector<Variant>& variants);


template<typename T>
void transfer_vector(std::vector<T>& dest, std::vector<T>& src)
{
    dest.insert(
        dest.end(),
        std::make_move_iterator(src.begin()),
        std::make_move_iterator(src.end())
    );

    src.clear();
    src.shrink_to_fit();
}


template<typename T>
void transfer_elem(std::vector<T>& dest, std::vector<T>& src, const size_t i)
{
    //src will NOT be resized
    dest.insert(
        dest.end(),
        std::make_move_iterator(src.begin() + i),
        std::make_move_iterator(src.begin() + i + 1)
    );
}

std::string to_tandem_rep(std::string_view seq);

std::vector<int> expand_coordinates(const Coord& coordinates);


std::string_view find_commonest_str(const std::vector<std::string_view>& v);


bool has_this(std::string_view query, std::string_view target_ptrn);


bool has_gaps(std::string_view cigar_string);


int count_repeats(std::string_view ptrn, std::string_view seq);


void diff_kmers(
    std::string_view query, 
    std::string_view subject, 
    const size_t k,
    Kmers& diff
);


int count_kmer_overlap(std::string_view seq, const Kmers& kmer_set);


int to_idx( 
    const int aln_start, const int target_pos, const CigarVec& cigar_vector
);

/*
struct MatchResult
{
    MatchResult(
        const int ss_start, const int ss_end,
        const Alignment& aln
    );
};*/
#endif 
