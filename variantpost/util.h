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
#include <unordered_set>

#include "fasta/Fasta.h"


typedef std::set<std::string_view> Kmers;
typedef std::vector<std::pair<int, int> > Coord;


struct UserParams{ 
    int mapq_thresh;
    char base_q_thresh;
    double lq_rate_thresh;
    int match_score;
    int mismatch_penal;
    int gap_open_penal;
    int gap_ext_penal;
    int kmer_size; 
    int local_thresh;

    UserParams(
        const int mapq_thresh,
        const int base_q_thresh,
        const double lq_rate_thresh,
        const int match_score,
        const int mismatch_penal,
        const int gap_open_penal,
        const int gap_ext_penal,
        const int kmer_size,
        const int local_thresh
    );
}; 


struct LocalReference{ 
    private:
        std::string _seq; 
    
    public:
        FastaReference fasta;
        std::string chrom;
        int start;
        int end;
        std::string_view seq;
        std::unordered_map<int, std::string_view> dict;

    LocalReference(
        const string& fastafile,
        const std::string& chrom,
        const int start,
        const int end
    );     
};


struct Variant
{
    int pos;
    std::string ref;
    std::string alt;
    
    bool has_n;
    bool is_clipped_segment;
    
    int ref_len;
    int alt_len;
    int variant_end_pos = pos + ref_len;
    int lpos = -1;
    int rpos = -1;
    
    bool is_substitute;
    bool is_ins;
    bool is_del;
    bool is_complex;
    bool is_shiftable;
    
    Variant(
        const int pos, 
        const std::string& ref, 
        const std::string& alt, 
        bool is_clipped_segment = false
    );
    
    void set_leftmost_pos(const LocalReference& loc_ref);
    
    void set_rightmost_pos(const LocalReference& loc_ref);

    bool is_equivalent(const Variant& v, const LocalReference& loc_ref) const;
    
    std::string minimal_repeat_unit() const;

    bool operator == (const Variant& rhs) const
    {
        return (pos == rhs.pos && ref == rhs.ref && alt == rhs.alt);
    }
};


std::vector<std::pair<char, int>> to_cigar_vector(std::string_view cigar_string);

void move_up_insertion(std::vector<std::pair<char, int>>& cigar_vector);

void splice_cigar(
    std::vector<std::pair<char, int>>& cigar_vector,
    const int start_offset,
    const std::vector<int>& genomic_positions,
    const Coord& coordinates
);

void parse_variants(
    const int aln_start, 
    std::string_view ref_seq, 
    std::string_view read_seq,
    std::string_view base_qualities, 
    const std::vector<std::pair<char, int>>& cigar_vector,
    const std::unordered_map<int, std::string_view>& ref_dict,
    std::vector<Variant>& variants, 
    std::string& non_ref_quals
);


template<typename A>
void transfer_vector(std::vector<A>& dest, std::vector<A>& src)
{
    dest.insert(
        dest.end(),
        std::make_move_iterator(src.begin()),
        std::make_move_iterator(src.end())
    );

    src.clear();
    src.shrink_to_fit();
}


template<typename A>
void transfer_elem(std::vector<A>& dest, std::vector<A>& src, const size_t i)
{
    //src will NOT be resized
    dest.insert(
        dest.end(),
        std::make_move_iterator(src.begin() + i),
        std::make_move_iterator(src.begin() + i + 1)
    );
}

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



/*
//typedef std::vector<std::pair<int, int> > Coord;
//typedef std::set<std::string> Kmers;

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
    src.shrink_to_fit();
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

std::string get_unspliced_ref_seq(const int loc_ref_start, const std::string& loc_ref_seq,
                                  const int aln_start, const int aln_end);

std::string get_spliced_ref_seq(const std::string& chrom, FastaReference& fr, 
                                const int aln_start,
                                const std::vector<std::pair<char, int>>& cigar_vector);

//obsolete 
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






bool has_gaps(std::string& cigar_string);


//new
struct UserParams{ 
    int mapq_thresh;
    char base_q_thresh;
    double lq_rate_thresh;
    int match_score;
    int mismatch_penal;
    int gap_open_penal;
    int gap_ext_penal;
    int kmer_size; 
    int local_thresh;
}; 

//new
struct LocalReference{ 
    FastaReference fasta;
    std::string chrom;
    int start;
    int end;
    std::string seq;
    std::unordered_map<int, char> dict;

    void set_up(const std::string& fastafile); 
};






struct RefSeq {
    std::string seq;
    int start;      //1-based first base pos
    int stop;       //1-based last base pos
};

struct Variant{
    int pos;
    std::string ref;
    std::string alt;
    //std::string chrom;
    
    bool is_clipped_segment;
    bool has_n;
    
    int ref_len;
    int alt_len;
    int variant_end_pos = pos + ref_len;
    int lpos = -1;
    int rpos = -1;
    
    bool is_substitute;
    bool is_ins;
    bool is_del;
    bool is_complex;
    bool is_shiftable;
    
    Variant(const int pos, const std::string& ref, const std::string& alt, bool is_clipped_segment = false);
    
    void left_aln(const int unspliced_local_reference_start, const std::unordered_map<int, char> & indexed_local_reference); 
    //bool is_shiftable(const std::unordered_map<int, char> & indexed_local_reference) const;
    int get_leftmost_pos(const int unspliced_local_reference_start, const std::unordered_map<int, char> & indexed_local_reference) const;
    
    void set_leftmost_pos(const LocalReference& loc_ref);
    
    
    int get_rightmost_pos(const int unspliced_local_reference_end, const std::unordered_map<int, char> & indexed_local_reference) const;
    
    void set_rightmost_pos(const LocalReference& loc_ref);

    bool is_equivalent(const Variant& v, const LocalReference& loc_ref) const;
    
    std::string minimal_repeat_unit() const;

    bool operator == (const Variant& rhs) const
    {
        return (pos == rhs.pos && ref == rhs.ref && alt == rhs.alt);
    }
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

/\*
struct Contig
{
    std::vector<PairwiseBaseAlignmnent> alignment = {};
    std::vector<int> skip_starts = {};
    std::vector<int> skip_ends = {};

    Contig() {}

    Contig(
        const std::vector<PairwiseBaseAlignmnent> & alignment,
        const std::vector<int> & skip_starts,
        const std::vector<int> & skip_ends
    ) : alignment(alignment), skip_starts(skip_starts), skip_ends(skip_ends) {}
};
*/

//void make_contig(const std::vector<RealignedGenomicSegment> & realns, Contig & contig);

/*
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
                                   std::string & non_ref_quals,
                                   bool include_clip_as_variant=false);

//new 
void parse_variants(const int aln_start, const std::string& ref_seq, const std::string& read_seq,
                    const std::string& base_qualities, const std::vector<std::pair<char, int>>& cigar_vector,
                    const std::unordered_map<int, char>& ref_dict,
                    std::vector<Variant>& variants, std::string& non_ref_quals, bool include_clips=false);



int count_repeats(const std::string & ptrn, const std::string & seq);



std::set<std::string> diff_kmers(const std::string & query,
                                 const std::string & subject,
                                 const size_t k);

int count_kmer_overlap(const string & seq, const std::set<std::string> & kmer_set);

/\*
std::set<std::string> make_kmers(const std::string & seq, const size_t k);

std::unordered_map<std::string, int> generate_kmer(const std::string & seq, 
                                                      const size_t k,
                                                      std::unordered_set<std::string> & kmers);
*/

/*
double euclidean_dist(const std::string & query,
                      const size_t k,
                      const std::unordered_map<std::string, int> & subject_kmer_cnt,
                      const std::unordered_set<std::string> & subject_kmers);

*/
#endif 
