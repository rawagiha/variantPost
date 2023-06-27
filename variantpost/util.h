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

    bool operator < (const Variant& rhs) const
    {
        if (pos < rhs.pos) 
        {    
            return true;
        }
        else if (rhs.pos < pos) 
        {    
            return false;
        }
        else //pos == rhs.pos
        {
            if (ref.size() < rhs.ref.size()) 
            {
                return true;
            }
            else if (rhs.ref.size() < ref.size())
            {
                return false;
            }
            else //same ref deleted
            {
                return (alt < rhs.alt);
            }
        }
    }
};


void find_shared_variants(
    std::vector<Variant>& shared,
    std::vector<Variant>& var_vec1,
    const std::vector<Variant>& var_vec2
);


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


int find_split_idx(
    const int read_start,
    const int target_pos,
    const std::vector<std::pair<char, int>>& cigar_vector
);


#endif 
