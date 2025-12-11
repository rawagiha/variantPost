#ifndef UTIL_H
#define UTIL_H

#include <set>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iterator>
#include <limits.h>
#include <algorithm>
#include <string_view>
#include <unordered_map>

#include "fasta/Fasta.h"
#include "ssw/ssw_cpp.h"

typedef std::vector<std::string> Strs;
typedef std::vector<int> Ints;
typedef std::vector<bool> Bools;
typedef std::set<std::string_view> Kmers;
typedef std::vector<std::pair<int, int>> Coord;
typedef std::vector<std::pair<char, int>> CigarVec;
typedef std::unordered_map<int, std::string_view> Dict;
typedef StripedSmithWaterman::Alignment Alignment;
typedef StripedSmithWaterman::Aligner Aligner;
typedef StripedSmithWaterman::Filter Filter;

struct Variant;
typedef std::vector<Variant> Vars;


/*********************************************************************************/
/*                               utility structs                                 */
/*********************************************************************************/

//------------------------------------------------------------------------------
// user parameters
struct UserParams { 
    UserParams(const int base_q_thresh, const double lq_rate_thresh, 
               const int match_score, const int mismatch_penal, 
               const int gap_open_penal, const int gap_ext_penal, 
               const int kmer_size, const int dimer_window, const int local_thresh);
    
    // thresholds
    char base_q_thresh; double lq_rate_thresh; int local_thresh;
    
    // local alignment params
    int match_score; int mismatch_penal; int gap_open_penal; int gap_ext_penal;
    int max_gap_open_penal = 10; int max_mismatch_penal = 10;
    
    // kmer related 
    int kmer_size; int dimer_window;
}; 

//------------------------------------------------------------------------------
// features derived from input reference genome (not from alignment data) 
struct LocalReference { 
    LocalReference(const std::string& fastafile, 
                   const std::string& chrom, const int start, const int end);
    
    void setFlankingBoundary(const Variant& target, const size_t window);
     
    FastaReference fasta;
    std::string chrom;
    
    int start; int end;
    
    // flanking region defined by 2-mer diversity
    int flanking_start = -1; int flanking_end = -1;
    bool has_flankings = false;
    int low_cplx_len = 0; // len of low compex seq within flanking
    
    //--------------------------------------------------------------------------
    // sequences 
    std::string_view seq; // read only
    std::string _seq; // data for string_view
    
    Dict dict; // dictitionary {pos, base}    
};

//------------------------------------------------------------------------------
struct Variant {
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
    
    void countRepeats(const LocalReference& loc_ref);

    //--------------------------------------------------------------------------
    // test variant identity after normalization
    bool isEquivalent(const Variant& v, const LocalReference& loc_ref) const;
    
    //--------------------------------------------------------------------------
    // inputs
    int pos; // 1-based genomic pos
    std::string ref; std::string alt; 
    std::string indel_seq; // only defined for simple indels
    std::string_view qual; // base qualities, NOT including padding (aln may begin with insertion) 
    
    //--------------------------------------------------------------------------
    // flanking sequences
    std::string_view lt_seq; // left flanking
    std::string_view mid_seq; // middle (inserted seq)
    std::string_view rt_seq; // right flanking

    //--------------------------------------------------------------------------
    // numeric data
    int ref_len = 0, alt_len = 0, event_len = 0; // allele len, event = the longer 
    int indel_len = 0; // len of inserted or deleted sequence
    int repeats = 0; // n of indel sequence repeats 
    int lpos = -1, rpos = -1; // left and right aligned positions
    int _end_pos = pos + ref_len; // end position before right-aligned
    int end_pos = rpos + ref_len; // end postion of event after right-aligned
    
    //--------------------------------------------------------------------------
    // boolean flags
    bool has_n = false; // true if "N" base in alleles
    bool has_flankings = false; // true if flanking sequences set
    bool is_snv = false, is_mnv = false, is_ins = false, is_del = false; 
    bool is_substitute = false;
    bool is_complex = false; // true if complex event
    bool is_shiftable = false; // true if left- and right-positions are different
};

//------------------------------------------------------------------------------
// make Variant hashable
template<> struct std::hash<Variant> {
    size_t operator () (const Variant& v) const noexcept {
        size_t h1 = std::hash<int>()(v.pos);
        size_t h2 = std::hash<std::string>()(v.ref);
        size_t h3 = std::hash<std::string>()(v.alt);
        return h1 ^ h2 ^ h3;
    }
};

//------------------------------------------------------------------------------
// overloaded operators
bool operator == (const Variant& lhs, const Variant& rhs); 
bool operator != (const Variant& lhs, const Variant& rhs);
bool operator < (const Variant& lhs, const Variant& rhs);

/*********************************************************************************/
/*                               utility functions                               */
/*********************************************************************************/

//------------------------------------------------------------------------------
// fill CigarVec from CIGAR string
void fill_cigar_vector(const std::string& cigar_str, CigarVec& cigar_vector);

//------------------------------------------------------------------------------
// overloaded version
void fill_cigar_vector(const std::vector<uint32_t>& cigar, CigarVec& cigar_vector);

//------------------------------------------------------------------------------
// fill Vars and index variants in the read by parsing CIGAR
void read2variants(const int aln_start, std::string_view ref_seq, 
                   std::string_view read_seq, std::string_view base_qualities, 
                   const CigarVec& cigar_vector, const Dict& ref_dict, 
                   Vars& variants, Ints& var_idx, Coord& idx2pos);


//void sw2nonrefs(const std::string& ref, const std::string& query, 
//                Alignment& aln, NonRefs& nrs);

//------------------------------------------------------------------------------
int find_target(LocalReference& loc_ref, const Variant& target, Vars& vars);

//------------------------------------------------------------------------------
// make seq with Vars and index it
void make_sequence(LocalReference& loc_ref, const Vars& variants, const int start, 
                   const int end, std::string& seq, Coord* p_idx2pos=nullptr);

//------------------------------------------------------------------------------
// introduce variant to existing seq/idx2pos 
void mutate_sequence(const Variant& variant, std::string& seq, Coord& idx2pos);

//------------------------------------------------------------------------------
// make a set of kmers from string 
void make_kmers(std::string_view seq, const size_t k, Kmers& kmers);

//------------------------------------------------------------------------------
// kmers specific to s1/s2 -> dkm1/dkm2
void differential_kmers(std::string_view s1, std::string_view s2,
                        const size_t k, Kmers& dkm1, Kmers& dkm2);

#endif 
