#include <set>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "util.h"
#include "fasta/Fasta.h"

// select most frequent str
// ------------------------------------------------------------
std::string find_commonest_str(const std::vector<std::string> & v)
{
    std::unordered_map<std::string, int> freq;
    for (const auto & elem : v) {
         freq[elem]++;
    }
   
    auto commonest = std::max_element(freq.begin(), freq.end(), [] (const auto & x, const auto & y) {return x.second < y.second;});

    return commonest->first;
}


// convert numeric base qual arr to FASTQ-style ASCII string
//----------------------------------------------------------------------------
inline char to_base_qual ( int q ) {
  return static_cast<char> ( q + 33 );
}

std::string to_fastq_qual ( const std::vector<int> &qvec ) {
  std::vector<char> res ( qvec.size() );
  std::transform ( qvec.begin(), qvec.end(), res.begin(), to_base_qual );
  std::string fq ( res.begin(), res.end() );
  return fq;
}


// parse cigar string: 10M4D3M2S -> {<'M', 10>, <'D', 4>, <'M', 3>, <'S', 2>}
//-----------------------------------------------------------------------------
inline std::pair<char, int> to_op_and_op_len ( const std::string &cigar ) {
  size_t last_idx = cigar.size() - 1;
  return std::make_pair ( cigar.substr ( last_idx, 1 ) [0],
                          std::stoi ( cigar.substr ( 0, last_idx ) ) );
}

std::vector<std::pair<char, int>> to_cigar_vector ( const std::string &
cigar_string ) {
  std::vector<std::pair<char, int>> cigarette;

  size_t pos = 0;
  size_t newpos;
  const size_t len = cigar_string.size();

  while ( pos < len ) {
    newpos = cigar_string.find_first_of ( "MIDNSHPX=", pos ) + 1;
    cigarette.emplace_back ( to_op_and_op_len ( cigar_string.substr ( pos,
                             newpos - pos ) ) );
    pos = newpos;
  }

  return cigarette;
}

int count_repeats(const std::string & ptrn,
                  const std::string & seq)
{
    int ptrn_len = ptrn.size();
    int seq_len = seq.size();
    
    int i = 0;
    int n = 0;
    while (ptrn_len * i < seq_len) {
        std::string sub_seq = seq.substr(ptrn_len * i, std::min(ptrn_len, seq_len - (ptrn_len * i)));
        if (ptrn == sub_seq) {
            ++n;
            ++i;
        }
        else return n;
    }

    return n;
}


void parse_splice_pattern(std::vector<std::pair<int, int>> & exons,
                          std::vector<std::pair<int, int>> & introns,
                          const std::vector<std::pair<char, int>> & cigar_vector,
                          const int start,
                          const int end)
{
    char op;
    int op_len;
    int curr_pos = start;
    for (const auto & c : cigar_vector) {
        op = c.first;
        op_len = c.second;
        
        switch (op) {
            case 'M':
            case 'S': //offset added
            case 'D':
            case 'X':
            case '=':
                curr_pos += op_len;
                break;
            case 'N':
                introns.emplace_back(curr_pos, (curr_pos + op_len - 1));
                curr_pos += op_len;
                break;
            default:
                break;
        }
    }

    curr_pos = start;
    for (const auto & i : introns) {
        exons.emplace_back(curr_pos, i.first - 1);
        curr_pos = i.second + 1;
    }
    exons.emplace_back(curr_pos, end);
}   


// fit local referenece to read alignment
// not considered del??
// do we need this??
//-----------------------------------------------------------------------------
std::string get_unspliced_ref_seq(const int aln_start, const int aln_end,
                                  const int unspliced_local_reference_start,
                                  const std::string &unspliced_local_reference) 
{
  int start_idx = aln_start - unspliced_local_reference_start;
  size_t expected_ref_len = aln_end - aln_start + 1;
   
  if (start_idx >= 0) {
    std::string fitted_ref = unspliced_local_reference.substr(start_idx, expected_ref_len);
    if (fitted_ref.size() == expected_ref_len) {
        return fitted_ref;
    }
  }

  // failed to fit
  return "";
}


// get spliced reference
// ----------------------------------------------------------------------------
std::string get_spliced_ref_seq(const std::string & chrom, const int aln_start, 
                                const std::vector<std::pair<char, int>> & cigar_vector,
                                FastaReference & fr)
{
    
    int curr_pos = aln_start - 1;
    
    char op;
    int op_len;
    std::string ref_seq;
    for (const auto & c : cigar_vector) {
        op = c.first;
        op_len = c.second;
        
        switch (op) {
            case 'M':
            case 'X':
            case 'D':
                ref_seq += fr.getSubSequence(chrom, curr_pos, op_len);
                curr_pos += op_len;
                break;
            case 'N':
                curr_pos += op_len;
                break;
            default:
                break;

        }

    }
    
    return ref_seq;     
}


// mapping genomic pos -> reference base
//-----------------------------------------------------------------------------
std::unordered_map<int, char> reference_by_position(const std::string & unspliced_local_reference, 
                                                    int unspliced_local_reference_start,
                                                    int unspliced_local_reference_end) 
{

  std::unordered_map<int, char> indexed_local_reference;

  int i = 0;
  int pos = unspliced_local_reference_start;
  while ( pos <= unspliced_local_reference_end ) {
    indexed_local_reference[pos] = unspliced_local_reference[i];
    ++i;
    ++pos;
  }

  return indexed_local_reference;
}

Variant::Variant 
(    int pos, 
     const std::string &ref, 
     const std::string &alt
) : pos (pos), ref (ref), alt (alt), 
    ref_len (ref.size()), alt_len (alt.size()),
    is_substitute ((alt_len == ref_len)),
    is_ins ((alt_len > ref_len)),
    is_del ((alt_len < ref_len))
{}

inline bool is_rotatable(const std::string & allele)
{
    return (allele[0] == allele[allele.size() - 1]);
}


inline void to_left(int & pos, 
                    int & variant_end_pos,
                    std::string & longer_allele,
                    std::string & shorter_allele,
                    const std::unordered_map<int, char> & indexed_local_reference)
{
    --pos;
    --variant_end_pos;
    char prev_base = indexed_local_reference.at(pos);
    longer_allele.pop_back();
    longer_allele.insert(0, 1, prev_base);
    shorter_allele = prev_base;
}


void left_align(int & pos,
                int & variant_end_pos, 
                std::string & ref, std::string & alt, 
                const bool is_ins,
                const int unspliced_local_reference_start,
                const std::unordered_map<int, char> & indexed_local_reference)
{
    std::string & longer_allele = (is_ins) ? alt : ref;
    std::string & shorter_allele = (is_ins) ? ref : alt;

    while (is_rotatable(longer_allele) && (unspliced_local_reference_start < pos)) {
        to_left(pos, variant_end_pos, longer_allele, shorter_allele, indexed_local_reference);
    }
}


void Variant::left_aln(const int unspliced_local_reference_start,
                       const std::unordered_map<int, char> & indexed_local_reference)
{
    left_align(pos, variant_end_pos, ref, alt, is_ins, unspliced_local_reference_start, indexed_local_reference);
}         


int Variant::get_leftmost_pos(const int unspliced_local_reference_start, 
                              const std::unordered_map<int, char> & indexed_local_reference) const
{
    if (is_substitute) return pos;   
    
    int pos_ = pos;
    int variant_end_pos_ = variant_end_pos;
    std::string ref_ = ref;
    std::string alt_ = alt;

    left_align(pos_, variant_end_pos_, ref_, alt_, is_ins, 
               unspliced_local_reference_start, indexed_local_reference);

    return pos_;
}


inline void to_right(int & variant_end_pos, 
                     std::string & longer_allele, std::string & shorter_allele,
                     const std::unordered_map<int, char> & indexed_local_reference)
{
    char next_base = indexed_local_reference.at(variant_end_pos);
    longer_allele.erase(0, 1);
    longer_allele += next_base;
    shorter_allele = next_base;
    ++variant_end_pos;
}


void right_align(int & pos, int & variant_end_pos, 
                 std::string & ref, std::string & alt, 
                 const bool is_ins,
                 const int unspliced_local_reference_end,
                 const std::unordered_map<int, char> & indexed_local_reference)
{
    std::string & longer_allele = (is_ins) ? alt : ref;
    std::string & shorter_allele = (is_ins) ? ref : alt;

    do {
        to_right(variant_end_pos, 
                 longer_allele, shorter_allele,
                 indexed_local_reference);
        ++pos;
    }
    while (is_rotatable(longer_allele) && (pos < unspliced_local_reference_end));
    
    // undo the last right shift
    to_left(pos, variant_end_pos, longer_allele, shorter_allele, indexed_local_reference);
}


int Variant::get_rightmost_pos(const int unspliced_local_reference_end, 
                               const std::unordered_map<int, char> & indexed_local_reference) const
{
    //NA: snv, mnv
    if (is_substitute) return pos;

    int pos_ = pos;
    int variant_end_pos_ = variant_end_pos;
    std::string ref_ = ref;
    std::string alt_ = alt;

    right_align(pos_, 
                variant_end_pos_, 
                ref_, alt_, is_ins,
                unspliced_local_reference_end, 
                indexed_local_reference);

    return pos_;
}


bool Variant::is_shiftable(const std::unordered_map<int, char> & indexed_local_reference) const
{

    if (is_substitute) {
        return false;
    }

    if (is_ins) {
        if (is_rotatable(alt)) {
            return true;
        }
        else {
            std::string longer = alt;
            std::string shorter = ref;
            int variant_end_pos_ = variant_end_pos;

            to_right(variant_end_pos_, longer, shorter, indexed_local_reference);

            return is_rotatable(longer);
        }
    }
    else {
        if (is_rotatable(ref)) {
            return true;
        }
        else {
            std::string longer = ref;
            std::string shorter = alt;
            int variant_end_pos_ = variant_end_pos;

            to_right(variant_end_pos_, longer, shorter, indexed_local_reference);

            return is_rotatable(longer);
        }
    }
}

bool Variant::is_equivalent(const Variant & other, 
                            const int unspliced_local_reference_start, 
                            const std::unordered_map<int, char> & indexed_local_reference) const
{
    std::vector<int> var_len_0{ref_len, alt_len};
    std::vector<int> var_len_1{other.ref_len, other.alt_len};

    if (var_len_0 == var_len_1){
        int pos_ = other.pos;
        int varient_end_pos_ = other.variant_end_pos;
        std::string ref_ = other.ref;
        std::string alt_ = other.alt;
        
        if (is_substitute) {
            return ((pos == pos_) & (ref == ref_) & (alt == alt_));
        } 
        else {
            left_align(pos_, varient_end_pos_, ref_, alt_, other.is_ins, 
                       unspliced_local_reference_start, 
                       indexed_local_reference);
            
            return ((pos == pos_) & (ref == ref_) & (alt == alt_));
        }               
    }
    else {
        return false;
    }
} 


std::string Variant::minimal_repeat_unit() const
{
    std::string seq = "";
    if (is_substitute) seq = alt;
    if (is_ins) seq = alt.substr(1, alt_len - 1);
    if (is_del) seq = ref.substr(1, ref_len - 1);
    
    int seq_size = seq.size();
    if (seq_size < 2) {
        return seq;
    }
    else {
        int mid_len = (int)seq_size / 2;

        int step = 1;
        while (step <= mid_len) {
            std::vector<std::string> tandems;
            for (int i = 0; step * i < seq_size; ++i) {
                tandems.push_back(seq.substr(step * i, std::min(step, seq_size - (step * i))));    
            }
            std::sort(tandems.begin(), tandems.end());
            tandems.erase(std::unique(tandems.begin(), tandems.end()),  tandems.end());
            if (tandems.size() == 1) return tandems[0];
            else ++step;    
        }
    }
    return seq;
}

std::vector<Variant> find_mapped_variants(const int aln_start, const int aln_end, 
                                          const std::string & ref_seq, 
                                          const std::string & read_seq,
                                          const std::string & base_qualities,
                                          const std::vector<std::pair<char, int>> & cigar_vector,
                                          std::string & non_ref_quals) 
{
    std::vector<Variant> variants;

    if (read_seq == ref_seq) return variants;
    
    char op;
    int op_len;
    int ref_idx = 0;
    int read_idx = 0;
    int pos = aln_start;

    for (std::vector<std::pair<char, int>>::const_iterator itr =
        cigar_vector.begin();
        itr != cigar_vector.end(); ++itr) {

        op = (*itr).first;
        op_len = (*itr).second;

        switch (op) {
            case 'M':
            case 'X':
            case '=':
                for (int i = 0; i < op_len; ++i) {
                    // snv
                    std::string ref = ref_seq.substr(ref_idx, 1);
                    std::string alt = read_seq.substr(read_idx, 1);
                    if (ref != alt) {
                        variants.emplace_back(pos, ref, alt); 
                        non_ref_quals += base_qualities.substr(read_idx, 1);
                    }
                    ++ref_idx;
                    ++read_idx;
                    ++pos;
                }
                break;
            case 'I':
                variants.emplace_back(pos - 1, 
                                      ref_seq.substr(ref_idx - 1, 1),
                                      ref_seq.substr(ref_idx - 1, 1) + read_seq.substr(read_idx, op_len));
                
                non_ref_quals += base_qualities.substr(read_idx, op_len);                 
                read_idx += op_len;
                break;
            case 'D':
                variants.emplace_back(pos - 1, 
                                      ref_seq.substr(ref_idx - 1, 1 + op_len), 
                                      ref_seq.substr(ref_idx - 1, 1));
                ref_idx += op_len;
                pos += op_len;
                break;
            case 'N':
                pos += op_len;
                break;
            case 'S':
                non_ref_quals += base_qualities.substr(read_idx, op_len);
                read_idx += (op_len);
                break;
            case 'H':
            case 'P':
                //
                break;
        }
    }
    return variants;
}

std::set<std::string> make_kmers(const std::string & seq, const size_t k)
{
    size_t n = seq.size();
    
    std::set<std::string> kmers;
    for (size_t i = 0; i < n - k; ++i) {
        kmers.insert(seq.substr(i, k));
    }

    return kmers;
}


std::unordered_map<std::string, int> generate_kmer(const std::string & seq, 
                                                   const size_t k,
                                                   std::unordered_set<std::string> & kmers)
{
    size_t n = seq.size();
    
    std::vector<std::string> v;
    std::unordered_map<std::string, int> kmer_cnt;
    for (size_t i = 0; i < n - k; ++i) {
        //++kmer_cnt[seq.substr(i, k)];
        v.push_back(seq.substr(i, k));
        kmers.insert(seq.substr(i, k));
    }
   
    return kmer_cnt;
}

inline int kmer_cnt_lookup(const std::string & kmer, const std::unordered_map<std::string, int> & kmer_cnt)
{
   std::unordered_map<std::string, int>::const_iterator it = kmer_cnt.find(kmer);
   if (it == kmer_cnt.end()) {
       return 0;
   }
   else {
       return (*it).second;
   }
}


double euclidean_dist(const std::string & query,
                      const size_t k,
                      const std::unordered_map<std::string, int> & subject_kmer_cnt,
                      const std::unordered_set<std::string> & subject_kmers)
{
    
    std::unordered_set<std::string> union_kmers;
    std::unordered_map<std::string, int> query_kmer_cnt = generate_kmer(query, k, union_kmers);
    
    /*
    for (const auto & kmer : subject_kmers) union_kmers.insert(kmer);
    
    double dist = 0.0;
    for (const auto & kmer : union_kmers) {
        
        double d = static_cast<double>(kmer_cnt_lookup(kmer, query_kmer_cnt) - kmer_cnt_lookup(kmer, subject_kmer_cnt));
        dist += (d * d);
    }
    
    return std::sqrt(dist);
    */
    return 0.1;
}
