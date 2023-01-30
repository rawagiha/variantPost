#ifndef LOCAL_ALIGNMENT_LIB_H
#define LOCAL_ALIGNMENT_LIB_H

#include <string>
#include <vector>
#include <utility>
#include <string.h>

#include "ssw/ssw_cpp.h"

typedef StripedSmithWaterman::Alignment Alignment;
typedef StripedSmithWaterman::Aligner Aligner;
typedef StripedSmithWaterman::Filter Filter;

/*
struct InputRead
{
    std::string read_seq;
    std::string base_qualities;
    int target_start;
    
    InputRead(
        const std::string & read_seq, 
        const std::string & base_qualities,
        int target_start = -1
    ) : read_seq(read_seq), 
        base_qualities(base_qualities), 
        target_start(target_start) {}

   InputRead() {}
};
*/

struct SimplifiedRead
{
    std::string seq;
    std::string base_qualities;
    int target_start;

    SimplifiedRead() {}

    SimplifiedRead(
        const std::string & seq,
        const std::string & base_qualities,
        int target_start = -1
    ) : seq(seq),
        base_qualities(base_qualities),
        target_start(target_start) {}
};

/*
struct MergedRead
{
    std::string merged_seq;
    std::string merged_qualities;
    int target_start;

    MergedRead(
        const std::string & merged_seq,
        const std::string & merged_qualities,
        const int target_start
    ) : merged_seq(merged_seq), 
        merged_qualities(merged_qualities), 
        target_start(target_start) {}

    MergedRead() {}
};
*/


SimplifiedRead merge_reads(const std::vector<SimplifiedRead> & inputs);

SimplifiedRead pairwise_stitch(const std::vector<SimplifiedRead> & inputs, bool lt_extention);


char match_to_contig
(
            const std::string & query, const bool is_dirty_query,
            const std::string & contig_seq, const std::string & ref_contig_seq, const std::vector<std::pair<int, int>> & decomposed_contig, const bool is_dirty_contig,
            const int n_tandem_repeats, const std::string & repeat_unit, const std::string & rv_repeat_unit, const bool is_complete_tandem_repeat,
            const std::pair<int, int> & repeat_boundary,
            const Filter & filter,
            const Aligner & aligner, 
            Alignment & alignment
);          

/*
struct Alignment {
    uint16_t alignment_score;
    int32_t  ref_begin;
    int32_t  ref_end;
    int32_t  query_begin;
    int32_t  query_end;
    std::string cigar_string;

    Alignment( uint16_t alignment_score, int32_t  ref_begin, int32_t  ref_end,
               int32_t  query_begin, int32_t  query_end, const std::string & cigar_string );
};

Alignment align( const std::string & ref,
                 const std::string & query,
                 const uint8_t & match_score,
                 const uint8_t & mismatch_penalty,
                 const uint8_t & gap_open_penalty,
                 const uint8_t & gap_extending_penalty );

struct ParsedVariant {
    bool is_indel;
    bool is_ins;
    bool is_del;
    uint16_t variant_len;
    std::string lt_ref;
    std::string lt_query;
    std::string lt_clipped_segment;
    std::string ins_seq;
    std::string del_seq;
    std::string ref_base;
    std::string alt_base;
    std::string rt_ref;
    std::string rt_query;
    std::string rt_clipped_segment;
    uint32_t genomic_pos;
};

std::vector<ParsedVariant> find_variants(
    const Alignment & alignment,
    const std::string & ref,
    const std::string & query,
    const uint32_t & genomic_ref_start = 0 );


std::pair<std::string, std::string>  flatten_reads(const std::pair<std::string, std::string> seed_read,
                                                   const std::vector<std::pair<std::string, std::string>> & reads);


char match_to_contig(const std::string & query, const bool is_dirty_query,
                     const std::string & contig_seq, const std::string & ref_contig_seq, const std::vector<std::pair<int, int>> & decomposed_contig, const bool is_dirty_contig,
                     const int n_tandem_repeats, const std::string & repeat_unit, const std::string & rv_repeat_unit, const bool is_complete_tandem_repeat,
                     const std::pair<int, int> & repeat_boundary);
*/

/*
char is_compatible(const std::string & contig,
                   const std::string & ref_contig,
                   const std::string & query,
                   const std::vector<std::pair<int, int>> & decomposed_contig,
                   const std::string & repeat_unit,
                   const std::string & reversed_repeat_unit,
                   const int expected_num_repeats,
                   const bool is_complete_tandem_repeat,
                   const std::pair<int, int> & boundary_indexes,
                   const bool is_dirty);
*/
/*
void merge_reads(std::vector<std::string> & seqs, std::vector<std::string> & seq_quals);
*/

#endif
