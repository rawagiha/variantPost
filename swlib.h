#ifndef SMITH_WATERMAN_LIB_H
#define SMITH_WATERMAN_LIB_H

#include <string>
#include <vector>
#include <string.h>

namespace sw {

struct Alignment {
    uint16_t alignment_score;
    int32_t  ref_begin;
    int32_t  ref_end;
    int32_t  query_begin;
    int32_t  query_end;
    std::string cigar_string;

    Alignment ( uint16_t alignment_score, int32_t  ref_begin, int32_t  ref_end,
                int32_t  query_begin, int32_t  query_end, const std::string & cigar_string );
};

Alignment align ( const std::string & ref,
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

std::vector<ParsedVariant> find_variants (
    const Alignment & alignment,
    const std::string & ref,
    const std::string & query,
    const uint32_t & genomic_ref_start = 0 );


std::string flatten_reads ( std::vector<std::vector<std::string>> & reads );
} //end of name sw

#endif
