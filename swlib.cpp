#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <regex>
#include <string.h>


#include "ssw_cpp.h"
#include "swlib.h"


sw::Alignment::Alignment ( uint16_t __alignment_score, int32_t  __ref_begin,
                           int32_t  __ref_end, int32_t  __query_begin, int32_t  __query_end,
                           const std::string & __cigar_string )
{
    alignment_score = __alignment_score;
    ref_begin = __ref_begin;
    ref_end = __ref_end;
    query_begin = __query_begin;
    query_end = __query_end;
    cigar_string = __cigar_string;
}


sw::Alignment sw::align ( const std::string & ref,
                          const std::string & query,
                          const uint8_t & match_score,
                          const uint8_t & mismatch_penalty,
                          const uint8_t & gap_open_penalty,
                          const uint8_t & gap_extending_penalty )
{

    int32_t mask_len = strlen ( query.c_str() ) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;

    StripedSmithWaterman::Alignment __alignment;

    StripedSmithWaterman::Aligner aligner ( match_score,
                                            mismatch_penalty,
                                            gap_open_penalty,
                                            gap_extending_penalty );

    StripedSmithWaterman::Filter filter;

    aligner.Align ( query.c_str(), ref.c_str(), ref.size(), filter, &__alignment,
                    mask_len );


    sw::Alignment alignment ( __alignment.sw_score, __alignment.ref_begin,
                              __alignment.ref_end, __alignment.query_begin, __alignment.query_end,
                              __alignment.cigar_string );

    return alignment;
}

//@Function:
//      Parse CIGAR to a list. 60=2D8I32=1S -> [60=, 2D, 8I, 32=, 1S]
std::vector<std::string> decompose_cigar_string ( const std::string &
        cigar_string )
{
    std::regex cigar_pattern ( R"(\d+[MIDNSHPX=])" );

    std::regex_iterator<std::string::const_iterator> cigar_itr (
        cigar_string.begin(),
        cigar_string.end(), cigar_pattern );
    std::regex_iterator<std::string::const_iterator> end_of_itr;

    std::vector<std::string> cigarette;
    while ( cigar_itr != end_of_itr ) {
        cigarette.push_back ( ( *cigar_itr ).str() );
        ++cigar_itr;
    }

    return cigarette;
}

bool is_gap ( std::string cigar_itr )
{
    return (
               ( cigar_itr.find ( "I" ) != std::string::npos )
               ||
               ( cigar_itr.find ( "D" ) != std::string::npos )
           );
}

//@Function:
//      Check if there exist consecutive gaps (I or D)
//      [60=, 2D, 8I, 32=, 1S] -> true
//      [60=, 6I, 2X, 32=, 1S] -> false
bool has_consecutive_gap ( const std::vector<std::string> & cigarette )
{
    bool prev_is_gap = false;
    for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {
        if ( ( prev_is_gap ) && is_gap ( ( *itr ) ) ) {
            return true;
        }
        prev_is_gap = is_gap ( ( *itr ) );
    }
    return false;
}


//@Function:
//      Merge gaps into a signle gap of specified type
//      [4I, 2I] -> "6I"
std::string concat_gaps ( const std::vector<std::string> & cigarette,
                          std::string  gap_type )
{
    std::string gaps;

    if ( !cigarette.empty() ) {
        uint16_t total_gap_len = 0;

        for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
                itr != cigarette.end(); ++itr ) {
            total_gap_len += std::stoi ( ( *itr ).substr ( 0, ( *itr ).length() - 1 ) );
        }

        gaps = std::to_string ( total_gap_len ) + gap_type;
    }

    return gaps;
}


//@Function:
//      Merge consecutive gaps and make insertion come first
//      [4=, 2I, 2D, 1I, 3=, 3D, 1I, 2D, 4I, 4=]
//      -> [4=, 3I, 2D, 3=, 5I, 5D, 4=]
void edit_cigar ( std::vector<std::string> & cigarette )
{
    std::vector<std::string> tmp, ins, del;

    bool prev_is_gap = false;
    for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {
        if ( is_gap ( ( *itr ) ) ) {
            if ( ( *itr ).find ( "I" ) != std::string::npos ) {
                ins.push_back ( *itr );
            }
            else {
                del.push_back ( *itr );
            }
            prev_is_gap = true;
        }
        else {
            if ( prev_is_gap ) {
                std::string merged_ins = concat_gaps ( ins, "I" );
                if ( !merged_ins.empty() ) {
                    tmp.push_back ( merged_ins );
                }
                std::string merged_del = concat_gaps ( del, "D" );
                if ( !merged_del.empty() ) {
                    tmp.push_back ( merged_del );
                }
                tmp.push_back ( *itr );
                ins.clear();
                del.clear();
            }
            else {
                tmp.push_back ( *itr );
            }

            prev_is_gap = false;
        }
    }
    std::swap ( cigarette, tmp );
}


std::vector<sw::ParsedVariant> sw::find_variants ( const sw::Alignment &
        alignment,
        const std::string & ref,
        const std::string & query,
        const uint32_t & genomic_ref_start )
{

    std::vector<std::string> cigarette = decompose_cigar_string (
            alignment.cigar_string );
    std::vector<sw::ParsedVariant> variants;

    if ( has_consecutive_gap ( cigarette ) ) {
        edit_cigar ( cigarette );
    }

    uint32_t genomic_pos = ( genomic_ref_start > 0 ) ? genomic_ref_start - 1 : 0;

    uint32_t ref_idx = alignment.ref_begin, query_idx = alignment.query_begin;

    for ( std::vector<std::string>::iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {

        char operation = ( *itr ).back(); // cigar operation
        uint16_t op_len = std::stoi ( ( *itr ).substr ( 0,
                                      ( *itr ).length() - 1 ) ); // operation length

        if ( operation == 'I' ) {
            sw::ParsedVariant ins;

            ins.is_indel = true;
            ins.is_ins = true;
            ins.is_del = false;

            ins.lt_ref = ref.substr ( 0, ref_idx );
            ins.lt_query = query.substr ( 0, query_idx );
            ins.ins_seq = query.substr ( query_idx, op_len );
            ins.variant_len = op_len;
            ins.rt_ref = ref.substr ( ref_idx );
            ins.rt_query = query.substr ( query_idx + op_len );

            ins.lt_clipped_segment = query.substr ( 0, alignment.query_begin );
            ins.rt_clipped_segment = query.substr ( alignment.query_end + 1,
                                                    query.length() - alignment.query_end );

            ins.genomic_pos = genomic_pos;

            variants.push_back ( ins );

            query_idx += op_len;

        }
        else if ( operation == 'D' ) {
            sw::ParsedVariant del;

            del.is_indel = true;
            del.is_ins = false;
            del.is_del = true;

            del.lt_ref = ref.substr ( 0, ref_idx );
            del.lt_query = query.substr ( 0, query_idx );
            del.del_seq = ref.substr ( ref_idx, op_len );
            del.variant_len = op_len;
            del.rt_ref = ref.substr ( ref_idx + op_len );
            del.rt_query = query.substr ( query_idx );

            del.lt_clipped_segment = query.substr ( 0, alignment.query_begin );
            del.rt_clipped_segment = query.substr ( alignment.query_end + 1,
                                                    query.length() - alignment.query_end );

            del.genomic_pos = genomic_pos;

            variants.push_back ( del );

            ref_idx += op_len;
            genomic_pos += op_len;
        }
        else if ( operation == 'X' ) {
            sw::ParsedVariant smv;      // single or multi-nucleotide variant

            smv.is_indel = false;
            smv.is_ins = false;
            smv.is_del = false;

            smv.lt_ref = ref.substr ( 0, ref_idx );
            smv.lt_query = query.substr ( 0, query_idx );
            smv.ref_base = ref.substr ( ref_idx, op_len );
            smv.alt_base = query.substr ( query_idx, op_len );
            smv.variant_len = op_len;
            smv.rt_ref = ref.substr ( ref_idx + op_len );
            smv.rt_query = query.substr ( query_idx + op_len );

            smv.lt_clipped_segment = query.substr ( 0, alignment.query_begin );
            smv.rt_clipped_segment = query.substr ( alignment.query_end + 1,
                                                    query.length() - alignment.query_end );

            smv.genomic_pos = genomic_pos;

            variants.push_back ( smv );

            ref_idx += op_len;
            query_idx += op_len;
            genomic_pos += op_len;
        }
        else if ( operation == 'S' ) {
            //pass for softclips
        }
        else {
            ref_idx += op_len;
            query_idx += op_len;
            genomic_pos += op_len;
        }

    }

    return variants;
}


void average_quals ( std::string & v1, std::string & v2 )
{
    auto v1_len = v1.size();
    auto v2_len = v2.size();

    //std::assert(v1_len == v2_len);

    for ( int i = 0; i < v1_len; i ++ ) {

        int q1 = static_cast<int> ( v1[i] );
        int q2 = static_cast<int> ( v2[i] );
        char q = static_cast<char> ( static_cast<int> ( ( q1 + q2 ) / 2 ) );
        //char q = 'l';
        v1[i] = q;
        v2[i] = q;
    }

}


//expect cleaneds read as input
std::string sw::stitch_two_reads ( const std::vector<std::string> & v1,
                                   const std::vector<std::string> & v2 )
{

    std::string read1 = v1[0];
    std::string qual1 = v1[1];
    std::string read2 = v2[0];
    std::string qual2 = v2[1];

    // gap-less alignment
    const int match_score = 2;
    const int mismatch_penalty = 0;
    const int gap_open_penalty = std::max ( read1.size(), read2.size() );
    const int gap_extention_penalty = 1;
    sw::Alignment aln = sw::align ( read1, read2, match_score, mismatch_penalty,
                                    gap_open_penalty, gap_extention_penalty );

    std::vector<std::string> cigar_vec = decompose_cigar_string (
            aln.cigar_string );
    const int read1_begin = aln.ref_begin;
    const int read1_end = aln.ref_end;
    const int read2_begin = aln.query_begin;
    const int read2_end = aln.query_end;

    // TODOtest stitchable
    // do by a separate func

    std::cout << aln.cigar_string << std::endl;
    std::cout << read1_begin << ", " << read1_end << std::endl;
    std::cout << read2_begin << ", " << read2_end << std::endl;

    std::string lt_ext, rt_ext, mread1, mqual1, mread2, mqual2;
    if ( read2_begin == 0 ) {
        lt_ext = read1.substr ( 0, read1_begin );
        rt_ext = read2.substr ( read2_end + 1 );

        mread1 = read1.substr ( read1_begin );
        mqual1 = qual1.substr ( read1_begin );
        mread2 = read2.substr ( 0, read2_end + 1 );
        mqual2 = qual2.substr ( 0, read2_end + 1 );
    }
    else if ( read1_begin == 0 ) {
        lt_ext = read2.substr ( 0, read2_begin );
        rt_ext = read1.substr ( read1_end + 1 );

        mread1 = read1.substr ( 0, read1_end + 1 );
        mqual1 = qual1.substr ( 0, read1_end + 1 );
        mread2 = read2.substr ( read2_begin );
        mqual2 = qual2.substr ( read2_begin );
    }
    else {
        // not stitchable
    }

    // do middle part
    std::string mid = "";
    const std::string::size_type mid_len = mread1.size();
    for ( int i = 0; i < mid_len; ++i ) {
        if ( mread1[i] == mread2[i] ) {
            mid += mread1[i];
        }
        else {
            mid += ( mqual1[i] >= mqual2[i] ) ? mread1[i] : mread2[i];
        }
    }

    average_quals ( mqual1, mqual2 );

    std::string stitched_read = lt_ext + mid + rt_ext;

    return stitched_read;
}



