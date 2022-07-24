#include <iostream>
#include <algorithm>
#include <string>
#include <deque>
#include <vector>
#include <iterator>
//#include <regex>
#include <random>
#include <string.h>


#include "ssw/ssw_cpp.h"
#include "swlib.h"


sw::Alignment::Alignment( uint16_t __alignment_score, int32_t  __ref_begin,
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


sw::Alignment sw::align( const std::string & ref,
                         const std::string & query,
                         const uint8_t & match_score,
                         const uint8_t & mismatch_penalty,
                         const uint8_t & gap_open_penalty,
                         const uint8_t & gap_extending_penalty )
{

    int32_t mask_len = strlen( query.c_str() ) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;

    StripedSmithWaterman::Alignment __alignment;

    StripedSmithWaterman::Aligner aligner( match_score,
                                           mismatch_penalty,
                                           gap_open_penalty,
                                           gap_extending_penalty );

    StripedSmithWaterman::Filter filter;

    aligner.Align( query.c_str(), ref.c_str(), ref.size(), filter, &__alignment,
                   mask_len );


    sw::Alignment alignment( __alignment.sw_score, __alignment.ref_begin,
                             __alignment.ref_end, __alignment.query_begin, __alignment.query_end,
                             __alignment.cigar_string );

    return alignment;
}

/*
//@Function:
//      Parse CIGAR to a list. 60=2D8I32=1S -> [60=, 2D, 8I, 32=, 1S]
std::vector<std::string> decompose_cigar_string( const std::string &
        cigar_string )
{
    std::regex cigar_pattern( R"(\d+[MIDNSHPX=])" );

    std::regex_iterator<std::string::const_iterator> cigar_itr(
        cigar_string.begin(),
        cigar_string.end(), cigar_pattern );
    std::regex_iterator<std::string::const_iterator> end_of_itr;

    std::vector<std::string> cigarette;
    while ( cigar_itr != end_of_itr ) {
        cigarette.push_back( ( *cigar_itr ).str() );
        ++cigar_itr;
    }

    return cigarette;
}
*/

inline std::vector<std::string> _decompose_cigar_string(const std::string & cigar_string)
{
    std::vector<std::string> cigarette;
    
    size_t pos = 0;
    size_t newpos = 0;
    const size_t len = cigar_string.size();
    
    while (pos < len) {
        newpos = cigar_string.find_first_of("MIDNSHPX=", pos) + 1;
        cigarette.push_back(cigar_string.substr(pos, newpos - pos));
        pos = newpos;
    }
    
    return cigarette; 
}


inline bool is_gap( std::string cigar_itr )
{
    return (
               ( cigar_itr.find( "I" ) != std::string::npos )
               ||
               ( cigar_itr.find( "D" ) != std::string::npos )
           );
}

//@Function:
//      Check if there exist consecutive gaps (I or D)
//      [60=, 2D, 8I, 32=, 1S] -> true
//      [60=, 6I, 2X, 32=, 1S] -> false
bool has_consecutive_gap(const std::vector<std::string> & cigarette)
{
    bool prev_is_gap = false;
    for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {
        if ( ( prev_is_gap ) && is_gap( ( *itr ) ) ) {
            return true;
        }
        prev_is_gap = is_gap( ( *itr ) );
    }
    return false;
}


//@Function:
//      Merge gaps into a signle gap of specified type
//      [4I, 2I] -> "6I"
std::string concat_gaps( const std::vector<std::string> & cigarette,
                         std::string  gap_type )
{
    std::string gaps;

    if ( !cigarette.empty() ) {
        uint16_t total_gap_len = 0;

        for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
                itr != cigarette.end(); ++itr ) {
            total_gap_len += std::stoi( ( *itr ).substr( 0, ( *itr ).length() - 1 ) );
        }

        gaps = std::to_string( total_gap_len ) + gap_type;
    }

    return gaps;
}


//@Function:
//      Merge consecutive gaps and make insertion come first
//      [4=, 2I, 2D, 1I, 3=, 3D, 1I, 2D, 4I, 4=]
//      -> [4=, 3I, 2D, 3=, 5I, 5D, 4=]
void edit_cigar( std::vector<std::string> & cigarette )
{
    std::vector<std::string> tmp, ins, del;

    bool prev_is_gap = false;
    for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {
        if ( is_gap( ( *itr ) ) ) {
            if ( ( *itr ).find( "I" ) != std::string::npos ) {
                ins.push_back( *itr );
            }
            else {
                del.push_back( *itr );
            }
            prev_is_gap = true;
        }
        else {
            if ( prev_is_gap ) {
                std::string merged_ins = concat_gaps( ins, "I" );
                if ( !merged_ins.empty() ) {
                    tmp.push_back( merged_ins );
                }
                std::string merged_del = concat_gaps( del, "D" );
                if ( !merged_del.empty() ) {
                    tmp.push_back( merged_del );
                }
                tmp.push_back( *itr );
                ins.clear();
                del.clear();
            }
            else {
                tmp.push_back( *itr );
            }

            prev_is_gap = false;
        }
    }
    std::swap( cigarette, tmp );
}


std::vector<sw::ParsedVariant> sw::find_variants( const sw::Alignment &
        alignment,
        const std::string & ref,
        const std::string & query,
        const uint32_t & genomic_ref_start )
{

    std::vector<std::string> cigarette = _decompose_cigar_string(
            alignment.cigar_string );
    std::vector<sw::ParsedVariant> variants;

    if ( has_consecutive_gap( cigarette ) ) {
        edit_cigar( cigarette );
    }

    uint32_t genomic_pos = ( genomic_ref_start > 0 ) ? genomic_ref_start - 1 : 0;

    uint32_t ref_idx = alignment.ref_begin, query_idx = alignment.query_begin;

    for ( std::vector<std::string>::iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {

        char operation = ( *itr ).back(); // cigar operation
        uint16_t op_len = std::stoi( ( *itr ).substr( 0,
                                     ( *itr ).length() - 1 ) ); // operation length

        if ( operation == 'I' ) {
            sw::ParsedVariant ins;

            ins.is_indel = true;
            ins.is_ins = true;
            ins.is_del = false;

            ins.lt_ref = ref.substr( 0, ref_idx );
            ins.lt_query = query.substr( 0, query_idx );
            ins.ins_seq = query.substr( query_idx, op_len );
            ins.variant_len = op_len;
            ins.rt_ref = ref.substr( ref_idx );
            ins.rt_query = query.substr( query_idx + op_len );

            ins.lt_clipped_segment = query.substr( 0, alignment.query_begin );
            ins.rt_clipped_segment = query.substr( alignment.query_end + 1,
                                                   query.length() - alignment.query_end );

            ins.genomic_pos = genomic_pos;

            variants.push_back( ins );

            query_idx += op_len;

        }
        else if ( operation == 'D' ) {
            sw::ParsedVariant del;

            del.is_indel = true;
            del.is_ins = false;
            del.is_del = true;

            del.lt_ref = ref.substr( 0, ref_idx );
            del.lt_query = query.substr( 0, query_idx );
            del.del_seq = ref.substr( ref_idx, op_len );
            del.variant_len = op_len;
            del.rt_ref = ref.substr( ref_idx + op_len );
            del.rt_query = query.substr( query_idx );

            del.lt_clipped_segment = query.substr( 0, alignment.query_begin );
            del.rt_clipped_segment = query.substr( alignment.query_end + 1,
                                                   query.length() - alignment.query_end );

            del.genomic_pos = genomic_pos;

            variants.push_back( del );

            ref_idx += op_len;
            genomic_pos += op_len;
        }
        else if ( operation == 'X' ) {
            sw::ParsedVariant smv;      // single or multi-nucleotide variant

            smv.is_indel = false;
            smv.is_ins = false;
            smv.is_del = false;

            smv.lt_ref = ref.substr( 0, ref_idx );
            smv.lt_query = query.substr( 0, query_idx );
            smv.ref_base = ref.substr( ref_idx, op_len );
            smv.alt_base = query.substr( query_idx, op_len );
            smv.variant_len = op_len;
            smv.rt_ref = ref.substr( ref_idx + op_len );
            smv.rt_query = query.substr( query_idx + op_len );

            smv.lt_clipped_segment = query.substr( 0, alignment.query_begin );
            smv.rt_clipped_segment = query.substr( alignment.query_end + 1,
                                                   query.length() - alignment.query_end );

            smv.genomic_pos = genomic_pos;

            variants.push_back( smv );

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


/*
namespace bases {
char bases[5] = {'A', 'C', 'G', 'T', 'N'};
}
*/

static char bases[5] = {'A', 'C', 'G', 'T', 'N'};

template <typename T>
std::vector<int> get_max_indices(const T & arr)
{
    std::vector<int> indices;

    auto it = std::max_element(std::begin(arr), std::end(arr));
    while (it != std::end(arr)) {
        indices.push_back(std::distance(std::begin(arr), it));
        it = std::find(std::next(it), std::end(arr), *it);
    }
    return indices;
}


struct BaseCount {

    double a = 0.0, c = 0.0, g = 0.0, t = 0.0, n = 0.0;
    double a_cnt = 0.0001, c_cnt = 0.0001, g_cnt = 0.0001, t_cnt = 0.0001,
           n_cnt = 0.0001;

    BaseCount() {}

    BaseCount(char base, char qual)
    {
        double _qual = static_cast<double>(qual);

        switch (base) {
        case 'A':
            a = _qual;
            ++a_cnt;
            break;
        case 'C':
            c = _qual;
            ++c_cnt;
            break;
        case 'G':
            g = _qual;
            ++g_cnt;
            break;
        case 'T':
            t = _qual;
            ++t_cnt;
            break;
        default:
            n = _qual;
            ++n_cnt;
            break;
        }
    }

    void add(char base, char qual)
    {
        int _qual = static_cast<int>(qual);

        switch (base) {
        case 'A':
            a += _qual;
            ++a_cnt;
            break;
        case 'C':
            c += _qual;
            ++c_cnt;
            break;
        case 'G':
            g += _qual;
            ++g_cnt;
            break;
        case 'T':
            t += _qual;
            ++t_cnt;
            break;
        default:
            n += _qual;
            ++n_cnt;
            break;
        }
    }

    std::vector<char> get_consensus()
    {
        double counts[] = {a_cnt, c_cnt, g_cnt, t_cnt, n_cnt};
        double quals[] = {a, c, g, t, n};

        std::vector<int> max_cnt_indices = get_max_indices(counts);

        int max_idx;
        if ( max_cnt_indices.size() > 1 ) {
            std::vector<int> max_qual_indices = get_max_indices(quals);
            max_idx = max_qual_indices[0];
        }
        else {
            max_idx = max_cnt_indices[0];
        }

        std::vector<char> ret = {bases[max_idx], static_cast<char>(quals[max_idx] / counts[max_idx])};

        return ret;
    }
};


void update( std::deque<BaseCount> & consensus,  const int & read2_begin,
             const std::string & lt_ext, const std::string & lt_qual,
             const std::string & mread2, const std::string & mqual2,
             const std::string & rt_ext, const std::string & rt_qual )
{
    if ( read2_begin == 0 ) {

        // update middle part
        const auto lt_len = lt_ext.size();
        for ( size_t i = lt_len; (i < consensus.size() && (i - lt_len) < mread2.size() - 1); ++i ) {
            consensus[i].add( mread2[i - lt_len], mqual2[i - lt_len] );
        }

        // rt extention
        for ( size_t i = 0; i < rt_ext.size(); ++i ) {
            consensus.emplace_back( BaseCount( rt_ext[i], rt_qual[i] ) );
        }
    }
    else {
        // update middle part
        for ( size_t i = 0; i < mread2.size(); ++i ) {
            consensus[i].add( mread2[i], mqual2[i] );
        }

        //lt extention
        for ( int i = lt_ext.size() - 1; i >= 0; --i ) {
            consensus.emplace_front( BaseCount( lt_ext[i], lt_qual[i] ) );
        }
    }
}


std::vector<std::string> get_consensus_contig( std::deque<BaseCount> &
        consensus )
{

    std::string consensus_seq;
    std::string consensus_qual;
    for ( std::deque<BaseCount>::iterator itr = consensus.begin();
            itr != consensus.end(); ++itr ) {
        std::vector<char> c = ( *itr ).get_consensus();
        consensus_seq += c[0];
        consensus_qual += c[1];
    }

    std::vector<std::string> ret {consensus_seq, consensus_qual};

    return ret;
}


void pairwise_stitch( std::deque<BaseCount> & consensus,
                      const std::vector<std::string> & v )
{

    std::vector<std::string> c = get_consensus_contig( consensus );
    
    std::string read1 = c[0];
    std::string qual1 = c[1];

    std::string read2 = v[0];
    std::string qual2 = v[1];

    // gap-less alignment
    const int match_score = 2;
    const int mismatch_penalty = 10;
    const int gap_open_penalty = std::max( read1.size(), read2.size() );
    const int gap_extention_penalty = 1;

    sw::Alignment aln = sw::align( read1, read2, match_score, mismatch_penalty,
                                   gap_open_penalty, gap_extention_penalty );

    std::vector<std::string> cigar_vec = _decompose_cigar_string(
            aln.cigar_string );
    
    const int read1_begin = aln.ref_begin;
    const int read1_end = aln.ref_end;
    const int read2_begin = aln.query_begin;
    const int read2_end = aln.query_end;

    //std::cout << aln.cigar_string << ", " << read1_begin << ", " << read1_end << ", " << read2_begin << ", " << read2_end << std::endl;
    
    std::string lt_ext, lt_qual, rt_ext, rt_qual, mread1, mqual1, mread2, mqual2;
    if ( read2_begin == 0 ) {
        lt_ext = read1.substr( 0, read1_begin );
        lt_qual = qual1.substr( 0, read1_begin );
        rt_ext = read2.substr( read2_end + 1 );
        rt_qual = qual2.substr( read2_end + 1 );

        mread2 = read2.substr( 0, read2_end + 1 );
        mqual2 = qual2.substr( 0, read2_end + 1 );
    }
    else if ( read1_begin == 0 ) {
        lt_ext = read2.substr( 0, read2_begin );
        lt_qual = qual2.substr( 0, read2_begin );
        rt_ext = read1.substr( read1_end + 1 );
        rt_qual = qual1.substr( read1_end + 1 );
        
        size_t m_len = read2_end - read2_begin + 1;
        mread2 = read2.substr( read2_begin, m_len );
        mqual2 = qual2.substr( read2_begin, m_len );
    }
    else {
        // not stitchable -> reconsider
    }


    // check for stichability needed
    update( consensus, read2_begin, lt_ext, lt_qual, mread2, mqual2, rt_ext,
            rt_qual );
   
}


std::string sw::flatten_reads( const std::vector<std::string> seed_read,
                               const std::vector<std::vector<std::string>> & reads )
{
    std::deque<BaseCount>  consensus;
    std::string seed_bases = seed_read[0];
    std::string seed_quals = seed_read[1];

    for ( size_t i = 0; i < seed_bases.size(); ++i ) {
        consensus.emplace_back( BaseCount( seed_bases[i], seed_quals[i] ) );
    }


    for ( std::vector<std::vector<std::string>>::const_iterator itr = reads.begin();
            itr != reads.end(); ++itr ) {
        pairwise_stitch( consensus, *itr );
    }


    std::string ans = "";
    std::string qual = "";
    for ( size_t i = 0; i < consensus.size(); ++i ) {
        BaseCount c = consensus[i];
        ans += c.get_consensus()[0];
        qual += c.get_consensus()[1];
        //std::cout << static_cast<int>( c.get_consensus()[1] ) << ", ";
    }

    std::cout << ans << std::endl;
    std::cout << qual << std::endl;

    return "str";
}

