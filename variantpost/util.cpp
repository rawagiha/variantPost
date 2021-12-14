#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "util.h"

inline char to_base_qual( int q )
{
    return static_cast<char>( q + 33 );
}

std::string to_fastq_qual( const std::vector<int> & qvec )
{
    std::vector<char> res( qvec.size() );
    std::transform( qvec.begin(), qvec.end(), res.begin(), to_base_qual );
    std::string fq( res.begin(), res.end() );
    return fq;
}


std::vector<std::string> to_cigar_vector( const std::string & cigar_string )
{
    std::vector<std::string> cigarette;

    size_t pos = 0;
    size_t newpos;
    size_t len = cigar_string.size();
    while ( pos < len ) {
        newpos = cigar_string.find_first_of( "MIDNSHPX=", pos ) + 1;
        cigarette.push_back( cigar_string.substr( pos, newpos - pos ) );
        pos = newpos;
    }

    return cigarette;
}

int get_clipped_len( const std::vector<std::string> & cigar_vector,
                     bool for_start )
{
    int clip_len = 0;

    if ( for_start ) {
        std::string first_op = cigar_vector.front();
        if ( first_op.find( 'S' ) != std::string::npos ) {
            clip_len = std::stoi( first_op.substr( 0, first_op.size() - 1 ) );
        }
    }
    else {
        std::string last_op = cigar_vector.back();
        if ( last_op.find( 'S' ) != std::string::npos ) {
            clip_len = std::stoi( last_op.substr( 0, last_op.size() - 1 ) );
        }
    }

    return clip_len;
}


std::string get_read_wise_ref_seq( int aln_start, int aln_end,
                                   int unspliced_local_reference_start,
                                   const std::string & unspliced_local_reference )
{
    int start_idx = aln_start - 1 - unspliced_local_reference_start;
    return unspliced_local_reference.substr( start_idx,
            ( aln_end - aln_start + 1 ) );
}

std::vector<std::string> split( const std::string & seq, int aln_start,
                                int aln_end, int pos, const std::vector<std::string> & cigar_vector )
{

}


// position indexed local reference
//-----------------------------------------------------------------------------
std::map<int, char> pos_index_reference( const std::string &
        unspliced_local_reference, int unspliced_local_reference_start,
        int unspliced_local_reference_end )
{
    std::map<char, int> indexed_local_reference;

    int i = 0;
    int pos = unspliced_local_reference_start;
    while ( pos <= unspliced_local_reference_end ) {
        indexed_local_reference[pos] = unspliced_local_reference[i];
        ++i;
        ++pos;
    }

    return indexed_local_reference
}

// simplified Variant object used in C++ codes
// ----------------------------------------------------------------------------
Variant::Variant( const std::string & chrom,
                  const int pos,
                  const std::string & ref,
                  const std::string & alt,
                  const int unspliced_local_reference_start,
                  const int unspliced_local_reference_end,
                  const std::map<int, char> & indexed_local_reference )
{
    chrom_ = chrom;
    pos_ = pos;
    ref_ = ref;
    alt_ = alt;
    unspliced_local_reference_start_ = unspliced_local_reference_start;
    unspliced_local_reference_end_ = unspliced_local_reference_end;
    indexed_local_reference_ = indexed_local_reference;

    ref_len_ = ref_.size();
    alt_len_ = alt_.size();

    is_substitute_ = ( alt_len_ == ref_len_ )
                     is_ins_ = ( alt_len_ > ref_len_ );
    is_del_ = ( alt_len_ < ref_len_ );

    variant_end_pos_ = pos_ + ref_len_
}

inline ends_with_same( const string & ref, int ref_len, const string & alt,
                       int alt_len )
{
    return ( ref[ref_len - 1] == alt[alt_len - 1] != 'N' )
}

inline void to_right( int & variant_end_pos, std::string & ref,
                      std::string & alt, const bool is_ins,
                      const std::map<int, char> & indexed_local_reference )
{
    char next_base = indexed_local_reference[variant_end_pos];
    if ( is_ins ) {
        alt.erase( 0, 1 );
        alt += next_base;
        ref = next_base;
    }
    else {
        ref.erase( 0, 1 );
        ref += next_base;
        alt = next_base;
    }
    ++variant_end_pos;
}

bool Variant::is_shiftable()
{

    if ( is_substitute_ ) {
        return false;
    }


    // true if left_alignable
    if ( ends_with_same( alt_, alt_len_, ref_, ref_len_ ) ) {
        return true;
    }

    // check if right_alignable
    int _variant_end_pos = variant_end_pos_;
    std::string _ref = ref_;
    std::string _alt_ = alt_;

    to_right( _variant_end_pos, _ref, _alt, is_ins_, indexed_local_reference_ );

    return ends_with_same( _ref, ref_len_, _alt, alt_len );
}

/*
bool Variant::operator == ( const Variant &rhs ) const
{
    std::string lhs_alt = alt_;
    std::string lhs_ref = ref_;

    std::string rhs_alt = rhs.alt_;
    std::string rhs_ref = rhs.ref_;
}
*/


//------------------------------------------------------------------------------
