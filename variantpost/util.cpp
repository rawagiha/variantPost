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

/*
std::vector<std::string> split( const std::string & seq, int aln_start,
                                int aln_end, int pos, const std::vector<std::string> & cigar_vector )
{

}
*/

// position indexed local reference
//-----------------------------------------------------------------------------
std::map<int, char> pos_index_reference( const std::string &
        unspliced_local_reference, int unspliced_local_reference_start,
        int unspliced_local_reference_end )
{
    std::map<int, char> indexed_local_reference;

    int i = 0;
    int pos = unspliced_local_reference_start;
    while ( pos <= unspliced_local_reference_end ) {
        indexed_local_reference[pos] = unspliced_local_reference[i];
        ++i;
        ++pos;
    }

    return indexed_local_reference;
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

    is_substitute_ = ( alt_len_ == ref_len_ );
    is_ins_ = ( alt_len_ > ref_len_ );
    is_del_ = ( alt_len_ < ref_len_ );

    variant_end_pos_ = pos_ + ref_len_;
}


inline bool is_rotatable( const std::string & allele )
{
    return ( allele[0] == allele[allele.size() - 1] );
}


inline void to_right( int & variant_end_pos, std::string & longer_allele,
                      std::string & shorter_allele,
                      const std::map<int, char> & indexed_local_reference )
{
    char next_base = indexed_local_reference.at( variant_end_pos );
    longer_allele.erase( 0, 1 );
    longer_allele += next_base;
    shorter_allele = next_base;
    ++variant_end_pos;
}

inline void to_left ( int & pos, std::string & longer_allele,
                      std::string & shorter_allele,
                      const std::map<int, char> & indexed_local_reference )
{
    --pos;
    char prev_base = indexed_local_reference.at( pos );
    std::cout << prev_base << " :: " << std::endl;
    longer_allele.pop_back();
    longer_allele.insert(0, 1, prev_base);
    shorter_allele = prev_base;
}

void left_align( int & pos, std::string & ref, std::string & alt, bool is_ins,
                 const int unspliced_local_reference_start,
                 const std::map<int, char> & indexed_local_reference )
{
    std::string & longer_allele = (is_ins) ? alt : ref;
    std::string & shorter_allele = (is_ins) ? ref : alt;

    std::cout << "prev " << longer_allele << " " << shorter_allele << std::endl;
    while ( ( is_rotatable( longer_allele ) ) & ( unspliced_local_reference_start <=
            pos ) ) {
        to_left( pos, longer_allele, shorter_allele, indexed_local_reference );
    }
    std::cout << "post " << longer_allele << " " << shorter_allele << std::endl;
}


bool Variant::is_shiftable()
{

    if ( is_substitute_ ) {
        return false;
    }

    if ( is_ins_ ) {
        if ( is_rotatable( alt_ ) ) {
            return true;
        }
        else {
            std::string longer = alt_;
            std::string shorter = ref_;
            int _variant_end_pos = variant_end_pos_;

            to_right( _variant_end_pos, longer, shorter, indexed_local_reference_ );

            return is_rotatable( longer );
        }
    }
    else {
        if ( is_rotatable( ref_ ) ) {
            return true;
        }
        else {
            std::string longer = ref_;
            std::string shorter = alt_;
            int _variant_end_pos = variant_end_pos_;

            to_right( _variant_end_pos, longer, shorter, indexed_local_reference_ );

            return is_rotatable( longer );
        }
    }
}



bool Variant::operator == ( const Variant & rhs ) const
{
    std::string lhs_ref = ref_;
    std::string lhs_alt = alt_;
    int lhs_pos = pos_;

    std::string rhs_alt = rhs.alt_;
    std::string rhs_ref = rhs.ref_;
    int rhs_pos = rhs.pos_;

    left_align( lhs_pos, lhs_ref, lhs_alt, is_ins_,
                unspliced_local_reference_start_, indexed_local_reference_ );
    left_align( rhs_pos, rhs_ref, rhs_alt, rhs.is_ins_,
                unspliced_local_reference_start_, indexed_local_reference_ );

    std::cout << lhs_pos << " : " << rhs_pos << std::endl;
    std::cout << lhs_ref << " : " << rhs_ref << std::endl;
    std::cout << lhs_alt << " : " << rhs_alt << std::endl;
    return ( ( lhs_pos == rhs_pos ) & ( lhs_ref == rhs_ref ) &
             ( lhs_alt == rhs_alt ) );
}



//------------------------------------------------------------------------------
