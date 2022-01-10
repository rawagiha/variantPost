#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>

#include "util.h"


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

// returns a list of intron start and stop positions (inclusive: [start, start])
std::vector<std::pair<int, int>> get_introns ( const
std::vector<std::pair<char, int>> &cigar_vector, const int aln_start ) {

  std::vector<std::pair<int, int>> introns;

  char operation;
  int operation_len;
  int current_pos = aln_start - 1;
  for ( std::vector<std::pair<char, int>>::const_iterator itr =
          cigar_vector.begin();
        itr != cigar_vector.end(); ++itr ) {

    operation = ( *itr ).first;
    operation_len = ( *itr ).second;

    switch ( operation ) {
      case 'M':
        current_pos += operation_len;
        break;
      case 'I':
        break;
      case 'D':
        current_pos += operation_len;
        break;
      case 'N':
        break;


    }
  }
}

// fit local referenece to read alignment
//-----------------------------------------------------------------------------
std::string get_read_wise_ref_seq ( int aln_start, int aln_end,
                                    int unspliced_local_reference_start,
                                    const std::string &unspliced_local_reference ) {
  int start_idx = aln_start - unspliced_local_reference_start;
  return unspliced_local_reference.substr ( start_idx,
         ( aln_end - aln_start + 1 ) );
}


// position indexed local reference
//-----------------------------------------------------------------------------
std::map<int, char> pos_index_reference ( const std::string &
    unspliced_local_reference, int unspliced_local_reference_start,
    int unspliced_local_reference_end ) {
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
/*
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
*/

Variant::Variant ( const std::string &chrom, int pos, const std::string &ref,
                   const std::string &alt ) : chrom_ ( chrom ), pos_ ( pos ), ref_ ( ref ),
  alt_ ( alt ), ref_len_ ( ref_.size() ), alt_len_ ( alt_.size() ),
  is_substitute_ ( ( alt_len_ == ref_len_ ) ),
  is_ins_ ( ( alt_len_ > ref_len_ ) ),
  is_del_ ( ( alt_len_ < ref_len_ ) )
{}


inline bool is_rotatable( const std::string & allele )
{
    return ( allele[0] == allele[allele.size() - 1] );
}


inline void to_left ( int & pos, std::string & longer_allele,
                      std::string & shorter_allele,
                      const std::map<int, char> & indexed_local_reference )
{
    --pos;
    char prev_base = indexed_local_reference.at( pos );
    longer_allele.pop_back();
    longer_allele.insert( 0, 1, prev_base );
    shorter_allele = prev_base;
}


void left_align( int & pos, std::string & ref, std::string & alt, bool is_ins,
                 const int unspliced_local_reference_start,
                 const std::map<int, char> & indexed_local_reference )
{
    std::string & longer_allele = ( is_ins ) ? alt : ref;
    std::string & shorter_allele = ( is_ins ) ? ref : alt;

    while ( ( is_rotatable( longer_allele ) ) & ( unspliced_local_reference_start <=
            pos ) ) {
        to_left( pos, longer_allele, shorter_allele, indexed_local_reference );
    }
}


int Variant::get_leftmost_pos(const int unspliced_local_reference_start, const std::map<int, char> & indexed_local_reference) const
{
    int pos = pos_;
    std::string ref = ref_ ;
    std::string alt = alt_;

    left_align( pos, ref, alt, is_ins_, unspliced_local_reference_start,
                indexed_local_reference );

    return pos;
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


void right_align( int & pos, int & variant_end_pos, std::string & ref,
                  std::string & alt, bool is_ins,
                  const int unspliced_local_reference_end,
                  const std::map<int, char> & indexed_local_reference )
{
    std::string & longer_allele = ( is_ins ) ? alt : ref;
    std::string & shorter_allele = ( is_ins ) ? ref : alt;

    do {
        to_right( variant_end_pos, longer_allele, shorter_allele,
                  indexed_local_reference );
        ++pos;
    }
    while ( ( is_rotatable ( longer_allele ) ) & ( pos <=
            unspliced_local_reference_end ) );

    to_left( pos, longer_allele, shorter_allele,
             indexed_local_reference ); // undo the last right shift
}


int Variant::get_rightmost_pos(const int unspliced_local_reference_end, const std::map<int, char> & indexed_local_reference) const
{
    int pos = pos_;
    int variant_end_pos = variant_end_pos_;
    std::string ref = ref_ ;
    std::string alt = alt_;

    right_align( pos, variant_end_pos, ref, alt, is_ins_,
                 unspliced_local_reference_end, indexed_local_reference );

    return pos;
}


bool Variant::is_shiftable(const std::map<int, char> & indexed_local_reference) const
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

            to_right( _variant_end_pos, longer, shorter, indexed_local_reference );

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

            to_right( _variant_end_pos, longer, shorter, indexed_local_reference );

            return is_rotatable( longer );
        }
    }
}


/*
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

    return ( ( lhs_pos == rhs_pos ) & ( lhs_ref == rhs_ref ) &
             ( lhs_alt == rhs_alt ) );
}
*/

// parse mapping to find Variant
//------------------------------------------------------------------------------

/*
inline void append_snv( std::vector<Variant> & variants,
                        const int ref_idx, const int read_idx, const int pos,
                        const std::string & ref_seq, const std::string & read_seq,
                        const std::string & chrom,
                        const int & unspliced_local_reference_start,
                        const int & unspliced_local_reference_end,
                        const std::map<int, char> & indexed_local_reference )
{
    std::string ref = ref_seq.substr( ref_idx, 1 );
    std::string alt = read_seq.substr( read_idx, 1 );

    if ( ref != alt ) {
        variants.emplace_back( chrom, pos, ref, alt,
                               unspliced_local_reference_start, unspliced_local_reference_end,
                               indexed_local_reference );
    }

}
*/

std::vector<Variant> find_mapped_variants ( const int aln_start,
    const int aln_end, const std::string &ref_seq, const std::string &read_seq,
    const std::vector<std::pair<char, int>> &cigar_vector,
    const std::string &chrom,
    const int &unspliced_local_reference_start,
    const int &unspliced_local_reference_end,
    const std::map<int, char> &indexed_local_reference ) {
  std::vector<Variant> variants;

  if ( read_seq == ref_seq ) {
    return variants;
  }

  char operation;
  int operation_len;
  int ref_idx=0;
  int read_idx=0;
  int pos = aln_start;

  for ( std::vector<std::pair<char, int>>::const_iterator itr =
          cigar_vector.begin();
        itr != cigar_vector.end(); ++itr ) {

    operation = ( *itr ).first;
    operation_len = ( *itr ).second;

    switch ( operation ) {
      case 'M':
        for ( int i = 0; i < operation_len; ++i ) {


          std::string ref = ref_seq.substr ( ref_idx, 1 );
          std::string alt = read_seq.substr ( read_idx, 1 );

          if ( ref != alt ) {
            variants.emplace_back ( chrom, pos, ref, alt );
            //std::cout << pos << "-" << ref << "-" << alt << std::endl;
          }
          ++ref_idx;
          ++read_idx;
          ++pos;
        }
        break;
      case 'I':
        variants.emplace_back ( chrom, pos - 1, ref_seq.substr ( ref_idx - 1, 1 ),
                                ref_seq.substr ( ref_idx - 1, 1 ) + read_seq.substr ( read_idx,
                                    operation_len ) );
        //std::cout << pos << "-" << ref_seq.substr( ref_idx - 1, 1 ) << "-" << read_seq.substr( read_idx, operation_len ) << std::endl;
        read_idx += operation_len;
        break;
      case 'D':
        variants.emplace_back ( chrom, pos - 1, ref_seq.substr ( ref_idx - 1,
                                1 + operation_len ), ref_seq.substr ( ref_idx - 1, 1 ) );
        ref_idx += operation_len;
        pos += operation_len;
        break;
      case 'N':
        pos += operation_len;
        break;
      case 'S':
        read_idx += ( operation_len );
        break;
      case 'H':
        //
        break;
      case 'P':
        //
        break;
      case 'X':
        /*
        append_snv( variants, ref_idx, read_idx, pos, ref_seq, read_seq, chrom,
                    unspliced_local_reference_start, unspliced_local_reference_end,
                    indexed_local_reference );

        */
        ++ref_idx;
        ++read_idx;
        ++pos;
        break;
      case '=':
        ++ref_idx;
        ++read_idx;
        ++pos;
        break;
    }

    //MIDNSHPX=

  }
  return variants;
}
