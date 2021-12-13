#include <vector>
#include <string>
#include <iostream>


#include "pileup_parser.h"
#include "util.h"


pileup::ParsedRead::ParsedRead( int unspliced_local_reference_start,
                                const std::string & unspliced_local_reference,
                                const std::string & read_name,
                                const bool is_reverse,
                                const std::string & cigar_string,
                                const int aln_start,
                                const int aln_end,
                                const std::string & read_seq,
                                const std::string & ref_seq,
                                const std::vector<int> & q,
                                const int mapq )
{
    read_name_ = read_name;
    is_reverse_ = is_reverse;

    cigar_string_ = cigar_string;
    cigar_vector_ = to_cigar_vector( cigar_string );

    aln_start_ = aln_start + 1; // 1-based
    aln_end_ = aln_end;  // pysam aln_end is 1-based
    
    read_seq_ = read_seq;
    if ( cigar_string.find( 'N' ) != std::string::npos ) {
        ref_seq_ = ref_seq;
        is_spliced_ = true;
    }
    else {
        ref_seq_ = get_read_wise_ref_seq( aln_start, aln_end,
                                         unspliced_local_reference_start, unspliced_local_reference );
        is_spliced_ = false;
    }

    base_qualities_ = to_fastq_qual( q );
    mapq_ = mapq;
}

void pileup::parse_pileup( int unspliced_local_reference_start,
                           const std::string & unspliced_local_reference,
                           const std::vector<std::string> & read_names,
                           const std::vector<bool> & are_reverse,
                           const std::vector<std::string> & cigar_strings,
                           const std::vector<int> & aln_starts,
                           const std::vector<int> & aln_ends,
                           const std::vector<std::string> & read_seqs,
                           const std::vector<std::string> & ref_seqs,
                           const std::vector<std::vector<int>> & quals,
                           const std::vector<int> & mapqs )
{
    size_t pileup_size = read_names.size();

    std::vector<pileup::ParsedRead> parsed;


    for ( size_t i = 0; i < pileup_size; ++i ) {
        parsed.emplace_back( pileup::ParsedRead(
                                 unspliced_local_reference_start,
                                 unspliced_local_reference,
                                 read_names[i],
                                 are_reverse[i],
                                 cigar_strings[i],
                                 aln_starts[i],
                                 aln_ends[i],
                                 read_seqs[i],
                                 ref_seqs[i],
                                 quals[i],
                                 mapqs[i]
                             ) 
        );
    }
}
