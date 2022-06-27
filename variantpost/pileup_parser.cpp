#include <map>
#include <random>
#include <vector>
#include <string>
#include <utility>
#include <iostream>


#include "pileup_parser.h"
#include "util.h"

/* read filter */


std::vector<pileup::ParsedRead> get_gapped_seed_reads(std::vector<pileup::ParsedRead> & parsed_reads, size_t n);


pileup::ParsedRead::ParsedRead( int unspliced_local_reference_start,
                                int unspliced_local_reference_end,
                                const std::string & unspliced_local_reference,
                                const std::string & read_name,
                                const bool is_reverse,
                                const std::string & cigar_string,
                                const int aln_start,
                                const int aln_end,
                                const std::string & read_seq,
                                const std::string & ref_seq,
                                const std::vector<int> & q,
                                const int mapq,
                                const Variant & target, 
                                const std::map<int, char> & indexed_local_reference)
{
    const std::string chrom = target.chrom_;
    read_name_ = read_name;
    is_reverse_ = is_reverse;

    cigar_string_ = cigar_string;   // may be removed later 
    cigar_vector_ = to_cigar_vector( cigar_string ); // may be removed later 

    aln_start_ = aln_start + 1; // 1-based
    aln_end_ = aln_end;  // pysam aln_end is 1-based
    
    /*may be removed later*/
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


    variants = find_mapped_variants( aln_start, aln_end, ref_seq_, read_seq_, cigar_vector_,
                                                 chrom, unspliced_local_reference_start, unspliced_local_reference_end,
                                                 indexed_local_reference );
    
    //std::string variant_str; 
    is_target = false;
    for (auto & v : variants) {
        variant_str += (std::to_string(v.pos_) + "_" + v.ref_ + "_" + v.alt_ + ";");
        if (!is_target) {
            is_target = target.is_equivalent(v, unspliced_local_reference_start, indexed_local_reference);      
        }
    }
}

void pileup::parse_pileup(
    const std::string & chrom,
    int pos,
    const std::string & ref,
    const std::string & alt,
    int unspliced_local_reference_start,
    int unspliced_local_reference_end,
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

    std::vector<pileup::ParsedRead> parsed_reads;

    std::map<int, char> aa = reference_by_position( unspliced_local_reference,
                             unspliced_local_reference_start, unspliced_local_reference_end );

    // target variant
    Variant trgt = Variant(chrom, pos, ref, alt);
    
    
    //Variant v = Variant( "1", 241661227,  "A",  "ATTT");
     //                    unspliced_local_reference_start, unspliced_local_reference_end, aa );
    //Variant vv = Variant( "1", 241661228,  "T",  "TTTT");
      //                    unspliced_local_reference_start, unspliced_local_reference_end, aa );

    //std::cout << v.is_equivalent(vv, unspliced_local_reference_start, aa) << std::endl;
    
    //std::cout << v.get_rightmost_pos(unspliced_local_reference_end, aa ) << " onaji " << vv.get_leftmost_pos( unspliced_local_reference_start, aa) << std::endl;
    //for (auto i : aa) std::cout << i.first << " : " << i.second  << std::endl;

    for ( size_t i = 0; i < pileup_size; ++i ) {
        parsed_reads.emplace_back( unspliced_local_reference_start,
                                 unspliced_local_reference_end,
                                 unspliced_local_reference,
                                 read_names[i],
                                 are_reverse[i],
                                 cigar_strings[i],
                                 aln_starts[i],
                                 aln_ends[i],
                                 read_seqs[i],
                                 ref_seqs[i],
                                 quals[i],
                                 mapqs[i],
                                 trgt,
                                 aa
                           );
 
    }
 
    int h = 0;
    for (auto & read: parsed_reads) {
        if (read.is_target) {
            h += 1;
        }   
    }
    std::cout << h << " num of target read" << std::endl;
    std::vector<pileup::ParsedRead> j = get_gapped_seed_reads(parsed_reads, 6);
    for (auto & read : j) {
        std::cout << read.aln_start_ << "  " << read.aln_end_ << std::endl;
    }
}

// select gapped_seed

std::vector<pileup::ParsedRead> get_gapped_seed_reads(std::vector<pileup::ParsedRead> & parsed_reads, size_t n=5)
{
    
    std::vector<pileup::ParsedRead> target_reads;
    std::vector<pileup::ParsedRead> seed_candidates;
    std::vector<std::string> variant_ptrns;
    
    for (auto & read : parsed_reads) {
        if ( read.is_target ) {
            target_reads.push_back(read);
            variant_ptrns.push_back(read.variant_str);       
        }
    } 
   
    std::string commonest_ptrn = find_commonest_str(variant_ptrns);
    
    
    for (auto & read : target_reads) {
        if ( read.variant_str == commonest_ptrn ){
            seed_candidates.push_back(read);
        }
    }

    if (seed_candidates.size() > n) {
        std::shuffle(seed_candidates.begin(), seed_candidates.end(), std::default_random_engine(123));
        std::vector<pileup::ParsedRead> tmp(seed_candidates.begin(), seed_candidates.begin() + n);
        return tmp;
    }
    
    return seed_candidates;   
}
