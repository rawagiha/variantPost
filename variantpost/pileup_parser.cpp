#include <map>
#include <random>
#include <vector>
#include <string>
#include <utility>
#include <iostream>

#include "util.h"
#include "swlib.h"
#include "pileup_parser.h"



/* read filter */


std::vector<pileup::ParsedRead> get_gapped_seed_reads(std::vector<pileup::ParsedRead> & parsed_reads, size_t n);


//check covering

// read classifer
char classify_read(const int aln_start,
                   const int aln_end,
                   const std::string & cigar_string,
                   const std::vector<std::pair<char, int>> & cigar_vec,
                   //const Variant & variant, 
                   const bool is_ins,
                   const bool is_del,
                   const int l_pos, 
                   const int r_pos,
                   const std::vector<Variant> & variants);

// read evaluate func (ref_seq, clipping, worth_realn)
// covering checker



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
    const std::string chrom = target.chrom;
    read_name_ = read_name;
    is_reverse_ = is_reverse;

    cigar_string_ = cigar_string;   // may be removed later 
    cigar_vector_ = to_cigar_vector( cigar_string ); // may be removed later 

    aln_start_ = aln_start + 1; // 1-based
    aln_end_ = aln_end;  // pysam aln_end is 1-based
    
    //std::pair<char, int> first_c, last_c = cigar_vector_[0], cigar_vector_.back();
    std::pair<char, int> first_c = cigar_vector_[0];
    std::pair<char, int> last_c = cigar_vector_.back();
    int start_offset = ( first_c.first == 'S' ) ? first_c.second : 0;
    int read_start = aln_start_ - start_offset;
    int end_offset = ( last_c.first == 'S' ) ? last_c.second : 0;   
    int read_end = aln_end + end_offset;
    
    
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

    std::vector<std::pair<int, int>> exons;
    std::vector<std::pair<int, int>> introns;
    parse_splice_pattern(exons, introns, cigar_vector_, read_start, read_end);
        
        /* 
        for (auto e : exons ) {
            std::cout << "(" << e.first << ", " << e.second << "), " ;
        }
        std::cout << std::endl;
        for (auto i : introns ) {
            std::cout << "(" << i.first << ", " << i.second << "), ";
        }
        std::cout << std::endl;
        */
    
    //remove check for empty str
    if ( !ref_seq_.empty() ) {
        variants = find_mapped_variants( aln_start, aln_end, ref_seq_, read_seq_, cigar_vector_,
                                         chrom, unspliced_local_reference_start, unspliced_local_reference_end,
                                         indexed_local_reference );
    }
    // read evaluations
    // test if clipped
    
    bool is_soft = false, is_hard = false;
    if ( cigar_string.find( 'S' ) == std::string::npos ) {
        is_soft = true;
    }

    if ( cigar_string.find( 'H' ) == std::string::npos ) {
        is_hard = true;
    }
    
    is_clipped = false;
    if ( is_soft || is_hard ) {
        is_clipped = true;
    }


    // test if read seq == reference
    is_ref_seq = false;
    if ( variants.empty() && ( !is_clipped ) ) {
        is_ref_seq = true;
    }
    
    if ( 1) {
    
        // test if target is aligned
        is_target = false;
        for (auto & v : variants) {
            
            // do for spl ptrn
            variant_str += (std::to_string(v.pos) + "_" + v.ref + "_" + v.alt + ";");
            if (!is_target) {
                is_target = target.is_equivalent(v, unspliced_local_reference_start, indexed_local_reference);      
            }
        }

        //
    
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
            //std::cout << read.read_name_ << " " << read.is_reverse_ << std::endl;
        } 
        else {
            //std::cout << read.read_name_ << " " << read.is_reverse_ << " " << read.cigar_string_ << std::endl;
        }  
    }
    std::cout << h << " num of target read" << std::endl;
    
    if ( h > 0 ) {
    std::vector<pileup::ParsedRead> j = get_gapped_seed_reads(parsed_reads, 6);
    std::vector<std::string> seed_read = {j[0].read_seq_, j[0].base_qualities_};
    
    std::cout << j[0].read_seq_ << "  " << j[0].base_qualities_ << std::endl;
    std::vector<std::vector<std::string>> _reads;
    for ( size_t i = 1; i < 6; ++i) {
        std::cout << j[i].read_seq_ << "  " << j[i].base_qualities_ << std::endl;
        std::vector<std::string> r = {j[i].read_seq_, j[i].base_qualities_};
        _reads.push_back(r);
    }
    
    //std::vector<std::string> lll = {"TAAATTTATGTAAATCACTTTGGACCCAGCATGTCCTTAGGTTTTACCCATTCGTCAAACTGCTCTGCTGTGAGATAGCCAAGTTCGATAGCAGTTTCCTTTAAGGTTGATCCATTTTTTTTGTGTGCTGTCTTAGCAATCTTTGCTGCCTTGTCATACCCTGAAGAAAAAATAAAAAGACGACATATGGGTTAGCAGTGA", "@E@JJJJJJJJJJJEHIIJJJJJJJJJJJJJJJJJJJJJJJJJJGHCECGEGJJJJJJIJJJIIJJJIJHHJJJJJJJJJJJHJJJJJJJJJJJJIJJJJHJJJJJJJJJJJJJJJJJJJJJIJIJIIJIHJJJJGHIJJJJJJHEIHJIJJJJIJJJJHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJE@"};
    //std::vector<std::vector<std::string>> _pp = {{"TTTCATTATAAATTTATGTAAATCACTTTGGACCCAGCATGTCCTTAGGTTTTACCCATTCGTCAAACTGCTCTGCTGTGAGATAGCCAAGTTCGATAGCAGTTTCCTTTAAGGTTGATCCATTTTTTTTGTGTGCTGTCT", "KFK<7KKKKFAKKKKKKKKKKKKF<KKKKKFAKKKFKKKKKKKKKKKKKKFKKFFKKFKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKA"}};
    
    sw::flatten_reads(seed_read, _reads);
    }
   
   // for (auto & read : j) {
   //     std::cout << read.aln_start_ << "  " << read.aln_end_ << std::endl;
   // }
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
