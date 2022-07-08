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
char classify_covering(const int lpos, 
                       const int pos,
                       const int rpos,
                       const bool is_shiftable,
                       const std::vector<std::pair<int, int>> & exons,
                       const std::vector<Variant> & variants);
 


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



pileup::ParsedRead::ParsedRead( const int & unspliced_local_reference_start,
                                const int & unspliced_local_reference_end,
                                const std::string & unspliced_local_reference,
                                const std::string & read_name,
                                const bool & is_reverse,
                                const std::string & cigar_string,
                                const int & aln_start,
                                const int & aln_end,
                                const std::string & read_seq,
                                const std::string & ref_seq_default,
                                const std::vector<int> & q,
                                const int & mapq,
                                const Variant & target, 
                                const int & lpos,
                                const int & pos,
                                const int & rpos,
                                const bool & is_shiftable,
                                const std::map<int, char> & indexed_local_reference
) : read_name(read_name), is_reverse(is_reverse), cigar_string(cigar_string), 
    aln_start(aln_start),  aln_end(aln_end), read_seq(read_seq), mapq(mapq) 
{
    //const std::string chrom = target.chrom;

    cigar_vector = to_cigar_vector( cigar_string ); // may be removed later 

    std::pair<char, int> first_cigar = cigar_vector[0];
    std::pair<char, int> last_cigar = cigar_vector.back();
    int start_offset = ( first_cigar.first == 'S' ) ? first_cigar.second : 0;
    int read_start = aln_start - start_offset;
    int end_offset = ( last_cigar.first == 'S' ) ? last_cigar.second : 0;   
    int read_end = aln_end + end_offset;
    
    
    /*may be removed later*/
    if ( cigar_string.find( 'N' ) != std::string::npos ) {
        ref_seq = ref_seq_default;
        is_spliced = true;
    }
    else {
        ref_seq = get_read_wise_ref_seq( aln_start, aln_end,
                                          unspliced_local_reference_start, unspliced_local_reference );
        is_spliced = false;
    }

    base_qualities = to_fastq_qual( q );

    std::vector<std::pair<int, int>> exons;
    std::vector<std::pair<int, int>> introns;
    parse_splice_pattern(exons, introns, cigar_vector, read_start, read_end);
        
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
    if ( !ref_seq.empty() ) {
        variants = find_mapped_variants( aln_start, aln_end, ref_seq, read_seq, cigar_vector,
                                         unspliced_local_reference_start, unspliced_local_reference_end,
                                         indexed_local_reference );
    }
    
    
    char cov = classify_covering(lpos, pos, rpos, is_shiftable, exons, variants);
    //std::cout << read_name << "  " << cov << " " << is_reverse << std::endl;
    
    // read evaluations
    // test if clipped
    
    if (read_name == "HISEQ:45:C2MKUACXX:1:1105:3122:96546") {
        std::cout << "yes!!" << std::endl;
        std::cout << variants.empty() << std::endl;
        for (auto & v : variants) {
            std::cout << v.pos << " " << v.ref << " " << v.alt << std::endl;
        }
    }
    
    bool is_soft = ( cigar_string.find( 'S' ) != std::string::npos ) ? true : false;
    bool is_hard = ( cigar_string.find( 'S' ) != std::string::npos ) ? true : false;
    bool is_clipped = ( is_soft || is_hard ) ? true : false;
    

    // test if read seq == reference
    is_ref_seq = (ref_seq == read_seq) ? true : false;
    
    if ( 1) {
    
        // test if target is aligned
        is_target = false;
        for (auto & v : variants) {
            
            // do for spl ptrn and clipping(pos only) -> rename to event_str
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
    int pos, //lpos?
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
    Variant trgt = Variant(pos, ref, alt);
    int lpos = trgt.get_leftmost_pos(unspliced_local_reference_start, aa);
    int rpos = trgt.get_rightmost_pos(unspliced_local_reference_end, aa);
    bool is_shiftable = (lpos != rpos) ? true : false;
    
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
                                 lpos,
                                 pos,
                                 rpos,
                                 is_shiftable,
                                 aa
                           );
 
    }
 
    int h = 0;
    for (auto & read: parsed_reads) {
        if (read.is_target) {
            h += 1;
            //std::cout << read.read_name << " " << read.is_reverse << std::endl;
        } 
        else {
            //std::cout << read.read_name_ << " " << read.is_reverse_ << " " << read.cigar_string_ << std::endl;
        }  
    }
    std::cout << h << " num of target read" << std::endl;
    
    if ( h > 0 ) {
    std::vector<pileup::ParsedRead> j = get_gapped_seed_reads(parsed_reads, 6);
    std::vector<std::string> seed_read = {j[0].read_seq, j[0].base_qualities};
    
    std::cout << j[0].read_seq << "  " << j[0].base_qualities << std::endl;
    std::vector<std::vector<std::string>> _reads;
    for ( size_t i = 1; i < 6; ++i) {
        std::cout << j[i].read_seq << "  " << j[i].base_qualities << std::endl;
        std::vector<std::string> r = {j[i].read_seq, j[i].base_qualities};
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

// check how the read covers the locus of interest
// 'A' : valid covering
// 'B' : possibly (may be a part of complex event) 
// 'C' : undetermined or not covering 
//-------------------------------------------------
char classify_covering(const int lpos, 
                       const int pos,
                       const int rpos, 
                       const bool is_shiftable,
                       const std::vector<std::pair<int, int>> & unspl_segments,
                       const std::vector<Variant> & variants)
{
    for (auto & segment : unspl_segments) {
        int start = segment.first;
        int stop = segment.second;  
        if (!is_shiftable) {
            if (ordered(start, pos, stop)) return 'A'; // simply covering
        }
        else {
            // repeats bounded 
            if (start <= lpos && rpos + 1 <= stop) return 'A'; // covering & repeat bounded
            // repeats left open  
            else if (ordered(lpos + 1, start, rpos)) {
                
                // variants exist between start and rpos
                for (auto & variant : variants) {
                    if (ordered(start, variant.pos, rpos)) {
                        return 'B'; // repeats broken
                    }
                }
                return 'L';    
            }           
            // repeats right open
            else if (ordered(lpos, stop, rpos)) {
                for (auto & variant : variants) {
                    if (ordered(lpos, variant.pos, stop)) {
                        return 'E'; // repeats broken
                    }
                }
                return 'R'; 
            }
        }
    }
    return 'D'; //uncovering  
}
