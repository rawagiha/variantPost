#include "util.h"
#include "reads.h"
#include "match.h"
#include "search.h"
#include "pileup.h"
#include "sequence_model.h"

SearchResult::SearchResult() {}

void _search_target(SearchResult& rslt,
                    
                    const std::string& fastafile,
    
                    const std::string& chrom, const int pos, 
                    const std::string& ref, const std::string& alt,
    
                    //const int mapq_thresh, 
                    const int base_q_thresh,
                    const double lq_base_rate_thresh,
                    const int match_score, const int mismatch_penal,
                    const int gap_open_penal, const int gap_ext_penal,
                    const int kmer_size, const int dimer_window,
                    const int local_thresh, const int ref_start, const int ref_end,
    
                    const Strs& read_names,
                    const Bools& are_reverse,
                    const Strs& cigar_strs,
                    const Ints& aln_starts,
                    const Ints& aln_ends,
                    const Strs& read_seqs,
                    const Strs& quals,
                    //const Ints& mapqs,                            //mapqs to be removed 
                    const Bools& are_from_first_bam,
                    const bool has_second)
{  
    // basic prep
    UserParams params(base_q_thresh, lq_base_rate_thresh,
                      match_score, mismatch_penal, gap_open_penal, gap_ext_penal, 
                      kmer_size, dimer_window, local_thresh);
    
    LocalReference loc_ref(fastafile, chrom, ref_start, ref_end);   
    
    Variant target(pos, ref, alt); target.setEndPos(loc_ref);
    
    // additional prep 
    loc_ref.setFlankingBoundary(target, params.dimer_window);
    // TODO how to return the result in this case??
    if (!loc_ref.has_flankings) return;
    
    // for repetitive indels and complex indel/MNV
    target.setFlankingSequences(loc_ref); 
    target.countRepeats(loc_ref);
    std::cout << "repeats " << target.repeats << std::endl;

    // pileup setup
    Pileup pileup(read_names, are_reverse, cigar_strs, 
                  aln_starts, aln_ends, read_seqs, quals, 
                  are_from_first_bam, has_second, 
                  params, loc_ref, target);
    
    pileup.gridSearch(params, loc_ref, target);
    if (pileup.u_cnt) {
        std::cout << "s cnt pre 0 "  << pileup.s_cnt << std::endl;
        pileup.setHaploTypes(loc_ref, target);
        
        std::cout << "s cnt pre 1 "  << pileup.s_cnt << std::endl;
        
        // less effective for repeats with long unit
        pileup.differentialKmerAnalysis(params, loc_ref, target);
        
        std::cout << "s cnt pre "  << pileup.s_cnt << std::endl;

        pileup.searchByRealignment(params, loc_ref, target);
  
        std::cout << "s cnt "  << pileup.s_cnt << "  " <<  pileup.n_cnt << std::endl;
        
        if (pileup.has_hiconf_support)
             match2haplotypes(pileup, read_seqs, params); 
        
    

        for (const auto& read : pileup.reads) {
            //std::cout << read.cigar_str << " " << read.smer << " " << read.nmer << " " << read.dist_to_non_target << std::endl;
            if (read.rank == 's') rslt.target_statuses.push_back(1);
            else if (read.rank == 'n' && !read.covered_in_clip) { 
                rslt.target_statuses.push_back(0); 
                //std::cout << read.name << " " << read.cigar_str << std::endl; 
                }
            else if (read.rank == 'u' || read.rank == 'y') rslt.target_statuses.push_back(-1);
            else rslt.target_statuses.push_back(-2);
        }
    } else {
    
       // std::cout << "here no u " << std::endl; 

        for (const auto& read : pileup.reads) {
            //std::cout << read.cigar_str << " " << read.smer << " " << read.nmer << " " << read.dist_to_non_target << std::endl;
            if (read.rank == 's') rslt.target_statuses.push_back(1);
            else if (read.rank == 'n' && !read.covered_in_clip) rslt.target_statuses.push_back(0);
            else if (read.rank == 'u' || read.rank == 'y') rslt.target_statuses.push_back(-1);
            else rslt.target_statuses.push_back(-2);
        }
    }
    
    
    /*
    if (pileup.has_hiconf_support) {
        if (pileup.u_cnt) {
            // up to 2nd most frequent haplotypes
            pileup.setHaploTypeByFrequency(); 
            pileup.setSequenceFromHaplotype(loc_ref);
            
            // w/ kmer specific to target hap & w/o specific to alt haps
            // -> 's'
            // w/ kmer specific to alt haps & w/o specific to target hap
            // -> 'n'
            pileup.reRankByKmer(params, loc_ref);
            
            //Here, the remaining reads are not resolved by kmer analysis
            
            // competitive alignment agains to all haplotypes
            // uniq best aln to target hap -> 's'
            // uniq best aln to (one of )alt haps -> 'n'
            // tie exists -> 'u' 
            match2haplotypes(pileup, read_seqs, params);
            
            for (const auto& read : pileup.reads) {
                if (read.rank == 's') rslt.target_statuses.push_back(1);
                else if (read.rank == 'n') rslt.target_statuses.push_back(0);
                else if (read.rank == 'u') rslt.target_statuses.push_back(-1);
                else rslt.target_statuses.push_back(-2);
            }
        } else {
            // best case 
            // just make a contig & report
            for (const auto& read : pileup.reads) {
                if (read.rank == 's') rslt.target_statuses.push_back(1);
                else if (read.rank == 'n') rslt.target_statuses.push_back(0);
                else if (read.rank == 'u') rslt.target_statuses.push_back(-1);
                else rslt.target_statuses.push_back(-2);
            }    
        }
    } else if (pileup.has_no_support) {
        //no result case
    } else {
        if (pileup.s_cnt) {
            
        }
        
        if (pileup.u_cnt) {
            SequenceModel seqm(pileup, loc_ref, target);        
            seqm.compareToRefByKmer(pileup, loc_ref, params);
            seqm.reRankByReAlignment(pileup, read_seqs, params);
        } else {
            //std::cout << "aho" << std::endl;        

        }
        //std::cout << seqm.flank_start << " " << seqm.target_start << " " << seqm.target_end << " " << seqm.flank_end << std::endl;
        
        //pileup.compareToRefByKmer(loc_ref, params, target);
        for (const auto& read : pileup.reads) {
            if (read.rank == 's') rslt.target_statuses.push_back(1);
            else if (read.rank == 'n') rslt.target_statuses.push_back(0);
            else if (read.rank == 'u') rslt.target_statuses.push_back(-1);
            else rslt.target_statuses.push_back(-2);
        }    
    }*/
}   
