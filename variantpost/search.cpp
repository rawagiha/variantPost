#include "util.h"
#include "reads.h"
#include "search.h"
#include "pileup.h"
#include "sequence_model.h"

SearchResult::SearchResult() {}

void _search_target(SearchResult& rslt,
                    
                    const std::string& fastafile,
    
                    const std::string& chrom, const int pos, 
                    const std::string& ref, const std::string& alt,
    
                    const int mapq_thresh, const int base_q_thresh,
                    const double lq_base_rate_thresh,
                    const int match_score, const int mismatch_penal,
                    const int gap_open_penal, const int gap_ext_penal,
                    const int kmer_size, const int local_thresh,
                    const int retarget_thresh, const int ref_start, const int ref_end,
    
                    const Strs& read_names,
                    const Bools& are_reverse,
                    const Strs& cigar_strs,
                    const Ints& aln_starts,
                    const Ints& aln_ends,
                    const Strs& read_seqs,
                    const Strs& quals,
                    const Ints& mapqs,                            //mapqs to be removed 
                    const Bools& are_from_first_bam)
{  
    // basic prep
    UserParams params(mapq_thresh, base_q_thresh, lq_base_rate_thresh,
                      match_score, mismatch_penal, gap_open_penal, gap_ext_penal, 
                      kmer_size, local_thresh, retarget_thresh);
    LocalReference loc_ref(fastafile, chrom, ref_start, ref_end);   
    
    /*should terminate if fail to set flanking*/ 
    
    Variant target(pos, ref, alt); target.setEndPos(loc_ref);
    
    // additional prep 
    loc_ref.setFlankingBoundary(target, params.min_dimer_cnt);
    target.setFlankingSequences(loc_ref); // for complex indel/MNV
    
    // pileup setup
    Pileup pileup(read_names, are_reverse, cigar_strs, 
                  aln_starts, aln_ends, read_seqs, quals, 
                  are_from_first_bam, params, loc_ref, target);
    
    if (pileup.has_hiconf_support) {
        if (pileup.u_cnt) {
            pileup.setHaploTypeByFrequency(); 
            pileup.setSequenceFromHaplotype(loc_ref);
            pileup.reRankByKmer(params, loc_ref);
            
            for (const auto& read : pileup.reads) {
                if (read.rank == 's') rslt.target_statuses.push_back(1);
                else if (read.rank == 'n') rslt.target_statuses.push_back(0);
                else if (read.rank == 'u') rslt.target_statuses.push_back(-1);
                else rslt.target_statuses.push_back(-2);
            }
        }
        else {
            for (const auto& read : pileup.reads) {
                if (read.rank == 's') rslt.target_statuses.push_back(1);
                else if (read.rank == 'n') rslt.target_statuses.push_back(0);
                else if (read.rank == 'u') rslt.target_statuses.push_back(-1);
                else rslt.target_statuses.push_back(-2);
            }    
        }
    }
    else if (pileup.has_no_support) {
        //no result case
    }
    else {
        SequenceModel seqm(pileup, loc_ref, target);        
        seqm.compareToRefByKmer(pileup, loc_ref, params);
        seqm.reRankByReAlignment(pileup, read_seqs, params);
        //std::cout << seqm.flank_start << " " << seqm.target_start << " " << seqm.target_end << " " << seqm.flank_end << std::endl;
        
        //pileup.compareToRefByKmer(loc_ref, params, target);
        for (const auto& read : pileup.reads) {
            if (read.rank == 's') rslt.target_statuses.push_back(1);
            else if (read.rank == 'n') rslt.target_statuses.push_back(0);
            else if (read.rank == 'u') rslt.target_statuses.push_back(-1);
            else rslt.target_statuses.push_back(-2);
        }    
    }
}   
