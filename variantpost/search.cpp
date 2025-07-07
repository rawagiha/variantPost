#include "util.h"
#include "reads.h"
#include "search.h"
#include "pileup.h"

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
    
                    const std::vector<std::string>& read_names,
                    const std::vector<bool>& are_reverse,
                    const std::vector<std::string>& cigar_strs,
                    const std::vector<int>& aln_starts,
                    const std::vector<int>& aln_ends,
                    const std::vector<std::string>& read_seqs,
                    const std::vector<std::string>& quals,
                    const std::vector<int>& mapqs,                            //mapqs to be removed 
                    const std::vector<bool>& are_from_first_bam)
{  
    // basic prep
    UserParams params(mapq_thresh, base_q_thresh, lq_base_rate_thresh,
                      match_score, mismatch_penal, gap_open_penal, gap_ext_penal, 
                      kmer_size, local_thresh, retarget_thresh);
    LocalReference loc_ref(fastafile, chrom, ref_start, ref_end);   
    Variant target(pos, ref, alt); target.setEndPos(loc_ref);
    
    // additional prep 
    loc_ref.setFlankingBoundary(target, params.min_dimer_cnt);
    target.setFlankingSequences(loc_ref); // for complex indel/MNV
    
    // pileup setup
    Pileup pileup(read_names, are_reverse, cigar_strs, 
                  aln_starts, aln_ends, read_seqs, quals, 
                  are_from_first_bam, params, loc_ref, target);
    
    if (pileup.s_cnt) {
        std::cout << pileup.s_cnt << " " << pileup.n_cnt << " " << pileup.u_cnt << std::endl;
        if (pileup.u_cnt) {
            pileup.setHaploTypeByFrequency(); 
            pileup.setSequenceFromHaplotype(loc_ref);
            pileup.differentialKmerAnalysis(params, loc_ref);
            //std::cout << pileup.seq0 << std::endl;
            //std::cout << pileup.seq1 << std::endl;
            //std::cout << pileup.seq2 << std::endl;
            //
            //std::cout << pileup.s_cnt << " " << pileup.n_cnt << " " << pileup.u_cnt << std::endl;
            for (const auto& read : pileup.reads) {
                if (read.rank == 's') rslt.target_statuses.push_back(1);
                else if (read.rank == 'n') rslt.target_statuses.push_back(0);
                else if (read.rank == 'u') rslt.target_statuses.push_back(-1);
                else rslt.target_statuses.push_back(-2);
            }
        }
        else {
        // best case
            for (const auto& read : pileup.reads) {
                if (read.rank == 's') rslt.target_statuses.push_back(1);
                else if (read.rank == 'n') rslt.target_statuses.push_back(0);
                else if (read.rank == 'u') rslt.target_statuses.push_back(-1);
                else rslt.target_statuses.push_back(-2);
            }    
        }
    }
    else if (pileup.u_cnt) {
        //kmer search -> asm
    }
    else {
        //no target no undetermined -> report resuls
    }
}   
