#include "util.h"
#include "reads.h"
#include "match.h"
#include "search.h"
#include "pileup.h"
#include "consensus.h"
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
    
    // pileup setup
    Pileup pileup(read_names, are_reverse, cigar_strs, 
                  aln_starts, aln_ends, read_seqs, quals, 
                  are_from_first_bam, has_second, 
                  params, loc_ref, target);
    if (!pileup.sig_s_hiconf.empty()) {
        Consensus con;
        for (const auto& pair : pileup.sig_s_hiconf) {
            std::cout << pair.first << " " << pair.second.size() << std::endl;
            for (const auto i : pair.second) {
                con._from_variants(pileup.reads[i].covering_start, pileup.reads[i].covering_end,  pileup.reads[i].variants, params, loc_ref);
            }
        }
        size_t kk = con.ref.size();
        for (size_t i = 0; i < kk; ++i){
            std::cout << con.pos[i] << " " << con.ref[i] << " " << con.alt[i] << std::endl;
        }
    }


    pileup.gridSearch(params, loc_ref, target);
    if (pileup.u_cnt) {
        pileup.setHaploTypes(loc_ref, target);
        pileup.differentialKmerAnalysis(params, loc_ref, target);
        pileup.searchByRealignment(params, loc_ref, target);
        
        if (pileup.has_hiconf_support)
             match2haplotypes(pileup, read_seqs, params); 
        
        for (const auto& read : pileup.reads) {
            if (read.rank == 's') rslt.target_statuses.push_back(1);
            else if (read.rank == 'n' && !read.covered_in_clip) { 
                rslt.target_statuses.push_back(0); 
            }
            else if (read.rank == 'u' || read.rank == 'y') rslt.target_statuses.push_back(-1);
            else rslt.target_statuses.push_back(-2);
        }
    } else {
        for (const auto& read : pileup.reads) {
            if (read.rank == 's') rslt.target_statuses.push_back(1);
            else if (read.rank == 'n' && !read.covered_in_clip) rslt.target_statuses.push_back(0);
            else if (read.rank == 'u' || read.rank == 'y') rslt.target_statuses.push_back(-1);
            else rslt.target_statuses.push_back(-2);
        }
    }

   
}   
