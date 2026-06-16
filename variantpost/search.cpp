#include "util.h"
#include "reads.h"
#include "match.h"
#include "search.h"
#include "pileup.h"
#include "consensus.h"
#include "variant_types.h"

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
    if (!loc_ref.has_flankings) { std::cout << "no flanking! " << std::endl; return; }
    
    // Annotation for low complex tags
    loc_ref.findLowComplexRegion();
    target.setFlankingSequences(loc_ref); 
    target.countRepeats(loc_ref);
    target.testForDeNovoRepeats(loc_ref);
    if (loc_ref.locplx_start <= target.pos
        && target.end_pos <= loc_ref.locplx_end)
    { target.in_homopolymer = true; }
    
    
    // pileup setup
    // 1. reads are check for target based on the original alignment
    // 2. non-target/non-ref variations located in read-middle, surrounded 
    //    with complex sequence (signature-u) are flagged for grid-search/haplotype comparison
    Pileup pileup(read_names, are_reverse, cigar_strs, 
                  aln_starts, aln_ends, read_seqs, quals, 
                  are_from_first_bam, has_second, 
                  params, loc_ref, target);
    
    if (pileup.has_second_bam) {
        pileup.inferGermlineHaplotype(params); 
        make_sequence(
            loc_ref, pileup.homo_vars, pileup.start, pileup.end, pileup.rseq, &pileup.i2pr);
    } else {
        make_sequence(loc_ref, {}, pileup.start, pileup.end, pileup.rseq, &pileup.i2pr);
    }
    
    pileup.gridSearch(params, loc_ref, target);
    if (pileup.u_cnt) {
        pileup.setHaploTypes(loc_ref, target);
        pileup.differentialKmerAnalysis(params, loc_ref, target);
        pileup.searchByRealignment(params, loc_ref, target);
        
        if (pileup.has_hiconf_support)
             match2haplotypes(pileup, read_seqs, params); 
        for (const auto& read : pileup.reads) {
            if (read.covering_ptrn == 'C') continue;
            bool is_first = (read.is_control) ? false : true;

            rslt.are_from_first_bam.push_back(is_first);
            if (read.rank == 's' ) { rslt.target_statuses.push_back(1); }
            else if (read.rank == 'n' && !read.covered_in_clip) { 
                rslt.target_statuses.push_back(0); 
            }
            //else if (read.rank == 'u' || read.rank == 'y') rslt.target_statuses.push_back(-1);
            else if (read.rank == 'y') rslt.target_statuses.push_back(-1);
            else rslt.target_statuses.push_back(-2);
        }
    } else {
        for (const auto& read : pileup.reads) {
            if (read.covering_ptrn == 'C') continue;
            bool is_first = (read.is_control) ? false : true;
            rslt.are_from_first_bam.push_back(is_first);
            if (read.rank == 's') rslt.target_statuses.push_back(1);
            else if (read.rank == 'n' && !read.covered_in_clip) rslt.target_statuses.push_back(0);
            else if (read.rank == 'u' || read.rank == 'y') rslt.target_statuses.push_back(-1);
            else rslt.target_statuses.push_back(-2);
        }
    }

    for (const auto& read : pileup.reads) {
        if (read.rank == 's') std::cout << read.name << " " << read.cigar_str <<  "  " << read.aln_start << std::endl;
    }

    
    // Realn against personalized genome
    Variant per(target.pos, target.ref, target.alt);
    personalize(pileup, loc_ref, params, target, per, rslt.pltseq, rslt.prtseq);    
    
    
    
    //if (per != target) {
    rslt.ppos = per.pos;
    rslt.pref = per.ref;
    rslt.palt = per.alt;
    //}    
    
    Consensus con;
    if (pileup.hiconf_read_idx > -1) {
        con._from_variants(pileup.start, pileup.end,  pileup.reads[pileup.hiconf_read_idx].variants, params, loc_ref);
        size_t kk = con.ref.size();
        rslt.positions = con.pos;
        rslt.ref_bases = con.ref;
        rslt.alt_bases = con.alt;
    } 
}   
