#include "util.h"
#include "match.h"
#include "search.h"
#include "pileup.h"
#include "consensus.h"
#include "variant_types.h"


// pasted from reads.cpp
inline auto append_num = [](std::string& target, auto val) {
    std::array<char, 20> buffer; 
    auto [ptr, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), val);
    if (ec == std::errc{}) {
        target.append(buffer.data(), ptr - buffer.data());
    }
};


SearchResult::SearchResult() {}

void SearchResult::setReadInfo(const Read& read, const char qual_thresh) {
    
    are_from_first_bam.push_back(!read.is_control);
    
    if (read.rank == Rank::Supporting) {
        target_statuses.push_back(1);
    } else if (read.rank == Rank::NotSupporting && !read.covered_in_clip) {
        target_statuses.push_back(0);
        //if (!read.is_quality_map || !read.qc_passed || !read.has_local_clip) return;
        
        // Collect variants from high-qual negatives
        for (const auto& v : read.variants) {
            if (v.mean_qual < qual_thresh) continue;
            
            hq_negative_cnts[VariantKey{v.pos, v.end_pos, v.ref, v.alt}]++;
        }    
    } else if (read.rank == Rank::Undetermined) {
        target_statuses.push_back(-1);
    } else {
        target_statuses.push_back(-2);
    }
} 

void SearchResult::finalize() {
    for (const auto& [key, cnt] : hq_negative_cnts) {
        std::cout << key.pos << " " << key.ref << " " << key.alt << " " << cnt << " " << std::endl;
        if (cnt < 2) continue;
        std::string v_str;
        v_str.reserve(256);
        append_num(v_str, key.pos);
        v_str.append("_").append(key.ref).append("_").append(key.alt);
        trans_vars.push_back(v_str);
        std::cout << trans_vars.size() << std::endl;
    }

}

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
    std::cout << loc_ref.flanking_start << " " << loc_ref.flanking_end << std::endl;
    target.setFlankingSequences(loc_ref); 
    std::cout << target.lt_seq << " has it >> " << target.rt_seq << std::endl;
    
    
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
        pileup.inferGermlineHaplotype(params, target.pos); 
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
            //if (read.covering_ptrn == 'C') continue;
            
            rslt.setReadInfo(read, params.base_q_thresh);
        }
    } else {
        for (const auto& read : pileup.reads) {
            //if (read.covering_ptrn == 'C') continue;
            
            rslt.setReadInfo(read, params.base_q_thresh);
        }
    }

    for (const auto& read : pileup.reads) {
        if (read.rank == Rank::Supporting) std::cout << read.name << " " << read.cigar_str <<  "  " << read.aln_start << std::endl;
    }
    
    // Realn against personalized genome
    if (pileup.has_second_bam && !target.is_substitute)
        personalize(pileup, loc_ref, params, target, rslt);    
    
    if (pileup.has_hiconf_support) {
        from_consensus_variant_list(rslt, loc_ref, pileup.start, pileup.end, pileup.hap0_vars);
    }
    
    rslt.finalize(); 
}   
