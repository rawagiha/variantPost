#include "pileup.h"
#include "match.h"

typedef std::vector<std::pair<std::string_view, size_t>> PatternCnt;

//------------------------------------------------------------------------------
inline auto less = [](const auto& x, const auto& y) { return x.second < y.second; };
inline auto more = [](const auto& x, const auto& y) { return x.second > y.second; };

inline void count_patterns(const Idx& ptrn_read_idx, PatternCnt& ptrn_cnts) {
    for (const auto& elem : ptrn_read_idx)
        ptrn_cnts.emplace_back(elem.first, elem.second.size()); 
    
    std::sort(ptrn_cnts.begin(), ptrn_cnts.end(), more);
}

//------------------------------------------------------------------------------
Pileup::Pileup(const Strs& names, const Bools& is_rv, const Strs& cigar_strs,
               const Ints& aln_starts, const Ints& aln_ends, const Strs& seqs,
               const Strs& quals, const Bools& is_first_bam, const bool has_scnd,
               const UserParams& params, LocalReference& loc_ref, Variant& target) {
     
    has_second = has_scnd; // may not need

    const int qc_start = (target.lpos - params.local_thresh < loc_ref.flanking_start) 
                       ? target.lpos - params.local_thresh : loc_ref.flanking_start;
    const int qc_end = (target.end_pos + params.local_thresh > loc_ref.flanking_end) 
                     ? target.end_pos + params.local_thresh : loc_ref.flanking_end;

    int ref_hap_n = 0; 
    sz = static_cast<int>(names.size()); reads.reserve(sz);
    for (int i = 0; i < sz; ++i) {
        reads.emplace_back(names[i], is_rv[i], cigar_strs[i], aln_starts[i], 
                           aln_ends[i], seqs[i], quals[i], is_first_bam[i]);
        
        auto& read = reads[i];
        read.setReference(loc_ref); if (read.is_na_ref) continue;
        
        read.setVariants(loc_ref); 
        read.parseCoveringPattern(loc_ref, target);
        if (read.covering_ptrn == 'C') continue;

        read.parseLocalPattern(loc_ref, target);
        read.qualityCheck(qc_start, qc_end, params.base_q_thresh, params.lq_rate_thresh);
        read.isStableNonReferenceAlignment(loc_ref); 

        if (read.covering_ptrn == 'A' && !read.idx2pos.empty())
            read.isCenterMapped(target);
        
        read.is_quality_map = (read.is_stable_non_ref && read.is_central_mapped);
         
        // search complex target by string match
        if (target.is_complex && read.covering_ptrn == 'A' && read.is_stable_non_ref) {
            if (target.is_substitute) read.hasTargetComplexSubstitute(target);
            else read.hasTargetComplexIndel(loc_ref, target);
        }
                
        // assign rank
        if (read.has_target) {
            read.rank = 's'; ++s_cnt;
            read.setSignatureStrings(params); sig_s[read.non_ref_sig].push_back(i);
            if (read.is_quality_map) sig_s_hiconf[read.non_ref_sig].push_back(i); 
        } else if (read.has_local_events) {
            read.rank = 'u'; ++u_cnt;
            read.setSignatureStrings(params);
            if (read.is_quality_map && read.qc_passed) sig_u[read.non_ref_sig].push_back(i);
            starts.push_back(read.covering_start); ends.push_back(read.covering_end);
        } else { 
            read.rank = 'n'; ++n_cnt;
            std::cout << read.dist_to_non_target << " " << read.cigar_str << std::endl;
            if (read.is_quality_map) ++ref_hap_n;                
        }       
    }
    
    if (!s_cnt && !u_cnt) { has_no_support = true; return; }
    
    kmer_sz = (loc_ref.low_cplx_len > params.kmer_size)
            ? loc_ref.low_cplx_len : params.kmer_size;

    int _strt = (loc_ref.flanking_start - kmer_sz > loc_ref.start) 
              ? loc_ref.flanking_start - kmer_sz : loc_ref.start;
    starts.push_back(_strt); 
    
    int _end = (loc_ref.flanking_end + kmer_sz < loc_ref.end)
             ? loc_ref.flanking_end + kmer_sz : loc_ref.end;  
    ends.push_back(_end);
    
    start = *std::min_element(starts.begin(), starts.end());
    end = *std::max_element(ends.begin(), ends.end());
        
    has_hiconf_support = (!sig_s_hiconf.empty()); has_ref_hap = (ref_hap_n); 
}

//------------------------------------------------------------------------------
void Pileup::gridSearch(const UserParams& params, 
                        LocalReference& loc_ref, const Variant& target) {
    if (s_cnt || sig_u.empty()) return;

    make_sequence(loc_ref, {}, start, end, rseq, &i2p_r);
    
    PatternCnt u_sig_cnts;
    count_patterns(sig_u, u_sig_cnts);
    
    //NOTE:
    // grid search is applied to reads 
    // 'u'-reads with flanking_start/end covered && central mapped 
    // all reads in sig_u are such quality-mapp reads
    std::string query;
    for (const auto& _sig : u_sig_cnts) {
        const auto& vars = reads[sig_u[_sig.first][0]].variants;
        make_sequence(loc_ref, vars, start, end, query);

        std::vector<Ints> grid; gap_grid(params, grid);
        bool res = search_over_grid(start, loc_ref, rseq, query, grid, target); 
        if (res) {
            const auto& read_idxes = sig_u[_sig.first];
            for (int _i : read_idxes) {
                reads[_i].rank = 's'; ++s_cnt; --u_cnt;
                sig_s[_sig.first].push_back(_i);
                sig_s_hiconf[_sig.first].push_back(_i);
            }
            std::cout << _sig.first << std::endl;
            sig_u.erase(_sig.first);
       }
        query.clear();
    }

    if (!s_cnt) return; // no change
    if(sig_s_hiconf.size()) has_hiconf_support = true;
     
    for (const auto& elem : sig_u) 
        for (int i : elem.second) { reads[i].rank = 'n'; ++n_cnt; --u_cnt; }
    
    // keep sig_u contents to set up non-supporting haplotypes
}

//------------------------------------------------------------------------------
void Pileup::setHaploTypes(LocalReference& loc_ref, const Variant& target) {
    PatternCnt s_sig_cnt; int idx0 = -1;
    if (!sig_s_hiconf.empty()) {
        count_patterns(sig_s_hiconf, s_sig_cnt);
        hap0 = s_sig_cnt[0].first; idx0 = sig_s_hiconf[hap0][0];
    } else if (!sig_s.empty()) {
        count_patterns(sig_s, s_sig_cnt);
         hap0 = s_sig_cnt[0].first; idx0 = sig_s[hap0][0];
    }
    
    if (idx0 > -1) {  
        make_sequence(loc_ref, reads[idx0].variants, start, end, seq0, &i2p_0);
    } else {
        make_sequence(loc_ref, {target}, start, end, seq0, &i2p_0);
    }
    
    PatternCnt u_sig_cnt; int idx1 = -1, idx2 = -1;
    if (!sig_u.empty()) {
        count_patterns(sig_u, u_sig_cnt);
        hap1 = u_sig_cnt[0].first; idx1 = sig_u[hap1][0];
        make_sequence(loc_ref, reads[idx1].variants, start, end, seq1);
        
        if (u_sig_cnt.size() > 1) {
            hap2 = u_sig_cnt[1].first; idx2 = sig_u[hap2][0];
            make_sequence(loc_ref, reads[idx2].variants, start, end, seq2);
        }
    }
    
    make_sequence(loc_ref, {}, start, end, rseq, &i2p_r);
    
    if (idx0 < 0 || (idx1 == -1 && idx2 == -1)) vs_ref_hap = true;    
}


//------------------------------------------------------------------------------
void Pileup::differentialKmerAnalysis(const UserParams& params, 
                                      LocalReference& loc_ref, const Variant& target) {
    
    Kmers km0, km1, km2, kmr, km12, km12r;   
    
    // target kmers
    make_kmers(seq0, kmer_sz, km0); 
    if (vs_ref_hap) {
        make_kmers(rseq, kmer_sz, km12r);
    } else {
        make_kmers(seq1, kmer_sz, km1); 
        make_kmers(seq2, kmer_sz, km2);
        if (has_ref_hap) make_kmers(rseq, kmer_sz, kmr);
        std::set_union(km1.begin(), km1.end(), km2.begin(), km2.end(),
                       std::inserter(km12, km12.begin()));
        std::set_union(km12.begin(), km12.end(), kmr.begin(), kmr.end(),
                       std::inserter(km12r, km12r.begin()));
    }
    
    std::set_difference(km0.begin(), km0.end(), km12r.begin(), km12r.end(),
                        std::inserter(kmers_t, kmers_t.end()));
    std::set_difference(km12r.begin(), km12r.end(), km0.begin(), km0.end(), 
                        std::inserter(kmers_nt, kmers_nt.end()));
    
    for (auto& read : reads) {
        if (read.rank != 'u') continue;
        
        for (const auto& kmer : kmers_nt) 
            if (read.seq.find(kmer) != std::string_view::npos) ++(read.nmer);
        for (const auto& kmer : kmers_t) 
            if (read.seq.find(kmer) != std::string_view::npos) ++(read.smer);
        
        if (read.smer && !read.nmer) { 
            if (has_hiconf_support || read.is_quality_map) { 
                read.rank = 's'; --u_cnt; ++s_cnt; 
            } 
            else { read.rank = 'y'; --u_cnt; ++y_cnt; } // likel'y' supporting 
        }
        // not do for complex cases? 
        if (read.nmer && !read.smer) { read.rank = 'n'; --u_cnt; ++n_cnt; }
    }
    has_likely_support = (y_cnt); 
}
