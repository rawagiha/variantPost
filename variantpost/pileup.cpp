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
             
    kmer_sz = (loc_ref.low_cplx_len > params.kmer_size)
            ? loc_ref.low_cplx_len : params.kmer_size;
    
    const int qc_start = (target.lpos - params.local_thresh < loc_ref.flanking_start) 
                       ? target.lpos - params.local_thresh : loc_ref.flanking_start;
    const int qc_end = (target.end_pos + params.local_thresh > loc_ref.flanking_end) 
                     ? target.end_pos + params.local_thresh : loc_ref.flanking_end;

    int ref_hap_n = 0, ineff_kmer = 0, overlapping = 0; 
    sz = static_cast<int>(names.size()); reads.reserve(sz);
    for (int i = 0; i < sz; ++i) {
        reads.emplace_back(names[i], is_rv[i], cigar_strs[i], aln_starts[i], 
                           aln_ends[i], seqs[i], quals[i], is_first_bam[i]);
         
        auto& read = reads[i];
        read.setReference(loc_ref); if (read.is_na_ref) continue;
        
        read.setVariants(loc_ref); 
        read.parseCoveringPattern(loc_ref, target);
        if (read.covering_ptrn == 'C') continue;

        read.parseLocalPattern(loc_ref, target, kmer_sz);
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
            read.checkByRepeatCount(target, has_excess_ins_hap);
            
            // sigle, but non-target mutation closed by complex seq (anti-pattern)
            // example:
            //      target    del(C)
            //                    *
            //        ref GCTAAAAACAAAATGC
            //      read  GCTAAAAAAAAAATGC  de-novo homopolymer
            if (read.is_central_mapped && read.has_anti_pattern) { 
                read.rank = 'n'; ++n_cnt; --u_cnt; 
            }

            if (read.rank != 'n') {
                starts.push_back(read.covering_start); ends.push_back(read.covering_end);
            }   
            
            // reads are qualified for undetermined signature profiling if
            // 1. capped with enough 2-mer diversity (is_stable_non_ref)
            // 2. mapped in 2nd/3rd readlen quartile (is_central_mapped)
            // 3. local freq of dirty base < thresh (qc_passed)
            if (read.is_quality_map && read.qc_passed) {
                sig_u[read.non_ref_sig].push_back(i); // register sig
                
                // annot for signature usage: (kmer, overlap, grid-search)
                //  kmer: ineffctive kmer if 1
                //  overlap: positionally overlap non-target if 1
                //  grid-search: do not use this sig for grid-search
                if (u_sig_annot.find(read.non_ref_sig) == u_sig_annot.end()) {
                    u_sig_annot[read.non_ref_sig] = {0, 0, 0};
                }

                if (read.ineffective_kmer) { ++ineff_kmer; u_sig_annot[read.non_ref_sig][0] = 1; }
                if (read.has_positional_overlap) { ++overlapping; u_sig_annot[read.non_ref_sig][1] = 1; }  
                if (read.rank == 'n') { ++overlapping; u_sig_annot[read.non_ref_sig][2] = 1; }
            }
        } else { 
            read.rank = 'n'; ++n_cnt;
            if (read.is_quality_map) ++ref_hap_n;                
        }       
    }
    
    if (!s_cnt && !u_cnt) { has_no_support = true; return; }
    
    int _strt = (loc_ref.flanking_start - kmer_sz > loc_ref.start) 
              ? loc_ref.flanking_start - kmer_sz : loc_ref.start;
    starts.push_back(_strt); 
    
    int _end = (loc_ref.flanking_end + kmer_sz < loc_ref.end)
             ? loc_ref.flanking_end + kmer_sz : loc_ref.end;  
    ends.push_back(_end);
    
    start = *std::min_element(starts.begin(), starts.end());
    // iso-seq reads may have covering start < loc_ref.start
    start = (start < loc_ref.start) ? loc_ref.start : start;
    end = *std::max_element(ends.begin(), ends.end());
    // iso 
    end = (loc_ref.end < end) ? loc_ref.end : end;
        
    has_hiconf_support = (!sig_s_hiconf.empty()); has_ref_hap = (ref_hap_n); 
    
    // rerank
    if (has_hiconf_support) { 
        for (const auto& elem : sig_u) {
            for (int i : elem.second) { 
                reads[i].rank = 'n'; ++n_cnt; --u_cnt; 
                u_sig_annot[reads[i].non_ref_sig][2] = 1;
            }
        }
    }
    //keep sig_u contents to set up non-supporting haplotypes
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
    // but check u_sig_annot 
    std::string query;
    for (const auto& _sig : u_sig_cnts) {

        if (u_sig_annot[_sig.first][2] == 1) continue;    
        
        std::cout << _sig.first << std::endl;

        const auto& vars = reads[sig_u[_sig.first][0]].variants;
        make_sequence(loc_ref, vars, start, end, query);

        std::vector<Ints> grid; gap_grid(params, grid);
        bool res = search_over_grid(start, loc_ref, rseq, query, grid, target); 
        if (res) {
            const auto& read_idxes = sig_u[_sig.first];
            for (int _i : read_idxes) {
                if ( reads[_i].rank == 'n') continue;
                reads[_i].rank = 's'; ++s_cnt; --u_cnt;
                sig_s[_sig.first].push_back(_i);
                sig_s_hiconf[_sig.first].push_back(_i);
            }
            sig_u.erase(_sig.first);
       }
        query.clear();
    }

    if (!s_cnt) return; // no change
    if(sig_s_hiconf.size()) has_hiconf_support = true;
    
    // rerank to 'n' only happens with has_hiconf_support 
    for (const auto& elem : sig_u) {
        for (int i : elem.second) { reads[i].rank = 'n'; ++n_cnt; --u_cnt; }
    }
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
        hiconf_read_idx = idx0; 
        make_sequence(loc_ref, reads[idx0].variants, start, end, seq0, &i2p_0);
    } else {
        make_sequence(loc_ref, {target}, start, end, seq0, &i2p_0);
    }
    
    PatternCnt u_sig_cnt; int idx1 = -1, idx2 = -1;
    
    std::cout << "empty " << sig_u.empty() << " " << has_excess_ins_hap << std::endl;
    if (!sig_u.empty()) {
        count_patterns(sig_u, u_sig_cnt);
        hap1 = u_sig_cnt[0].first; idx1 = sig_u[hap1][0];
        make_sequence(loc_ref, reads[idx1].variants, start, end, seq1);
        if (u_sig_cnt.size() > 1) {
            hap2 = u_sig_cnt[1].first; idx2 = sig_u[hap2][0];
            make_sequence(loc_ref, reads[idx2].variants, start, end, seq2);
        }
    } else if (has_excess_ins_hap) {
    // haplotype with additional ins-repeats may be clipped and may not be captured    
        std::string alt_added = target.alt + target.alt.substr(1);
        Vars vlst = {Variant(target.pos, target.ref, alt_added)};
        make_sequence(loc_ref, vlst, start, end, seq1); idx1 = 0; //pseudo index
    }
    
    if (rseq.empty())
        make_sequence(loc_ref, {}, start, end, rseq, &i2p_r);
    
    if (idx0 < 0) vs_ref_hap = true; //non-target non-ref hap may exist
    if (idx1 == -1 && idx2 == -1) {
        vs_ref_hap = true; no_non_target_haps = true; 
    }   
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
        std::cout << read.rank << " here " << s_cnt << " " << u_cnt << std::endl;
        for (const auto& kmer : kmers_nt) 
            if (read.seq.find(kmer) != std::string_view::npos) ++(read.nmer);
        for (const auto& kmer : kmers_t) 
            if (read.seq.find(kmer) != std::string_view::npos) ++(read.smer);
        
        if (kmers_nt.size() && read.smer && !read.nmer) { 
            if (has_hiconf_support || no_non_target_haps || read.is_quality_map ) { 
                read.rank = 's'; --u_cnt; ++s_cnt; 
            } 
            else { read.rank = 'y'; --u_cnt; ++y_cnt; } // likel'y' supporting 
        }
        std::cout << read.rank << " " << s_cnt << " " << u_cnt << std::endl;
        // exclude complex cases? 
        if (kmers_t.size() && read.nmer && !read.smer) { read.rank = 'n'; --u_cnt; ++n_cnt; }
    }
    has_likely_support = (y_cnt); 
}

//------------------------------------------------------------------------------
void Pileup::searchByRealignment(const UserParams& params,
                                 LocalReference& loc_ref, const Variant& target) {
    //if (s_cnt || !u_cnt || seq0.empty()) return;
    if (s_cnt || seq0.empty()) return;
    
    int fss = -1, fse = -1, fes = -1, fee = -1, ts = -1, te = -1;  
    
    //std::cout << "y cnt " << y_cnt << std::endl;
     
    for (const auto& elem : i2p_0) {
        if (elem.second == loc_ref.flanking_start) fss = elem.first;
        if (ts < 0 && elem.second == target.pos) ts = elem.first;
        if (elem.second == target.end_pos) te = elem.first;
        if (elem.second == loc_ref.flanking_end) fee = elem.first;
    }
    
    if (fss == -1 || ts == -1 || te == -1 || fee == -1) return;
    fse = fss + params.dimer_window; fes = fee - params.dimer_window;
      
    Alignment aln; Filter filter;
    Aligner aligner(params.match_score, params.mismatch_penal,
                    params.gap_open_penal, params.gap_ext_penal);
    
    //std::cout << seq0 << std::endl;
    aligner.SetReferenceSequence(seq0.c_str(), seq0.size());
    //std::cout << fss << " " << fse << " " << ts << " " << te << " " << fes << " " << fee << std::endl; 
    std::string u_seq; std::bitset<3> check_points;
    for (const auto& elem : sig_u) {
        
        if (!u_sig_annot[elem.first][0]) continue; // realn only if kmer analysis wont work
        
        make_sequence(loc_ref, reads[elem.second[0]].variants, start, end, u_seq);
        int32_t mask_len = strlen(u_seq.c_str()) < 30 ? 15 : strlen(u_seq.c_str()) / 2;
        aligner.Align(u_seq.c_str(), filter, &aln, mask_len);
        //std::cout << u_seq << " " << aln.cigar_string << std::endl;
        check_match_pattern(aln, check_points, fss, fse, ts, te, fes, fee);
        if (check_points.count() == 3) {
            for (auto n : elem.second) reads[n].rank = 's';            
        } else if (check_points.test(1) && u_sig_annot[elem.first][1]) {
            for (auto n : elem.second) reads[n].rank = 's';
        } else if (!check_points.test(1)) {
            for (auto n : elem.second) reads[n].rank = 'n';
        } 
        u_seq.clear(); check_points.reset();
    }
    
    // loop over u (make u_sig_cnt as pileup member)
    //
    // next loop over remaining 'u'/'y'-reads

    //aligner.SetReferenceSequence(seq0.c_str(), seq0.size());
}
