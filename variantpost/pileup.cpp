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

    int qc_start = loc_ref.flanking_start;
    if (loc_ref.start <= target.pos - target.event_radius 
        && target.pos - target.event_radius < loc_ref.flanking_start) {
        qc_start = target.pos - target.event_radius; 
    }
    int qc_end = loc_ref.flanking_end;
    if (target.pos + target.event_radius <= loc_ref.end 
        && loc_ref.flanking_end < target.pos + target.event_radius) {
        qc_end = target.pos + target.event_radius;
    }
    
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

        read.parseLocalPattern(loc_ref, target, params);

        read.qualityCheck(qc_start, qc_end, params.base_q_thresh, params.lq_rate_thresh);
        read.trimLowQualBases(params.base_q_thresh);    
        read.isStableNonReferenceAlignment(loc_ref); 

        if (read.covering_ptrn == 'A' && !read.idx2pos.empty())
            read.isCenterMapped(target);

        read.is_quality_map = (read.is_stable_non_ref && read.is_central_mapped);
         
        //consider multipl vars in flnk
        if (target.is_complex 
            && read.covering_ptrn == 'A' && read.qc_passed && read.flnk_v_cnt > 1) {
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
            
            //if (target.is_substitute && read.covering_ptrn == 'A') {
            //    read.rank == 'n'; +n_cnt; --u_cnt;
            //    continue;
            //}


            read.checkByRepeatCount(target, has_excess_ins_hap);
            
            // Rank "anti-pattern" & "smaller change" as 'n'
            // Anti-pattern:
            // Single, non-target mutation closed with high-cplx sequence (anti-pattern)
            // example:
            //      target    del(C)
            //                    *
            //        ref GCTAAAAACAAAATGC
            //      read  GCTAAAAAAAAAATGC  de-novo homopolymer
            // Smaller change (non-complex):
            // Fewer base changes than target (e.g. a SNV for indel target)
            // 
            // Store such variants in ns_vars 
            if (read.is_stable_non_ref && (read.has_anti_pattern || read.has_smaller_change)) { 
                if (read.has_anti_pattern) anti_sigs.push_back(read.non_ref_sig);
                
                std::cout << " n reads " << read.name << " " << read.cigar_str << std::endl;
                read.rank = 'n'; ++n_cnt; --u_cnt;
                for (auto& v : read.variants) { 
                    if (v.in_target_flnk) { 
                        v.testForDeNovoRepeats(loc_ref); ns_vars[v]++; 
                    } 
                }
            }

            if (read.rank != 'n') {
                starts.push_back(read.covering_start); ends.push_back(read.covering_end);
            }   
            
            // Reads are qualified for undetermined signature profiling if
            // 1. capped with enough 2-mer diversity (is_stable_non_ref)
            // 2. mapped in 2nd/3rd readlen quartile (is_central_mapped)
            // 3. local freq of dirty base < thresh (qc_passed)
            if (read.is_stable_non_ref && read.qc_passed && !read.has_local_clip) {
                sig_u[read.non_ref_sig].push_back(i); // register sig
                
                // annot for signature usage: (kmer, overlap, grid-search)
                // kmer: ineffctive kmer if 1
                // overlap: positionally overlap non-target if 1
                // grid-search: do not use this sig for grid-search
                if (u_sig_annot.find(read.non_ref_sig) == u_sig_annot.end()) {
                    u_sig_annot[read.non_ref_sig] = {0, 0, 0};
                }

                if (read.ineffective_kmer) { ++ineff_kmer; u_sig_annot[read.non_ref_sig][0] = 1; }
                if (read.has_positional_overlap) { ++overlapping; u_sig_annot[read.non_ref_sig][1] = 1; }  
                if (read.rank == 'n') { u_sig_annot[read.non_ref_sig][2] = 1; }
            }
        } else { 
            if (read.has_local_clip || read.covered_in_clip) { read.rank = 'u'; ++u_cnt; }
            else { read.rank = 'n'; ++n_cnt; if (read.is_quality_map) ++ref_hap_n; }                
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
                if (reads[i].rank != 'n') {
                    reads[i].rank = 'n'; ++n_cnt; --u_cnt; 
                    u_sig_annot[reads[i].non_ref_sig][2] = 1;
                }
            }
        }
    }
    
    if (!sig_s_hiconf.empty()) {
        PatternCnt s_sig_cnt; 
        count_patterns(sig_s_hiconf, s_sig_cnt);
        hap0 = s_sig_cnt[0].first; 
        hiconf_read_idx = sig_s_hiconf[hap0][0]; 
    }
    
    std::cout << s_cnt << " " << n_cnt << " " << u_cnt << std::endl;
}


//-----------------------------------------------------------------------------
// This method tends to overanalyze 
// -> tighter conditions to apply
// Should be run when targe is suspected to be a part of complex indels 
void Pileup::gridSearch(const UserParams& params, 
                        LocalReference& loc_ref, const Variant& target) {
    
    // Condition 0
    // Skip if target is homopolymer extending 
    if (s_cnt || sig_u.empty() || target.denovo_rep > 2) return;

    make_sequence(loc_ref, {}, start, end, rseq, &i2p_r);
    
    PatternCnt u_sig_cnts;
    count_patterns(sig_u, u_sig_cnts);
    
    std::cout << ns_vars.size() << " check size " << std::endl;
    const int tot = s_cnt + n_cnt + u_cnt; 
    // Condition 1
    // Skip if background hap involves homoplymer variations
    //   Step1: check for homopolymer involvment (denovo_rep)
    //   Step2: define as background hap if freq is > 0.1
    for (const auto& elem : ns_vars) {
        std::cout << elem.first.pos << " " << elem.first.ref << " " << elem.first.alt << " " << elem.first.denovo_rep << " " << elem.second << std::endl;
        if (elem.first.denovo_rep > 1 && elem.second * 10 > tot) return;
    }

    std::string query;
    for (const auto& _sig : u_sig_cnts) {
        // Condition 3
        // Skip with anti-pattern/smaller than target changes
        if (u_sig_annot[_sig.first][2] == 1) continue;    
        
        Vars vars;
        bool skip_this = false;
        for (auto& v : reads[sig_u[_sig.first][0]].variants) {
            // Condition 4
            // Skip with variations recurrently found 
            //           in reads already ranked as 'n' 
            //           also in target flank start/end                                   
            if (ns_vars.count(v) && ns_vars[v] > 1) {
                skip_this = true; break;
            }
            
            // Condition 5
            // Do not realn with homopolymer extending vars
            v.testForDeNovoRepeats(loc_ref);
            if (v.in_target_flnk && v.denovo_rep > 2) {
                skip_this = true; break;
            }
            
            // Condition 6
            // Do not use low quality variations
            if (v.mean_qual > params.base_q_thresh) vars.push_back(v);
        }
        if (skip_this ||vars.empty()) continue;
        
        make_sequence(loc_ref, vars, start, end, query);

        // grid-search
        std::vector<Ints> grid; gap_grid(params, grid);
        const bool res = search_over_grid(start, loc_ref, params, vars.size(), rseq, query, grid, target);
        
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
    
    // Exit if no supporting newly found
    if (!s_cnt) return;
    
    has_hiconf_support = true;
    
    //LOGIC: if target is aligned in a sigature, the remaining ones
    //are likely non-supproting
    for (const auto& elem : sig_u) {
        for (int i : elem.second) { reads[i].rank = 'n'; ++n_cnt; --u_cnt; }
    }
    // keep sig_u contents to set up non-supporting haplotypes
}


// LOGIC: target is already found (s_cnt > 0) 
//        -> other variations are unlikely supporting
//           (assuming target is supported by a single hap)
//        
//        variations found as anti pattern are non-supporting
bool is_settable(const int s_cnt, 
                 const std::string_view hap, 
                 const std::vector<std::string_view>& anti_sigs) {
    if (s_cnt) return true;
    if (std::find(anti_sigs.begin(), anti_sigs.end(), hap) 
        != anti_sigs.end()) return true;
    return false;
}

//------------------------------------------------------------------------------
void Pileup::setHaploTypes(LocalReference& loc_ref, const Variant& target) {
    PatternCnt s_sig_cnt; 
    
    if (hiconf_read_idx > -1) {
        make_sequence(
            loc_ref, reads[hiconf_read_idx].variants, start, end, seq0, &i2p_0);
    } else if (!sig_s.empty()) {
        PatternCnt s_sig_cnt;
        count_patterns(sig_s, s_sig_cnt); hap0 = s_sig_cnt[0].first;
        make_sequence(
            loc_ref, reads[sig_s[hap0][0]].variants, start, end, seq0, &i2p_0);
    }

    /*if (!sig_s_hiconf.empty()) {
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
    }*/
    
    PatternCnt u_sig_cnt; int idx1 = -1, idx2 = -1;
    
    if (!sig_u.empty()) {
        count_patterns(sig_u, u_sig_cnt);
        
        hap1 = u_sig_cnt[0].first; 
        std::cout << hap1 << std::endl;
        if (is_settable(s_cnt, hap1, anti_sigs)) {
            idx1 = sig_u[hap1][0];  
            make_sequence(loc_ref, reads[idx1].variants, start, end, seq1);
            if (seq0 == seq1) seq1.clear();
        }
        
        if (u_sig_cnt.size() > 1) {
            hap2 = u_sig_cnt[1].first;
            std::cout << hap2 << std::endl;
            if (is_settable(s_cnt, hap2, anti_sigs)) {
                idx2 = sig_u[hap2][0];
                make_sequence(loc_ref, reads[idx2].variants, start, end, seq2);
                if (seq0 == seq2 || seq1 == seq2 ) seq2.clear();
            }
        }
    } else if (has_excess_ins_hap) {
    // Haplotype with additional ins-repeats may be clipped and may not be captured    
        std::string alt_added = target.alt + target.alt.substr(1);
        Vars vlst = {Variant(target.pos, target.ref, alt_added)};
        make_sequence(loc_ref, vlst, start, end, seq1); idx1 = 0; //pseudo index
    }
    
    if (seq1.empty() && seq2.empty()) no_non_target_haps = true;
         
    // Always set reference hap
    make_sequence(loc_ref, {}, start, end, rseq, &i2p_r);
}


// Only count kmers overlapping target radius 
bool is_target_covering(const int start, const int end,
                        const Read& read, const std::string_view& kmer, 
                        const std::unordered_map<size_t, int>& dict) {
    const size_t i = read.seq.find(kmer); 
    if (i == std::string_view::npos) return false;
    
    const size_t j = i + kmer.size() - 1;   
    if (dict.find(i) == dict.end() || dict.find(j) == dict.end()) {
        // LOGIC lipped segments have no positional info
        if (read.covered_in_clip || read.has_local_clip) return true;
        else return false;
    }
    if (dict.at(j) < start || end < dict.at(i)) return false;

    return true;
}

//------------------------------------------------------------------------------
void Pileup::differentialKmerAnalysis(const UserParams& params, 
                                      LocalReference& loc_ref, const Variant& target) {
    
    // Target connecting homopolymers -> skip
    if (target.denovo_rep == 4) return;
    
    Kmers km0, km1, km2, kmr, km12, km12r;   
    
    // Target kmers
    make_kmers(seq0, kmer_sz, km0); 
    
    // Comparators
    make_kmers(seq1, kmer_sz, km1);
    make_kmers(seq2, kmer_sz, km2);
    std::set_union(km1.begin(), km1.end(), km2.begin(), km2.end(),
                   std::inserter(km12, km12.begin()));
    make_kmers(rseq, kmer_sz, kmr);
    std::set_union(km12.begin(), km12.end(), kmr.begin(), kmr.end(),
                   std::inserter(km12r, km12r.begin()));
    
    // Kmer set diff
    std::set_difference(km0.begin(), km0.end(), km12r.begin(), km12r.end(),
                        std::inserter(kmers_t, kmers_t.end()));
    std::set_difference(km12r.begin(), km12r.end(), km0.begin(), km0.end(),
                        std::inserter(kmers_nt, kmers_nt.end()));

    std::cout << "hap0 " << seq0 <<  " " << kmers_t.size() << " " << kmers_nt.size() << std::endl;
    std::cout << "hap1 " << seq1 << std::endl;
    std::cout << "hap2 " << seq2 << std::endl;
    std::cout << "hapr " << rseq << std::endl; 
    
    const int rad_start = target.pos - target.event_radius;
    const int rad_end = target.pos + target.event_radius;
    for (auto& read : reads) {
        std::cout << read.name << " " << read.cigar_str << " " << read.rank << std::endl;
        if (!read.qc_passed || read.rank != 'u') continue;
         
        std::unordered_map<size_t, int> dict;
        for (auto& elem : read.idx2pos) 
            dict[static_cast<size_t>(elem.first)] = elem.second;
        
        for (const auto& kmer : kmers_nt) {
            if (is_target_covering(rad_start, rad_end, read, kmer, dict)) 
                ++(read.nmer);
        }
        
        for (const auto& kmer : kmers_t) {
            if (!is_target_covering(rad_start, rad_end, read, kmer, dict)) 
                continue;
            
            // Positive count further must be in quality region
            int _i = static_cast<int>(read.seq.find(kmer));
            if (read.qs <= _i && _i + kmer_sz <= read.qe) ++(read.smer);
        }

        
        if (kmers_nt.size() && read.smer > read.nmer) { 
            read.rank = 'y'; --u_cnt; ++y_cnt; // likely supporting
            if (!read.nmer && read.smer > 1) {
                if (has_hiconf_support 
                    || no_non_target_haps 
                    || read.is_quality_map) {
                   read.rank = 's'; --y_cnt; ++s_cnt;
                }
            }
        } 
        
        if (kmers_t.size() && read.nmer && !read.smer) { 
            read.rank = 'n'; --u_cnt; ++n_cnt; 
        }
    }
    
    has_likely_support = (y_cnt);
    
    std::cout << " y cnt ! " << y_cnt << std::endl; 
     
}

//------------------------------------------------------------------------------
//NOTE THIS STEP DOES NOT MAKE SENSE
//WHY DO THIS HERE? 
//REALN FOR READS, NOT RE_CONSTRUSTED SEQ FROM U_SIG, MAY BE OK.
void Pileup::searchByRealignment(const UserParams& params,
                                 LocalReference& loc_ref, const Variant& target) {
    //if (s_cnt || !u_cnt || seq0.empty()) return;
    if (s_cnt || seq0.empty()) return;
    
    int fss = -1, fse = -1, fes = -1, fee = -1, ts = -1, te = -1;  
    
    std::cout << "y cnt " << y_cnt << std::endl;
     
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
    //std::string u_seq; i
    
    std::bitset<3> check_points;
    for (const auto& read : reads) {
        if (!read.qc_passed) continue;
        if (read.rank != 'y' && !read.covered_in_clip) continue;
        std::string ss{read.seq};
        int32_t mask_len = strlen(ss.c_str()) < 30 ? 15 : strlen(ss.c_str()) / 2;
        aligner.Align(ss.c_str(), filter, &aln, mask_len);
        std::cout << ss << " " << aln.cigar_string << std::endl;
        check_match_pattern(aln, check_points, fss, fse, ts, te, fes, fee); 
        std::cout <<  check_points.count() << std::endl;
        check_points.reset();
    }
    
    
    /*
    for (const auto& elem : sig_u) {
        
        // realn only if ineffctive kmer and no anti-pattern
        if (!u_sig_annot[elem.first][0] || u_sig_annot[elem.first][2]) continue; // realn only if kmer analysis wont work
        
        make_sequence(loc_ref, reads[elem.second[0]].variants, start, end, u_seq);
        int32_t mask_len = strlen(u_seq.c_str()) < 30 ? 15 : strlen(u_seq.c_str()) / 2;
        aligner.Align(u_seq.c_str(), filter, &aln, mask_len);
        std::cout << u_seq << " " << aln.cigar_string << std::endl;
        check_match_pattern(aln, check_points, fss, fse, ts, te, fes, fee);
        if (check_points.count() == 3) {
            for (auto n : elem.second) { if (reads[n].rank != 'n') reads[n].rank = 's'; }           
        } else if (check_points.test(1) && u_sig_annot[elem.first][1]) {
            for (auto n : elem.second) { if (reads[n].rank != 'n') reads[n].rank = 's'; }
        } else if (!check_points.test(1)) { 
            for (auto n : elem.second) reads[n].rank = 'n';
        } 
        u_seq.clear(); check_points.reset();
    }*/
    
    // loop over u (make u_sig_cnt as pileup member)
    //
    // next loop over remaining 'u'/'y'-reads

    //aligner.SetReferenceSequence(seq0.c_str(), seq0.size());
}
