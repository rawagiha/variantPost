#include "pileup.h"
#include "match.h"

typedef std::vector<std::pair<std::string_view, size_t>> PatternCnt;

//------------------------------------------------------------------------------
inline void count_patterns(const Idx& ptrn_read_idx, PatternCnt& ptrn_cnts) noexcept {
    ptrn_cnts.reserve(ptrn_cnts.size() + ptrn_read_idx.size()); // no-realloc
    for (const auto& [pattern, indices] : ptrn_read_idx) {
        ptrn_cnts.emplace_back(pattern, indices.size()); 
    }
    std::sort(ptrn_cnts.begin(), ptrn_cnts.end(), [](const auto& a, const auto& b) noexcept {
        return a.second > b.second; // more
    });
}

//------------------------------------------------------------------------------
// Set Pileup start. Long spliced reads (e.g., Iso-Seq) can have covering start < loc_ref.start
// -> Set pileup start to be >= loc_ref.start
[[nodiscard]] inline int set_start(Ints& starts, const LocalReference& loc_ref, const int kmer_sz) {
    const int s = (loc_ref.flanking_start - kmer_sz > loc_ref.start)
                ? loc_ref.flanking_start - kmer_sz : loc_ref.start;
    starts.push_back(s);
    
    int start = *std::min_element(starts.cbegin(), starts.cend());
    return std::max(start, loc_ref.start);     
}


//------------------------------------------------------------------------------
// Set Pileup end
[[nodiscard]] inline int set_end(Ints& ends, const LocalReference& loc_ref, const int kmer_sz) {
    const int e = (loc_ref.flanking_end + kmer_sz < loc_ref.end)
                ? loc_ref.flanking_end + kmer_sz : loc_ref.end;
    ends.push_back(e);

    int end = *std::max_element(ends.cbegin(), ends.cend());
    return std::min(end, loc_ref.end);
}


//------------------------------------------------------------------------------
//inline auto less = [](const auto& x, const auto& y) { return x.second < y.second; };
//inline auto more = [](const auto& x, const auto& y) { return x.second > y.second; };

// ORIG
//inline void count_patterns(const Idx& ptrn_read_idx, PatternCnt& ptrn_cnts) {
//    for (const auto& elem : ptrn_read_idx)
//        ptrn_cnts.emplace_back(elem.first, elem.second.size()); 
//    
//    std::sort(ptrn_cnts.begin(), ptrn_cnts.end(), more);
//}

//-----------------------------------------------------------------------------
// Long spliced reads (e.g., Iso-Seq) can have covering start < loc_ref.start
// -> Set pileup start to be >= loc_ref.start
//int set_start(Ints& starts, LocalReference& loc_ref, const int kmer_sz) {
//    int s = (loc_ref.flanking_start - kmer_sz > loc_ref.start)
//         ? loc_ref.flanking_start - kmer_sz : loc_ref.start;
//    starts.push_back(s);
    
    // Make sure loc_ref.start <= start
//    int start = *std::min_element(starts.begin(), starts.end());
//    start = (start < loc_ref.start) ? loc_ref.start : start; 
//    return start;     
//}

//-----------------------------------------------------------------------------
// Pileup end counterpart
//int set_end(Ints& ends, LocalReference& loc_ref, const int kmer_sz) {
//    int e = (loc_ref.flanking_end + kmer_sz < loc_ref.end)
//          ? loc_ref.flanking_end + kmer_sz : loc_ref.end;
//    ends.push_back(e);
//
//    int end = *std::max_element(ends.begin(), ends.end());
//    end = (loc_ref.end < end) ? loc_ref.end : end;
//    return end;
//}

//------------------------------------------------------------------------------
Pileup::Pileup(const Strs& names, const Bools& is_rv, const Strs& cigar_strs,
               const Ints& aln_starts, const Ints& aln_ends, const Strs& seqs,
               const Strs& quals, const Bools& is_first_bam, const bool has_scnd,
               const UserParams& params, LocalReference& loc_ref, Variant& target) {
    
    kmer_sz = std::max(loc_ref.low_cplx_len, params.kmer_size);

    int qc_start = loc_ref.flanking_start;
    if (loc_ref.start <= target.pos - target.event_radius && 
        target.pos - target.event_radius < loc_ref.flanking_start) {
        qc_start = target.pos - target.event_radius; 
    }

    int qc_end = loc_ref.flanking_end;
    if (target.pos + target.event_radius <= loc_ref.end && 
        loc_ref.flanking_end < target.pos + target.event_radius) {
        qc_end = target.pos + target.event_radius;
    }
    
    int ref_hap_n = 0, ineff_kmer = 0, overlapping = 0; 
    
    sz = static_cast<int>(names.size()); 
    Ints starts, ends; 
    reads.reserve(sz); 
    starts.reserve(sz); 
    ends.reserve(sz);
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

        if (read.covering_ptrn == 'A' && !read.idx2pos.empty()) {
            read.isCenterMapped(target);
        }

        read.is_quality_map = (read.is_stable_non_ref && read.is_central_mapped);

        // Search target complex variant in reads with multiple events in flanking
        // region (flnk_v_cnt > 1)
        if (target.is_complex 
            && read.covering_ptrn == 'A' && read.qc_passed && read.flnk_v_cnt > 1) {
            if (target.is_substitute) read.hasTargetComplexSubstitute(target);
            else read.hasTargetComplexIndel(loc_ref, target);
        }
                
        // Rank reads
        if (read.has_target) {
            read.rank = 's'; ++s_cnt;
            read.setSignatureStrings(params); sig_s[read.non_ref_sig].push_back(i);
            if (read.is_quality_map) sig_s_hiconf[read.non_ref_sig].push_back(i); 
        } else if (read.has_local_events) {
            read.rank = 'u'; ++u_cnt;
            read.setSignatureStrings(params);             
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
                read.rank = 'n'; ++n_cnt; --u_cnt;
                for (auto& v : read.variants) { 
                    if (v.in_target_flnk) { v.testForDeNovoRepeats(loc_ref); ns_vars[v]++; } 
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
                sig_u[read.non_ref_sig].push_back(i);
                
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
            else { read.rank = 'n'; ++n_cnt; } 
            
            // Catch reference haplotype in control BAM here
            if (read.is_control && read.is_quality_map) ++ref_hap_n;              
        }       
    }
    // No supporting and no undetermined reads (e.g., all non-supporting or no coverage)
    if (!s_cnt && !u_cnt) { has_no_support = true; return; }
    
    // Local reference start <= Pileup start
    start = set_start(starts, loc_ref, kmer_sz); 
    
    // Pileup end <= Local reference end
    end = set_end(ends, loc_ref, kmer_sz);
     
    // Set reference sequence to Pileup
    make_sequence(loc_ref, {}, start, end, rseq, &i2p_r);
    
    // Set if reference haplotype detected in control BAM
    // (control BAM may not be supplied)
    has_ref_hap = (ref_hap_n);
    
    // Set if high confidence target has been found
    has_hiconf_support = (!sig_s_hiconf.empty());
    if (!has_hiconf_support) return;
    
    // Sort high conf. target pattern by occurrence
    PatternCnt s_sig_cnt; count_patterns(sig_s_hiconf, s_sig_cnt); 
    // Set the most freq. pattern as target haplotype (hap0)
    hap0 = s_sig_cnt[0].first; hiconf_read_idx = sig_s_hiconf[hap0][0];

    // Heuristic: Target alignment would be unique 
    // -> may not be true if mapping algorithms are difference (e.g., DNA vs. RNA)
    if (!has_scnd) return; 

    // If the target is found in a unambiguously-mapped read,
    // patters in other unambigously mapped reads (sig_u) would be non-supporting. 
    for (const auto& elem : sig_u) {
        const auto& read_indexes =  elem.second; 
        for (const auto i : read_indexes) {
            if (reads[i].rank == 'n') continue;
            reads[i].rank = 'n'; ++n_cnt; --u_cnt;
            u_sig_annot[reads[i].non_ref_sig][2] = 1; // Don't use for grid-search
        }
    } 
}


//-----------------------------------------------------------------------------
// This method tends to overanalyze 
// -> tighter conditions to apply
// Should be run when targe is suspected to be a part of complex indels 
void Pileup::gridSearch(const UserParams& params, 
                        LocalReference& loc_ref, const Variant& target) {
    
    // Condition 0
    // Skip if target is homopolymer extending 
    if (s_cnt || sig_u.empty() || target.denovo_rep > 2) { std::cout << "grid exit 1 target denov "  << target.denovo_rep << std::endl; return;}

    PatternCnt u_sig_cnts;
    count_patterns(sig_u, u_sig_cnts);
    
    std::cout << ns_vars.size() << " check size " << std::endl;
    const int tot = s_cnt + n_cnt + u_cnt; 
    // Condition 1
    // Skip if background hap involves homoplymer variations
    //   Step1: check for homopolymer involvment (denovo_rep)
    //   Step2: define as background hap if freq is > 0.1
    for (const auto& elem : ns_vars) {
    //    std::cout << elem.first.pos << " " << elem.first.ref << " " << elem.first.alt << " " << elem.first.denovo_rep << " " << elem.second << " " << tot << std::endl;
        if (elem.first.denovo_rep > 1 && elem.second * 10 > tot) { std::cout << "grid exit 2 denov " << elem.first.denovo_rep << " freq "  << elem.second << " " << tot << std::endl; return;}
    }
    
    std::string query;
    std::vector<Ints> grid;
    gap_grid(params, grid);
    for (const auto& [sig, _] : u_sig_cnts) {
        // Condition 3
        // Skip with anti-pattern/smaller than target changes
        //if (u_sig_annot[_sig.first][2] == 1) continue;    
        
        Vars vars;
        bool skip_this = false;
        for (auto& v : reads[sig_u[sig][0]].variants) {
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
        std::cout << "grid skip happedn " << skip_this << " " << sig << std::endl;
        if (skip_this ||vars.empty()) continue;
        
        make_sequence(loc_ref, vars, start, end, query);

        // grid-search
        const bool res = search_over_grid(start, loc_ref, params, vars.size(), rseq, query, grid, target); 
        if (res) {
            const auto& read_idxes = sig_u[sig];
            for (int _i : read_idxes) {
                if ( reads[_i].rank == 'n') continue;
                reads[_i].rank = 's'; ++s_cnt; --u_cnt;
                sig_s[sig].push_back(_i);
                sig_s_hiconf[sig].push_back(_i);
            }
            sig_u.erase(sig);
       } 
       query.clear();
    }
    
    // Exit if no supporting newly found
    if (!s_cnt) return;
    
    has_hiconf_support = true;
    
    //LOGIC: if target is aligned in a sigature, the remaining ones
    //are likely non-supproting
    for (const auto& [_, indexes] : sig_u) {
        for (int i : indexes) { reads[i].rank = 'n'; ++n_cnt; --u_cnt; }
    }
    // keep sig_u contents to set up non-supporting haplotypes
}

//-----------------------------------------------------------------------------
// LOGIC: target is already found (s_cnt > 0) 
//        -> other variations are unlikely supporting
//           (assuming target is supported by a single hap)
//        
//        variations found as anti pattern are non-supporting
[[nodiscard]] inline bool hap_settable(const int s_cnt, 
                         const std::string_view hap, 
                         const std::vector<std::string_view>& anti_sigs) noexcept {
    if (s_cnt) return true;
    return std::find(anti_sigs.cbegin(), anti_sigs.cend(), hap) != anti_sigs.cend();
}

//ORIG
//inline bool hap_settable(const int s_cnt, 
//                         const std::string_view hap, 
//                         const std::vector<std::string_view>& anti_sigs) {
//    if (s_cnt) return true;
//    if (std::find(anti_sigs.begin(), anti_sigs.end(), hap) 
//        != anti_sigs.end()) return true;
//    return false;
//}

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
    } else {
        make_sequence(loc_ref, {target}, start, end, seq0, &i2p_0);
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
    
    int cntl_cov = 0;
    for (const auto& read : reads) {
        if (!read.fail_to_cover_flankings && read.covering_ptrn == 'A' && read.is_control) ++cntl_cov;
    }
    std::cout << cntl_cov << " control coverage " << std::endl; 
    
    
    if (!sig_u.empty()) {
        count_patterns(sig_u, u_sig_cnt);
        
        //TODO this is test purpose
        const int hap1_cnt = static_cast<int>(u_sig_cnt[0].second);
        const int hap2_cnt = static_cast<int>(u_sig_cnt[1].second);
        bool set_hap1 = false, set_hap2 = false;
        if (hap1_cnt > 4 && hap1_cnt * 5 >= cntl_cov) set_hap1 = true;
        if (hap2_cnt > 4 && hap2_cnt * 5 >= cntl_cov) set_hap2 = true; 

        hap1 = u_sig_cnt[0].first; 
        std::cout << hap1 << " " << u_sig_cnt[0].second << std::endl;
        if (hap_settable(s_cnt, hap1, anti_sigs) && set_hap1) {
            idx1 = sig_u[hap1][0];  
            make_sequence(loc_ref, reads[idx1].variants, start, end, seq1, &i2p_1);
            if (seq0 == seq1) seq1.clear();
        }
        
        if (u_sig_cnt.size() > 1) {
            hap2 = u_sig_cnt[1].first;
            std::cout << hap2 << " " << u_sig_cnt[1].second << std::endl;
            if (hap_settable(s_cnt, hap2, anti_sigs) && set_hap2) {
                idx2 = sig_u[hap2][0];
                make_sequence(loc_ref, reads[idx2].variants, start, end, seq2, &i2p_2);
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
}

[[nodiscard]] bool is_target_covering(const int start, const int end,
                        const Read& read, const std::string_view& kmer, 
                        const std::unordered_map<size_t, int>& dict) noexcept {
    const size_t i = read.seq.find(kmer); 
    if (i == std::string_view::npos) return false;
    
    const size_t j = i + kmer.size() - 1;   
    
    auto it_i = dict.find(i);
    auto it_j = dict.find(j);

    if (it_i == dict.end() || it_j == dict.end()) {
        return (read.covered_in_clip || read.has_local_clip);
    }
    
    if (it_j->second < start || end < it_i->second) return false;
    return true;
}

[[nodiscard]] inline bool has_the_sig(const std::string_view sig, const Strs& anti_sig) noexcept {
    return std::find(anti_sig.cbegin(), anti_sig.cend(), sig) != anti_sig.cend();
}


/*
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

inline bool has_the_sig(const std::string& sig, const Strs& anti_sig) {
    auto it = std::find(anti_sig.begin(), anti_sig.end(), sig);
    return (it != anti_sig.end());
}*/

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
                   std::inserter(km12, km12.end()));
    make_kmers(rseq, kmer_sz, kmr);
    std::set_union(km12.begin(), km12.end(), kmr.begin(), kmr.end(),
                   std::inserter(km12r, km12r.end()));
    
    // Kmer set diff
    std::set_difference(km0.begin(), km0.end(), km12r.begin(), km12r.end(),
                        std::inserter(kmers_t, kmers_t.end()));
    std::set_difference(km12r.begin(), km12r.end(), km0.begin(), km0.end(),
                        std::inserter(kmers_nt, kmers_nt.end()));

    std::cout << "hap0 " << seq0 <<  " " << kmers_t.size() << " " << kmers_nt.size() << std::endl;
    std::cout << "hap1 " << seq1 << std::endl;
    std::cout << "hap2 " << seq2 << std::endl;
    std::cout << "hapr " << rseq << std::endl; 
    
    const size_t target_kmer_sz = kmers_t.size();
    const size_t non_target_kmer_sz = kmers_nt.size(); 
    const int rad_start = target.pos - target.event_radius;
    const int rad_end = target.pos + target.event_radius;
    Strs anti_sig;
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

        
        if (non_target_kmer_sz && read.smer > read.nmer) { 
            std::cout << " this is y " << read.name << " " << read.cigar_str << " " << read.smer << " " << read.nmer << " " << read.non_ref_sig << std::endl;
            read.rank = 'y'; --u_cnt; ++y_cnt; //Likely supporting
            if (!read.nmer && read.smer > 1) {
                if (has_hiconf_support 
                    || no_non_target_haps 
                    || read.is_quality_map) {
                   read.rank = 's'; --y_cnt; ++s_cnt;
                }
            }
        } else if (target_kmer_sz && read.smer <= read.nmer) { 
            if (!read.smer && read.nmer) {
                read.rank = 'n'; --u_cnt; ++n_cnt; 
                anti_sig.push_back(read.non_ref_sig);
            } else if (read.smer && read.nmer) {
                read.rank = 'z'; --u_cnt; ++z_cnt; //Unlikely supporting
            }
        } 
    }
    
    // reclassify
    for (auto& read : reads) {
        const char rank_ = read.rank;
        if (rank_ == 's' || rank_ == 'n') continue; 
        const bool has_it = has_the_sig(read.non_ref_sig, anti_sig);
        if (has_it){
            read.rank = 'n'; ++n_cnt;
            if (rank_ == 'u') --u_cnt;
            else if (rank_ == 'y') --y_cnt;
            else --z_cnt;
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
    for (auto& read : reads) {
        if (!read.qc_passed) continue;
        if (read.rank != 'y' && !read.covered_in_clip) continue;
        
        std::string ss{read.seq}; // may be improved 

        int32_t mask_len = strlen(ss.c_str()) < 30 ? 15 : strlen(ss.c_str()) / 2;
        aligner.Align(ss.c_str(), filter, &aln, mask_len);
        std::cout << ss << " " << aln.cigar_string << std::endl;
        check_match_pattern(aln, check_points, fss, fse, ts, te, fes, fee); 
        std::cout <<  check_points.count() << std::endl;
        
        // perferct match
        if (check_points.count() == 3) {
            read.rank = 's'; --y_cnt; ++s_cnt;
        }
        
        check_points.reset();
    }
    
    if (s_cnt) {
        for (auto& read : reads) { 
            if (read.rank == 'y') { read.rank = 's'; --y_cnt; ++s_cnt; }
        }
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
