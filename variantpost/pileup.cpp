#include "pileup.h"
#include "match.h"
#include <unordered_set>
#include <queue>
#include <cstring>

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
Pileup::Pileup(const Strs& names, const Bools& is_rv, const Strs& cigar_strs,
               const Ints& aln_starts, const Ints& aln_ends, const Strs& seqs,
               const Strs& quals, const Bools& is_first_bam, bool has_scnd,
               const UserParams& params, LocalReference& loc_ref, Variant& target) 
    : has_second_bam(has_scnd) {
    if (has_second_bam) control_bam_cov = 0;
    
    sz = static_cast<int>(names.size());
    kmer_sz = std::max(loc_ref.low_cplx_len, params.kmer_size);

    const int target_left_lim = target.pos - target.event_radius;
    const int qc_start = (loc_ref.start <= target_left_lim && target_left_lim < loc_ref.flanking_start)
                         ? target_left_lim : loc_ref.flanking_start;
    const int target_right_lim = target.pos + target.event_radius;
    const int qc_end = (target_right_lim <= loc_ref.end && loc_ref.flanking_end < target_right_lim)
                       ? target_right_lim: loc_ref.flanking_end;
    
    //--- Search the target complex variant if the read satisfies
    //  1. is completely covered the locus ('A') 
    //  2. has multiple variants within flanking region (flnk_v_cnt > 1)
    auto search_complex_target = [&](Read& read) {
        if (target.is_complex && read.covering_ptrn == 'A' && read.qc_passed && read.flnk_v_cnt > 1) {
            target.is_substitute ? read.hasTargetComplexSubstitute(target)
                                 : read.hasTargetComplexIndel(loc_ref, target);
        }
    };
    
    int ref_hap_n = 0, ineff_kmer = 0, overlapping = 0;
    
    //--- Annotate non-reference varitions needing further analysis to determine if
    //     supporting the target as Undetermined Signature. The read must satisfy 
    //     1. has non-ambigously-mapped, non-reference pattern (is_stable_non_ref)
    //     2. clean bases (qc_passed)
    //     3. read is not clipped near the locus of interest
    auto annotate_undetermined_signatures = [&](Read& read, int read_idx) {
        if (read.is_stable_non_ref && read.qc_passed && !read.has_local_clip) {
            const auto& sig = read.non_ref_sig;
            sig_u[sig].push_back(read_idx);
            
            // TODO explain the following part later
            auto& annot = u_sig_annot[sig];
            if (annot.empty()) annot = {0, 0, 0};

            if (read.ineffective_kmer)       { ++ineff_kmer;  annot[0] = 1; }
            if (read.has_positional_overlap) { ++overlapping; annot[1] = 1; }
            if (read.rank == 'n')            {                annot[2] = 1; }
        }
    };
    
    //--- Demote reads with anti-patterns to non-supporting ('n') 
    //    
    auto process_anti_patterns = [&](Read& read) {
        if (read.is_stable_non_ref && (read.has_anti_pattern || read.has_smaller_change)) { 
            read.rank = 'n'; ++n_cnt; --u_cnt;
            
            if (read.has_anti_pattern) anti_sigs.push_back(read.non_ref_sig);
            for (auto& v : read.variants) { 
                if (v.in_target_flnk) { v.testForDeNovoRepeats(loc_ref); ns_vars[v]++; } 
            }
        }
    };
   
    Ints starts, ends;  
    reads.reserve(sz); starts.reserve(sz); ends.reserve(sz);
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
            if (has_second_bam && read.is_control 
                &&!read.fail_to_cover_flankings) ++control_bam_cov; 
        }
        
        search_complex_target(read);
                
        //--- Read classification
        if (read.has_target) {
            read.rank = 's'; ++s_cnt;
            read.setSignatureStrings(params); sig_s[read.non_ref_sig].push_back(i);
            if (read.is_quality_map) { sig_s_hiconf[read.non_ref_sig].push_back(i); has_hiconf_support = true; }
        } else if (read.has_local_events) {
            read.rank = 'u'; ++u_cnt;
            read.setSignatureStrings(params);             
            read.checkByRepeatCount(target, has_excess_ins_hap);
            
            process_anti_patterns(read);
            annotate_undetermined_signatures(read, i);

            if (read.rank != 'n') {
                starts.push_back(read.covering_start); ends.push_back(read.covering_end);
            }   
        } else { 
            if (read.has_local_clip || read.covered_in_clip) { read.rank = 'u'; ++u_cnt; }
            else { read.rank = 'n'; ++n_cnt; } 
            
            // Catch reference haplotype in control BAM here
            if (read.is_control && read.is_quality_map) { ++ref_hap_n; has_ref_hap = true; }              
        }       
    }
    
    if (!s_cnt && !u_cnt) { has_no_support = true; return; }
    
    //--- Define pileup start/end within local reference start/end 
    start = set_start(starts, loc_ref, kmer_sz); 
    end = set_end(ends, loc_ref, kmer_sz);
    make_sequence(loc_ref, {}, start, end, rseq, &i2pr);
    
    //--- Focus on high confidence case (supporting reads with no mapping ambiguity)
    //    --> pick up most frequent pattern 
    if (!has_hiconf_support) return;
    else {
        PatternCnt s_sig_cnt; count_patterns(sig_s_hiconf, s_sig_cnt); 
        hap0 = s_sig_cnt[0].first; hiconf_read_idx = sig_s_hiconf[hap0][0];
    }
     
    //--- Focus on a single BAM file && high-confidence already found case (exit earlier)
    //    Undetermined singatures are also well-mapped and are unlikely to be supporting.
    //    (target is unique.). With second BAM, this may not be the case (DNA vs. RNA)
    if (!has_second_bam) return; 
    else { 
        for (const auto& [_, read_indexes] : sig_u) {
            for (const auto i : read_indexes) {
                if (reads[i].rank == 'n') continue;
                reads[i].rank = 'n'; ++n_cnt; --u_cnt;
                u_sig_annot[reads[i].non_ref_sig][2] = 1; // Don't use for grid-search
            }
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

//------------------------------------------------------------------------------
void Pileup::setHaploTypes(LocalReference& loc_ref, const Variant& target) {
    PatternCnt s_sig_cnt; 
    
    if (hiconf_read_idx > -1) {
        //if (seq0.empty()) {
            make_sequence(
                loc_ref, reads[hiconf_read_idx].variants, start, end, seq0, &i2p0);
        //}
    } else if (!sig_s.empty()) {
        PatternCnt s_sig_cnt;
        count_patterns(sig_s, s_sig_cnt); 
        hap0 = s_sig_cnt[0].first;
        if (seq0.empty()) {
            hap0 = s_sig_cnt[0].first;
            make_sequence(
                loc_ref, reads[sig_s[hap0][0]].variants, start, end, seq0, &i2p0);
        } else {
            auto& read = reads[sig_s[hap0][0]];
            if (read.is_stable_non_ref) {
                make_sequence(loc_ref, read.variants, start, end, seq0, &i2p0); 
            } 
        }

    } else {
        make_sequence(loc_ref, {target}, start, end, seq0, &i2p0);
    }

    PatternCnt u_sig_cnt; int idx1 = -1, idx2 = -1;
    
    //--- Paired BAM file case
    if (has_second_bam && control_bam_cov) {
        // pass for now
    }
    
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
            make_sequence(loc_ref, reads[idx1].variants, start, end, seq1, &i2p1);
            if (seq0 == seq1) seq1.clear();
        }
        
        if (u_sig_cnt.size() > 1) {
            hap2 = u_sig_cnt[1].first;
            std::cout << hap2 << " " << u_sig_cnt[1].second << std::endl;
            if (hap_settable(s_cnt, hap2, anti_sigs) && set_hap2) {
                idx2 = sig_u[hap2][0];
                make_sequence(loc_ref, reads[idx2].variants, start, end, seq2, &i2p2);
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


[[nodiscard]] inline bool is_target_covering_fast(
    const int rad_start, const int rad_end,
    const Read& read, const std::string_view kmer, const int kmer_sz)  noexcept
{
    size_t pos_in_read = read.seq.find(kmer);
    if (pos_in_read == std::string_view::npos) return false;

    size_t end_idx = pos_in_read + kmer_sz - 1;

    if (read.idx2pos.empty() || end_idx >= read.idx2pos.size()) {
        return false;
    }

    int genome_start = read.idx2pos[pos_in_read];
    int genome_end   = read.idx2pos[end_idx];

    if (genome_start == -1 || genome_end == -1) {
        return (read.covered_in_clip || read.has_local_clip);
    }

    return !(genome_end < rad_start || genome_start > rad_end);
}


[[nodiscard]] inline bool has_the_sig(const std::string_view sig, const Strs& anti_sig) noexcept {
    return std::find(anti_sig.cbegin(), anti_sig.cend(), sig) != anti_sig.cend();
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
        //for (auto& elem : read.idx2pos) 
        //    dict[static_cast<size_t>(elem.first)] = elem.second;
        //for (size_t j = 0; j < read.idx2pos.size(); ++k) {
        //    dict[j] = read.idx2pos[j];
        //} 
        
        for (const auto& kmer : kmers_nt) {
           //if (is_target_covering(rad_start, rad_end, read, kmer, dict)) 
           //     ++(read.nmer);
           if (is_target_covering_fast(rad_start, rad_end, read, kmer, kmer_sz)) {
                    ++read.nmer;
           }  
        }
        
        for (const auto& kmer : kmers_t) {
            //if (!is_target_covering(rad_start, rad_end, read, kmer, dict)) 
            //    continue;
            if (!is_target_covering_fast(rad_start, rad_end, read, kmer, kmer_sz))
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
    
    for (int i = 0; i < static_cast<int>(i2p0.size()); ++i) {
        if (i2p0[i] == loc_ref.flanking_start) fss = i;
        if (ts < 0 && i2p0[i] == target.pos) ts = i;
        if (i2p0[i] == target.end_pos) te = i;
        if (i2p0[i] == loc_ref.flanking_end) fee = i;
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

struct DBGNode {
    std::string kmer;
    int coverage = 0;
    std::vector<std::string> out_edges;
};

inline void dfs_find_longest_path(
    const std::string& current,
    const std::unordered_map<std::string, DBGNode>& graph,
    std::unordered_set<std::string>& visited,
    std::string& current_path,
    std::string& best_path,
    const size_t max_depth = 2000
) noexcept {
    if (current_path.size() > max_depth) return;

    auto it = graph.find(current);
    if (it == graph.end() || it->second.out_edges.empty()) {
        if (current_path.size() > best_path.size()) {
            best_path = current_path;
        }
        return;
    }

    // 高速化のため、カバレッジ（リード支持数）が高いエッジを優先して走査
    auto edges = it->second.out_edges;
    std::sort(edges.begin(), edges.end(), [&](const std::string& a, const std::string& b) noexcept {
        auto it_a = graph.find(a);
        auto it_b = graph.find(b);
        int cov_a = (it_a != graph.end()) ? it_a->second.coverage : 0;
        int cov_b = (it_b != graph.end()) ? it_b->second.coverage : 0;
        return cov_a > cov_b;
    });

    for (const auto& next_kmer : edges) {
        if (visited.find(next_kmer) == visited.end()) {
            visited.insert(next_kmer);
            size_t orig_sz = current_path.size();
            current_path.push_back(next_kmer.back()); // K-merの最後の塩基を追加

            dfs_find_longest_path(next_kmer, graph, visited, current_path, best_path, max_depth);

            current_path.resize(orig_sz); // バックトラック
            visited.erase(next_kmer);
        }
    }

    // どの分岐も行けなかった場合はそこまでのパスを評価
    if (current_path.size() > best_path.size()) {
        best_path = current_path;
    }
}

//==============================================================================
// Pileupクラスのメソッド: searchByDeBruijnGraph（完全実装）
//==============================================================================
void Pileup::searchByDeBruijnGraph(const UserParams& params,
                                   LocalReference& loc_ref,
                                   const Variant& target)
{
    // 早期リターン条件: 既に十分な高信頼度支持（s_cnt）が得られているか、
    // あるいは検証対象の 'u' リードがそもそも存在しない場合はスキップ
    //if (s_cnt || u_cnt == 0) return;
    if (s_cnt || y_cnt == 0) return;

    // 1. K-merサイズの決定（paramsの指定またはデフォルト値）
    // DBGの分岐耐性を上げるため、通常はリード長やリピート長に応じて動的に決定
    const size_t k = (kmer_sz > 0) ? static_cast<size_t>(kmer_sz) : 31;
    if (k < 1) return;

    // 2. De Bruijn Graphの構築
    // 要件3: "read.rank == 'u' のみからつくる"
    std::unordered_map<std::string, DBGNode> dbg;
    std::unordered_map<std::string, int> in_degrees;

    
    std::unordered_map<std::string, int> kmer_counts;
    for (const auto& read : reads) {
        if (!read.qc_passed || read.rank != 'y') continue;
        //if (read.has_smaller_change || read.has_anti_pattern) continue;
        //if (!read.has_local_clip && !read.is_stable_non_ref) continue;
        if (read.seq.size() < k) continue;
        
        for (size_t i = 0; i <= read.seq.size() - k; ++i) {
            kmer_counts[std::string(read.seq.substr(i, k))]++;
        }
    }
    
    const int min_coverage = 2; 
    
    for (const auto& read : reads) {
        if (!read.qc_passed || read.rank != 'y') continue;
        //if (read.has_smaller_change || read.has_anti_pattern) continue;
        //if (!read.has_local_clip && !read.is_stable_non_ref) continue;
        if (read.seq.size() < k) continue;

        // クリップ部分やインデル部分を網羅するため、リード全体のウィンドウをスライド
        for (size_t i = 0; i <= read.seq.size() - k; ++i) {
            std::string kmer{read.seq.substr(i, k)};
            
            if (kmer_counts[kmer] < min_coverage) continue;

            auto& node = dbg[kmer];
            if (node.kmer.empty()) {
                node.kmer = kmer;
            }
            node.coverage = kmer_counts[kmer];

            // 次のK-merへのエッジを張る
            if (i < read.seq.size() - k) {
                std::string next_kmer{read.seq.substr(i + 1, k)};
                // 重複エッジの登録を防ぐ
                if (kmer_counts[next_kmer] >= min_coverage){
                    if (std::find(node.out_edges.begin(), node.out_edges.end(), next_kmer) == node.out_edges.end()) {
                        node.out_edges.push_back(next_kmer);
                        in_degrees[next_kmer]++;
                    }
                }
            }
        }
    }

    if (dbg.empty()) return;

    // 3. グラフのソースノード（入次数0、またはカバレッジ最大のノード）を起点に決定
    std::string start_kmer;
    int max_src_cov = -1;
    for (const auto& [kmer, node] : dbg) {
        if (in_degrees[kmer] == 0 && node.coverage > max_src_cov) {
            max_src_cov = node.coverage;
            start_kmer = kmer;
        }
    }
    // 入次数0が見つからない（ループ構造など）場合は、最大カバレッジのノードを起点にする
    if (start_kmer.empty()) {
        for (const auto& [kmer, node] : dbg) {
            if (node.coverage > max_src_cov) {
                max_src_cov = node.coverage;
                start_kmer = kmer;
            }
        }
    }

    // 4. アセンブリパスの探索（最有力コンティグの生成）
    std::unordered_set<std::string> visited_kmers;
    std::string current_path = start_kmer;
    std::string assembled_seq = start_kmer;

    visited_kmers.insert(start_kmer);
    dfs_find_longest_path(start_kmer, dbg, visited_kmers, current_path, assembled_seq);

    // 長すぎる、または短すぎるアセンブリ配列は異常値として除外
    if (assembled_seq.size() < k + 10) return;

    // 要件4: "assembly したtarget を含む配列をseq0に代入"
    seq0 = assembled_seq;
    i2p0.clear();
    i2p0.resize(seq0.size(), -1); // 初期値としてゲノム座標未定義（-1）を設定

    // 5. アセンブリ配列（seq0）をリファレンス配列（rseq）に対してアラインメントし、座標変換を計算
    // 要件4: "genomic positionはi2pをしてかえす"
    // ここでは、Smith-Waterman/Needleman-Wunsch に基づく内部アライナー（Aligner）を利用
    Aligner dgb_aligner(params.match_score, params.mismatch_penal,
                        params.gap_open_penal, params.gap_ext_penal);

    // リファレンス配列（rseq: クラス内で定義・抽出済み）をセット
    dgb_aligner.SetReferenceSequence(rseq.c_str(), rseq.size());

    Alignment aln_res;
    Filter filter_dummy;
    int32_t mask_len = (seq0.size() < 30) ? 15 : static_cast<int32_t>(seq0.size() / 2);

    dgb_aligner.Align(seq0.c_str(), filter_dummy, &aln_res, mask_len);

    // CIGAR文字列およびアラインメント開始位置をパースして i2p0 マップを構築
    // start は Pileup::start (ゲノム上のリファレンス座標)
    int ref_offset = aln_res.ref_begin; // rseq内でのローカルオフセット
    int query_idx = aln_res.query_begin;

    // CIGAR解析用の簡易ステートマシン
    // 例: "10M2I5D20M" などのオペレーションを処理
    const std::string& cigar = aln_res.cigar_string;
    size_t cigar_len = cigar.size();
    size_t c_idx = 0;

    while (c_idx < cigar_len) {
        int length = 0;
        while (c_idx < cigar_len && std::isdigit(cigar[c_idx])) {
            length = length * 10 + (cigar[c_idx] - '0');
            c_idx++;
        }
        if (c_idx >= cigar_len) break;
        char op = cigar[c_idx++];

        switch (op) {
            case 'M': // Match or Mismatch
            case 'X':
            case 'A':
                for (int l = 0; l < length; ++l) {
                    if (query_idx < static_cast<int>(i2p0.size())) {
                        // ゲノム上の絶対座標 = クラス全体のベースゲノム座標(start) + リファレンス内相対オフセット
                        i2p0[query_idx] = start + ref_offset;
                    }
                    query_idx++;
                    ref_offset++;
                }
                break;
            case 'I': // Insertion to Reference (Queryには塩基があり、Refにはない)
                for (int l = 0; l < length; ++l) {
                    if (query_idx < static_cast<int>(i2p0.size())) {
                        // 挿入セグメントの内部座標。周辺のゲノム境界を割り振るか、
                        // インデル検出用に特定のアンカー座標（直前のリファレンス座標等）をマッピング
                        i2p0[query_idx] = start + ref_offset;
                    }
                    query_idx++;
                }
                break;
            case 'D': // Deletion from Reference (Queryには塩基がなく、Refにある)
                ref_offset += length;
                break;
            case 'S': // Soft-clip
                query_idx += length;
                break;
            default:
                break;
        }
    }

    // 6. アセンブリされたターゲット候補配列（seq0）に対して、'u' リードがサポートしているか再評価
    // 要件1, 2の解決: Complex Indel や長いソフトクリップを持つリードを、DBGによる共通コンティグを介して救い出す
    int target_left_lim = target.pos - target.event_radius;
    int target_right_lim = target.pos + target.event_radius;
    
    /*
    bool contig_contains_target = false;
    
    // コンティグの各塩基のゲノム座標（i2p0）から、target.pos に対応するインデックスを探す
    auto it_pos = std::find(i2p0.begin(), i2p0.end(), target.pos);
    
    if (it_pos != i2p0.end()) {
        size_t local_idx = std::distance(i2p0.begin(), it_pos);
        
        // 変異（Indel含む）の長さを考慮してコンティグから実際の配列を切り出す
        // SNVであれば1文字、Indelであればその長さ分
        size_t check_len = target.alt.size();
        if (local_idx + check_len <= seq0.size()) {
            std::string assembled_allele = seq0.substr(local_idx, check_len);
            
            // コンティグの該当座標の配列が、target.alt または target.ref と一致するかを厳密に判定
            if (assembled_allele == target.alt || assembled_allele == target.ref) {
                contig_contains_target = true;
            }
        }
    }*/

    //if (!contig_contains_target) { seq0.clear(); i2p0.clear(); return; }

    // コンティグ内にターゲットバリアントが存在しているかを検出する（簡易シグネチャ探索、またはアラインメント境界確認）
    // ここでは、アセンブリ配列自身がtarget indelをカバーしているかを検証
   // bool contig_contains_target = false;

    // targetのref/alt配列パターンがseq0に見出せるか、座標系がtarget.posに交差しているかを判定
    //if (seq0.find(target.alt) != std::string::npos || seq0.find(target.ref) != std::string::npos) {
    //    contig_contains_target = true;
    //}

    //if (!contig_contains_target) {
        // もし文字列一致で不十分な場合、座標空間（i2p0）にtarget.pos周辺が含まれているかでフォールバック
     //   auto it_pos = std::find(i2p0.begin(), i2p0.end(), target.pos);
      //  if (it_pos != i2p0.end()) {
       //     contig_contains_target = true;
        //}
    //}

    // コンティグがターゲットを正しく回収できていれば、このコンティグ（パス）を構成した 'u' リードを 's' へ昇格
    //if (contig_contains_target) {
    Aligner read_aligner(params.match_score, params.mismatch_penal,
                         params.gap_open_penal, params.gap_ext_penal);
        // 今回アセンブリされた高信頼性配列（seq0）を新たなリファレンスとしてセット
    read_aligner.SetReferenceSequence(seq0.c_str(), seq0.size());

    std::bitset<3> check_points;
    int fss = -1, fee = -1, ts = -1, te = -1;

        // i2p0 マップから各アプローチウィンドウ（Flanking / Target）のローカル座標を決定
    for (int i = 0; i < static_cast<int>(i2p0.size()); ++i) {
        if (i2p0[i] == loc_ref.flanking_start) fss = i;
        if (ts < 0 && i2p0[i] == target.pos) ts = i;
        if (i2p0[i] == target.end_pos) te = i;
        if (i2p0[i] == loc_ref.flanking_end) fee = i;
    }

        // 座標定義が有効な場合のみ、各リードのリアラインメント・チェックをおこなう
    if (fss != -1 && ts != -1 && te != -1 && fee != -1) {
        int fse = fss + params.dimer_window;
        int fes = fee - params.dimer_window;

        for (auto& read : reads) {
            if (!read.qc_passed || read.rank != 'y') continue;
            //if (read.has_smaller_change || read.has_anti_pattern) continue;
            //if (!read.has_local_clip && !read.is_stable_non_ref) continue;
            if (read.seq.size() < k) continue;
            
            Alignment read_aln;
            int32_t r_mask_len = (read.seq.size() < 30) ? 15 : static_cast<int32_t>(read.seq.size() / 2);

            // リードを seq0 コンティグに対してアラインメント（ソフトクリップの解消）
            read_aligner.Align(read.seq.data(), filter_dummy, &read_aln, r_mask_len);

            // match.h の外部ヘルパー関数でターゲットに合致するかチェック
            check_match_pattern(read_aln, check_points, fss, fse, ts, te, fes, fee);

            // 完全一致、または複雑なIndelの一部として高スコアでマッチした場合
            //if (check_points.count() == 3 || (check_points.test(1) && read.has_positional_overlap)) {
            if (check_points.count() == 3 ) {
                read.rank = 's'; // 確定支持リードへ昇格
                ++s_cnt; --y_cnt;

                // シグネチャマップの更新
                sig_s[read.non_ref_sig].push_back(&read - &reads[0]);
                sig_s_hiconf[read.non_ref_sig].push_back(&read - &reads[0]);
                has_hiconf_support = true;
            }
            check_points.reset();
        }
   }
   
   //if (!s_cnt) { seq0.clear(); i2p0.clear(); } 
}


