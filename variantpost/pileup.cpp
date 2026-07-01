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
        
        read.is_quality_map = (read.is_stable_non_ref && read.is_central_mapped);
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


void Pileup::inferGermlineHaplotype(const UserParams& params) {
    if (this->sz == 0) return;

    homo_vars.reserve(this->sz / 10); 
    hap1_vars.reserve(this->sz / 10);
    hap2_vars.reserve(this->sz / 10);
            
    struct CountTrack {
        int total_reads = 0;
        int total_depth = 0;
    };
    std::unordered_map<VariantKey, CountTrack, VariantKeyHash> global_counts;

    for (int i = 0; i < this->sz; ++i) {
        const auto& read = this->reads[i];
        if (!read.is_control) continue;
        if (read.covering_ptrn != 'A' || !read.qc_passed || !read.is_central_mapped) continue;
        
        for (const auto& v : read.variants) {
            if (v.mean_qual < params.base_q_thresh) continue;
            global_counts[VariantKey{v.pos, v.ref, v.alt}].total_reads++;
        }
    }

    if (global_counts.empty()) return;

    //--- Sort variant postions to calculate VAF
    std::vector<long long> unique_positions;
    for (const auto& [key, _] : global_counts) {
        unique_positions.push_back(key.pos);
    }
    std::sort(unique_positions.begin(), unique_positions.end());
    unique_positions.erase(std::unique(unique_positions.begin(), unique_positions.end()), unique_positions.end());
    
    //--- Calculate depths at sorted variant positions
    std::unordered_map<long long, int> pos_depth;
    
    for (int i = 0; i < this->sz; ++i) {
        const auto& read = this->reads[i];
        
        if (!read.is_control) continue;
         
        //--- To include ref reads. This is because that ref reads are not annotated for most flags
        if (!read.is_ref) {
            if (read.covering_ptrn != 'A' || !read.qc_passed || !read.is_central_mapped) continue;
        } 
        
        //--- Count only aligned segments (exclude spliced skips)
        for (const auto& seg : read.aligned_segments) {
            int seg_start = seg.first;
            int seg_end = seg.second;
            auto it = std::lower_bound(unique_positions.begin(), unique_positions.end(), seg_start);
            while (it != unique_positions.end() && *it <= seg_end) {
                pos_depth[*it]++;
                ++it;
            }
        }
    }

    //--- Map the depth to global count
    for (auto& [key, track] : global_counts) {
        track.total_depth = pos_depth[key.pos];
    }

    //--- Profile positions for the presence of variants (any)
    std::unordered_map<long long, std::vector<std::pair<VariantKey, int>>> pos_to_vars;
    for (const auto& [key, track] : global_counts) {
        if (track.total_depth == 0 || track.total_reads < 2) continue;
        pos_to_vars[key.pos].push_back({key, track.total_reads});
    }
    
    //--- Het/Hom classification
    std::vector<VariantKey> backbone_keys;
    backbone_keys.reserve(global_counts.size());
    
    for (auto& [pos, vars] : pos_to_vars) {
        // Coverage at the locus (pos)
        int total_depth = global_counts[vars[0].first].total_depth;
        
        // Sort variants with N of supporting reads (total_reads)
        std::sort(vars.begin(), vars.end(), [](const auto& a, const auto& b) {
            return a.second > b.second; // 降順
        });

        const auto& major_key = vars[0].first;
        int major_reads = vars[0].second;

        // Total count of NON-reference (alt) reads
        int total_alt_reads = 0;
        for (const auto& v : vars) total_alt_reads += v.second;

        // Get the reference read counts as total_depth(coverage) - total_alt-reads
        int ref_reads = total_depth - total_alt_reads;
        if (ref_reads < 0) ref_reads = 0; //zero-clip for safety

        double ref_vaf = static_cast<double>(ref_reads) / total_depth;
        double major_vaf = static_cast<double>(major_reads) / total_depth;

        bool is_snv = (major_key.ref.length() == 1 && major_key.alt.length() == 1);
        std::cout << major_key.pos << " " << major_key.ref << " " << major_key.alt 
                  << " Major_VAF:" << major_vaf << " Ref_VAF:" << ref_vaf << std::endl;

        // Homozygous if ref reads is so few (< 5%) 
        // or the locus has a dominant variant (SNV: vaf > 0.9 or indel: vaf > 0.7                    
        if (ref_vaf < 0.05 || ((is_snv && major_vaf > 0.9) || (!is_snv && major_vaf > 0.7))) {
            Variant v(static_cast<int>(major_key.pos), major_key.ref, major_key.alt, "");
            homo_vars.push_back(v); 
        }
        // Othewise Het, but require a vaf > 0.2 to remove technical artifacts
        else if (major_vaf >= 0.2) {
            backbone_keys.push_back(major_key);
        }
    }

    std::sort(homo_vars.begin(), homo_vars.end());
       
    std::sort(backbone_keys.begin(), backbone_keys.end(), [](const VariantKey& a, const VariantKey& b) {
        return a.pos < b.pos;
    });

    int N = backbone_keys.size();
    if (N == 0) return; // no het case 

    struct Edge { int same = 0; int diff = 0; };
    std::vector<std::unordered_map<int, Edge>> adj(N);

    std::vector<std::pair<int, bool>> read_calls;
    read_calls.reserve(64);

    for (int i = 0; i < this->sz; ++i) {
        const auto& read = this->reads[i];
        if (!read.is_control || read.covering_ptrn != 'A' || !read.qc_passed) continue;

        read_calls.clear();

        for (int j = 0; j < N; ++j) {
            const auto& b_key = backbone_keys[j];

            if (read.aln_end < b_key.pos) break;

            if (b_key.pos < read.aln_start) continue;


            bool is_mechanically_covered = false;
            for (const auto& seg : read.mapped_segments) {
                if (seg.first <= b_key.pos && b_key.pos <= seg.second) {
                    is_mechanically_covered = true;
                    break;
                }
            }

            if (!is_mechanically_covered) continue;

            bool has_alt = false;
            for (const auto& v : read.variants) {
                if (v.pos == b_key.pos && v.ref == b_key.ref && v.alt == b_key.alt) {
                    has_alt = true; // ALT
                    break;
                }
            }

            read_calls.push_back({j, has_alt});
        }

        // --- Within the same read, testing the linkage
        if (read_calls.size() < 2) continue;

        for (size_t a = 0; a < read_calls.size(); ++a) {
            for (size_t b = a + 1; b < read_calls.size(); ++b) {
                int idxA = read_calls[a].first;
                int idxB = read_calls[b].first;


                if (read_calls[a].second == read_calls[b].second) {
                    adj[idxA][idxB].same++;
                    adj[idxB][idxA].same++;
                } else {
                    adj[idxA][idxB].diff++;
                    adj[idxB][idxA].diff++;
                }
            }
       }
    }

    // Pick up the most reliable HET variant as the starting point
    // (Maximize coverage while penalizing deviation from VAF 0.5)
    int best_seed_idx = 0;
    double best_seed_score = -1.0;
    for (int i = 0; i < N; ++i) {
        int coverage = global_counts[backbone_keys[i]].total_reads;
        int depth = global_counts[backbone_keys[i]].total_depth;
        if (depth == 0) continue;

        double vaf = static_cast<double>(coverage) / depth;
        double score = coverage * (0.5 - std::abs(vaf - 0.5));

        if (score > best_seed_score) {
            best_seed_score = score;
            best_seed_idx = i;
        }
    }

    std::vector<Phase> phased_array(N, Phase::Unphased);
    phased_array[best_seed_idx] = Phase::Hap1; // designate the start as Hap1

    std::vector<int> queue;
    queue.reserve(N);
    queue.push_back(best_seed_idx);
    size_t q_head = 0;

    while (q_head < queue.size()) {
        int curr = queue[q_head++];
        Phase curr_phase = phased_array[curr];

        for (int neighbor = 0; neighbor < N; ++neighbor) {
            if (phased_array[neighbor] != Phase::Unphased) continue;

            const Edge& edge = adj[curr][neighbor];
            int total_links = edge.same + edge.diff;

            if (total_links >= 2) {
                double purity = static_cast<double>(std::max(edge.same, edge.diff)) / total_links;

                // same /(same + diff) > 0.75 or < 0.25
                if (purity >= 0.75) {
                    if (edge.same > edge.diff) {
                        phased_array[neighbor] = curr_phase;
                        queue.push_back(neighbor);
                    } else {
                        phased_array[neighbor] = (curr_phase == Phase::Hap1) ? Phase::Hap2 : Phase::Hap1;
                        queue.push_back(neighbor);
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N; ++i) {
        auto& key = backbone_keys[i];
        auto& phase = phased_array[i];        
        if (phase == Phase::Hap1) {
                hap1_vars.emplace_back(static_cast<int>(key.pos), key.ref, key.alt, "");
        } else if (phase == Phase::Hap2) {
                hap2_vars.emplace_back(static_cast<int>(key.pos), key.ref, key.alt, "");
        }
    }
    std::sort(hap1_vars.begin(), hap1_vars.end());
    std::sort(hap2_vars.begin(), hap2_vars.end());
}

//-----------------------------------------------------------------------------
// This method tends to overanalyze 
// -> tighter conditions to apply
// Should be run when targe is suspected to be a part of complex indels 
void Pileup::gridSearch(const UserParams& params, 
                        LocalReference& loc_ref, const Variant& target) {
    
    // Condition 0
    // Skip if target is homopolymer extending 
    if (s_cnt || sig_u.empty() || target.denovo_rep > 2) { return;}

    PatternCnt u_sig_cnts;
    count_patterns(sig_u, u_sig_cnts);
    
    const int tot = s_cnt + n_cnt + u_cnt; 
    // Condition 1
    // Skip if background hap involves homoplymer variations
    //   Step1: check for homopolymer involvment (denovo_rep)
    //   Step2: define as background hap if freq is > 0.1
    for (const auto& elem : ns_vars) {
        if (elem.first.denovo_rep > 1 && elem.second * 10 > tot) { return;}
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


inline void collect_target_hap_variants(const Reads& reads, Vars& vars) {
    size_t total_size = 0;
    for (const auto& read : reads) total_size += read.variants.size();
        vars.reserve(total_size);
    
    for (const auto& read : reads) {
        if (read.rank == 's') 
            vars.insert(vars.end(), read.variants.begin(), read.variants.end());
    }

    std::sort(vars.begin(), vars.end());

    auto last = std::unique(vars.begin(), vars.end());
    vars.erase(last, vars.end());
}

inline void prep_vars(Vars& dest, const Vars& source) {
    for (const auto& v : source) dest.push_back(v);
    std::sort(dest.begin(), dest.end());
}


//------------------------------------------------------------------------------
void Pileup::setHaploTypes(LocalReference& loc_ref, const Variant& target) {
    PatternCnt s_sig_cnt; 
    
    if (hiconf_read_idx > -1) {
        std::cout << " hi gonf " << hiconf_read_idx << " " << reads[hiconf_read_idx].cigar_str << std::endl;
        make_sequence(
            loc_ref, reads[hiconf_read_idx].variants, start, end, seq0, &i2p0);
        std::cout << " hap zero " << seq0 << std::endl;
    } else if (!sig_s.empty()) {
        PatternCnt s_sig_cnt;
        count_patterns(sig_s, s_sig_cnt); 
        hap0 = s_sig_cnt[0].first;
        make_sequence(
            loc_ref, reads[sig_s[hap0][0]].variants, start, end, seq0, &i2p0);
    } else {
        if (homo_vars.empty()) // NOTE: homo_vars empty if has_second_bam is false
            make_sequence(loc_ref, {target}, start, end, seq0, &i2p0);
    }
     
    if (has_second_bam) {
        // keep variants on target haplotype
        if (hiconf_read_idx > -1)  
            collect_target_hap_variants(reads, hap0_vars); 

        // personalize target with homozygous if unset
        if (seq0.empty()) {
            Vars hap0_full = {target};
            prep_vars(hap0_full, homo_vars); 
            make_sequence(loc_ref, hap0_full, start, end, seq0, &i2p0);
        }
        for (const auto& v : homo_vars) {
            std::cout << v.pos << " " << v.ref << " " << v.alt << " homo var " << std::endl;
        }


        Vars hap1_full = homo_vars;
        prep_vars(hap1_full, hap1_vars); 
        for (const auto& v : hap1_vars)
            std::cout << v.pos << " " << v.ref << " " << v.alt << " hap1 var " << std::endl;
         
        Vars hap2_full = homo_vars;
        prep_vars(hap2_full, hap2_vars);
        for (const auto& v : hap2_vars)
            std::cout << v.pos << " " << v.ref << " " << v.alt << " hap2 var " << std::endl;
        
        // REF/REF case
        if (hap1_full.empty() && hap2_full.empty()) {
            is_ref_hom = true; return;
        }

        // NOTE: hap1_vars empty but hap2_vars Non-empty never occurs
        if (!hap1_full.empty()) {
            make_sequence(loc_ref, hap1_full, start, end, seq1, &i2p1);    
            std::cout << seq1 << " hap 111 " << std::endl;
        }
        
        // REF/non_REF case
        if (has_ref_hap) return; 
        
        // homozygous for non_REF
        if (hap1_full == hap2_full) {
            is_alt_hom = true; return;
        }

        // heterozygous for non_REF
        // this may include case hap1_vars = {homo_snp1, het_snp1} hap2_vars = {homo_snp1} 
        if (!hap2_full.empty())
            make_sequence(loc_ref, hap2_full, start, end, seq2, &i2p2); 
        
        // Exit here.
        // We only set the background haplotypes if control alignments provided  
        return;
    }
    
    PatternCnt u_sig_cnt; int idx1 = -1, idx2 = -1;

    if (!sig_u.empty()) {
        count_patterns(sig_u, u_sig_cnt);
        hap1 = u_sig_cnt[0].first; 
        if (hap_settable(s_cnt, hap1, anti_sigs)) {
            idx1 = sig_u[hap1][0]; // Not -1 
            make_sequence(loc_ref, reads[idx1].variants, start, end, seq1, &i2p1);
            if (seq0 == seq1) seq1.clear();
        }
        
        if (u_sig_cnt.size() > 1) {
            hap2 = u_sig_cnt[1].first;
            if (hap_settable(s_cnt, hap2, anti_sigs)) {
                idx2 = sig_u[hap2][0]; // Not -1
                make_sequence(loc_ref, reads[idx2].variants, start, end, seq2, &i2p2);
                if (seq0 == seq2 || seq1 == seq2 ) seq2.clear();
            }
        }
    } else if (has_excess_ins_hap) {
    // Haplotype with additional ins-repeats may be clipped and may not be captured    
    // TODO: add example or remove this part -> too complex look so ad-hoc
        std::string alt_added = target.alt + target.alt.substr(1);
        Vars vlst = {Variant(target.pos, target.ref, alt_added)};
        make_sequence(loc_ref, vlst, start, end, seq1, &i2p1); idx1 = 0; //pseudo index
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

        std::cout << read.smer << " " << read.nmer << " " << read.cigar_str << " " << read.rank << std::endl; 
        if (non_target_kmer_sz && read.smer > read.nmer) { 
            read.rank = 'y'; --u_cnt; ++y_cnt; //Likely supporting
            std::cout << read.name <<  " y nct " << std::endl;
            if (!read.nmer && read.smer > 1) {
                if (has_hiconf_support 
                    || no_non_target_haps 
                    || read.is_quality_map
                    || target.indel_len > 19) {
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
        if (rank_ == 's' || rank_ == 'y' || rank_ == 'n') continue; 
        const bool has_it = has_the_sig(read.non_ref_sig, anti_sig);
        if (has_it){
            read.rank = 'n'; ++n_cnt;
            if (rank_ == 'u') --u_cnt;
            //else if (rank_ == 'y') --y_cnt;
            else --z_cnt;
        }
    }

    has_likely_support = (y_cnt);
}

//------------------------------------------------------------------------------
//NOTE THIS STEP DOES NOT MAKE SENSE
//WHY DO THIS HERE? 
//REALN FOR READS, NOT RE_CONSTRUSTED SEQ FROM U_SIG, MAY BE OK.
void Pileup::searchByRealignment(const UserParams& params,
                                 LocalReference& loc_ref, const Variant& target) {
    //if (s_cnt || !u_cnt || seq0.empty()) return;
    //if (s_cnt || seq0.empty()) return;
    
    int fss = -1, fse = -1, fes = -1, fee = -1, ts = -1, te = -1;  
    
    
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
    
    aligner.SetReferenceSequence(seq0.c_str(), seq0.size());
    
    std::bitset<3> check_points;
    for (auto& read : reads) {
        if (read.rank == 's' || read.rank == 'n' || !read.qc_passed) continue;
        if (read.rank != 'y' && !read.covered_in_clip && !read.has_local_clip) continue;
        //if (read.rank != 'y' && !read.covered_in_clip) continue;
        
        std::string ss{read.seq}; // may be improved 

        int32_t mask_len = strlen(ss.c_str()) < 30 ? 15 : strlen(ss.c_str()) / 2;
        aligner.Align(ss.c_str(), filter, &aln, mask_len);
        check_match_pattern(aln, check_points, fss, fse, ts, te, fes, fee); 
        
        std::cout << read.cigar_str << " " << read.smer << " " << read.nmer << " " << check_points.count() << std::endl;
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
}

