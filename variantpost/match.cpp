#include "match.h"
#include "search.h"
#include "hap_likelihood.h"

void aln2variants(Alignment& aln, Vars& vars, const int start,
                  const std::string& ref, const std::string& query) {
    if (!aln.cigar_string.size()) return;
    std::cout << aln.cigar_string << " " << aln.ref_begin << std::endl;
    std::cout << ref << std::endl;
    std::cout << query << std::endl;
    CigarVec cigar_vec; fill_cigar_vector(aln.cigar, cigar_vec);
    
    char op = '\0'; int pos = start;
    int op_len = 0, ri = aln.ref_begin, qi = aln.query_begin;
    for (const auto& cigar : cigar_vec) {
        op = cigar.first; op_len = cigar.second;
        switch (op) {
            case '=':
                ri += op_len; qi += op_len; pos += op_len;
                break;
            case 'X':
                vars.emplace_back(pos, 
                                  ref.substr(ri, op_len),
                                  query.substr(qi, op_len));
                ri += op_len; qi += op_len; pos += op_len;
                break;
            case 'D':
                if (ri) {
                    vars.emplace_back(pos - 1,
                                     ref.substr(ri - 1, op_len + 1),
                                     ref.substr(ri - 1, 1));
                    ri += op_len; pos += op_len;
                } break;
            case 'I':
                if (ri) {
                    std::string base = ref.substr(ri - 1, 1);
                    vars.emplace_back(pos - 1,
                                     ref.substr(ri - 1, 1),
                                     base.append(query.substr(qi, op_len)));
                    qi += op_len;
                } break;
            case 'S':
                break;
        }
    }
} 

//-----------------------------------------------------------------------------
// fss: flanking start
// fse: flanking end: flanking start + linguistic cplx search window
// ts: target start
// te: target end after right alignment
// fes: flanking end start: flanking end -  linguistic cplx search window
// fee: flankine end end
void check_match_pattern(Alignment& aln, 
                         std::bitset<3>& check_points,
                         const int fss, const int fse, 
                         const int ts, const int te,
                         const int fes, const int fee) {
    if (!aln.cigar_string.size()) return;

    CigarVec cigar_vec; fill_cigar_vector(aln.cigar, cigar_vec);
    char op = '\0'; int op_len = 0, idx = aln.ref_begin;
    for (const auto& cigar : cigar_vec) {
        op = cigar.first; op_len = cigar.second;
        switch (op) {
            case '=':
                if (idx <= fss && fse <= idx + op_len) check_points.set(0);
                if (idx <= ts && te <= idx + op_len) check_points.set(1); 
                if (idx <= fes && fee <= idx + op_len) check_points.set(2);
                idx += op_len; break;
            case 'X': case 'D':
                idx += op_len; break;
            case 'I': case 'S':
                break;
        }
    }
    //check_points fully set -> all matched. 
}

//------------------------------------------------------------------------------
// uint8_t (penal) -> 255 max
void gap_grid(const UserParams& params, std::vector<Ints>& grid) {
    int max_go = (params.gap_open_penal * 2 <= 255) ? params.gap_open_penal * 2 : 255;
    int max_ge = (params.gap_ext_penal * 2 <= 255) ? params.gap_ext_penal * 2 : 255; 
    for (int go = 1; go <= max_go; go +=2) 
        for (int ge = 0; ge <= max_ge; ++ge) {
            grid.push_back({params.match_score, params.mismatch_penal, go, ge}); 
            grid.push_back({params.match_score, params.mismatch_penal * 3, go, ge});
        } 
}

bool search_over_grid(const int start, LocalReference& loc_ref,
                      const UserParams& params, const size_t n_vars,
                      const std::string& refseq, const std::string& query, 
                      const std::vector<Ints>& grid, const Variant& target) {
    Alignment aln; Filter filter;
    auto qseq = query.c_str();
    Vars vars;
    int32_t mask_len = strlen(qseq) < 30 ? 15 : strlen(qseq) / 2;
    for (const auto& g : grid) {
        Aligner aligner(g[0], g[1], g[2], g[3]);
        aligner.SetReferenceSequence(refseq.c_str(), refseq.size());
        aligner.Align(qseq, filter, &aln, mask_len);
        aln2variants(aln, vars, start, refseq, query);
        std::cout << std::endl;
        for(auto& v : vars) {
            std::cout << v.pos << " " << v.ref << " " << v.alt << "- ";
        }
        if (vars.size() > n_vars) return false;
        std::cout << " total " << vars.size() << std::endl;
        if (find_target(loc_ref, params, target, vars) > -1) return true;
        vars.clear(); 
    }
    return false;    
}


//------------------------------------------------------------------------------
void match2haplotypes(Pileup& pileup, const Strs& read_seqs, const UserParams& params) {
    std::bitset<4> valid_idx;
    if (pileup.seq0.size()) { valid_idx.set(0); } if (pileup.seq1.size()) { valid_idx.set(1); }
    if (pileup.seq2.size()) { valid_idx.set(2); } if (pileup.rseq.size()) { valid_idx.set(3); }
    
    Alignment aln; Filter filter;
    Aligner aligner(params.match_score, params.mismatch_penal,
                    params.gap_open_penal, params.gap_ext_penal);
     
    for (int i = 0; i < pileup.sz; ++i) {
        if (!pileup.reads[i].qc_passed || pileup.reads[i].rank != Rank::Undetermined || pileup.reads[i].smer == 0) continue; 
        
        const auto& query = read_seqs[i].c_str();
        int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;  
        std::vector<int> scores = {0, 0, 0, 0};
        for (int j = 0; j < 4; ++j) {
            if (valid_idx.test(j)) {
                if (j == 0) {
                    aligner.SetReferenceSequence(pileup.seq0.c_str(), pileup.seq0.size());
                    std::cout << j << " " << pileup.seq0 << std::endl;
                } else if (j == 1) {
                    aligner.SetReferenceSequence(pileup.seq1.c_str(), pileup.seq1.size());
                    std::cout << j << " " << pileup.seq1 << std::endl;
                } else if (j == 2) {
                    aligner.SetReferenceSequence(pileup.seq2.c_str(), pileup.seq2.size());
                    std::cout << j << " " << pileup.seq2 << std::endl;
                } else {
                    aligner.SetReferenceSequence(pileup.rseq.c_str(), pileup.rseq.size());
                    std::cout << j << " " << pileup.rseq << std::endl;
                }
                aligner.Align(query, filter, &aln, mask_len); scores[j] = aln.sw_score;
            } 
        }
         
        // aln against target hap is solely highest
        if (scores[0] > scores[1] && scores[0] > scores[2] && scores[0] > scores[3]) {
            pileup.reads[i].rank = Rank::Supporting; ++pileup.s_cnt; --pileup.u_cnt; 
        } else if (scores[0] == scores[1] || scores[0] == scores[2] || scores[0] == scores[3]) {
            // tie score to one of non-target haplotypes -> remain as 'u'
        } else {
            pileup.reads[i].rank = Rank::NotSupporting; ++pileup.n_cnt; --pileup.u_cnt;
        }
        std::cout <<  pileup.reads[i].rank << std::endl;
    }
}


inline size_t count_overlap(const Vars& v1, const Vars& v2) {
    Vars intersection;
    intersection.reserve(std::min(v1.size(), v2.size()));
    std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(intersection));

    return intersection.size();
}

inline void fill_matched_sides(size_t curr, 
                               const size_t last, 
                               const CigarVec& cigar_vec,
                               std::vector<std::pair<int, int>>& matched_sides) {
    if (curr == 0 && last == 0) return;

    int lt_matched_len = 0, rt_matched_len = 0;
    
    char prev_op = '\0', next_op = '\0';
    if (curr > 0) {
        prev_op = cigar_vec[curr - 1].first;
        if (prev_op == '=')
            lt_matched_len = cigar_vec[curr - 1].second;
    }
    
    if (curr < last) {
        next_op = cigar_vec[curr + 1].first;
        if (next_op == '=')
            rt_matched_len = cigar_vec[curr + 1].second;
    } 
    
    matched_sides.emplace_back(lt_matched_len, rt_matched_len); 
}

inline bool is_in_flankings(const int pos, const int flanking_start, const int flanking_end) {
    return (flanking_start <= pos && pos <= flanking_end); 
}

struct ReAlnQc {
     bool alignable_as_snvs = false;
     bool too_many_events = false;
     bool out_of_flanking = false;
     bool is_likely_simple = false;
     int start = 0, end = 0; // flanking start/end pos
     int expected_event_num = 0;
     int local_thresh = 0;

     ReAlnQc(int a, int b, int c, int d) 
        : start(a), end(b), expected_event_num(c), local_thresh(d) { }

     bool pass() const {
        return (!too_many_events && !out_of_flanking);
     }
};


void realn_to_perfonalized_genome(
    Aligner& alngr, Filter& fltr, Alignment& aln, Variant& pv, ReAlnQc& qc, 
    const std::string& mut_seq, const std::string& hap_seq, const Ints& i2p) {
    
    const char* query = mut_seq.c_str();
    int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;
    alngr.SetReferenceSequence(hap_seq.c_str(), hap_seq.size());
    alngr.Align(query, fltr, &aln, mask_len);
    
    std::cout << aln.cigar_string << " realn " << std::endl; 
    CigarVec cigar_vec;
    fill_cigar_vector(aln.cigar_string, cigar_vec);
    
    Vars vars; 
    std::vector<std::pair<int, int>> matched_sides; // To test for complex indels  
    
    int i = aln.ref_begin, j = aln.query_begin; 
    
    int cnt = 0; // count variants in flankings
    size_t curr = 0;
    const size_t last = cigar_vec.size() - 1; 
    bool has_indels_in_flankings = false;
    for (const auto& [op, op_len] : cigar_vec) {
        switch (op) {
            case '=':
                i += op_len; j += op_len; break;
            case 'X':
                vars.emplace_back(i2p[i], hap_seq.substr(i, op_len), mut_seq.substr(j, op_len));
                
                if (is_in_flankings(i2p[i], qc.start, qc.end)) ++cnt;
                
                fill_matched_sides(curr, last, cigar_vec, matched_sides);
                i += op_len; j+= op_len;    
                break;
            case 'I': {
                Variant v(i2p[i - 1], hap_seq.substr(i - 1, 1), mut_seq.substr(j - 1 , op_len +1));
                
                if (is_in_flankings(i2p[i - 1], qc.start, qc.end)) {
                     ++cnt; has_indels_in_flankings = true;
                }

                int lt_lim = std::max(i - 1 - 10, 0);
                v.sample_lt_seq = hap_seq.substr(lt_lim, i - lt_lim);
                v.sample_rt_seq = hap_seq.substr(i);
                vars.push_back(v);
                fill_matched_sides(curr, last, cigar_vec, matched_sides);
                j += op_len;
                break; 
            }
            case 'D': {
                Variant v(i2p[i - 1], hap_seq.substr(i - 1, op_len + 1), mut_seq.substr(j - 1 , 1));
                
                if (is_in_flankings(i2p[i - 1], qc.start, qc.end)) {
                    ++cnt; has_indels_in_flankings = true;
                }
                
                int lt_lim = std::max(i - 1 - 10, 0);
                v.sample_lt_seq = hap_seq.substr(lt_lim, i - lt_lim);
                v.sample_rt_seq = hap_seq.substr(i + op_len);
                vars.push_back(v);
                fill_matched_sides(curr, last, cigar_vec, matched_sides);
                i += op_len;
                break;
           }
            case 'S':
                //j += op_len;
                break;
        }
        ++curr;
    } 
    
    // Early exit if the alignment on personalized genome is possibly SNVs
    // QC will PASS
    if (!has_indels_in_flankings) {
       qc.alignable_as_snvs = true; return;
    }

    // Early eixt if the alingment on personalized genome is with indels but more complicated 
    // -> complicated variations are unlikely to occur 
    if (qc.expected_event_num == 1 && cnt > qc.expected_event_num) {
        qc.too_many_events = true; return;
    }

    int closest = -1, min_dist = INT_MAX, dist = INT_MAX;
    for (int i = 0; i < static_cast<int>(vars.size()); ++i) {
        const auto& v = vars[i]; 
        if (is_in_flankings(v.pos, qc.start, qc.end)) {
            // Skip if it is SNV, while there are indels in the flankings
            if (has_indels_in_flankings && (v.ref.size() == v.alt.size())) continue;

            dist = std::abs(v.pos - pv.pos); // target.pos is supplied by pv (initialized by target)
            if (dist < min_dist) {closest = i; min_dist = dist; } 
        }
    }
    
    // Early exit if no indels in the flanking 
    // -> Realignment on personalzed genome is too different from the original 
    if (closest == -1) {
        qc.out_of_flanking = true;  return; 
    }

    pv = std::move(vars[closest]); 
    std::cout << closest <<  " closet " << " " << pv.pos << " " << pv.ref << " " << pv.alt << " " << cnt << " v cnt on ref " << qc.expected_event_num  << std::endl; 
       
    // Complex indels on the personalized genome
    if (matched_sides[closest].first >= qc.local_thresh && matched_sides[closest].second >= qc.local_thresh) {
        // exclude MNVs 
        if (!pv.is_complex) qc.is_likely_simple = true;
    }
}

inline bool expect_to_share_hets(const int start, const int end, Vars& vars, LocalReference& loc_ref, const Variant& target) {
    const int target_start = target.lpos;
    const int target_end = target.end_pos; 
    
    for (auto& v : vars) {
        v.setEndPos(loc_ref);
        // v is on the target contig
        if (v.end_pos < start || end < v.lpos) continue;
        // v does not interact with target -> should be found if on the same hap
        if (v.end_pos < target_start || target_end < v.lpos) return true;
    }
    return false;
}

inline void fill_result(const Variant& v, const ReAlnQc& qc, SearchResult& rslt) {
    rslt.ppos = v.pos;
    rslt.pref = v.ref;
    rslt.palt = v.alt;
    rslt.pltseq = v.sample_lt_seq;
    rslt.prtseq = v.sample_rt_seq;    
    
    if (!qc.pass()) return;
    
    rslt.personalized = true; // Successful analysis meant here NOT allele changed.   
    rslt.possible_snv = qc.alignable_as_snvs;
    rslt.likely_simple_on_personalized_genome = qc.is_likely_simple; 
}

void personalize(Pileup& pileup, LocalReference& loc_ref, const UserParams& params, const Variant& target, SearchResult& rslt) {
        
    if (!pileup.has_hiconf_support || pileup.seq0.empty()) return;
    
    // REF/REF case -> no personalization
    if (pileup.is_ref_hom){ 
        rslt.personalized = true;  // background haplotype successfully identified as REF/REF
        return;                    // -> safe to phase to a complex event 
    }

    std::string hiconf_seq = pileup.seq0;
    
    // Heuristic for longer (>20bp) indels to align it as a signle gap
    int gap_open_penal = params.gap_open_penal;
    int gap_ext_penal =  params.gap_ext_penal;
    if (target.indel_len > 19) {
        gap_open_penal = 10; gap_ext_penal = 0;
    }
    Alignment aln; Filter filter;
    Aligner aligner(params.match_score, params.mismatch_penal,
                    gap_open_penal, gap_ext_penal);
    
    size_t with_hap1 = count_overlap(pileup.hap0_vars, pileup.hap1_vars);
     
    Variant pv1(target.pos, target.ref, target.alt);
    Variant pv2(target.pos, target.ref, target.alt);
    
    std::cout << pv1.ref << " " << pv1.alt << " this is target " << std::endl;
    if (pileup.is_alt_het)
        std::cout << "This is Het/Het" << std::endl;
        
    HapLL::RepeatInfo ri1, ri2, rir;
    ReAlnQc qc_hap1(loc_ref.flanking_start, loc_ref.flanking_end, pileup.v_cnt, params.local_thresh);
    ReAlnQc qc_hap2 = qc_hap1;
    
    bool hap1_inferred = false, hap2_inferred = false;
    if (pileup.is_alt_het) {    
        size_t with_hap2 = count_overlap(pileup.hap0_vars, pileup.hap2_vars);
        
        std::cout << "over lap this hap1 " << with_hap1 << " " <<  " over lap this hap2 " << with_hap2 << std::endl;
        if (with_hap1 > with_hap2) {
            realn_to_perfonalized_genome(aligner, filter, aln, pv1, qc_hap1, hiconf_seq, pileup.seq1, pileup.i2p1);
            if (qc_hap1.pass()) hap1_inferred = true;
        } else if (with_hap1 < with_hap2) {
            realn_to_perfonalized_genome(aligner, filter, aln, pv2, qc_hap2, hiconf_seq, pileup.seq2, pileup.i2p2);
            if (qc_hap2.pass()) hap2_inferred = true;
        } else {
            // Evaluate hap1 vs. hap2 by likelihood 
            realn_to_perfonalized_genome(aligner, filter, aln, pv1, qc_hap1, hiconf_seq, pileup.seq1, pileup.i2p1);
            realn_to_perfonalized_genome(aligner, filter, aln, pv2, qc_hap2, hiconf_seq, pileup.seq2, pileup.i2p2);
            
            if (qc_hap1.pass() && qc_hap2.pass()) {
                 if (!qc_hap1.alignable_as_snvs && !qc_hap2.alignable_as_snvs) {
                    double ll_hap1 = HapLL::evaluate_variant(pv1, ri1);
                    double ll_hap2 = HapLL::evaluate_variant(pv2, ri2);
                    if (ll_hap1 >= ll_hap2){ hap1_inferred = true; } else { hap2_inferred = true; } 
                 } else {
                    hap1_inferred = qc_hap1.alignable_as_snvs; hap2_inferred = qc_hap2.alignable_as_snvs;
                 }  
            } else {
                hap1_inferred = qc_hap1.pass(); hap2_inferred = qc_hap2.pass();     
            }
        }
    } else {
        // Test if the target hap would share variants if it is on hap 1
        bool would_share = expect_to_share_hets(pileup.start, pileup.end, pileup.hap1_vars, loc_ref, target);
        
        std::cout << " this is compareison against hap1 vs ref. something shared with hap1?? " << with_hap1 << std::endl;
        std::cout << " it would share, it target is on hap1? " << would_share << std::endl;
        
        // No sharing while it should -> ref hap
        if (would_share && with_hap1 == 0) return; 
        
        realn_to_perfonalized_genome(aligner, filter, aln, pv1, qc_hap1, hiconf_seq, pileup.seq1, pileup.i2p1);
        std::cout << qc_hap1.pass() << " <- qc passed for hap1?? " << std::endl;
        if (!qc_hap1.pass()) return;

        if (with_hap1 || qc_hap1.alignable_as_snvs) {
            hap1_inferred = true;
        } else {           
            // Evaluate hap1 vs. ref hap by likelihood
            Variant pv_ref = target;
            pv_ref.sample_lt_seq = std::string{target.lt_seq}; pv_ref.sample_rt_seq = std::string{target.rt_seq};
            double ll_hap1 = HapLL::evaluate_variant(pv1, ri1);
            double ll_ref  = HapLL::evaluate_variant(pv_ref, rir);
            if (ll_hap1 >= ll_ref) hap1_inferred = true;
        }

        std::cout << pv1.ref << " " << pv1.alt << " realn " << std::endl; 
    }
    
    if (hap1_inferred) {
        fill_result(pv1, qc_hap1, rslt);
    } else if (hap2_inferred) {
        fill_result(pv2, qc_hap2, rslt);
    } else {
        // inferrence failed  
    }
    
    std::cout << rslt.ppos << " " <<  rslt.pref<< " " << rslt.palt << " " << rslt.pltseq << " " << rslt.prtseq << std::endl;
}
