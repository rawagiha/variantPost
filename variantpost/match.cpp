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
        if (!pileup.reads[i].qc_passed || pileup.reads[i].rank != 'u' || pileup.reads[i].smer == 0) continue; 
        
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
            pileup.reads[i].rank = 's'; ++pileup.s_cnt; --pileup.u_cnt; 
        } else if (scores[0] == scores[1] || scores[0] == scores[2] || scores[0] == scores[3]) {
            // tie score to one of non-target haplotypes -> remain as 'u'
        } else {
            pileup.reads[i].rank = 'n'; ++pileup.n_cnt; --pileup.u_cnt;
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

void realn_to_perfonalized_genome(
    LocalReference& loc_ref,
    Aligner& alngr, Filter& fltr, Alignment& aln,
    Variant& pv, 
    bool& is_personalized, 
    bool& is_likely_simple, 
    const std::string& mut_seq, 
    const std::string& hap_seq, const Ints& i2p) {
    
    const char* query = mut_seq.c_str();
    int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;
    alngr.SetReferenceSequence(hap_seq.c_str(), hap_seq.size());
    alngr.Align(query, fltr, &aln, mask_len);
    
    std::cout << aln.cigar_string << " realn " << std::endl; 
    CigarVec cigar_vec;
    fill_cigar_vector(aln.cigar_string, cigar_vec);
    
    Vars vars; 
    std::vector<std::pair<int, int>> matched_sides;
    int i = aln.ref_begin, j = aln.query_begin;
    size_t curr = 0;
    const size_t last = cigar_vec.size() - 1; 
    for (const auto& [op, op_len] : cigar_vec) {
        switch (op) {
            case '=':
                i += op_len; j += op_len; break;
            case 'X':
                vars.emplace_back(i2p[i], hap_seq.substr(i, op_len), mut_seq.substr(j, op_len));
                fill_matched_sides(curr, last, cigar_vec, matched_sides);
                i += op_len; j+= op_len;    
                break;
            case 'I': {
                Variant v(i2p[i - 1], hap_seq.substr(i - 1, 1), mut_seq.substr(j - 1 , op_len +1));
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
    
    int closest = -1, min_dist = INT_MAX, dist = INT_MAX;
    for (int i = 0; i < static_cast<int>(vars.size()); ++i) {
        const auto& v = vars[i]; 
        if (loc_ref.flanking_start <= v.pos && v.pos <= loc_ref.flanking_end) {
            dist = std::abs(v.pos - pv.pos); // pv.pos == target.pos
            if (dist < min_dist) {closest = i; min_dist = dist; } 
        }
    }
    std::cout << closest << " closet " << " " << loc_ref.locplx_start << " " << loc_ref.locplx_end << std::endl; 
    if (closest == -1) return; 
    pv = std::move(vars[closest]); 
    
    // Checking for complex indels on the personalized genome
    if (matched_sides[closest].first > 19 && matched_sides[closest].second > 19) {
        // exclude MNVs 
        if (!pv.is_complex) is_likely_simple = true;
    }

    is_personalized = true;

}

inline void fill_result(const Variant& v, SearchResult& rslt) {
    rslt.ppos = v.pos;
    rslt.pref = v.ref;
    rslt.palt = v.alt;
    rslt.pltseq = v.sample_lt_seq;
    rslt.prtseq = v.sample_rt_seq;    
}

void personalize(const Pileup& pileup, LocalReference& loc_ref, const UserParams& params, const Variant& target, SearchResult& rslt) {
        
    if (pileup.hiconf_read_idx < 0 || pileup.seq0.empty()) return;
    
    // REF/REF case -> no personalization
    if (pileup.is_ref_hom) return;
    
    std::string hiconf_seq = pileup.seq0;
    
    Alignment aln; Filter filter;
    Aligner aligner(params.match_score, params.mismatch_penal,
                    params.gap_open_penal, params.gap_ext_penal);
    
    size_t with_hap1 = count_overlap(pileup.hap0_vars, pileup.hap1_vars);
     
    Variant pv1(target.pos, target.ref, target.alt);
    Variant pv2(target.pos, target.ref, target.alt);
    
    std::cout << pv1.ref << " " << pv1.alt << " snv? " << std::endl;

        
    HapLL::RepeatInfo ri1, ri2, rir;
    bool is_personalized_hap1 = false, is_personalized_hap2 = false;
    bool is_likely_simple_hap1 = false, is_likely_simple_hap2 = false;
    if (pileup.is_alt_het) {    
        size_t with_hap2 = count_overlap(pileup.hap0_vars, pileup.hap2_vars);
        
        if (with_hap1 > with_hap2) {
            // assigned to hap1. 
            // is_personalized_hap1 will remainfalse if the realignment fails
            realn_to_perfonalized_genome(loc_ref, aligner, filter, aln, pv1, 
                                         is_personalized_hap1, is_likely_simple_hap1, hiconf_seq, pileup.seq1, pileup.i2p1);
        } else if (with_hap1 < with_hap2) {
            // assigned to hap2
            // is_personalized_hap2 will remainfalse if the realignment fails
            realn_to_perfonalized_genome(loc_ref, aligner, filter, aln, pv2, 
                                         is_personalized_hap2, is_likely_simple_hap2, hiconf_seq, pileup.seq2, pileup.i2p2);
        } else {
            // inference hap1 vs. hap2
            realn_to_perfonalized_genome(loc_ref, aligner, filter, aln, pv1, 
                                         is_personalized_hap1, is_likely_simple_hap1, hiconf_seq, pileup.seq1, pileup.i2p1);
            realn_to_perfonalized_genome(loc_ref, aligner, filter, aln, pv2, 
                                         is_personalized_hap2, is_likely_simple_hap2, hiconf_seq, pileup.seq2, pileup.i2p2);

            double ll_hap1 = is_personalized_hap1 ? HapLL::evaluate_variant(pv1, ri1) : -std::numeric_limits<double>::infinity();
            double ll_hap2 = is_personalized_hap2 ? HapLL::evaluate_variant(pv2, ri2) : -std::numeric_limits<double>::infinity();

            if (ll_hap1 >= ll_hap2) {
                is_personalized_hap2 = false; // discard hap2, favor hap1 on exact ties
            } else {
                is_personalized_hap1 = false; // discard hap1
            }
        }
    } else {
        realn_to_perfonalized_genome(loc_ref, aligner, filter, aln, pv1, 
                                     is_personalized_hap1, is_likely_simple_hap1, hiconf_seq, pileup.seq1, pileup.i2p1);
        
        std::cout << pv1.ref << " " << pv1.alt << " snv? " << std::endl; 
        // inference hap1 vs. hap_ref
        if (!with_hap1) {
            Variant pv_ref = target;
            std::cout << target.lt_seq << " <- tar lt, tar rt -> " << target.rt_seq << std::endl; 
            pv_ref.sample_lt_seq = std::string{target.lt_seq};
            pv_ref.sample_rt_seq = std::string{target.rt_seq};
            std::cout << pv_ref.sample_lt_seq << " " << pv_ref.sample_rt_seq << std::endl; 
            double ll_hap1 = HapLL::evaluate_variant(pv1, ri1);
            double ll_ref  = HapLL::evaluate_variant(pv_ref, rir);
            
            std::cout << ll_hap1 << " " << ll_ref << std::endl;
             
            if (ll_hap1 < ll_ref) {
                is_personalized_hap1 = false;
            } 
        }   
    }
    
    if (is_personalized_hap1) {
        fill_result(pv1, rslt);
        rslt.likely_simple_on_personalized_genome = is_likely_simple_hap1;
    } else if (is_personalized_hap2) {
        fill_result(pv2, rslt);
        rslt.likely_simple_on_personalized_genome = is_likely_simple_hap2;
    } 
    
    std::cout << rslt.ppos << " " <<  rslt.pref<< " " << rslt.palt << " " << rslt.pltseq << " " << rslt.prtseq << std::endl;
}
