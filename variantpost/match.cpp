#include "match.h"

//------------------------------------------------------------------------------
NonMatch::NonMatch(const int i_, const string& ref_, const string& alt_) 
    : i(i_), ref(ref_), alt(alt_) { 
    
    if (ref.size() == alt.size()){
        if (ref.size() == 1) is_snv = true;
        else is_mnv = true;
    } else if (ref.size() < alt.size()) {
        is_ins = true;
    } else is_del = true;  
}

//------------------------------------------------------------------------------
void aln2nms(Alignment& aln, NMS& nms,
             const std::string& ref, const std::string& query) {
    if (!aln.cigar_string.size()) return;

    CigarVec cigar_vec; fill_cigar_vector(aln.cigar, cigar_vec);
    
    char op = '\0';
    int op_len = 0, ri = aln.ref_begin, qi = aln.query_begin;
    for (const auto& cigar : cigar_vec) {
        op = cigar.first; op_len = cigar.second;
        switch (op) {
            case '=':
                ri += op_len; qi += op_len;
                break;
            case 'X':
                nms.emplace_back(ri, 
                                 ref.substr(ri, op_len),
                                 query.substr(qi, op_len));
                ri += op_len; qi += op_len;
                break;
            case 'D':
                if (ri) {
                    nms.emplace_back(ri - 1,
                                     ref.substr(ri - 1, op_len + 1),
                                     ref.substr(ri - 1, 1));
                    ri += op_len;
                } break;
            case 'I':
                if (ri) {
                    std::string base = ref.substr(ri - 1, 1);
                    nms.emplace_back(ri - 1,
                                     ref.substr(ri - 1, 1),
                                     base.append(query.substr(qi, op_len)));
                    qi += op_len;
                } break;
            case 'S':
                break;
        }
    }
} 


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


/*
//------------------------------------------------------------------------------
int find_target(const int start, const NMS& nms, 
                LocalReference& loc_ref, const Variant& target) {
    for (size_t j = 0; j < nms.size(); ++j) {
        if (target.pos == start + nms[j].i
            && target.ref == nms[j].ref && target.alt == nms[j].alt) return nms[j].i;
         
        //swappable
        if (j + 1 < nms.size()) {
            if (target.is_ins && nms[j].is_del && nms[j + 1].is_ins  
        }

    
    }
    return -1;
}*/

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



// crop_function 

void personalize(const Pileup& pileup, LocalReference& loc_ref, const UserParams& params, const Variant& target, Variant& per) {
    std::cout << " get here hiconf idx " <<  pileup.hiconf_read_idx << std::endl;
    
    
    if (pileup.hiconf_read_idx < 0) return;
    if (pileup.hap1.empty() && pileup.hap2.empty()) return;
    
    Vars hiconfvars;
    for (const auto& v : pileup.reads[pileup.hiconf_read_idx].variants) 
        if (v.mean_qual >= params.base_q_thresh) hiconfvars.push_back(v);
    
    std::string hiconf_seq;
    make_sequence(loc_ref, hiconfvars, pileup.start, pileup.end, hiconf_seq);
    
    auto query = hiconf_seq.c_str();
    int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;
    Alignment aln; Filter filter;    
    Aligner aligner(params.match_score, params.mismatch_penal,
                    params.gap_open_penal, params.gap_ext_penal);
    
    int aln_score = INT_MAX, which_hap = 1;
    std::string cigar_str = ""; 
    if (!pileup.seq1.empty()) {
        aligner.SetReferenceSequence(pileup.seq1.c_str(), pileup.seq1.size());
        aligner.Align(query, filter, &aln, mask_len);
        aln_score = aln.sw_score; cigar_str = aln.cigar_string;
    }    
    
    if (!pileup.seq2.empty()) {
        aligner.SetReferenceSequence(pileup.seq2.c_str(), pileup.seq1.size());
        aligner.Align(query, filter, &aln, mask_len);
        std::cout << "second " << aln.sw_score << std::endl;
        if (aln.sw_score > aln_score) {
            which_hap = 2; cigar_str = aln.cigar_string;
        }
    }

    CigarVec cigar_vec;
    fill_cigar_vector(cigar_str, cigar_vec); 
    const std::string& background = (which_hap == 1) ? pileup.seq1 : pileup.seq2;
    const Coord& i2p = (which_hap == 1) ? pileup.i2p_1 : pileup.i2p_2;
    char op = '\0'; int op_len = 0, i = 0, j = 0;
    Vars vars;
    for (const auto& c : cigar_vec) {
        op = c.first; op_len = c.second;
        switch (op) {
            case '=':
                i += op_len; j += op_len; break;
            case 'X':
                vars.emplace_back(i2p[i].second, background.substr(i, op_len), hiconf_seq.substr(j, op_len));
                i += op_len; j+= op_len;   
                break;
            case 'I':
                vars.emplace_back(i2p[i - 1].second, background.substr(i - 1, 1), hiconf_seq.substr(j - 1 , op_len +1));
                j += op_len;
                std::cout << i << " " << i2p[i].first << " " << i2p[i].second << std::endl;
                std::cout << i2p[i - 2].first << " " << i2p[i - 2].second << std::endl;
                std::cout << i2p[i - 1].first << " " << i2p[i - 1].second << std::endl;
                std::cout << i2p[i - 1 + op_len].first << " " << i2p[i - 1 + op_len ].second << std::endl;
                std::cout << i2p[i + op_len].first << " " << i2p[i + op_len ].second << std::endl;
                std::cout << i2p[i + op_len + 1].first << " " << i2p[i + op_len + 1].second << std::endl;
                break;
            case 'D':
                vars.emplace_back(i2p[i - 1].second, background.substr(i - 1, op_len + 1), hiconf_seq.substr(j - 1 , 1));

                std::cout << i << " " << i2p[i].first << " " << i2p[i].second << std::endl;
                std::cout << i2p[i - 1].first << " " << i2p[i - 1].second << std::endl;
                std::cout << i2p[i - 1 + op_len].first << " " << i2p[i - 1 + op_len ].second << std::endl;
                std::cout << i2p[i + op_len].first << " " << i2p[i + op_len ].second << std::endl;
                std::cout << i2p[i + op_len + 1].first << " " << i2p[i + op_len + 1].second << std::endl;
                i += op_len;
                break;
            case 'S':
                j += op_len;
                break;
        }
        
    } 
    int closest = -1, min_dist = INT_MAX, dist = INT_MAX;
    for (int i = 0; i < static_cast<int>(vars.size()); ++i) {
        auto& v = vars[i]; 
        std::cout << v.pos << " " << v.ref << " " << v.alt << std::endl; 
        if (!v.is_substitute) {
            dist = std::abs(v.pos - target.pos);
            if (dist < min_dist) {closest = i; min_dist = dist; } 
        }
    } 
    if (closest == -1) { per = target; return; }
    bool is_same = (vars[closest] == target);
    std::cout << vars[closest].pos << " " << vars[closest].ref << " " << vars[closest].alt << " " << is_same << std::endl;
    if (is_same) { per = target; } else { per = vars[closest]; }
}
