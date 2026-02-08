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

void check_match_pattern(Alignment& aln, 
                         std::bitset<3>& check_points,
                         const int fss, const int fse, 
                         const int ts, const int te,
                         const int fes, const int fee) {
    if (!aln.cigar_string.size()) return;

    CigarVec cigar_vec; fill_cigar_vector(aln.cigar, cigar_vec);
    //std::cout << aln.cigar_string << std::endl;
    char op = '\0'; int op_len = 0, idx = aln.ref_begin;
    for (const auto& cigar : cigar_vec) {
        op = cigar.first; op_len = cigar.second;
        switch (op) {
            case '=':
                //std::cout << idx << " " << idx + op_len << std::endl;
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
        std::cout << std::endl;
        if (find_target(loc_ref, target, vars) > -1) return true;
        vars.clear(); 
    }
    return false;    
}

/*
void index_non_matches(const CigarVec& cigar_vector, int i, NMS& nms) {
    char op = '\0'; int op_len = 0;  
    for (const auto& cigar : cigar_vector) {
        op = cigar.first; op_len = cigar.second;
        switch (op) {
            case 'S':
                nms.emplace_back(i, 'S', 0); break;  
            case 'X':
                nms.emplace_back(i, 'X', op_len); i += op_len; break; 
            case 'I':
                nms.emplace_back(i - 1, 'I', op_len); break;
            case 'D':
                nms.emplace_back(i - 1, 'D', op_len); i += op_len; break;
            default:
                i += op_len; break;
        }
    }
} 

void update_consensus(std::vector<Consensus>& cons, Alignment& aln) {
    CigarVec cigar_vector; fill_cigar_vector(aln.cigar, cigar_vector);
    char op = '\0'; int op_len = 0, i = aln.ref_begin;
    for (const auto& cigar : cigar_vector) {
        op = cigar.first; op_len = cigar.second;
        switch (op) {
            case '=':
                for (int j = 0; j < op_len; ++i, ++j) 
                    cons[i].cnt++; 
                break;
            case 'X':
                for (int j = 0; j < op_len; ++i, ++j) {
                    cons[i].cnt++; cons[i].non_matches.emplace_back('X', 1);
                } break;
            case 'I':
                cons[i - 1].non_matches.emplace_back('I', op_len); break;
            case 'D':
                cons[i - 1].non_matches.emplace_back('D', op_len); i += op_len; break;
        }
    }
}; 



//------------------------------------------------------------------------------
void match_flag(Alignment& aln, const int fl_start, const int event_start,
                const int event_end, const int fl_end, CF& cf, MF& mf) {
    CigarVec cigar_vector; fill_cigar_vector(aln.cigar, cigar_vector);
    
    // no overlap with event region
    if (aln.ref_end <= event_start || event_end <= aln.ref_begin) return; 
    
    NMS nms; index_non_matches(cigar_vector, aln.ref_begin, nms); 
    
    if (aln.ref_begin <= fl_start) cf.set(0); //flanking start covered
    if (aln.ref_begin <= event_start) cf.set(1); // event start covered
    if (event_end <= aln.ref_end) cf.set(2); // event end covered 
    if (fl_end <= aln.ref_end) cf.set(3); //flanking end covered 
   
    // flag if non-matches exist
    // flanking start to event start, event region, event end to flanking end
    for (const auto& nm : nms) {
        if (fl_start <= nm.i && nm.i < event_start) mf.set(0);
        if (event_start <=  nm.i && nm.i < event_end) mf.set(1);
        if (nm.i == event_end) { 
            if (nm.op == 'I' || nm.op == 'D') mf.set(2);
            else mf.set(1);
        }
        if (event_end < nm.i && nm.i <= fl_end) mf.set(2);
    }
}

//------------------------------------------------------------------------------
void set_up_comparison(const std::string& seq_t, const std::string& seq_nt0, 
                       const std::string& seq_nt1, const std::string& seq_nt2,
                       std::vector<Aligner>& aligners, std::bitset<4>& valid_idx) {
    if (seq_t.size()) {
        aligners[0].SetReferenceSequence(seq_t.c_str(), seq_t.size()); valid_idx.set(0);
    }
    if (seq_nt0.size()) {
        aligners[1].SetReferenceSequence(seq_nt0.c_str(), seq_nt0.size()); valid_idx.set(1);
    }
    if (seq_nt1.size()) {
        aligners[2].SetReferenceSequence(seq_nt1.c_str(), seq_nt1.size()); valid_idx.set(2);
    }
    if (seq_nt2.size()) {
        aligners[3].SetReferenceSequence(seq_nt2.c_str(), seq_nt2.size()); valid_idx.set(3);
    }
} 
*/

/*
void sw2nrs(const std::string& ref, const std::string& query, Alignment& aln, NRS& nrs) {
    CigarVec cigar_vector; fill_cigar_vector(aln.cigar, cigar_vector);
    char op = '\0'; int op_len = 0, ri = aln.ref_begin, qi = aln.query_begin;
    for (const auto& cigar : cigar_vector) {
        op = cigar.first; op_len = cigar.second;
        switch (op) {
            case '=':
                ri += op_len; qi += op_len; break;
            case 'X':
                nrs.emplace_back(ri, ref.substr(ri, op_len), query.substr(qi, op_len));
                ri += op_len; qi += op_len; break;
            case 'D':
                nrs.emplace_back(ri - 1, ref.substr(ri - 1 , op_len + 1), ref.substr(ri - 1 , 1));
                ri += op_len; break;
            case 'I':
                nrs.emplace_back(ri - 1, ref.substr(ri - 1 , 1), 
                                 ref.substr(ri - 1 , 1).append(query.substr(qi, op_len)));
                qi += op_len; break;
        }     
    }     
} 

//-----------------------------------------------------------------------------
// penalties for fewer & longer gaps + 
void penalty_grid(const UserParams& params, std::vector<std::vector<int>>& grid) {
    for (int i = param.mismatch_penal; i <= params.max_mismatch_penal; i += 2) {
        for (int j = param.gap_open_penal; j <= params.max_gap_open_penal; j += 1) {
            grid.emplace_back(i, j, 0); grid.emplace_back(i, j, 1); 
        }
    }
     
}
*/
/*
void gridsearch() {
    
    Alignment aln; Filter filter;
    
    Aligner aligner(params.match_score, params.mismatch_penal,
                    params.gap_open_penal, params.gap_ext_penal);
   // aligner.SetReferenceSequence(//c_str(), size()//);
    int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;
    aligner.Align(query, filter, &aln, mask_len);
    
    NSR nsr;
    sw2nrs; //test for target

    for (auto g : grid) {
        aligner.SetGapPenalty(g_, g__);
        aligner.Align(query, filter, &aln, mask_len);
        //test nsr
    }

}
*/     

/*
int expected_increase(const UserParams& params, Variant& target) {
    if (target.is_substitute) 
        return params.mismatch_penal * target.alt_len;
    else {
        int fill_gap = params.gap_open_penal 
                     + (target.indel_len - 1) * params.gap_ext_penal;
        
        int newly_mapped = 0;
        if (target.is_ins) 
            newly_mapped = params.match_score * target.indel_len;  
        return fill_gap + newly_mapped;
    }
}

//-----------------------------------------------------------------------------
//
void rerank_by_realn(Pileup& pileup, const Strs& read_seqs, 
                     const UserParams& params, Variant& target) {
    Aligner aligner_hap0(params.match_score, params.mismatch_penal, 
                         params.gap_open_penal, params.gap_ext_penal);
    aligner_hap0.SetReferenceSequence(pileup.seq0.c_str(), 
                                      pileup.seq0.size());
    
    Aligner aligner_ref(params.match_score, params.mismatch_penal,
                        params.gap_open_penal, params.gap_ext_penal);
    aligner_ref.SetReferenceSequence(pileup.rseq.c_str(), 
                                     pileup.rseq.size());
    //std::cout << pileup.seq0 << std::endl;
    //std::cout << pileup.rseq << std::endl;
    //std::cout << pileup.es << " " << pileup.ee  << std::endl; 
    Alignment aln_hap0, aln_ref; Filter filter_hap0, filter_ref;
    for (int i = 0; i < pileup.sz; ++i) {
        if (pileup.reads[i].rank != 'y') continue;
        
        const auto& query = read_seqs[i].c_str(); 
        int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;
        aligner_hap0.Align(query, filter_hap0, &aln_hap0, mask_len);
        //std::cout << query << " " << aln_hap0.cigar_string << std::endl;
        CF cf; MF mf;
        match_flag(aln_hap0, pileup.fs, pileup.es, pileup.ee, pileup.fe, cf, mf);
        
        // no event region overlap
        if (!cf.count()) { pileup.reads[i].rank = '\0'; continue; } 
        // perfect match with smer > 0 
        if (cf.count() == 4) { 
            if (!mf.count()) { pileup.reads[i].rank = 's'; continue; }
            if (mf.test(1)) { std::cout << "mismatch! " << pileup.reads[i].rank << " " << pileup.reads[i].smer << std::endl; continue; } 
            
            // potentially complex events
            if (pileup.reads[i].has_positional_overlap) {
                aligner_ref.Align(query, filter_ref, &aln_ref, mask_len);
                if (aln_hap0.sw_score <= aln_ref.sw_score) { 
                    pileup.reads[i].rank = 'n'; continue; 
                }
                //std::cout << aln_ref.cigar_string << std::endl;
                if (aln_ref.sw_score < aln_hap0.sw_score) { pileup.reads[i].rank = 's'; continue; }
                //else { std::cout << "unexpected" << aln_hap0.sw_score << " " << aln_ref.sw_score << " " << expected_increase(params, target) << std::endl;; }
                
            }
        }
    }
}

*/

//------------------------------------------------------------------------------
void match2haplotypes(Pileup& pileup, const Strs& read_seqs, const UserParams& params) {
    std::bitset<4> valid_idx;
    if (pileup.seq0.size()) { valid_idx.set(0); } if (pileup.seq1.size()) { valid_idx.set(1); }
    if (pileup.seq2.size()) { valid_idx.set(2); } if (pileup.rseq.size()) { valid_idx.set(3); }
    
    Alignment aln; Filter filter;
    Aligner aligner(params.match_score, params.mismatch_penal,
                    params.gap_open_penal, params.gap_ext_penal);
     
    for (int i = 0; i < pileup.sz; ++i) {
        if (pileup.reads[i].rank != 'u') continue; 
        
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
