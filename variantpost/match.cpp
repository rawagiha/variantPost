#include "match.h"

typedef std::vector<NonMatch> NMS;

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
        if (fl_start <= nm.i && nm.i< event_start) mf.set(0);
        if (event_start <= nm.i && nm.i <= event_end) mf.set(1);
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

//------------------------------------------------------------------------------
void match2haplotypes(Pileup& pileup, const std::vector<std::string>& read_seqs,
                      const UserParams& params) {
    std::bitset<4> valid_idx;
    if (pileup.seq0.size()) { valid_idx.set(0); } if (pileup.seq1.size()) { valid_idx.set(1); }
    if (pileup.seq2.size()) { valid_idx.set(2); } if (pileup.rseq.size()) { valid_idx.set(3); }
    
    Alignment aln; Filter filter;
    Aligner aligner(params.match_score, params.mismatch_penal,
                    params.gap_open_penal, params.gap-ext_penal);
     
    for (int i = 0; i < pileup.sz; ++i) {
        if (pileup.reads[i].rank != 'u') continue; 
        
        const auto& query = read_seqs[i].c_str();
        int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;  
        
        std::vector<int> scores = {0, 0, 0, 0};
        for (int j = 0; j < 4; ++j) {
            if (valid_idx.test(j)) {
                if (j == 0) {
                    aligner.SetReferenceSequence(pileup.seq0.c_str(), pileup.seq0.size());
                } else if (j == 1) {
                    aligner.SetReferenceSequence(pileup.seq1.c_str(), pileup.seq1.size());
                } else if (j == 2) {
                    aligner.SetReferenceSequence(pileup.seq2.c_str(), pileup.seq2.size());
                } else {
                    aligner.SetReferenceSequence(pileup.rseq.c_str(), pileup.rseq.size());
                }
                aligner.Align(query, filter, &aln, mask_len); scores[j] = aln.sw_score;
                std::cout << aln.cigar_string << " " << aln.sw_score << " " << pileup.reads[i].is_control << " " << pileup.reads[i].name << std::endl;
            } 
        }
         
        // aln against target hap is solely highest
        if (scores[0] > scores[1] && scores[0] > scores[2] && scores[0] > scores[3]) {
            pileup.reads[i].rank = 's'; ++pileup.s_cnt; --pileup.u_cnt; 
            std::cout << scores[0] << " "  << scores[1] << " " << scores[2] << " " <<  scores[3] << std::endl;
        } else if (scores[0] == scores[1] || scores[0] == scores[2] || scores[0] == scores[3]) {
            /* tie score to one of non-target haplotypes -> remain as 'u'*/
        } else {
            pileup.reads[i].rank = 'n'; ++pileup.n_cnt; --pileup.u_cnt;
        }
    }
}
