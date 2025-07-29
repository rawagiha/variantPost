#include "sequence_model.h"
#include "match.h"

// 1. sort 's' read by start
// 2. make coverage stat
// 3. cound variants
// 4. take consensus 

//------------------------------------------------------------------------------
/*
void make_sequence_2(const int start, const int end, 
             const Variants& variants, std::string& seq, Coord& idx2pos) {
        
}*/ 


//------------------------------------------------------------------------------
SequenceModel::SequenceModel(Pileup& pileup, LocalReference& loc_ref, Variant& target) {
    std::vector<int> starts, ends;
    for (const auto& read : pileup.reads) {
        if (read.rank == 'u' && read.qc_passed) {
            starts.push_back(read.covering_start); 
            ends.push_back(read.covering_end);
        }
    }

    if (!starts.empty()) {
        start = *std::min_element(starts.begin(), starts.end());
        end = *std::max_element(ends.begin(), ends.end());
    }
    
    if (target.pos <= start) start = loc_ref.start;
    if (end <= target._end_pos) end = loc_ref.end;

    Vars variants = {target};
    // from pileup.h
    make_sequence(loc_ref, variants, start, end, seq, &idx2pos);
    
         
    // set key index 
    target_start = target.pos - start;
    for (auto elem : idx2pos) {
        if (flank_start < 0 && loc_ref.flanking_start <= elem.second) flank_start = elem.first;
        if (flank_end < 0 && loc_ref.flanking_end <= elem.second) flank_end = elem.first;
        if (target_end < 0 && target.end_pos <= elem.second) target_end = elem.first; 
    }
}

//------------------------------------------------------------------------------
void SequenceModel::compareToRefByKmer(Pileup& pileup, LocalReference& loc_ref, UserParams& params) {
    auto rseq = loc_ref.seq.substr(start - loc_ref.start, end - start);

    Kmers km0, km1, kmt, kmr;
    int sz = (loc_ref.low_cplx_len > params.kmer_size) 
           ? loc_ref.low_cplx_len  : params.kmer_size;
   
    make_kmers(seq, sz, km0); make_kmers(rseq, sz, km1);
    std::set_difference(km0.begin(), km0.end(), km1.begin(), km1.end(),
                        std::inserter(kmt, kmt.end()));
    has_target_spec_kmers = (!kmt.empty());

    std::set_difference(km1.begin(), km1.end(), km0.begin(), km0.end(),
                        std::inserter(kmr, kmr.end()));
    has_ref_spec_kmers = (!kmr.empty());
    for (auto& read : pileup.reads) {
        if (read.rank == 'u') {
            for (const auto& kmer : kmt)
                if (read.seq.find(kmer) != std::string_view::npos) read.smer++;
            for (const auto& kmer : kmr)
                if (read.seq.find(kmer) != std::string_view::npos) read.nmer++;
        }
    }
} 

//------------------------------------------------------------------------------
void SequenceModel::reRankByReAlignment(Pileup& pileup, const std::vector<std::string>& read_seqs, UserParams& params) {
    Aligner aligner(params.match_score, params.mismatch_penal, params.gap_open_penal, params.gap_ext_penal);
    aligner.SetReferenceSequence(seq.c_str(), seq.size()); 
    Alignment aln; Filter filter;
    for (int i = 0; i < pileup.sz; ++i) {
        
        if (!pileup.reads[i].smer || pileup.reads[i].nmer) continue;
        
        // hereafter smers > 0 with no nmers
        const auto& query = read_seqs[i].c_str(); 
        int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;
        aligner.Align(query, filter, &aln, mask_len);
         
        CF cf; MF mf;
        match_flag(aln, flank_start, target_start, target_end, flank_end, cf, mf);
        // no event region overlap
        if (!cf.count()) { pileup.reads[i].rank = '\0'; continue; } 
        // perfect match with smer > 0 
        if (cf.count() == 4 && !mf.count()) { pileup.reads[i].rank = 's'; continue; }
        // perfect coverage with perfct event match (allow flanking non-matches)
        //if (cf.count() == 4 && !mf.test(1)) { pileup.reads[i].rank = 's'; continue; }  
    }
}
