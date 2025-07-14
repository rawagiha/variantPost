#include "sequence_model.h"

//------------------------------------------------------------------------------
SequenceModel::SequenceModel(Pileup& pileup, LocalReference& loc_ref, Variant& target) {
    std::vector<int> starts, ends;
    for (const auto& read : pileup.reads) {
        if (read.rank == 'u' && read.qc_passed) {
            starts.push_back(read.covering_start); 
            ends.push_back(read.covering_end);
        }
    }

    start = *std::min_element(starts.begin(), starts.end());
    end = *std::max_element(ends.begin(), ends.end());
    
    int idx = 0, pos = start;
    seq.append(loc_ref._seq.substr(start - loc_ref.start, target.pos - start));
    for (/* */; pos < target.pos; pos++) idx2pos.emplace_back(idx++, pos);

    seq.append(target.alt); target_start = target.pos - start -1;
    for(size_t i = 0; i < target.alt.size(); ++i) idx2pos.emplace_back(idx++, pos);

    pos = target._end_pos;
    seq.append(loc_ref._seq.substr(pos - loc_ref.start, end - pos));
    for (/* */; pos < end; ++pos) idx2pos.emplace_back(idx++, pos);
   
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
    int flen = loc_ref.flanking_end - loc_ref.flanking_start;
    int sz = (flen / 2 > params.kmer_size) ? flen / 2  : params.kmer_size; //may be revised
    
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
    Aligner aligner(params.match_score, params.mismatch_penal, 255, 255);
    aligner.SetReferenceSequence(seq.c_str(), seq.size()); 
    Alignment aln; Filter filter;
    std::cout << seq << std::endl;
    std::cout << flank_start << " " << target_start << " " << target_end << " " << flank_end << std::endl;
    for (int i = 0; i < pileup.sz; ++i) {
        if (pileup.reads[i].smer > 0 && !pileup.reads[i].nmer) {
            const auto& query = read_seqs[i].c_str(); 
            int32_t mask_len = strlen(query) < 30 ? 15 : strlen(query) / 2;
            
            aligner.Align(query, filter, &aln, mask_len);
            std::cout << read_seqs[i] << std::endl;
            std::cout << aln.ref_begin << " " << aln.ref_end << " " << aln.query_begin << " " << aln.query_end << " " << aln.cigar_string << std::endl;  
            
        }
    }


}
