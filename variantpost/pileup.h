#ifndef PILEUP_H
#define PILEUP_H

#include "util.h"
#include "reads.h"

typedef std::vector<Read> Reads;
typedef std::unordered_map<std::string_view, int> Freq;

//------------------------------------------------------------------------------
// features in a collection of reads (pileup)
struct Pileup {   
    //--------------------------------------------------------------------------    
    Pileup(const Strs& read_names, const Bools& are_reverse, const Strs& cigar_strs,
           const Ints& aln_starts, const Ints& aln_end, const Strs& read_seqs,
           const Strs& quals, const Bools& are_from_first_bam,
           const UserParams& params, LocalReference& loc_ref, Variant& target);
    
    //--------------------------------------------------------------------------
    void setHaploTypeByFrequency();
    
    //--------------------------------------------------------------------------
    void setSequenceFromHaplotype(LocalReference& loc_ref);

    void reRankByKmer(UserParams& params, LocalReference& loc_ref);
    void compareToRefByKmer(LocalReference& loc_ref, UserParams& params, const Variant& t);

    //--------------------------------------------------------------------------    
    Reads reads; 
    
    //--------------------------------------------------------------------------
    Freq freq_s_h, freq_s, freq_u; // s_h: supporiting hiconf, u: undetermined
     
    //--------------------------------------------------------------------------
    // metrics
    int sz = -1; // pileup size
    int s_cnt = 0, n_cnt = 0, u_cnt = 0; //cnt for supporting, non, undetermined
    
    //--------------------------------------------------------------------------
    int start = -1, end = -1; //haplotype sequence start/end
    
    //--------------------------------------------------------------------------
    // haplotype = any non-reference patterns in a read within target_locus 
    //             +/- target_event_len in a read.
    // see setHaploTypeByFrequency in "pileup.cpp" for LOGIC
    std::string_view hap0; // pattern supporting target
    std::string_view hap1; // inferred (from rank 'u') as non_supporting 1
    std::string_view hap2; // as non_supporting 2.
    
    //--------------------------------------------------------------------------
    // set by Pileup by inferring from the actual reads
    std::string seq0; // seq with hap0 (target)
    std::string seq1; // seq with hap1 (non-target 1)
    std::string seq2; // seq with hap2 (non-target 2). This may be refseq
    std::string rseq; // reference hap. this will be set if necessary

    //--------------------------------------------------------------------------
    // set by SequenceModel by inserting the target into reference seq
    //std::string tseq; // refseq with target 

    //--------------------------------------------------------------------------
    // set by either of Pileup or SequenceModel if necessary
    //std::string rseq; // refseq
    
    //--------------------------------------------------------------------------
    Kmers kmers_t; // kmers specific to target
    Kmers kmers_nt; // kmers specific to non_target  
    
    //--------------------------------------------------------------------------
    // flags
    bool has_hiconf_support = false; // center-aligned + surrounded by complex seq
    bool has_no_support = false; // no support no undetermined
    bool has_ref_hap = false; // reference haplotype likely exists in background

};
#endif
