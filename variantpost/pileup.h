#ifndef PILEUP_H
#define PILEUP_H

#include "util.h"
#include "reads.h"

typedef std::vector<Read> Reads;
typedef std::unordered_map<std::string_view, std::vector<int>> Idx;

//------------------------------------------------------------------------------
// features in a collection of reads (pileup)
struct Pileup {   
    //--------------------------------------------------------------------------    
    Pileup(const Strs& read_names, const Bools& are_reverse, const Strs& cigar_strs,
           const Ints& aln_starts, const Ints& aln_end, const Strs& read_seqs,
           const Strs& quals, const Bools& are_from_first_bam, const bool has_second,
           const UserParams& params, LocalReference& loc_ref, Variant& target);
    
    void gridSearch(const UserParams& params, LocalReference& loc_ref, const Variant& target); 
    
    void setHaploTypes(LocalReference& loc_ref, const Variant& target);  
    
    void differentialKmerAnalysis(const UserParams& params,
                                  LocalReference& loc_ref, const Variant& target); 
    
    void searchByRealignment(const UserParams& params,
                             LocalReference& loc_ref, const Variant& target);
    
    //--------------------------------------------------------------------------    
    Reads reads; 
    
    //--------------------------------------------------------------------------
    Idx sig_s_hiconf, sig_s, sig_u, u_sig_annot;
    int hiconf_read_idx = -1; // for index-case read in Pileup.reads
      
    //--------------------------------------------------------------------------
    // metrics
    int sz = -1; // pileup size
    int s_cnt = 0, n_cnt = 0, u_cnt = 0, y_cnt = 0; //cnt for supporting, non, undetermined
    
    //--------------------------------------------------------------------------
    Ints starts, ends; // storing coverings starts/ends
    int start = -1, end = -1; //typically, start = min(starts), end = max(ends)
                              
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
    Coord i2p_0; // index/pos pair for hap0
    Coord i2p_r; // index/pos pair for reference hap
    int es = -1, ee = -1, fs = -1, fe = -1;

    //--------------------------------------------------------------------------
    int kmer_sz = 0;
    Kmers kmers_t; // kmers specific to target
    Kmers kmers_nt; // kmers specific to non_target  
    
    //--------------------------------------------------------------------------
    // flags
    bool has_second = false; // second BAM input
    bool has_hiconf_support = false; // center-aligned + surrounded by complex seq
    bool has_no_support = false; // no support no undetermined
    bool has_likely_support = false; // reads likely supporting the target
    bool has_ref_hap = false; // reference haplotype likely exists in background
    bool has_excess_ins_hap = false; // haplotype with ins-repeat > expected
    bool vs_ref_hap = false; // true to perform kmer analysis vs ref hap
    bool no_non_target_haps = false; // non-target haplotypes not inferred
    bool has_valid_boundary = false; // valid flanking start/end event start/end   
};

#endif
