#ifndef PILEUP_H
#define PILEUP_H

#include "util.h"
#include "reads.h"
#include "variant_types.h"

// Fw declaration for doc purposes (they are in util.h)
struct UserParams;
struct LocalReference;
struct Variant;

using Reads = std::vector<Read>;
using Idx   = std::unordered_map<std::string_view, std::vector<int>>;

//------------------------------------------------------------------------------
struct Pileup {   
    //--------------------------------------------------------------------------    
    Pileup(const Strs& read_names, const Bools& are_reverse, const Strs& cigar_strs,
           const Ints& aln_starts, const Ints& aln_end, const Strs& read_seqs,
           const Strs& quals, const Bools& are_from_first_bam, const bool has_second,
           const UserParams& params, LocalReference& loc_ref, Variant& target);
    
    void inferGermlineHaplotype(const UserParams& param, const int target_pos);    
    void gridSearch(const UserParams& params, LocalReference& loc_ref, const Variant& target); 
    void setHaploTypes(LocalReference& loc_ref, const Variant& target);  
    void differentialKmerAnalysis(const UserParams& params, LocalReference& loc_ref, const Variant& target); 
    void searchByRealignment(const UserParams& params, LocalReference& loc_ref, const Variant& target);
    
    //--------------------------------------------------------------------------    
    Reads reads; 
    
    //--------------------------------------------------------------------------
    Vars hap0_vars;
    Vars hap1_vars;
    Vars hap2_vars;
    Vars homo_vars;

    //--------------------------------------------------------------------------
    // variants between flank start/end in non-supporting reads
    std::unordered_map<Variant, int> ns_vars; 
    
    Ints hiconf_s_read_idx;
    Ints s_read_idx;
    
    //--------------------------------------------------------------------------
    // signature related 
    Idx sig_u, u_sig_annot; // {u_sig : anno flags} 
    std::vector<std::string_view> anti_sigs; // sigs unlikely supporting target
      
    //--------------------------------------------------------------------------
    // metrics
    int sz = -1;
    int s_cnt = 0, n_cnt = 0, u_cnt = 0, y_cnt = 0, z_cnt = 0; //cnt for supporting, non, undetermined
    int v_cnt = 0; // variants in flanking regions
     
    //--------------------------------------------------------------------------
    //Ints starts, ends; // storing coverings starts/ends
    int start = -1, end = -1; //typically, start = min(starts), end = max(ends)
    int low_cplx_start = -1, low_cplx_end = -1;
    
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
    Ints i2p0, i2p1, i2p2, i2pr;
    int es = -1, ee = -1, fs = -1, fe = -1;

    //--------------------------------------------------------------------------
    int kmer_sz = 0;
    Kmers kmers_t; // kmers specific to target
    Kmers kmers_nt; // kmers specific to non_target  
    
    //--------------------------------------------------------------------------
    // flags
    bool has_second_bam = false; // second BAM input
    bool has_hiconf_support = false; // center-aligned + surrounded by complex seq
    bool has_no_support = false; // no support no undetermined
    bool has_likely_support = false; // reads likely supporting the target
    bool has_ref_hap = false; // reference haplotype likely exists in background
    bool has_excess_ins_hap = false; // haplotype with ins-repeat > expected
    bool no_non_target_haps = false; // non-target haplotypes not inferred
    bool is_ref_hom = false; 
    bool is_ref_het = false;
    bool is_alt_hom = false;
    bool is_alt_het = false;
};

#endif
