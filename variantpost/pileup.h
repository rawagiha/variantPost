#ifndef PILEUP_H
#define PILEUP_H

#include "util.h"
#include "reads.h"

typedef std::vector<Read> Reads;
typedef std::unordered_map<std::string_view, int> Freq;

//------------------------------------------------------------------------------
// features at pileup level
struct Pileup {   
    //--------------------------------------------------------------------------    
    Pileup(const std::vector<std::string>& read_names,
           const std::vector<bool>& are_reverse, 
           const std::vector<std::string>& cigar_strs,
           const std::vector<int>& aln_starts,
           const std::vector<int>& aln_ends,
           const std::vector<std::string>& read_seqs,
           const std::vector<std::string>& quals,
           const std::vector<bool>& are_from_first_bam,
           const UserParams& params,
           LocalReference& loc_ref, Variant& target);
    
    //--------------------------------------------------------------------------
    void setHaploTypeByFrequency();
    
    //--------------------------------------------------------------------------
    void setSequenceFromHaplotype(LocalReference& loc_ref);

    void differentialKmerAnalysis(UserParams& params, LocalReference& loc_ref);
    //--------------------------------------------------------------------------    
    Reads reads; 
    
    //--------------------------------------------------------------------------
    Freq freq_s_h, freq_s, freq_u; // s_h: supporiting hiconf, u: undetermined
     
    //--------------------------------------------------------------------------
    // metrics
    int sz = -1; // pileup size
    int s_cnt = 0, n_cnt = 0, u_cnt = 0; //cnt for supporting, non, undetermined
    
    //--------------------------------------------------------------------------
    // haplotype = any non-reference patterns in a read within target_locus 
    //             +/- target_event_len in a read.
    //             reference sequence is not considered. 
    std::string_view hap0; // pattern supporting target
    std::string_view hap1; // inferred (from undetermined) as non_supporting 1
    std::string_view hap2; // as non_supporting 2

    //--------------------------------------------------------------------------
    std::string seq0;
    std::string seq1;
    std::string seq2;
    
    //--------------------------------------------------------------------------
    Kmers kmers0; 
    Kmers kmers12;
    
    //--------------------------------------------------------------------------
    // flags
    //bool has_hiconf_target = false; // center-aligned + surrounded by complex seq


};
#endif
