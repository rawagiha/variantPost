#ifndef SEQUENCE_MODEL_H
#define SEQUENCE_MODEL_H

#include "pileup.h"


//------------------------------------------------------------------------------
struct SequenceModel {
    //--------------------------------------------------------------------------
    SequenceModel(Pileup& pileup, LocalReference& loc_ref, Variant& target);

    void compareToRefByKmer(Pileup& pileup, LocalReference& loc_ref, UserParams& params); 
    void reRankByReAlignment(Pileup& pileup, const Strs& read_seqs, UserParams& params);
     
    //--------------------------------------------------------------------------
    std::string seq; int start = -1, end = -1; Coord idx2pos;
    
    //--------------------------------------------------------------------------
    int target_start = -1, target_end = -1, flank_start = -1, flank_end = -1;

    //--------------------------------------------------------------------------
    bool has_target_spec_kmers = false;
    bool has_ref_spec_kmers = false;
};



#endif
