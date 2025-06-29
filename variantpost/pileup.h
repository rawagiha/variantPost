#ifndef PILEUP_H
#define PILEUP_H

#include "util.h"
#include "reads.h"

typedef std::vector<Read> Reads;

//------------------------------------------------------------------------------
// features at pileup level
struct Pileup
{   
    Pileup(const std::vector<std::string>& read_names,
           const std::vector<bool>& are_reverse, 
           const std::vector<std::string>& cigar_strs,
           const std::vector<int>& aln_starts,
           const std::vector<int>& aln_ends,
           const std::vector<std::string>& read_seqs,
           const std::vector<std::vector<int>>& quals,
           const std::vector<int>& mapqs,
           const std::vector<bool>& are_from_first_bam,
           const UserParams& params,
           LocalReference& loc_ref, Variant& target);
    
    void setLocalHaploTypes(LocalReference& loc_ref);
    
    
    Reads reads; 

};
#endif
