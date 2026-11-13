#ifndef SEARCH_H
#define SEARCH_H

#include "reads.h"
#include "variant_types.h"

struct SearchResult 
{
    //-------------------------------------------------------------------------
    SearchResult();
    
    //-------------------------------------------------------------------------
    void setReadInfo(const Read& read, const char qual_thresh);

    //-------------------------------------------------------------------------
    // feed trans_vars with recurrent HQ variants in non-target reads
    void finalize();
    
    // 8x-byte members 
    Ints positions;
    Strs ref_bases;
    Strs alt_bases;
    
    Ints target_statuses;
    Bools are_from_first_bam;
    Strs trans_vars;
    
    // stores high quality variant-key and its count in reads not carrying 
    // the target variant (negative). 
    std::unordered_map<VariantKey, int, VariantKeyHash> hq_negative_cnts;
     
    std::string ref;
    std::string alt;
    std::string pref;
    std::string palt;
    std::string prtseq;
    std::string pltseq;
    
    // stores homo vars (common to hap1/2)
    std::string homos;

    // stores het vars in hap parental to hap0
    std::string parental_hap_hets;  
    
    // local backgound genotype
    std::string loc_gt = "ref_hom";
    
    // inferrence method
    std::string infer_by;
    
    //4x-byte members
    int ppos = -1;
    int parental_hap = -1; // the origin of hap0
    
    //1x-byte members
    bool likely_simple_on_personalized_genome = false; 
    bool personalized = false;
    bool possible_snv = false;
};


//interface with python
void _search_target(
    SearchResult&,
    const std::string&,
    const std::string&,
    const int,
    const std::string&,
    const std::string&,
    //const int,
    const int, 
    const double,
    const int,
    const int,
    const int,
    const int,
    const int,
    const int,
    const int,
    const int,
    const int,
    const std::vector<std::string>&,
    const std::vector<bool>&, 
    const std::vector<std::string>&,
    const std::vector<int>&, 
    const std::vector<int>&,
    const std::vector<std::string>&, 
    const std::vector<std::string>&, 
    //const std::vector<int>&,
    const std::vector<bool>&,
    const bool
);

#endif
