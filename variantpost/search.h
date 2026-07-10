#ifndef SEARCH_H
#define SEARCH_H

#include "reads.h"
#include "variant_types.h"

struct SearchResult 
{
    bool likely_simple_on_personalized_genome = false; 
    bool personalized = false;
    bool possible_snv = false;
     
    std::vector<int> positions;
    std::vector<std::string> ref_bases;
    std::vector<std::string> alt_bases;
    //std::vector<std::string> base_quals;
    //std::vector<int> skip_starts;
    //std::vector<int> skip_ends;
    int ppos = -1;
    std::string ref;
    std::string alt;
    std::string pref;
    std::string palt;
    std::string prtseq;
    std::string pltseq;
    //std::vector<std::string> read_names;
    //std::vector<bool> are_reverse;
    std::vector<int> target_statuses;
    std::vector<bool> are_from_first_bam;
    std::vector<std::string> trans_vars;

    // Variants in high-quality non-supporting reads
    // -> to be used for complex variant reconstruction
    std::unordered_map<VariantKey, int, VariantKeyHash> hq_negative_cnts;
    
    //bool is_retargeted;

    SearchResult();
                                                
    void setReadInfo(const Read& read, const char qual_thresh);
    void finalize(); 
    
    //void fill_read_info(const Reads& reads, const int target_status);
    
    /*
    void report(
        const Contig& contig,
        Reads& targets,
        Reads& non_targets,
        Reads& undetermined,
        const bool is_retargeted = false
    );*/
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
