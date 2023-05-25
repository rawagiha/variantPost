#ifndef PILEUP_PROCESSOR_H
#define PILEUP_PROCESSOR_H

#include <string>
#include <vector>

struct ProcessedPileup 
{
    std::vector<int> positions;
    std::vector<std::string> ref_bases;
    std::vector<std::string> alt_bases;
    std::vector<std::string> base_quals;
    std::vector<int> skip_starts;
    std::vector<int> skip_ends;
    int target_pos;
    std::string ref;
    std::string alt;
    std::vector<std::string> read_names;
    std::vector<bool> are_reverse;
    std::vector<int> target_statuses;
    std::vector<bool> are_from_first_bam;

    ProcessedPileup();
                                                
    ProcessedPileup(
        const std::vector<int> & positions,
        const std::vector<std::string> & ref_bases,
        const std::vector<std::string> & alt_bases,
        const std::vector<std::string> & base_quals,
        const std::vector<int> & skip_starts,
        const std::vector<int> & skip_ends,
        const int target_pos,
        const std::string ref,
        const std::string alt,
        std::vector<std::string> & read_names,
        std::vector<bool> & are_reverse,
        std::vector<int> & target_statuses,
        std::vector<bool> & are_from_first_bam);
};

    ProcessedPileup  cpp_process_pileup(
                        const std::string &,
                        const std::string &,
                        const int,
                        const std::string &,
                        const std::string &,
                        const int,
                        const int, 
                        const double,
                        const int,
                        const int,
                        const int,
                        const int,
                        const int,
                        const int,
                        const int,
                        const std::vector<std::string> &,
                        const std::vector<bool> &, 
                        const std::vector<std::string> &,
                        const std::vector<int> &, 
                        const std::vector<int> &,
                        const std::vector<std::string> &, 
                        const std::vector<std::vector<int>> &, 
                        const std::vector<int> &,
                        const std::vector<bool> &
    );

#endif
