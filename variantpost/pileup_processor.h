#ifndef PILEUP_PROCESSOR_H
#define PILEUP_PROCESSOR_H

#include <string>
#include <vector>

namespace pp {

    struct ProcessedPileup 
    {
        std::string contig;
        int target_pos;
        std::string ref;
        std::string alt;
        std::vector<std::string> read_names;
        std::vector<bool> are_reverse;
        std::vector<bool> are_target;
        std::vector<bool> are_from_first_bam;

        ProcessedPileup();
                                                
        ProcessedPileup(const std::string & contig,
                        const int target_pos,
                        const std::string ref,
                        const std::string alt,
                        std::vector<std::string> & read_names,
                        std::vector<bool> & are_reverse,
                        std::vector<bool> & are_target,
                        std::vector<bool> & are_from_first_bam);
    };

    ProcessedPileup  process_pileup(const std::string &,
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
                                    //const std::string &, 
                                    const std::vector<std::string> &,
                                    const std::vector<bool> &, 
                                    const std::vector<std::string> &,
                                    const std::vector<int> &, 
                                    const std::vector<int> &,
                                    const std::vector<std::string> &, 
                                    //const std::vector<std::string> &,
                                    const std::vector<std::vector<int>> &, 
                                    const std::vector<int> &,
                                    const std::vector<bool> &
    );
} // end of namespace "pp (process pileup)"
#endif
