#ifndef PILEUP_PROCESSOR_H
#define PILEUP_PROCESSOR_H

#include <string>
#include <vector>


namespace pp {


std::string  process_pileup(const std::string &,
                   const std::string &,
                   int,
                   const std::string &,
                   const std::string &,
                   int, 
                   int,
                   int,
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
