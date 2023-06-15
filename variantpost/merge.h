#ifndef MERGER_H
#define MERGER_H


#include <string>
#include "ssw/ssw_cpp.h"


typedef StripedSmithWaterman::Alignment Alignment;
typedef StripedSmithWaterman::Aligner Aligner;
typedef StripedSmithWaterman::Filter Filter;


struct Seq
{
    std::string seq = "";
    std::string base_quals = "";
    int target_start = -1;

    Seq();

    Seq(
        const std::string& seq,
        const std::string& base_quals,
        int target_start = -1
    );

    bool empty() const;
};

Seq merge_reads(const std::vector<Seq>& inputs);

#endif
