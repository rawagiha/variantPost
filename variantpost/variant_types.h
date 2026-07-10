#ifndef VARIANT_TYPES_H
#define VARIANT_TYPES_H

#include "util.h"
#include <string>
#include <functional> 


// For haplotype inference
enum class Phase {
    Unphased,
    Hap1,
    Hap2
};

struct VariantKey {
    long long pos;
    long long end_pos;
    std::string ref;
    std::string alt;

    bool operator==(const VariantKey& o) const {
        return pos == o.pos && end_pos == o.end_pos && ref == o.ref && alt == o.alt;
    }
};

struct VarStat {
    int occurrence = 0; 
    int quality_occurrence = 0; 
    int read_depth = 0;
    int ref_depth = 0;
    double vaf = 0.0;
};

struct VariantKeyHash {
    std::size_t operator()(const VariantKey& k) const {
        std::size_t h1 = std::hash<long long>{}(k.pos);
        std::size_t h2 = std::hash<long long>{}(k.end_pos);
        std::size_t h3 = std::hash<std::string>{}(k.ref);
        std::size_t h4 = std::hash<std::string>{}(k.alt);
        
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};

#endif 
