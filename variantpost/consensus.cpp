#include "search.h"
#include "consensus.h"

void from_consensus_variant_list(SearchResult& rslt, LocalReference& loc_ref,
                                 const int start, const int end, const Vars& consensus) {
    int curr = start;
    auto it = consensus.begin();
    while (it != consensus.end()) {
        while (curr < it->pos) {
            rslt.positions.push_back(curr);
            std::string ref_base = std::string(1, loc_ref.base_at(curr));
            rslt.ref_bases.push_back(ref_base);
            rslt.alt_bases.push_back(ref_base);
            ++curr;
        }

        if (curr > it->pos) {
            ++it;
            continue;
        }
        
        //Hereafter, curr == it->pos cases
        rslt.positions.push_back(it->pos);
        rslt.ref_bases.push_back(it->ref);
        rslt.alt_bases.push_back(it->alt);
        
        auto next_it = std::next(it);
        
        // This is the last variant
        if (next_it == consensus.end()) {
            curr = it->_end_pos;
            it = next_it; //exit loop as next_it is the ender
            continue;
        }
        
        // Current variant do not pass the next
        // Note: _end_pos points to curr + ref length
        //   atCGAtcg
        //   atGAGtcg
        //
        //   Or,complex indels followed by substitutes
        //
        //   atCGAtgc
        //   at--Gtgc
        if (it->_end_pos <= next_it->pos) {
            curr = it->_end_pos;
            it = next_it; // proceed to next (not the ender) 
            continue;
        } 

        // Tricky case where current variant passes the next (i,e,. pos(next) < _end_pos (curr))
        //  Complex indels with preceding mutations
        //    atCGAtcg   atC  GAtcg    atCGA  tcg
        //    atG--tcg   atCAT--tcg    atC--ATtcg
        
        int max_end_pos = it->_end_pos;
                
        while (next_it != consensus.end() && next_it->pos < max_end_pos) {
            rslt.positions.push_back(next_it->pos);
            rslt.ref_bases.push_back(next_it->ref);
            rslt.alt_bases.push_back(next_it->alt);
                    
            if (next_it->_end_pos > max_end_pos) {
                max_end_pos = next_it->_end_pos;
            }
            ++next_it;
        }
                
        curr = max_end_pos;
        it = next_it;
    }
    
    for (/**/; curr <= end; ++curr) {
        rslt.positions.push_back(curr);
        std::string ref_base = std::string(1, loc_ref.base_at(curr));
        rslt.ref_bases.push_back(ref_base);
        rslt.alt_bases.push_back(ref_base);
    }
}


void variant_consensus(const std::vector<Read>& reads, const Ints& idx, Vars& consensus) {
    
    consensus.clear();
    
    std::unordered_map<VariantKey, VarStat, VariantKeyHash> variant_map;
    std::map<int, int> coverage_diff;

    for (const int i : idx) {
        coverage_diff[reads[i].covering_start] += 1;
        coverage_diff[reads[i].covering_end + 1] -= 1;

        for (const auto& v : reads[i].variants) {
            VariantKey vk{v.pos, v.end_pos, v.ref, v.alt};
            
            auto [it, inserted] = variant_map.emplace(vk, VarStat{});
            if (inserted) {
                it->second.occurrence = 1;
            } else {
                it->second.occurrence++;
            }
        }
    }

    std::map<int, int> depth_map;
    int current_depth = 0;
    for (const auto& event : coverage_diff) {
        current_depth += event.second;
        depth_map[event.first] = current_depth;
    }

    auto get_coverage = [&depth_map](int pos) {
        auto it = depth_map.upper_bound(pos);
        if (it == depth_map.begin()) return 0;
        --it;
        return it->second;
    };

    for (auto& [vk, vs] : variant_map) {
        vs.read_depth = get_coverage(vk.pos);
    }

    std::map<int, std::vector<std::pair<VariantKey, VarStat>>> variants_by_locus;
    for (const auto& [vk, vs] : variant_map) {
        variants_by_locus[vk.pos].push_back({vk, vs});
    }

    for (const auto& [locus, key_stat_vec] : variants_by_locus) {
        int max_read_depth = 0;
        int total_alt_count = 0;
        
        const VariantKey* most_frequent_key = nullptr;
        const VarStat* most_frequent_stat = nullptr;

        for (const auto& [vk, vs] : key_stat_vec) {
            max_read_depth = std::max(max_read_depth, vs.read_depth);
            total_alt_count += vs.occurrence;

            if (most_frequent_stat == nullptr || vs.occurrence > most_frequent_stat->occurrence) {
                most_frequent_key = &vk;
                most_frequent_stat = &vs;
            }
        }

        if (max_read_depth < 1) continue;

        int ref_count = std::max(0, max_read_depth - total_alt_count);

        if (most_frequent_stat != nullptr && most_frequent_stat->occurrence > ref_count) {
            consensus.push_back(Variant{
                static_cast<int>(most_frequent_key->pos),
                most_frequent_key->ref,
                most_frequent_key->alt
            });
        }
    }
}
