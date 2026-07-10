#include "consensus.h"

const std::string REF = ".";
const std::string D = ">";

void fill_ref_alt(const std::string& refalt, Strs& ref, Strs& alt, 
                  const int pos, int& del_ref_end) {
    std::string ref_allele = refalt.substr(0, refalt.find(D));
    std::string alt_allele = refalt.substr(refalt.find(D) + 1);
    ref.push_back(refalt.substr(0, refalt.find(D)));
    alt.push_back(refalt.substr(refalt.find(D) + 1));
    del_ref_end = pos + static_cast<int>(ref_allele.size());
}

std::string find_commenest_change(const Strs& refalts) {
    std::unordered_map<std::string, int> cnts;
    for (const auto& elem : refalts) cnts[elem]++;
   
    auto max_iter = std::max_element(cnts.begin(), cnts.end(),
         [](const auto& a, const auto& b) { return a.second < b.second; });
    
    return max_iter->first;  
}

void Consensus::_from_variants(const int start_, const int end_, const Vars& vars, 
                               UserParams& params, LocalReference& loc_ref) {
    if (start_ < start) { start = start_; } if (end < end_) { end = end_; }
    
    //fill reference segments
    int curr = start_;
    for (const auto& var : vars) {
        while (curr < var.pos) {
            aln[curr++].push_back(REF);
        }
        
        if (var.mean_qual < params.base_q_thresh) aln[curr].push_back(REF);
        else aln[curr].push_back(var.ref + ">" + var.alt);
        curr = var._end_pos;
    }
    while (curr <= end_) { aln[curr++].push_back(REF); }
   
   int _cov = 0, ref_cnt = 0, del_ref_end = INT_MIN;
   for (const auto& pair : aln) {
       if (pair.first < loc_ref.start || loc_ref.end < pair.first) continue;    
       
       if (pair.first < del_ref_end) continue;

       pos.push_back(pair.first);
       cov.push_back(static_cast<int>(pair.second.size())); //cov is no-empty
       
       //ref.push_back(loc_ref.dict[pair.first]);
       ref_cnt = std::count(pair.second.begin(), pair.second.end(), REF);
       if (ref_cnt * 2 >= cov.back()) {
           std::string s = std::string(1, loc_ref.base_at(pair.first));
           ref.push_back(s); alt.push_back(ref.back());
       } else {
           std::string refaln = find_commenest_change(pair.second);
           fill_ref_alt(refaln, ref, alt, pair.first, del_ref_end);
       } 

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
