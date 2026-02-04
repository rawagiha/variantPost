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
           std::string s = std::string(loc_ref.dict[pair.first]);
           ref.push_back(s); alt.push_back(ref.back());
       } else {
           std::string refaln = find_commenest_change(pair.second);
           fill_ref_alt(refaln, ref, alt, pair.first, del_ref_end);
       } 

   }
}
