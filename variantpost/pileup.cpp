#include "pileup.h"

//------------------------------------------------------------------------------
inline auto rarer = [](const auto& x, const auto& y) { return x.second < y.second; };

//------------------------------------------------------------------------------
Pileup::Pileup(const std::vector<std::string>& read_names,
               const std::vector<bool>& are_reverse, 
               const std::vector<std::string>& cigar_strs,
               const std::vector<int>& aln_starts,
               const std::vector<int>& aln_ends,
               const std::vector<std::string>& read_seqs,
               const std::vector<std::string>& quals,
               const std::vector<bool>& are_from_first_bam,
               const UserParams& params,
               LocalReference& loc_ref,
               Variant& target) {
    const int qc_start = target.lpos - params.local_thresh;
    const int qc_end = target.end_pos + params.local_thresh;
    
    sz = static_cast<int>(read_names.size()); reads.reserve(sz);
    for (int i = 0; i < sz; ++i) {
        reads.emplace_back(read_names[i], are_reverse[i], cigar_strs[i],
                           aln_starts[i], aln_ends[i], read_seqs[i],
                           quals[i], are_from_first_bam[i]);
        auto& read = reads[i];
        read.setReference(loc_ref); if (read.is_na_ref) continue;
        
        read.setVariants(loc_ref); 
        read.parseCoveringPattern(loc_ref, target);
        if (read.covering_ptrn == 'C') continue;
        
        read.parseLocalPattern(loc_ref, target);
        read.qualityCheck(qc_start, qc_end, params);
        
        if (loc_ref.has_flankings) 
            read.isStableNonReferenceAlignment(loc_ref); 

        // covered and has non_ref pattern 
        if (read.covering_ptrn == 'A' && !read.idx2pos.empty())
            read.isCenterMapped(target);
        
        // search complex target by string match
        if (target.is_complex         && 
            read.covering_ptrn == 'A' && read.is_stable_non_ref)
                read.hasTargetComplexVariant(loc_ref, target);
        
        // assign rank
        if (read.has_target) {
            read.setSignatureStrings(params); read.rank = 's'; ++s_cnt;
            freq_s[read.non_ref_sig]++;
            if (read.is_stable_non_ref) freq_s_h[read.non_ref_sig]++;
        }
        else if (read.has_local_events) {
            read.setSignatureStrings(params); read.rank = 'u'; ++u_cnt;
            if (read.is_central_mapped && read.is_stable_non_ref && 
                read.is_control) freq_u[read.non_ref_sig]++;
        }
        else { read.rank = 'n'; ++n_cnt; }      
    }
}
    
//------------------------------------------------------------------------------
// LOGIC:
// Given that target is already identified...
// 1. recurrently aligned as non-target
// 2. located around the middle of read
// 3. surrounded by complex sequence
// 4. found in second BAM (if applicable)
// -> Such patterns are inferred as non_target haplotypes. Up to 2 assuming diploid
void Pileup::setHaploTypeByFrequency() {
    if (!freq_s_h.empty())
        hap0 = std::max_element(freq_s_h.begin(), freq_s_h.end(), rarer)->first;
    else // freq_s must be non empty
        hap0 = std::max_element(freq_s.begin(), freq_s.end(), rarer)->first;

    std::string_view ptrn_1, ptrn_2; int cnt_1 = 0, cnt_2 = 0;
    for (const auto& elem : freq_u) {
        if (elem.second > cnt_1) {
            cnt_2 = cnt_1; ptrn_2 = ptrn_1;
            cnt_1 = elem.second; ptrn_1 = elem.first; 
        }
        else if (elem.second > cnt_2 && elem.second < cnt_1) {
            cnt_2 = elem.second; ptrn_2 = elem.first;
        }
    }
    if (cnt_1 > 1) { hap1 = ptrn_1; } if (cnt_2 > 1) { hap2 = ptrn_2; }
}

//------------------------------------------------------------------------------
void make_sequence(LocalReference& loc_ref, 
                   const std::vector<Variant>& variants, 
                   const int start, const int end, std::string& seq) {
    if (variants.empty()) return;

    seq.append(loc_ref._seq.substr(start - loc_ref.start, variants.front().pos - start));
    for (size_t i = 0; i < variants.size(); ++i) {
        seq.append(variants[i].alt);
        if (i < variants.size() - 1) {
            seq.append(loc_ref._seq.substr(variants[i]._end_pos - loc_ref.start, 
                                           variants[i + 1].pos - variants[i]._end_pos));
        }
    }
    seq.append(loc_ref._seq.substr(variants.back()._end_pos - loc_ref.start,
                                   end - variants.back()._end_pos));
}

//------------------------------------------------------------------------------
void Pileup::setSequenceFromHaplotype(LocalReference& loc_ref)
{
    int idx0 = -1, idx1 = -1, idx2 = -1;
    std::vector<int> starts, ends;
    for (int i = 0; i < sz; ++i) {
        auto& read = reads[i];
        if (read.rank == 's' && read.non_ref_sig == hap0) { 
            if (idx0 < 0) idx0 = i;
            starts.push_back(read.covering_start); ends.push_back(read.covering_end); 
        }
        else if (read.rank == 'u') {
            if (!hap1.empty()) {
                if (idx1 < 0 && read.non_ref_sig == hap1) {
                    idx1 = i; read.rank = 'n'; ++n_cnt; --u_cnt;
                }
            }
            if (!hap2.empty()) {
                if (idx2 < 0 && read.non_ref_sig == hap2) {
                    idx2 = i; read.rank = 'n'; ++n_cnt; --u_cnt;
                }
            }
        }
    }
    
    int start = *std::min_element(starts.begin(), starts.end()); 
    int end = *std::max_element(ends.begin(), ends.end());

    if (idx0 > -1)
        make_sequence(loc_ref, reads[idx0].variants, start, end, seq0);
    if (idx1 > -1) 
        make_sequence(loc_ref, reads[idx1].variants, start, end, seq1);
    if (idx2 > -1)
        make_sequence(loc_ref, reads[idx2].variants, start, end, seq2);
}

//------------------------------------------------------------------------------
void Pileup::differentialKmerAnalysis(UserParams& params, LocalReference& loc_ref) {
    Kmers km0, km1, km2, km12;
    
    int  flen = loc_ref.flanking_end - loc_ref.flanking_start;
    int sz = (flen / 2 > params.kmer_size) ? flen / 2  : params.kmer_size;
       
    make_kmers(seq0, sz, km0); make_kmers(seq1, sz, km1); make_kmers(seq2, sz, km2);
    
    std::set_union(km1.begin(), km1.end(), km2.begin(), km2.end(),
                   std::inserter(km12, km12.begin()));
     
    std::set_difference(km0.begin(), km0.end(), km12.begin(), km12.end(),
                        std::inserter(kmers0, kmers0.end()));
    std::set_difference(km12.begin(), km12.end(), km0.begin(), km0.end(), 
                        std::inserter(kmers12, kmers12.end()));
    
    for (auto& read : reads) {
        int neg = 0, pos = 0;
        if (read.rank == 'u') {
            for (const auto& kmer : kmers12) {
                if (read.seq.find(kmer) != std::string_view::npos) ++neg;
            }
            for (const auto& kmer : kmers0) {
                if (read.seq.find(kmer) != std::string_view::npos) ++pos;
            }
        }
        if (pos > 0 && !neg) { read.rank = 's'; --u_cnt; ++s_cnt; }
        if (neg > 0 && !pos) { read.rank = 'n'; --u_cnt; ++n_cnt; }
    }
}

/*
//------------------------------------------------------------------------------
void setSequenceFromPattern(LocalReference& loc_ref)
{
}

void parse_sigstr_to_seq(std::string_view ptrn, const Reads& reads) {
    for (const auto& read : reads) {
        if (read.non_ref_sig == ptrn) {
            for (const auto& v : read.variants)
                std::cout << v.pos << " " << v.ref << " " << v.alt << std::endl;
        }
    }
}

//------------------------------------------------------------------------------
void Pileup::setSupportingPattern() {
    std::unordered_map<std::string_view, int> freq;
    
    for (const auto& read : reads) {
        if (read.rank == 's') freq[read.non_ref_sig]++;
        if (!has_hiconf_target && read.is_central_mapped && read.is_stable_non_ref)
            has_hiconf_target = true;
    }
    
    supporting_ptrn = std::max_element(freq.begin(), freq.end(), smaller)->first;
}

//------------------------------------------------------------------------------
void Pileup::setLocalHaploTypes() {
    std::unordered_map<std::string_view, int> freq;
    
    for (const auto& read : reads)
        if (read.rank == 'u' && read.is_stable_non_ref &&      
            read.is_control  && read.is_central_mapped ) 
            freq[read.non_ref_sig]++;
    
    std::string_view ptrn_1, ptrn_2; int cnt_1 = 0, cnt_2 = 0;
    
    for (const auto& elem : freq) {
        if (elem.second > cnt_1) {
            cnt_2 = cnt_1; ptrn_2 = ptrn_1;
            cnt_1 = elem.second; ptrn_1 = elem.first; 
        }
        else if (elem.second > cnt_2 && elem.second < cnt_1) {
            cnt_2 = elem.second; ptrn_2 = elem.first;
        }
    }
    
    // LOGIC:
    // Given that target is already identified...
    //
    // 1. recurrently aligned as non-target
    // 2. located around the middle of read
    // 3. surrounded by complex sequence
    // 4. found in second BAM (if applicable)
    // 
    // -> Such patterns to be used to recover the background haplotype
    //    Up to 2 assuming diploid
    if (cnt_1 > 1) hap1 = ptrn_1; 
    if (cnt_2 > 2) hap2 = ptrn_2;

    // demote if haplotypes are inferred with high confident target
    if ((hap1.empty() && hap2.empty()) || !has_hiconf_target) return;
    
    for (auto& read : reads) {
        if (read.rank == 'u') {
            if (read.non_ref_sig == hap1 || read.non_ref_sig == hap2) {
                read.rank = 'n'; --und; ++non_suppr;
            }
        }
    }    
}

//------------------------------------------------------------------------------
void setSequenceFromPattern(LocalReference& loc_ref)
{
    std::vector<int> starts, ends;
    for (const auto& read : reads) {
            
    }
}
*/

/*
//helper
void find_start_end_idx(const std::pair<int, int>& idx2pos,
                        const int start_pos, const int end_pos, 
                        int& start_idx, int& end_idx)
{
    for (const auto& elem : idx2pos)
    {
        if (!start_idx && elem.second == start_pos) start_idx = elem.first;
        if (!end_idx && elem.second == end_pos) end_idx = elem.first;
        if (start_idx && end_idx) return;
    }
} 

void compareStableSegmentLen(LocalReference& loc_ref)
{
    size_t expected_len = 0;
    for (const auto& read : reads)
    {
        if (read.rank == 's' && read.is_stable_non_ref
            read.non_ref_sig == supporting_ptrn)
        {   
            const auto& i2p = read.idx2pos; size_t i = 0, j = 0;
            for (const auto& e : i2p)
            {
                if (!i && i2p.second == loc_ref.flanking_start) i = e.first;
                if (!j && i2p.second == loc_ref.flanking_end) j = e.first;
            }
            if (i && j)
            {
                expected_len = j - i + 1; break;
            }
        }
    }
    
    if (!expected_len) return;
    
    for (auto& read : reads)
    { 
        if (read.rank == 'u' && read.is_stable_non_ref)
        {
            const auto& i2p = read.idx2pos; size_t i = 0, j = 0;
        }
    } 
}
*/

