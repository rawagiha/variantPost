#include "util.h"
#include "reads.h"

//------------------------------------------------------------------------------
// NOTE: initilizer list reordered to suppress warning
Read::Read(std::string_view name_, const bool is_rv, 
           std::string_view cigar_str_, const int start, const int end, 
           std::string_view seq_, const std::string_view quals_, const bool is_frm_frst) 
    : name(name_), base_quals(quals_), aln_start(start), aln_end(end), seq(seq_),
      cigar_str(cigar_str_), is_reverse(is_rv), is_control(is_frm_frst) {
    cigar_vector = to_cigar_vector(cigar_str_);
    start_offset 
        = cigar_vector.front().first == 'S' ? cigar_vector.front().second : 0;
    end_offset 
        = cigar_vector.back().first == 'S' ? cigar_vector.back().second : 0;

    read_start = aln_start - start_offset; read_end = aln_end + end_offset;
}

//------------------------------------------------------------------------------
void Read::setReference(LocalReference& loc_ref) {
    int loc_ref_len = static_cast<int>(loc_ref.seq.size());
    int pos_2 = read_start; // for coordinate setup    
  
    if (cigar_str.find('N') == std::string_view::npos) {
        int _idx = aln_start - loc_ref.start; // idx on local reference 
        int _ref_len = aln_end - aln_start + 1; // expected refseq len
        if (_idx >= 0 && _idx + _ref_len <= loc_ref_len)
            ref_seq = loc_ref.seq.substr(_idx, _ref_len);
    }          
    else {   
        int pos_1 = aln_start - 1; // for refseq setup
        char op = '\0'; int op_len = 0;
        for (const auto& c : cigar_vector) {
            op = c.first; op_len = c.second;

            switch (op) {
                case 'M': case '=': case 'X': case 'D':
                    spliced_ref_seq += 
                        loc_ref.fasta.getSubSequence(loc_ref.chrom, pos_1, op_len);
                    pos_1 += op_len; pos_2 += op_len;
                    break;
                case 'S':
                    pos_2 += op_len;
                    break;
                case 'N':
                    skipped_segments.emplace_back(pos_2, (pos_2 + op_len -1));
                    pos_1 += op_len; pos_2 += op_len;
                    break;
                default:
                    break;
            }        
        }
        ref_seq = spliced_ref_seq; // covert to string_view
    }

    // reference seq flags
    is_na_ref = (ref_seq.empty()); if (is_na_ref) return;
    is_ref = (seq == ref_seq);
    
    //aligned segment
    pos_2 = read_start; // reset to read_start
    for (const auto& seg : skipped_segments) {
        aligned_segments.emplace_back(pos_2, seg.first - 1);
        pos_2 = seg.second + 1;
    }
    aligned_segments.emplace_back(pos_2, read_end);
}

//------------------------------------------------------------------------------
// NOTE: use only after setReference
void Read::setVariants(LocalReference& loc_ref) {
    if (is_ref) return;
    
    // from util.h
    read2variants(aln_start, ref_seq, seq, base_quals, 
                   cigar_vector, loc_ref.dict, variants, idx2pos);
}

//------------------------------------------------------------------------------
// 'A': target locus + shiftable region is completely covered
// 'B': partially
// 'C': not covered or skipped by splicing
void Read::parseCoveringPattern(LocalReference& loc_ref, const Variant& target) {
    // skipped
    for (const auto& seg : skipped_segments)
        if (seg.first < target.lpos && target.rpos < seg.second) return;
        
    // completly covered
    for (const auto& seg : aligned_segments) {
        if (seg.first <= target.lpos && target.rpos <= seg.second) {
            covering_start = seg.first; covering_end = seg.second; 
            covering_ptrn = 'A'; return;
        }
    }
             
    // partially covered (left- and right-partial)
    bool is_partial = false;
    covering_start = aligned_segments.front().first;
    covering_end = aligned_segments.back().second;

    if (target.lpos <= covering_start && covering_start <= target.end_pos) {
        covering_end = aligned_segments.front().second;
        // starts with softclip
        if (start_offset) is_partial = true;
        // has variants in shiftable
        if (!variants.empty() && variants.front().pos <= target.end_pos)
            is_partial = true;
    } 
    else if (target.lpos <= covering_end && covering_end <= target.end_pos) {
        covering_start = aligned_segments.back().first;
        // ends with softclip
        if (end_offset) is_partial = true;
        if (!variants.empty() && target.lpos <= variants.back().pos)
            is_partial = true;
    }
    if (is_partial) covering_ptrn = 'B';
}

//------------------------------------------------------------------------------
void Read::qualityCheck(const int start, const int end, const UserParams& params) {   
    
    const int last_ = static_cast<int>(base_quals.size()) - 1;
    
    int i = 0, j = last_; //index
    for (const auto& elem : idx2pos) {    
        if (!i && start <= elem.second) i = elem.first;
        if (j == last_ && end > elem.second) j = elem.first;  
    }
    if (i == j) return;
    
    int cnt = 0;
    for (int k = i; k <= j; ++k)
        if (base_quals[k] < params.base_q_thresh ) ++cnt;
    qc_passed = (static_cast<double>(cnt) / (j - i) <= params.lq_rate_thresh); 
}

//------------------------------------------------------------------------------
void Read::parseLocalPattern(LocalReference& loc_ref, const Variant& target) {
    int idx = 0; std::vector<int> dists;
    for (auto& v : variants) {
        if (target.isEquivalent(v, loc_ref)) {
            has_target = true; target_idx = idx;
            target_pos = v.pos; target_ref = v.ref; target_alt = v.alt;
            return;
        }
        
        // distances from target region to non-target variants
        v.setEndPos(loc_ref);
        if (target.end_pos < v.lpos || v.rpos < target.lpos) {
            int d1 = std::abs(target.end_pos - v.lpos);
            int d2 = std::abs(v.rpos - target.lpos);
            if (d1 < d2) dists.push_back(d1); else dists.push_back(d2);    
        } 
        else dists.push_back(0); // 0 if overlapped
        ++idx;
    }
    
    //clipping
    if (start_offset)
        dists.push_back(std::abs(target.pos - aln_start));
    if (end_offset)
        dists.push_back(std::abs(target.pos - aln_end));
    if (!dists.empty()) 
        dist_to_non_target = *std::min_element(dists.begin(), dists.end());
    if (dist_to_non_target <= target.event_len + target.rpos - target.lpos) 
        has_local_events = true;
}

inline bool is_dirty(const Variant& v, const char thresh) {
    for (const auto& q : v.qual) if (q < thresh) return true;
    return false; 
}

//------------------------------------------------------------------------------
// NOTE: qualityCheck before use
void Read::setSignatureStrings(const UserParams& params) {
    if(!qc_passed) return;

    for (const auto& v : variants) {
        if (is_dirty (v, params.base_q_thresh)) continue;
        non_ref_sig += (std::to_string(v.pos) + ":" + v.ref + ">" + v.alt + ";");
    }

    std::string s1, s2;
    for (const auto& s : skipped_segments) {
        s1 = std::to_string(s.first); s2 = std::to_string(s.second);
        non_ref_sig += ("spl=" + s1 + "-" + s2 + "|");
        splice_sig += ("spl=" + s1 + "-" + s2 + "|");
    }
    
    if (start_offset) 
        non_ref_sig += ("lt=" + std::to_string(aln_start - 1) + "$");
    
    if (end_offset) 
        non_ref_sig += ("rt=" + std::to_string(aln_end + 1));
}

//------------------------------------------------------------------------------
// non-ref alignments bounded by steep increase in dimer diversity (flanking)
void Read::isStableNonReferenceAlignment(LocalReference& loc_ref) {
    // fail to cove flanking boundaries defined by maxinum dimer-diverisity increase
    if (loc_ref.flanking_start < covering_start || covering_end < loc_ref.flanking_end) {
        fail_to_cover_flankings = true; return;
    } 

    // must have variants
    if (variants.empty() || !qc_passed) return;
    is_stable_non_ref = true;  
}

//------------------------------------------------------------------------------
// test if target locus is mapped in the middle of mapped read len tertiles
void Read::isCenterMapped(const Variant& target) {
    int d = INT_MAX, i = -1;
    for (const auto& elem : idx2pos) {
        if (std::abs(elem.second - target.pos) < d) {
            i = elem.first; d = std::abs(elem.second - target.pos);
        }
    }
   
    int tertile = static_cast<int>(idx2pos.size() / 3);
    is_central_mapped = (tertile <= i && i <= tertile * 2);             
}

//------------------------------------------------------------------------------
// use if input is complex indel
void Read::hasTargetComplexIndel(LocalReference& loc_ref, const Variant& target) {
    int i = -1;
    for (const auto& elem : idx2pos) {
        if (elem.second == loc_ref.flanking_start) { i = elem.first; break; }
    }
    if (i < 0) return;
    
    int lt_len = static_cast<int>(target.lt_seq.size());
    int rt_len = static_cast<int>(target.rt_seq.size());
    if (seq.substr(i, lt_len) != target.lt_seq) return;
    if (seq.substr(i + lt_len, target.alt_len) != target.mid_seq) return;
    if (seq.substr(i + lt_len + target.alt_len, rt_len) != target.rt_seq) return;
    has_target = true; 
}

//------------------------------------------------------------------------------
// use if input is complex substitute 
void Read::hasTargetComplexSubstitute(const Variant& target) {
    int i = -1;
    for (const auto& elem : idx2pos) {
        if (elem.second == target.pos) { i = elem.first; break; }
    }
    if (i < 0) return;
    if (seq.substr(i, target.alt_len) == target.alt) has_target = true;
}
