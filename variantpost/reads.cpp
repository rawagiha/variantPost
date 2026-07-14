#include "util.h"
#include "reads.h"


std::ostream& operator<<(std::ostream& os, Rank rank) {
    switch (rank) {
        case Rank::Supporting:        return os << "Supporting";
        case Rank::LikelySupporting:  return os << "LikelySupporting";
        case Rank::UnlikelySupporting: return os <<"UnlikelySupporting";
        case Rank::NotSupporting:     return os << "NotSupporting";
        case Rank::Undetermined:      return os << "Undetermined";
        case Rank::NotAnalyzed:  return os << "NotAnalyzed";  
    }   
    return os << "NotAnalyzed";
}

inline int find_first_read_idx(const Ints& idx2pos, int target_pos, int start_offset, int end_offset) {
    if (idx2pos.empty()) return -1;
    
    auto begin_it = idx2pos.begin() + start_offset;
    auto end_it = idx2pos.end() - end_offset;
    
    auto it = std::lower_bound(begin_it, end_it, target_pos);
    
    if (it != end_it && *it == target_pos) {
        return static_cast<int>(it - idx2pos.begin());
    }
    return -1;
}

inline bool is_variant(const int idx, const Ints& var_idx) {
    return std::binary_search(var_idx.begin(), var_idx.end(), idx);
}

//------------------------------------------------------------------------------
// NOTE: initilizer list reordered to suppress warning
Read::Read(std::string_view name_, bool is_rv, 
           const std::string& cigar_str_, int start, int end, 
           std::string_view seq_,  std::string_view quals_, bool is_frm_frst) 
    : name(name_), base_quals(quals_), seq(seq_), cigar_str(cigar_str_), aln_start(start), aln_end(end),
      is_reverse(is_rv), is_control(is_frm_frst) {
    
    cigar_vector.reserve(8); // heuristic hard cord. 8 is sufficient in most cases 
    fill_cigar_vector(cigar_str_, cigar_vector);
    start_offset 
        = cigar_vector.front().first == 'S' ? cigar_vector.front().second : 0;
    end_offset 
        = cigar_vector.back().first == 'S' ? cigar_vector.back().second : 0;

    read_start = aln_start - start_offset; read_end = aln_end + end_offset;
}

//------------------------------------------------------------------------------
void Read::setReference(LocalReference& loc_ref) {
    const int loc_ref_len = static_cast<int>(loc_ref.seq.size());
    //int pos_2 = read_start; // for coordinate setup    
    
    Coord unmapped_segments;
    if (cigar_str.find('N') == std::string_view::npos 
       && cigar_str.find('D') == std::string_view::npos) {
        const int _idx = aln_start - loc_ref.start; // idx on local reference 
        const int _ref_len = aln_end - aln_start + 1; // expected refseq len
        if (_idx >= 0 && _idx + _ref_len <= loc_ref_len)
            ref_seq = loc_ref.seq.substr(_idx, _ref_len);
    } else {   
        int pos_1 = aln_start - 1; // for refseq setup
        int pos_2 = read_start; 
        
        size_t required_capacity = 0;
        for (const auto& [op, op_len] : cigar_vector) {
            if (op == 'M' || op == '=' || op == 'X' || op == 'D') required_capacity += op_len;
        }
        
        spliced_ref_seq.reserve(required_capacity);
        for (const auto& [op, op_len] : cigar_vector) {
            switch (op) {
                case 'M': case '=': case 'X': case 'D': {
                    if (op == 'D') {
                        unmapped_segments.emplace_back(pos_1, (pos_1 + op_len));
                    }
                    spliced_ref_seq += loc_ref.fasta.getSubSequence(loc_ref.chrom, pos_1, op_len);
                    pos_1 += op_len; pos_2 += op_len;
                    break;
                }
                case 'S':
                    pos_2 += op_len;
                    break;
                case 'N':
                    unmapped_segments.emplace_back(pos_2, (pos_2 + op_len -1));
                    skipped_segments.emplace_back(pos_2, (pos_2 + op_len -1));
                    pos_1 += op_len; pos_2 += op_len;
                    break;
                default: break;
            }        
        }
        ref_seq = spliced_ref_seq; // covert to string_view
    }

    // reference seq flags
    is_na_ref = (ref_seq.empty()); if (is_na_ref) return;
    is_ref = (seq == ref_seq);
    
    // aligned segments
    int current_pos = read_start; //<- I don't member why reset to read_start (keep for now, definitely investigate laater)
                                  //july 3. I think that this was to include cases here the target pos is only covered by softclipped read
                                  //in long indels. aligned_segments is misleading name. Reconder
    current_pos = aln_start;
    for (const auto& [skip_start, skip_end] : skipped_segments) {
        aligned_segments.emplace_back(current_pos, skip_start - 1);
        current_pos = skip_end + 1;
    }
    aligned_segments.emplace_back(current_pos, read_end);

    current_pos = aln_start;
    for (const auto& [un_start, un_end] : unmapped_segments) {
        mapped_segments.emplace_back(current_pos, un_start - 1);
        current_pos = un_end + 1;
    }
    mapped_segments.emplace_back(current_pos, aln_end);
}

//------------------------------------------------------------------------------
// NOTE: use only after setReference
void Read::setVariants(LocalReference& loc_ref) {
    if (is_ref) return;
    
    // from util.h
    read2variants(aln_start, ref_seq, seq, base_quals, 
                  cigar_vector, loc_ref, variants, var_idx, idx2pos);

    for (auto& v : variants) { 
        var_pos.push_back(v.pos);

        v.isInFlanking(loc_ref); 
        if (v.in_target_flnk) ++flnk_v_cnt;
    }
    
    // sorted -> var_idx will be changed.
    // This is okay. var_idx is only used for quality check
    if (variants.size() > 1) {
        std::sort(variants.begin(), variants.end()); 
        std::sort(var_pos.begin(), var_pos.end()); 
    }
}

bool Read::IsRefAt(const int pos) const {
    bool res = false;
    if (pos < covering_start || covering_end < pos) return res; 
    for (const auto& [start, end] : mapped_segments) {
        if (start <= pos && pos <= end) {
            res = (std::binary_search(var_pos.begin(), var_pos.end(), pos)) ? false : true;
            break;
        }
    }
    
    return res; 
}

//------------------------------------------------------------------------------
void Read::parseCoveringPattern(LocalReference& loc_ref, const Variant& target) {
    // skipped
    for (const auto& [seg_start, seg_end] : skipped_segments)
        if (seg_start < target.lpos && target.rpos < seg_end) return;
         
    if ((read_start <= target.lpos && target.lpos <= aln_start) ||
        (aln_end <= target.rpos && target.rpos <= read_end)) {
        covered_in_clip = true;
    }
       
    // completly covered
    for (const auto& [seg_start, seg_end] : aligned_segments) {
        if (seg_start <= target.lpos && target.rpos < seg_end) {
            covering_start = seg_start; covering_end = seg_end; 
            covering_pattern = CoveringPattern::Full;
            return;
        }
    }
             
    // partially covered (left- and right-partial)
    bool is_partial = false;
    covering_start = aligned_segments.front().first;
    covering_end = aligned_segments.back().second;

    if (target.lpos <= covering_start && covering_start <= target.end_pos) {
        covering_end = aligned_segments.front().second;
        // starts with softclip or has variants in shiftable region
        if (start_offset || (!variants.empty() && variants.front().pos <= target.end_pos)) is_partial = true; 
    } 
    else if (target.lpos <= covering_end && covering_end <= target.end_pos) {
        covering_start = aligned_segments.back().first;
        // ends with softclip or has variants in shiftable region
        if (end_offset || (!variants.empty() && target.lpos <= variants.back().pos)) is_partial = true;
    }
    if (is_partial) covering_pattern = CoveringPattern::Partial;
}

//------------------------------------------------------------------------------
//Find consecutive low-qaul bases from ends
//qs: quality segment start 
//qe: quality segment end
void Read::trimLowQualBases(const char qual_thresh) {
    const int n = static_cast<int>(base_quals.size());
    int i = 0, j = n - 1;
    while (i < n && base_quals[i] < qual_thresh) ++i;
    qs = i; 
    while (j >= i && base_quals[j] < qual_thresh) --j;
    qe = j;
}

//------------------------------------------------------------------------------
//Calc. freq of low quality non-ref bases between start/end
void Read::qualityCheck(const int start, const int end, 
                        const char qual_thresh, const double freq_thresh) {   
    if (is_ref) qc_passed = true;
    
    if (idx2pos.empty()) return;

    auto it_i = std::find_if(idx2pos.begin(), idx2pos.end(), [start](int pos){ return pos >= start; });
    auto it_j = std::find_if(idx2pos.rbegin(), idx2pos.rend(), [end](int pos){ return pos < end; });

    if (it_i == idx2pos.end() || it_j == idx2pos.rend()) return;
    
    int i = static_cast<int>(std::distance(idx2pos.begin(), it_i));
    int j = static_cast<int>(std::distance(idx2pos.begin(), it_j.base()) - 1);

    if (i >= j) return;
     
    int non_ref_cnt = 0, dirty_non_ref_cnt = 0;
    for (int k = i; k <= j; ++k) {
        if (is_variant(k, var_idx)) {
            ++non_ref_cnt;
            if (base_quals[k] < qual_thresh) ++dirty_non_ref_cnt; 
        }
    }
    
    qc_passed = non_ref_cnt ? (static_cast<double>(dirty_non_ref_cnt) / non_ref_cnt <= freq_thresh) : true;
}

//------------------------------------------------------------------------------
void Read::parseLocalPattern(LocalReference& loc_ref, 
                             const Variant& target, const UserParams& params) {
    if (!target.is_complex) {
        int idx = find_target(loc_ref, params, target, variants); 
        if (idx != -1) {
            has_target = true; target_idx = idx;
            auto& v = variants[idx];
            target_pos = v.pos; target_ref = v.ref; target_alt = v.alt;
            return;
        }
    }
     
    int kmer_size = params.kmer_size; // TODO consider if this is needed.
    
    int tmp_d = INT_MAX, dis_kmer = 0, anti_ptrn = 0, 
        total_base_changes = 0, sub_cnt = 0, indel_cnt = 0;
    for (size_t i = 0; i < variants.size(); ++i) {
        auto& v  = variants[i];

        if (v.is_substitute) ++sub_cnt; else ++indel_cnt;

        if (target.end_pos < v.lpos || v.end_pos < target.lpos) {
            tmp_d = std::abs(target.end_pos - v.lpos);
            if (dist_to_non_target > tmp_d) dist_to_non_target = tmp_d;
            tmp_d = std::abs(v.rpos - target.lpos);
            if (dist_to_non_target > tmp_d) dist_to_non_target = tmp_d;    
        } else { 
            tmp_d = 0;
        }
        if (dist_to_non_target > tmp_d) dist_to_non_target = tmp_d;
        
        if (v.in_target_flnk) {
            ++anti_ptrn; total_base_changes += v.event_len;
            
            // TODO Add logic here
            //if (!target.is_complex && v.is_shiftable) has_anti_pattern = true;
        }
        
        // ineffective kmer
        // ....var1..target...var2...
        // if targer is the middle of complex event cluster, kmers from target + ref
        // may not be matched
        // TODO: AS OF 2026.feb.2, not used. keep this? 
        // UPDATE: realn aginst ref with target included only for complex indels 
        //         with ineffective kmer?? (Feb.23.)
        if (i + 1 < variants.size() 
            && v.pos <= target.pos && target.pos <= variants[i + 1].pos) {
            if (target.pos - v.pos <= kmer_size 
                && variants[i + 1].pos - target.pos <= kmer_size) ++dis_kmer; 
        }
    }
    has_positional_overlap = (!dist_to_non_target); ineffective_kmer = (dis_kmer); 

    //clipping
    int dist_to_clip = INT_MAX;
    if (start_offset) {
        if (read_start <= target.lpos && target.lpos <= aln_start) {
            dist_to_clip = 0;
        } else {
            tmp_d = std::abs(target.pos - aln_start);
            if (dist_to_clip > tmp_d) dist_to_clip = tmp_d;
        }
    }
    if (end_offset) {
        if (aln_end <= target.rpos && target.rpos <= read_end) {  
            dist_to_clip = 0;
        } else {
            tmp_d = std::abs(target.rpos - aln_end); 
            if (dist_to_clip > tmp_d) dist_to_clip = tmp_d;
        }
    }

    if (dist_to_non_target <= target.event_radius) 
        has_local_events = true;
    if (dist_to_clip <= target.event_radius)
        has_local_clip = true;

    if (!has_local_clip) {
        // only 1 non-target event in flnk start/end (>1 may be complex-target) 
        if (anti_ptrn == 1) has_anti_pattern = true; 
        
        // smaller N of bases changed than target
        if (!target.is_complex 
            && target.event_len > total_base_changes) has_smaller_change = true;

        // indel target but no indel in flnk start/end
        if (!target.is_substitute && !indel_cnt && (anti_ptrn == sub_cnt)) 
            has_anti_pattern = true;
    }
}

void Read::checkByRepeatCount(const Variant& target, bool& has_excess_ins_hap) {
    if (!target.repeats || has_target) return;
    
    auto it = std::find(idx2pos.rbegin(), idx2pos.rend(), target.pos + 1);
    if (it == idx2pos.rend()) return;

    int idx = static_cast<int>(std::distance(idx2pos.begin(), it.base()) - 1);

    std::string_view rt_side = seq.substr(idx);
    std::string_view lt_side = seq.substr(0, idx);
    std::string_view indel_seq_sv = target.indel_seq;

    const int rt_len = static_cast<int>(rt_side.size());
    const int indel_len = target.indel_len;

    int rep = 0;
    for (int i = 0; rt_len - i >= indel_len; i += indel_len) {
        if (rt_side.substr(i, indel_len) == indel_seq_sv) ++rep; else break;
    }

    for (int i = idx - indel_len; i >= 0; i -= indel_len) {
        if (lt_side.substr(i, indel_len) == indel_seq_sv) ++rep; else break;
    }

    if (rep != target.repeats) {
        has_anti_pattern = true;
        covered_in_clip = false;
        if (target.is_ins && rep > target.repeats) {
            is_central_mapped = true;
            has_excess_ins_hap = true;
        }
    }
}

inline auto append_num = [](std::string& target, auto val) {
    std::array<char, 20> buffer; 
    auto [ptr, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), val);
    if (ec == std::errc{}) {
        target.append(buffer.data(), ptr - buffer.data());
    }
};


//------------------------------------------------------------------------------
// NOTE: qualityCheck before use
void Read::setSignatureStrings(const UserParams& params) {
    if (!qc_passed) return;
    
    non_ref_sig.clear();
    non_ref_sig.reserve(256); 

    for (const auto& v : variants) {
        bool dirty = false;
        for (const auto q : v.qual) {
            if (q < params.base_q_thresh) { dirty = true; break; }
        }
        if (dirty) continue;

        append_num(non_ref_sig, v.pos);
        non_ref_sig.append(":").append(v.ref).append(">").append(v.alt).append(";");        
    }

    splice_sig.clear();
    splice_sig.reserve(128);
    
    for (const auto& [s_start, s_end] : skipped_segments) {
        non_ref_sig.append("spl=");
        splice_sig.append("spl=");
    
        append_num(non_ref_sig, s_start);
        append_num(splice_sig, s_start);

        non_ref_sig.append("-");
        splice_sig.append("-");

        append_num(non_ref_sig, s_end);
        append_num(splice_sig, s_end);

        non_ref_sig.append("|");
        splice_sig.append("|");
    }

    if (start_offset) {
        non_ref_sig.append("lt=");
        append_num(non_ref_sig, aln_start - 1);
        non_ref_sig.append("$");
    }
    if (end_offset) {
        non_ref_sig.append("rt=");
        append_num(non_ref_sig, aln_end + 1);
    }    
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
// True if flanking start/end is mapped within
// Test if target locus is mapped in extreme ends (1/8 or 8/8)
// TODO use this for long-read cases?
void Read::isCenterMapped(const Variant& target) {
    
    int d = INT_MAX, i = -1, thresh = 0;
    if (is_ref) {
        i = target.pos - aln_start;
        if (i > 0) {
            thresh = (aln_end - aln_start + 1) / 8;
            is_central_mapped = (thresh <= i && i <= thresh * 7);
            return;
        }
    }
    
    if (idx2pos.empty()) return;
    thresh = static_cast<int>(idx2pos.size() / 8);
     
    for (int j = 0; j < static_cast<int>(idx2pos.size()); ++j) {
        if (std::abs(idx2pos[j] - target.pos) < d) {
            i = j; d = std::abs(idx2pos[j] - target.pos); 
        }
    }

    if (i < qs || qe < i) has_local_events = false; 
    
    is_central_mapped = (thresh <= i && i <= thresh * 7);             
}

//------------------------------------------------------------------------------
// use if input is complex indel
void Read::hasTargetComplexIndel(LocalReference& loc_ref, const Variant& target) {
    
    const size_t mid_len = target.mid_seq.size();
    
    bool left_match = false;
    auto left_it = seq.find(target.lt_seq);
    if (left_it != std::string_view::npos) {
        left_it += target.lt_seq.size();
        left_match = (seq.substr(left_it, mid_len) == target.mid_seq);
    }

    bool right_match = false;
    auto right_it = seq.find(target.rt_seq);
    if (right_it != std::string_view::npos && right_it >= mid_len) {
        right_match = (seq.substr(right_it - mid_len,  mid_len) == target.mid_seq);
    }
   
    if (left_match && right_match) { has_target = true; return; }
    return;

    int i = -1;
    //for (const auto& elem : idx2pos) {
    //    if (elem.second == loc_ref.flanking_start) { i = elem.first; break; }
    //}
    for (int j = 0; j < static_cast<int>(idx2pos.size()); ++j) {
        if (idx2pos[j] == loc_ref.flanking_start) i = j;
    }

    if (i < 0) return;
    
    int lt_len = static_cast<int>(target.lt_seq.size());
    int rt_len = static_cast<int>(target.rt_seq.size());
    
    const int seq_len = static_cast<int>(seq.size());
    
    if (seq_len < i + lt_len + target.alt_len) return;
    
    if (seq.substr(i, lt_len) != target.lt_seq) return;
    if (seq.substr(i + lt_len, target.alt_len) != target.mid_seq) return;
    if (seq.substr(i + lt_len + target.alt_len, rt_len) != target.rt_seq) return;
    has_target = true; 
}

//------------------------------------------------------------------------------
// use if input is complex substitute 
void Read::hasTargetComplexSubstitute(const Variant& target) {
    int i = -1;
    //for (const auto& elem : idx2pos) {
    //    if (elem.second == target.pos) { i = elem.first; break; }
    //}
    for (int j = 0; j < static_cast<int>(idx2pos.size()); ++j) {
        if (idx2pos[j] == target.pos) i = j;
    }
    if (i < 0) return;
    if (seq.substr(i, target.alt_len) == target.alt) has_target = true;
}
