#include "util.h"
#include "reads.h"

//------------------------------------------------------------------------------
// NOTE: initilizer list reordered to suppress warning
Read::Read(std::string_view name_, const bool is_rv, 
           const std::string& cigar_str_, const int start, const int end, 
           std::string_view seq_, const std::string_view quals_, const bool is_frm_frst) 
    : name(name_), base_quals(quals_), aln_start(start), aln_end(end), seq(seq_),
      cigar_str(cigar_str_), is_reverse(is_rv), is_control(is_frm_frst) {
    fill_cigar_vector(cigar_str_, cigar_vector);
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
    } else {   
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
                   cigar_vector, loc_ref.dict, variants, var_idx, idx2pos);

    if (variants.size() > 1)
        std::sort(variants.begin(), variants.end()); 
}

//------------------------------------------------------------------------------
// 'A': target locus + shiftable region is completely covered
// 'B': partially
// 'C': not covered or skipped by splicing
void Read::parseCoveringPattern(LocalReference& loc_ref, const Variant& target) {
    // skipped
    for (const auto& seg : skipped_segments)
        if (seg.first < target.lpos && target.rpos < seg.second) return;
         
    if (read_start <= target.lpos && target.lpos <= aln_start)
        covered_in_clip = true;
    else if (aln_end <= target.rpos && target.rpos <= read_end)
        covered_in_clip = true;
       
    // completly covered
    for (const auto& seg : aligned_segments) {
        if (seg.first <= target.lpos && target.rpos < seg.second) {
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

inline bool is_variant(const int idx, const Ints& var_idx) {
    auto it = std::find(var_idx.begin(), var_idx.end(), idx);
    return (it != var_idx.end()); 
}

//------------------------------------------------------------------------------
void Read::trimLowQualBases(const char qual_thresh) {
    int i = 0, n = static_cast<int>(base_quals.size()), j = n - 1;
    bool is_low_qual = true;
    while (is_low_qual && i < n) {
        is_low_qual = (base_quals[i] < qual_thresh); ++i;
    } qs = (--i);
    
    is_low_qual = true;
    while (is_low_qual && j >= 0) {
        is_low_qual = (base_quals[j] < qual_thresh); --j;
    } qe = (++j);
}

//------------------------------------------------------------------------------
void Read::qualityCheck(const int start, const int end, 
                        const char qual_thresh, const double freq_thresh) {   
    int i = 0, j = 0; //index
    for (const auto& elem : idx2pos) {    
        if (!i && start <= elem.second) i = elem.first;
        if (elem.second < end) j = elem.first;  
    }
    if (i == j) return;
     
    int non_ref_cnt = 0, dirty_non_ref_cnt = 0;
    for (int k = i; k <= j; ++k) {
        if (is_variant(k, var_idx)) {
            ++non_ref_cnt;
            if (base_quals[k] < qual_thresh) ++dirty_non_ref_cnt; 
        }
    }
    
    if (dirty_non_ref_cnt) 
        qc_passed = (static_cast<double>(dirty_non_ref_cnt) / non_ref_cnt <= freq_thresh);
    else qc_passed = true;
}


//------------------------------------------------------------------------------
void Read::parseLocalPattern(LocalReference& loc_ref, 
                             const Variant& target, const int kmer_size) {
    if (!target.is_complex) {
        int idx = find_target(loc_ref, target, variants); 
        if (idx != -1) {
            has_target = true; target_idx = idx;
            auto& v = variants[idx];
            target_pos = v.pos; target_ref = v.ref; target_alt = v.alt;
            return;
        }
    }
     
    int tmp_d = INT_MAX, dis_kmer = 0, anti_ptrn = 0, total_base_changes = 0;
    for (size_t i = 0; i < variants.size(); ++i) {
        // distances from target region to non-target variants
        // note this may not be found by pos comparison alone 
        // for example for long deletions
        auto& v = variants[i]; v.setEndPos(loc_ref);
        if (target.end_pos < v.lpos || v.end_pos < target.lpos) {
            tmp_d = std::abs(target.end_pos - v.lpos);
            if (dist_to_non_target > tmp_d) dist_to_non_target = tmp_d;
            tmp_d = std::abs(v.rpos - target.lpos);
            if (dist_to_non_target > tmp_d) dist_to_non_target = tmp_d;    
        } else { 
            tmp_d = 0;
        }
        if (dist_to_non_target > tmp_d) dist_to_non_target = tmp_d;
        
        if (loc_ref.flanking_start <= v.pos && v._end_pos <= loc_ref.flanking_end) {
            ++anti_ptrn; total_base_changes += v.event_len;
        }
        
        // ineffective kmer
        // ....var1..target...var2...
        // if targer is the middle of complex event cluster, kmers from target + ref
        // may not be matched
        // TODO: AS OF 2026.feb.2, not used. keep this? 
        if (i + 1 < variants.size() 
            && v.pos <= target.pos && target.pos <= variants[i + 1].pos) {
            if (target.pos - v.pos <= kmer_size 
                && variants[i + 1].pos - target.pos <= kmer_size) ++dis_kmer; 
        }
    }
    has_positional_overlap = (!dist_to_non_target); ineffective_kmer = (dis_kmer); 
    has_anti_pattern = (anti_ptrn == 1); // > 1 may be complex 
    has_smaller_change = (target.event_len > total_base_changes); 

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

    const int event_radius = target.event_len + target.rpos - target.lpos;
    if (dist_to_non_target <= event_radius) 
        has_local_events = true;
    if (dist_to_clip <= event_radius)
        has_local_clip = true;
}

void Read::checkByRepeatCount(const Variant& target, bool& has_excess_ins_hap) {
    if (!target.repeats || has_target) return;
    
    int idx = -1;
    for (const auto& elem : idx2pos)  
        if (elem.second == target.pos + 1) { idx = elem.first; break; }
    if (idx < 0) return;

    std::string_view rt_side = seq.substr(idx), 
                     lt_side = seq.substr(0, idx),
                     indel_seq_sv = target.indel_seq;
    
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
        //repeats are counted regardless of clipping -> deactivated
        covered_in_clip = false; 
        //more ins repeats -> relative read loc. does not matter  
        if (target.is_ins && rep > target.repeats) {
            is_central_mapped = true; has_excess_ins_hap = true; 
        }
    }
}

inline bool is_dirty(const Variant& v, const char thresh) {
    for (const auto& q : v.qual) { if (q < thresh) return true; }
    return false; 
}

//------------------------------------------------------------------------------
// NOTE: qualityCheck before use
void Read::setSignatureStrings(const UserParams& params) {
    if(!qc_passed) return;

    for (const auto& v : variants) {
        // use only clean variants for signature
        if (is_dirty(v, params.base_q_thresh)) continue;
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
// RELAXED: test if target locus is mapped in 2nd/3rd quartile
void Read::isCenterMapped(const Variant& target) {
    int d = INT_MAX, i = -1;
    for (const auto& elem : idx2pos) {
        if (std::abs(elem.second - target.pos) < d) {
            i = elem.first; d = std::abs(elem.second - target.pos);
        }
    }
    
    if (i < qs || qe < i) has_local_events = false; 
    
    int quantile = static_cast<int>(idx2pos.size() / 4);
    is_central_mapped = (quantile <= i && i <= quantile * 3);             
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
