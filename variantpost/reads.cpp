#include "util.h"
#include "reads.h"

//------------------------------------------------------------------------------
// helper function used in Read::parseCoveringPattern
inline bool is_ascending(const int a, const int b, const int c)
{
    return (a <= b && b <= c);
}

//------------------------------------------------------------------------------
// NOTE: initilizer list reordered to suppress warning
Read::Read(std::string_view name_, const bool is_rv, 
           std::string_view cigar_str_, const int start, const int end, 
           std::string_view seq_, const std::vector<int>& quals_, const int mapq_, 
           const bool is_frm_frst) 
    : name(name_), mapq(mapq_), aln_start(start), aln_end(end), seq(seq_), 
      cigar_str(cigar_str_), is_reverse(is_rv), is_from_first_bam(is_frm_frst)
{
    for (const auto& q : quals_)
        base_quals.push_back(static_cast<char>(q + 33));
    
    cigar_vector = to_cigar_vector(cigar_str_);
    
    start_offset 
        = cigar_vector.front().first == 'S' ? cigar_vector.front().second : 0;
    end_offset 
        = cigar_vector.back().first == 'S' ? cigar_vector.back().second : 0;

    read_start = aln_start - start_offset; read_end = aln_end + end_offset;
}

//------------------------------------------------------------------------------
void Read::setReference(LocalReference& loc_ref)
{
    int loc_ref_len = static_cast<int>(loc_ref.seq.size());
    int pos_2 = read_start; // for coordinate setup    
  
    if (cigar_str.find('N') == std::string::npos)
    {
        int _idx = aln_start - loc_ref.start; // idx on local reference 
        int _ref_len = aln_end - aln_start + 1; // expected refseq len
        if (_idx >= 0 && _idx + _ref_len <= loc_ref_len)
            ref_seq = loc_ref.seq.substr(_idx, _ref_len);
    }          
    else//spliced reference 
    {   
        int pos_1 = aln_start - 1; // for refseq setup
        char op = '\0'; int op_len = 0;
        for (const auto& c : cigar_vector)
        {
            op = c.first; op_len = c.second;

            switch (op)
            {
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
    for (const auto& seg : skipped_segments)
    {
        aligned_segments.emplace_back(pos_2, seg.first - 1);
        pos_2 = seg.second + 1;
    }
    aligned_segments.emplace_back(pos_2, read_end);
}

//------------------------------------------------------------------------------
// NOTE: use only after setReference
void Read::setVariants(LocalReference& loc_ref)
{
    if (is_ref) return;
    
    // from util.h
    read2variants(aln_start, ref_seq, seq, base_quals, 
                   cigar_vector, loc_ref.dict, variants, idx2pos);
}

//------------------------------------------------------------------------------
// 'A': target locus + shiftable region is completely covered
// 'B': partially
// 'C': not covered or skipped by splicing
void Read::parseCoveringPattern(LocalReference& loc_ref, const Variant& target)
{
    // skipped
    for (const auto& seg : skipped_segments)
        if (seg.first < target.lpos && target.rpos < seg.second) return;
        
    // completly covered
    for (const auto& seg : aligned_segments)
    {
        if (seg.first <= target.lpos && target.rpos <= seg.second)
        {
            covering_start = seg.first; covering_end = seg.second; 
            covering_ptrn = 'A'; return;
        }
    }
             
    // partially covered (left- and right-partial)
    bool is_partial = false;
    covering_start = aligned_segments.front().first;
    covering_end = aligned_segments.back().second;

    if (is_ascending(target.lpos, covering_start, target.end_pos)) 
    {
        covering_end = aligned_segments.front().second;
        // starts with softclip
        if (start_offset) is_partial = true;
        // has variants in shiftable
        if (!variants.empty() && variants.front().pos <= target.end_pos)
            is_partial = true;
    } 
    else if (is_ascending(target.lpos, covering_end, target.end_pos))
    {
        covering_start = aligned_segments.back().first;
        // ends with softclip
        if (end_offset) is_partial = true;
        if (!variants.empty() && target.lpos <= variants.back().pos)
            is_partial = true;
    }

    if (is_partial) covering_ptrn = 'B';
}

//------------------------------------------------------------------------------
void Read::qualityCheck(const int start, const int end, const UserParams& params)
{   
    
    const int last_ = static_cast<int>(base_quals.size()) - 1;
    
    int i = 0, j = last_; //index
    for (const auto& elem : idx2pos)
    {    
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
void Read::parseLocalPattern(LocalReference& loc_ref, const Variant& target)
{
    int idx = 0; std::vector<int> dists;
    for (auto& v : variants)
    {
        if (target.isEquivalent(v, loc_ref))
        {
            has_target = true; target_idx = idx;
            target_pos = v.pos; target_ref = v.ref; target_alt = v.alt;
            return;
        }
        
        // distances from target region to non-target variants
        v.setEndPos(loc_ref);
        if (v.lpos < target.lpos || target.end_pos < v.rpos)
        {
            int lt_d = target.lpos - v.lpos;
            int rt_d = v.rpos - target.end_pos;
            if (lt_d < rt_d) dists.push_back(lt_d);
            else dists.push_back(rt_d);    
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
}

//------------------------------------------------------------------------------
// NOTE: qualityCheck before use
void Read::setSignatureStrings()
{
    if(!qc_passed) return;

    for (const auto& v : variants) 
    {
        non_ref_sig += (std::to_string(v.pos) + "_" + v.ref + "_" + v.alt + ";");
    }

    std::string s1, s2;
    for (const auto& s : skipped_segments) 
    {
        s1 = std::to_string(s.first); s2 = std::to_string(s.second);
        non_ref_sig += ("spl=" + s1 + "-" + s2 + ";");
        splice_sig += ("spl=" + s1 + "-" + s2 + ";");
    }
    
    if (start_offset) 
        non_ref_sig += ("lt_clip_end=" + std::to_string(aln_start - 1) + ";");
    
    if (end_offset) 
        non_ref_sig += ("rt_clip_start=" + std::to_string(aln_end + 1));
}

//------------------------------------------------------------------------------
// non-ref alignments bounded by steep increase in dimer diversity (flanking)
void Read::isStableNonReferenceAlignment(LocalReference& loc_ref)
{
    
    // fail to cove flanking boundaries defined by maxinum dimer-diverisity increase
    if (loc_ref.flanking_start < covering_start || covering_end < loc_ref.flanking_end) 
        fail_to_cover_flankings = true;
     
    // must have variants 
    if (is_ref      || fail_to_cover_flankings || 
        !qc_passed  || start_offset            || end_offset) return;
    
    // no variants out of flankings
    if (variants.front().pos < loc_ref.flanking_start ||
        loc_ref.flanking_end < variants.back().pos) return;
    
    is_stable_non_ref = true;  
}


//------------------------------------------------------------------------------
// only applicable if input is complex
void Read::hasTargetComplexVariant(LocalReference& loc_ref, const Variant& target)
{
    
    int i = -1;
    for (const auto& elem : idx2pos)
    {
        if (elem.second == loc_ref.flanking_start) 
        {    
            i = elem.first; break;
        }
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

//------------------------------------------------------------------------------
void prep_reads(const std::vector<std::string>& read_names,
                const std::vector<bool>& are_reverse, 
                const std::vector<std::string>& cigar_strs,
                const std::vector<int>& aln_starts,
                const std::vector<int>& aln_ends,
                const std::vector<std::string>& read_seqs,
                const std::vector<std::vector<int>>& quals,
                const std::vector<int>& mapqs,
                const std::vector<bool>& are_from_first_bam,
                LocalReference& loc_ref,
                Variant& target,
                const UserParams& params,
                Reads& reads)
{
    
    const int qc_start = target.lpos - params.local_thresh;
    const int qc_end = target.end_pos + params.local_thresh;
    size_t n_reads = read_names.size(); reads.reserve(n_reads);
    for (size_t i = 0; i < n_reads; ++i)
    {
        if (mapqs[i] < params.mapq_thresh) continue; //can do this at Python?
        
        reads.emplace_back(read_names[i], are_reverse[i], cigar_strs[i],
                           aln_starts[i], aln_ends[i], read_seqs[i],
                           quals[i], mapqs[i], are_from_first_bam[i]);

        reads[i].setReference(loc_ref); 
        
        if (reads[i].is_na_ref) continue;
        reads[i].setVariants(loc_ref); 

        reads[i].parseCoveringPattern(loc_ref, target);
        if (reads[i].covering_ptrn == 'C') continue;
        
        reads[i].parseLocalPattern(loc_ref, target);
        reads[i].qualityCheck(qc_start, qc_end, params); 
        
        if (loc_ref.has_flankings)
            reads[i].isStableNonReferenceAlignment(loc_ref); 
        
        reads[i].is_analyzable = true;
        
        //
        if (target.is_complex && reads[i].covering_ptrn == 'A' && reads[i].is_stable_non_ref)
                reads[i].hasTargetComplexVariant(loc_ref, target);
                
        reads[i].is_analyzable = true;
    }
    reads.shrink_to_fit();
}

//------------------------------------------------------------------------------
void triage_reads(Reads& reads, Reads& supportings, Reads& candidates, 
                  Reads& non_supportings, UserParams& params)
{
    size_t buff = reads.size(); 
    supportings.reserve(buff); candidates.reserve(buff); non_supportings.reserve(buff);

    for (size_t i = 0; i < buff; ++i)
    {
        if (!reads[i].is_analyzable) continue;

        if(reads[i].has_target)
        {
            reads[i].setSignatureStrings(); transfer_elem(supportings, reads, i);
        }
        else if (reads[i].dist_to_non_target < params.local_thresh)
        {   
            reads[i].setSignatureStrings(); transfer_elem(candidates, reads, i);
        }
        else 
            transfer_elem(non_supportings, reads, i);
    }

    reads.clear(); 
    supportings.shrink_to_fit(); candidates.shrink_to_fit(); non_supportings.shrink_to_fit();

    /*std::cout << supportings.size() << " " << candidates.size() << " " << non_supportings.size() << " " << std::endl;

    for (auto r : non_supportings)
        std::cout << r.is_reverse << " " << r.name << " " << r.covering_ptrn << std::endl;
    */
}



//------------------------------------------------------------------------------
// QC
//void Read::qualityCheck(const Variant& target, const UserParams& params)
//{
    /**/
//}

/*
std::string_view get_unspliced_ref_seq(
    const int loc_ref_start, 
    std::string_view loc_ref_seq,
    const int aln_start, 
    const int aln_end
)
{
  int start_idx = aln_start - loc_ref_start;
  size_t expected_ref_len = aln_end - aln_start + 1;
  
  if (0 <= start_idx && start_idx <= int(loc_ref_seq.size())) 
  {
    std::string_view fitted_ref = loc_ref_seq.substr(start_idx, expected_ref_len);
    if (fitted_ref.size() == expected_ref_len) //may be short near rt end
    {
        return fitted_ref;
    }
  }

  // failed to fit
  return "";
}


std::string_view get_spliced_ref_seq(
    std::string& ref_seq,
    const std::string& chrom, 
    FastaReference & fr,
    const int aln_start, 
    const std::vector<std::pair<char, int>> & cigar_vector
)
{
    int curr_pos = aln_start - 1;
    
    char op;
    int op_len;
    for (const auto & c : cigar_vector) {
        op = c.first;
        op_len = c.second;
        
        switch (op) 
        {
            case 'M':
            case 'X':
            case 'D':
                ref_seq += fr.getSubSequence(chrom, curr_pos, op_len);
                curr_pos += op_len;
                break;
            case 'N':
                curr_pos += op_len;
                break;
            default:
                break;
        }

    }
    
    return static_cast<std::string_view>(ref_seq);     
}


void annot_ref_seq(Read& read, LocalReference& loc_ref)
{ 
    if (read.cigar_str.find('N') != std::string::npos)
    {
        read.ref_seq = get_spliced_ref_seq(
            read.spliced_ref_seq,
            loc_ref.chrom, loc_ref.fasta, 
            read.aln_start, read.cigar_vector
        );
    }
    else
    {
        read.ref_seq = get_unspliced_ref_seq(
            loc_ref.start, loc_ref.seq,
            read.aln_start, read.aln_end
        );
    }
    
    read.is_ref = (read.ref_seq == read.seq);
    read.is_na_ref = (read.ref_seq.empty());
}


void annot_splice_pattern(Read& read)
{
    if (read.is_na_ref) return;
    
    char op;
    int op_len;
    int curr_pos = read.read_start;
    for (const auto& c : read.cigar_vector) 
    {
        op = c.first;
        op_len = c.second;
        
        switch (op) 
        {
            case 'M':
            case 'S': //read_start used
            case 'D':
            case 'X':
            case '=':
                curr_pos += op_len;
                break;
            case 'N':
                read.skipped_segments.emplace_back(
                    curr_pos, (curr_pos + op_len - 1)
                );
                curr_pos += op_len;
                break;
            default:
                break;
        }
    }

    curr_pos = read.read_start;
    for (const auto& segment : read.skipped_segments) 
    {
        read.aligned_segments.emplace_back(curr_pos, segment.first - 1);
        curr_pos = segment.second + 1;
    }
    read.aligned_segments.emplace_back(curr_pos, read.read_end);
}   

inline bool is_ascending(const int a, const int b, const int c)
{
    return (a <= b && b <= c);
}
*/

/*                                                            
void covering_patterns(Read& read, const Variant& target)   
{
    if (!read.skipped_segments.empty()) 
    {
        for (const auto& segment : read.skipped_segments) 
        {
            if (segment.first < target.lpos && target.rpos < segment.second) 
            {
                // this read skips the target locus 
                read.covering_ptrn = 'X';
                read.is_na_ref = true;
                return;
            }
        }
    }
    
    for (const auto& segment : read.aligned_segments) 
    {
        read.covering_start = segment.first;
        read.covering_end = segment.second;
 
        if (!target.is_shiftable)
        {
            if (is_ascending(read.covering_start, target.pos, read.covering_end)) 
            {
                //this read covers the target locus
                read.covering_ptrn = 'A';
                return;
            }
        } 
        //shiftable indel 
        else  
        {
            // both bound
            if (read.covering_start <= target.lpos 
                && target.rpos + 1 <= read.covering_end) 
            {
                read.covering_ptrn = 'A';
                return;
            }
            // left unbound
            else if (is_ascending(target.lpos + 1, read.covering_start, target.rpos)) 
            {
                // first & lt-clipped
                if (&segment == &read.aligned_segments.front() && read.start_offset) 
                {
                     //out of shift but needs further check 
                     read.covering_ptrn = 'B';
                     return;
                }

                // variants exist between start and rpos
                for (const auto& variant : read.variants) 
                {
                    if (is_ascending(read.covering_start, variant.pos, target.rpos)) 
                    {
                        read.covering_ptrn = 'B'; 
                        return;
                    }
                }
                
                // unbound with reference bases
                // out of shift no further check needed
                read.covering_ptrn = 'C';
                return;
            }
            // right unbound
            else if (is_ascending(target.lpos, read.covering_end, target.rpos)) 
            {
                // last & rt-clipped
                if (&segment == &read.aligned_segments.back() && read.end_offset) 
                {
                    read.covering_ptrn = 'B';
                    return;
                }
                
                // variants exist between lpos and stop
                for (const auto& variant : read.variants) 
                {
                    if (is_ascending(target.lpos, variant.pos, read.covering_end)) 
                    {
                        read.covering_ptrn = 'B';
                        return; 
                    }
                }
                
                // right reference bases only -> undetermined
                read.covering_ptrn = 'C';
                return;
            }
        }
    }
    
    // this read dose not cover the target locus
    // but may contain target (eg, large del)
    read.covering_ptrn = 'D';
}


void annot_covering_ptrn(
    Read& read, 
    const Variant& target, 
    LocalReference& loc_ref,
    bool is_retargeted
)
{    
    if (read.is_na_ref)
    {
        read.covering_ptrn = 'X';
        return;
    }
    else
    {
        if (!read.is_ref)
        {
            if (!is_retargeted)
            {
                parse_variants(
                    read.aln_start, read.ref_seq, 
                    read.seq, read.base_quals, read.cigar_vector, 
                    loc_ref.dict, read.variants, read.non_ref_quals
                );
            }
         }
    }

    covering_patterns(read, target);

    if (read.covering_ptrn == 'A')
    {
        if (read.aln_start <= target.pos 
            && target.pos <= read.aln_end) read.is_tight_covering = true;
        
        // experimental central score 
        if (target.pos - read.covering_start > read.covering_end - target.pos)
        {
            read.central_score = (read.covering_end - target.pos) 
                / static_cast<double>(read.covering_end - read.covering_start);
        }
        else
        {
            read.central_score = (target.pos - read.covering_start)
                / static_cast<double>(read.covering_end - read.covering_start);
        }
    }     
}


void annot_target_info(
    Read& read, 
    const Variant& target, 
    const LocalReference& loc_ref
)  
{    
    if (read.variants.empty() || read.is_na_ref) return;
       
    //int target_end = target.variant_end_pos - 1;
    int target_start = target.lpos;
    int target_end = target.rpos + target.ref_len - 1;
    //int target_end = target.rpos - 1;
    bool is_target = false, is_already_found = false;

    int idx = 0;
    std::vector<int> dist;
    for (auto& variant : read.variants) 
    {
        if (is_already_found) is_target = false;
        else is_target = target.is_equivalent(variant, loc_ref);
        
        if (is_target) 
        {
            read.has_target = true;   
            read.target_pos = variant.pos;
            read.target_ref = variant.ref;
            read.target_alt = variant.alt;
            read.variants_target_idx = idx;
        }
        else 
        {
            // non-targets within target's shiftable region
            if (target_start <= variant.pos && variant.pos <= target_end) 
            {      
                variant.is_overlapping = true;
                
                // magic numbers 
                // should abolish? 
                if (target.indel_len < 4 && read.central_score > 0.15)
                {
                    read.incomplete_shift = true; 
                }
                dist.push_back(0);
            }
            
            variant.set_leftmost_pos(loc_ref);
            variant.rightAlign(loc_ref);

            int event_region_start = variant.lpos;
            int event_region_end = variant.rpos + (variant.ref_len - 1);
            
            //shiftable-non-targets overlapping with target 
            //(target may not be shiftable)
            if (!(event_region_end < target_start
                || target_end < event_region_start))
            {
                variant.is_overlapping = true;
                dist.push_back(0);
            }

            if (event_region_end < target_start
                || target_end < event_region_start)
            {
                //non-overlap -> pass
                dist.push_back(std::abs(target_start - event_region_end));
            }
            else
            {
                //overlap cases
                int lt_d = std::abs(target.pos - event_region_start);
                int rt_d = std::abs(event_region_end - target.pos);

                if (lt_d < rt_d) dist.push_back(lt_d);
                else dist.push_back(rt_d);

            }
        }
        ++idx;
    }
   
    if (!dist.empty()) 
    {    
        read.dist_to_non_target = *std::min_element(dist.begin(), dist.end());
    }
}   


void annot_clip_pattern(Read& read, const Variant& target)
{
    if (read.is_na_ref) return;
    
    if (!read.start_offset && !read.end_offset) 
    {    
        read.clip_ptrn = 'U';    
        return;
    }
    
    // spliced side check
    bool is_left_spliced = (read.aln_start < read.covering_start) ? true : false;
    bool is_right_spliced = (read.covering_end < read.aln_end) ? true : false;
    
    //int lt_dist = std::abs(target.lpos - read.aln_start);
    int lt_dist = std::max(0, target.lpos - read.aln_start);
    //int rt_dist = std::abs(read.aln_end - (target.rpos + target.ref_len - 1));
    int rt_dist = std::max(0, read.aln_end - (target.rpos + target.ref_len - 1));
    bool is_lefty = (lt_dist < rt_dist) ? true : false;

    // both clipped
    if (read.start_offset && read.end_offset)
    {
        // target is middle exon
        if (is_left_spliced && is_right_spliced) 
        {
            read.clip_ptrn = 'U';
            return;
        }
        else if (is_lefty) 
        {
            read.dist_to_clip = lt_dist;
            read.clip_ptrn = 'L';
            return;            
        }
        else 
        {
            read.dist_to_clip = rt_dist;
            read.clip_ptrn = 'R';
            return;
        }
    }
    //left-only clipped
    else if (read.start_offset) 
    {
        if (is_left_spliced) 
        {
            read.clip_ptrn = 'U';
            return;
        }
        else 
        {
            read.dist_to_clip = lt_dist;
            read.clip_ptrn = 'L';
            return;
        }
     }
    //right-only clipped
    else  
    {
        if (is_right_spliced) 
        {
            read.clip_ptrn = 'U';
            return;
        }
        else 
        {
            read.dist_to_clip = rt_dist;
            read.clip_ptrn = 'R';
            return;
        }
    }
}


void annot_non_ref_signature(Read& read)
{
    if (read.is_ref || read.is_na_ref) return;   
    
    std::string sig = "non_ref:", spl_sig = "spl:";
    
    for (const auto& variant : read.variants) 
    {
        sig += (
            std::to_string(variant.pos) 
            + "_" 
            + variant.ref 
            + "_" 
            + variant.alt 
            + ";"
        );
    } 
    
    for (const auto& segment : read.skipped_segments) 
    {
        sig += ("spl=" + std::to_string(segment.first) 
                 + "-" + std::to_string(segment.second) + ";");
            
        spl_sig += ("spl=" + std::to_string(segment.first) 
                    + "-" + std::to_string(segment.second) + ";");
    }
    
    if (read.start_offset) 
    {
        sig += ("lt_clip_end=" + std::to_string(read.aln_start - 1) + ";");
    }
    
    if (read.end_offset) 
    {    
        sig += ("rt_clip_start=" + std::to_string(read.aln_end + 1));
    }

    read.non_ref_signature = sig;
    read.splice_signature = spl_sig;
}
                                                                           

int eval_loc_uniq(
    Read& read,
    const UserParams& user_params,
    const LocalReference& loc_ref
)
{
    
    if (read.clip_ptrn != 'U' 
        || read.variants.empty() 
        || read.nonref_lq_rate >= user_params.lq_rate_thresh) return 0;
    
    //spliced reads may have very large lt/rt_len
    int lt_len = read.variants.front().pos - read.aln_start;
    int rt_len = read.aln_end - read.variants.back().variant_end_pos + 1;  
    int read_len = static_cast<int>(read.seq.size());
    if (lt_len >= read_len || rt_len >= read_len) return 0;    

    // lt fragment
    bool lt_uniq = false;
    int lt_start_pos = -1;
    const size_t lt_most = loc_ref.seq.find(read.seq.substr(0, lt_len));
    if (lt_most != std::string_view::npos)
    {
        const size_t _lt_most = loc_ref.seq.rfind(read.seq.substr(0, lt_len));
        if (lt_most == _lt_most)
        { 
            lt_start_pos = loc_ref.start + static_cast<int>(lt_most);
            if (lt_start_pos ==  read.aln_start) lt_uniq = true;
            //std::cout << "lt pos:" << lt_start_pos << " aln start:" << read.aln_start << std::endl;      
        }
    }
    
    // rt fragment
    bool rt_uniq = false;
    int rt_start_pos = -1;
    const size_t rt_most = loc_ref.seq.rfind(read.seq.substr(read.seq.size() - rt_len));
    if (rt_most != std::string_view::npos)
    {    
        const size_t _rt_most = loc_ref.seq.find(read.seq.substr(read.seq.size() - rt_len));
        if (rt_most == _rt_most)
        {
            rt_start_pos = loc_ref.start + static_cast<int>(rt_most);
            if (rt_start_pos == (read.aln_end - rt_len + 1)) rt_uniq = true;     
            //std::cout << "rt pos:" << rt_start_pos << " actual:" <<  read.aln_end - rt_len + 1<< std::endl;
        }
    } 
    
    if (lt_uniq && rt_uniq) 
    {
        return 1;
    }
    else 
    {     
        return -1;
    }
}


void annot_local_ptrn(
    Read& read, 
    const Variant& target,
    const UserParams& user_params, 
    const LocalReference& loc_ref
)
{
    //trivial cases
    if (read.is_ref || read.covering_ptrn == 'X') 
    {    
        // 'N': no futher check
        read.local_ptrn = 'N';
        return;
    }
     
    if (read.has_target)
    {     
        read.local_ptrn = 'A';
        return;
    }
      
    read.local_uniqueness = eval_loc_uniq(read, user_params, loc_ref);
    // filter likey non-gapped read
    if (read.local_uniqueness > 0 && !has_gaps(read.cigar_str))
    {
        read.local_ptrn = 'N';
        return;
    }
    
    switch (read.covering_ptrn)
    {
        case 'X':
        case 'C':
            read.local_ptrn = 'N';
            return;
        case 'B':
            // 'B': possible target
            read.local_ptrn = 'B';
            return;
        case 'D':
            if (target.rpos + target.ref_len < read.covering_start 
                || read.covering_end < target.lpos)
            {
                read.local_ptrn = 'N';
                return;
            }
            break;
        default:
            break;
    }
    
    //int read_len = read.seq.length();
    //int variant_free_dist = target.ref_len + target.alt_len + 10; 
    //double clip_free_dist = read_len / 2;
    
    if (read.dist_to_non_target <= user_params.local_thresh 
        || read.dist_to_clip <= 2 * user_params.local_thresh) 
    {    
        read.may_be_complex = true;
    }
                      
    //partial match to target for cplx or multiallelic
    if (!read.dist_to_non_target) 
    {    
        read.local_ptrn = 'B';
        return;
    }

    if (read.covering_ptrn == 'D') 
    {
        if (read.covering_end <= target.pos) 
        {
            if (read.clip_ptrn == 'R' 
                //|| read.dist_to_non_target <= user_params.local_thresh) 
                || read.dist_to_non_target <= target.pos - target.lpos) 
            {
                read.local_ptrn = 'B';
                return;
            }
        }
        else if (target.pos <= read.covering_start) 
        {
            if (read.clip_ptrn == 'L' 
                //|| read.dist_to_non_target <= user_params.local_thresh) 
                || read.dist_to_non_target <= target.variant_end_pos - target.pos ) 
            {
                read.local_ptrn = 'B';
                return;
            }
        }
        
        read.local_ptrn = 'N';
        return;
    }
    else 
    {
        if (read.clip_ptrn == 'U') 
        {
            if (user_params.local_thresh < read.dist_to_non_target) 
            {
                read.local_ptrn = 'N';  
            }
            else read.local_ptrn = 'B';
        }
        else 
        {
            if (user_params.local_thresh < read.dist_to_non_target 
                && 2 * user_params.local_thresh < read.dist_to_clip) 
            {    
                read.local_ptrn = 'N';
            }
            else read.local_ptrn = 'B';
        }
    }
}


void count_matched_ends(Read& read)
{
    if (read.is_ref 
        || read.covering_ptrn == 'X' 
        || read.covering_ptrn == 'C') return;
    
    const bool has_variants = !read.variants.empty();
    
    if (read.cigar_vector.front().first != 'S')
    {
        if (has_variants)
        {
            read.lt_end_matches = read.variants.cbegin()->pos - read.aln_start;
        }
        else
        {
            read.lt_end_matches = read.aln_end - read.aln_start;
        }           
    }       
    else
    {
         read.lt_end_matches = 0;
    }
   
    
    if (read.cigar_vector.back().first != 'S')
    {
        if (has_variants)
        {
            read.rt_end_matches = read.aln_end 
                                  - read.variants.crbegin()->variant_end_pos;
        }
        else
        {
            read.rt_end_matches = read.aln_end - read.aln_start;
        }   
    }
    else
    {
        read.rt_end_matches = 0;
    }
}


void eval_local_quality(Read& read, const Variant& target, const UserParams& user_params)
{   
    const int eval_start = target.lpos - 3, eval_end = target.rpos + int(target.ref.size()) + 2;
    
    char op = '\0';
    double tot_cnt = 0.0;
    int op_len = 0, read_idx = 0, lq_cnt = 0, curr_pos = read.read_start;
    for (const auto& c : read.cigar_vector)
    {
        op = c.first;
        op_len = c.second;
        switch (op)
        {
            case 'M':
            case 'X':
            case '=':
            case 'S':
                for (int i = 0; i < op_len; ++i)
                {
                    if (eval_start <= curr_pos && curr_pos <= eval_end) 
                    {
                        ++tot_cnt;
                        if (read.base_quals[read_idx] < user_params.base_q_thresh) ++ lq_cnt;
                    }
                    ++read_idx;
                    ++curr_pos;                              
                }
                break;                   
            case 'D':
            case 'N':
                curr_pos += op_len;
                break;
            case 'I':
                if (eval_start <= curr_pos && curr_pos <= eval_end)
                {
                    for (int i = 0; i < op_len; ++i)
                    {
                         ++tot_cnt;
                         if (read.base_quals[read_idx] < user_params.base_q_thresh) ++ lq_cnt;   
                    }
                    ++read_idx;
                }
                else read_idx += op_len;
                break;
            default:
                break;
        }    
    }
    
    if (tot_cnt)
    {
        read.local_lq_rate = lq_cnt / tot_cnt; 
    }
}


void eval_read_quality(Read& read, const UserParams& user_params)
{
    if (!read.non_ref_quals.empty())
    {
        double non_ref_cnt = 0.0, overall_cnt = 0.0;
        for (const char q : read.non_ref_quals)
        {
            if (q <= user_params.base_q_thresh) non_ref_cnt += 1.0;
        }
        
        for (const char q : read.base_quals)
        {
            if (q <= user_params.base_q_thresh) overall_cnt += 1.0;
        }

        read.nonref_lq_rate = non_ref_cnt / read.non_ref_quals.length();
        read.overall_lq_rate = non_ref_cnt / read.seq.length();
    }
}


void annotate_reads(
    Reads& reads, 
    const Variant& target, 
    const UserParams& user_params, 
    LocalReference& loc_ref,
    bool is_retargeted
)
{
    for (auto& read : reads)
    {
        if (!is_retargeted)
        {
            //annot_ref_seq(read, loc_ref);
            read.setReference(loc_ref);
            read.setVariants(loc_ref);
            //annot_splice_pattern(read);
        }
        
        read.parseCoveringPattern(loc_ref, target);
        annot_covering_ptrn(read, target, loc_ref, is_retargeted);
        annot_target_info(read, target, loc_ref);
        annot_clip_pattern(read, target);
        annot_non_ref_signature(read);  
        annot_local_ptrn(read, target, user_params, loc_ref);
        count_matched_ends(read);
        eval_read_quality(read, user_params);
        //eval_local_quality(read, target, user_params);
    }
}


void classify_reads(
    Reads& reads, 
    Reads& targets, 
    Reads& candidates, 
    Reads& non_targets, 
    const UserParams& user_params
)
{
    size_t max_size = reads.size();
    targets.reserve(max_size);
    candidates.reserve(max_size);
    non_targets.reserve(max_size);    
     
    for (size_t i = 0; i < max_size; ++i)
    {
        
        if (reads[i].name == "WIGTC-HISEQ:6:1307:10350:85982")
        {
            std::cout << reads[i].covering_ptrn << " " << reads[i].local_ptrn << " " << reads[i].dist_to_non_target << " " << reads[i].variants.size() << std::endl;
        }
        
        
        if (reads[i].local_ptrn == 'A') 
        {
            transfer_elem(targets, reads, i);
        }
        else
        {
            if (reads[i].mapq >= user_params.mapq_thresh
               && reads[i].overall_lq_rate < user_params.lq_rate_thresh)
            {
                if (reads[i].local_ptrn == 'B') 
                {    
                    transfer_elem(candidates, reads, i);
                }
                else if (reads[i].local_ptrn == 'N' 
                    && reads[i].is_tight_covering)
                {
                    transfer_elem(non_targets, reads, i);
                }        
             }
         }
    }  
    
    reads.clear();
    targets.shrink_to_fit();
    candidates.shrink_to_fit();
    non_targets.shrink_to_fit();

    for (auto j : candidates)
    {
        if (! j.is_reverse) std::cout << j.name << std::endl;
    }
}*/ 
