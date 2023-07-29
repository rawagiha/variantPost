#include <vector>
#include <string>
#include <algorithm>
#include <string_view>

#include "reads.h"
#include "util.h"


inline char offset_qual(int q) 
{
    return static_cast<char>(q + 33);
}


std::string to_qual_str (const std::vector<int>& qual_vec) 
{
    std::vector<char> res (qual_vec.size());
    std::transform (qual_vec.begin(), qual_vec.end(), res.begin(), offset_qual);
    std::string q_str (res.begin(), res.end());
    return q_str;
}


Read::Read(
    std::string_view name, 
    const bool is_reverse,
    std::string_view cigar_str,
    const int aln_start,
    const int aln_end,
    std::string_view seq,
    const std::vector<int>& quals,
    const int mapq,
    const bool is_from_first_bam
) : name(name),
    is_reverse(is_reverse), 
    is_from_first_bam(is_from_first_bam),
    mapq(mapq), 
    aln_start(aln_start), 
    aln_end(aln_end), 
    seq(seq),
    cigar_str(cigar_str) 
{
    base_quals = to_qual_str(quals);
    
    cigar_vector = to_cigar_vector(cigar_str);
    start_offset 
        = cigar_vector.front().first == 'S' ? cigar_vector.front().second : 0;
    end_offset 
        = cigar_vector.back().first == 'S' ? cigar_vector.back().second : 0;

    read_start = aln_start - start_offset;
    read_end = aln_end + end_offset;
}


void sort_by_start(Reads& reads)
{
    std::sort(
        reads.begin(), reads.end(), [](const Read& a, const Read& b)
        {
            return a.read_start <  b.read_start ? true : false;
        }
    );

    std::sort(
        reads.begin(), reads.end(), [](const Read& a, const Read& b)
        {
            return a.aln_start <  b.aln_start ? true : false;
        }
    );
}


void sort_by_kmer(Reads& reads)
{
    std::sort(
        reads.begin(), reads.end(), [](const Read & a, const Read & b)
        {
            return a.kmer_score >  b.kmer_score ? true : false;
        }
    ); 
}


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
            // non-ref exists within target's shiftable region
            if (target_start <= variant.pos && variant.pos <= target_end) 
            {      
                if (target.indel_len < 4 && read.central_score > 0.15)
                {
                    read.incomplete_shift = true; 
                }
                dist.push_back(0);
            }

            variant.set_leftmost_pos(loc_ref);
            variant.set_rightmost_pos(loc_ref);

            int event_region_start = variant.lpos;
            int event_region_end = variant.rpos + (variant.ref_len - 1);

            if (event_region_end <= target_start) 
            {
                dist.push_back(target_start - event_region_end);
            }
            else if (target_end <= event_region_start) 
            {
                dist.push_back(event_region_start - target_end);
            }
            else 
            {
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
                                                                            

bool is_locally_unique(
    Read& read,
    const LocalReference& loc_ref
)
{
    if (read.clip_ptrn != 'U' || read.variants.empty()) return false;
    
    int lt_len = read.variants.front().pos - read.aln_start;
    int rt_len = read.aln_end - read.variants.back().variant_end_pos + 1;  
    int read_len = static_cast<int>(read.seq.size());
    if (lt_len  >= read_len || rt_len >= read_len) return false;    
 
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

    if (lt_uniq && rt_uniq) return true;

    return false;
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
      
    read.is_loc_uniq = is_locally_unique(read, loc_ref);
    if (read.is_loc_uniq && !has_gaps(read.cigar_str))
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
        if (read.covering_end < target.pos) 
        {
            if (read.clip_ptrn == 'R' 
                || read.dist_to_non_target <= user_params.local_thresh) 
            {
                read.local_ptrn = 'B';
                return;
            }
        }
        else if (target.pos < read.covering_start) 
        {
            if (read.clip_ptrn == 'L' 
                || read.dist_to_non_target <= user_params.local_thresh) 
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
            annot_ref_seq(read, loc_ref);
            annot_splice_pattern(read);
        }
        
        annot_covering_ptrn(read, target, loc_ref, is_retargeted);
        

        annot_target_info(read, target, loc_ref);
        annot_clip_pattern(read, target);
        annot_non_ref_signature(read);  
        annot_local_ptrn(read, target, user_params, loc_ref);
        eval_read_quality(read, user_params);
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
} 
