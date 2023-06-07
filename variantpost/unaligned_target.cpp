#include <set>
#include <string>
#include <vector>
#include <cstdlib>
#include <climits>
#include <utility>
#include <cstring>
#include <iterator>
#include <algorithm>

#include "util.h"
#include "localn.h"
#include "aligned_target.h"
#include "unaligned_target.h"
#include "read_classifier.h"
#include "fasta/Fasta.h"



//TODO spliced cases
void mock_target_seq(
        std::string & mocked_seq,
        std::vector<std::pair<int, int>> & layout,
        std::string & ref_seq,
        const Variant & target,
        FastaReference & fr,
        const size_t kmer_size)
{
    int lt_frag_end = target.pos;
    int rt_frag_start = lt_frag_end + (target.ref.size() - 1);
    
    size_t k = 3;
    
    int offset = 0;
    if (target.is_complex) offset = 1; 
    std::string lt_frag = fr.getSubSequence(target.chrom, target.pos - offset - kmer_size * k, kmer_size * k);
    
    layout.emplace_back(0, lt_frag.size());
    
    std::string mid_frag = "";
    if (target.is_complex)
    {
        mid_frag = target.alt;
    }
    else
    {
        if (target.is_ins)
        {
            mid_frag = target.alt.substr(1);
        }
    }

    layout.emplace_back(lt_frag.size(), mid_frag.size());
    std::string rt_frag = fr.getSubSequence(target.chrom,  rt_frag_start, kmer_size * k);
    layout.emplace_back(lt_frag.size() + mid_frag.size() + 1, rt_frag.size()); 
       
    mocked_seq = (lt_frag + mid_frag + rt_frag);
    ref_seq = fr.getSubSequence(target.chrom,  
                                target.pos - kmer_size * k,  
                                target.ref.size() - 1 + kmer_size * k * 2);
}

void kmer_annot_candidates(const Variant & target, 
                           FastaReference & fr,
                           std::vector<Read> & candidates,
                           std::string & mocked_seq,
                           std::vector<std::pair<int, int>> & mocked_seq_layout,
                           std::string & ref_seq,
                           const size_t kmer_size)
                           
{
    mock_target_seq(mocked_seq, mocked_seq_layout, ref_seq, target, fr, kmer_size);

    std::set<std::string> target_kmers = diff_kmers(mocked_seq, ref_seq, kmer_size);
    
    for (auto & read : candidates)
    {
        read.kmer_score = count_kmer_overlap(read.seq, target_kmers);
    }
}


void sift_by_kmer(Variant & target, FastaReference & fr,
                  std::vector<Read> & candidates,
                  std::vector<Read> & non_targets)
{
    //std::sort(candidates.begin(), candidates.end(),
    //          [](const Read & a, const Read & b){return a.kmer_score > b.kmer_score ? true : false;});
    
    sort_by_kmer(candidates);
    
    for (std::vector<Read>::reverse_iterator i = candidates.rbegin();
            i != candidates.rend(); ++i)
    {
        if ((*i).kmer_score)
        {
            break;
        }
        else 
        {
            non_targets.insert(non_targets.end(), 
                               std::make_move_iterator(candidates.rbegin()), 
                               std::make_move_iterator(candidates.rbegin() + 1));
            candidates.pop_back();
        }
    }
}


int dist_to_closest(const Variant & target, 
                    int & min_idx, 
                    std::vector<Variant> & variants,
                    bool indel_only = false)
{
    const int target_pos = target.pos;
    int min_dist = INT_MAX, tmp_dist = INT_MAX;
    for (size_t i = 0; i < variants.size(); ++i)
    {
        if (variants[i].ref.find("N") != std::string::npos
            ||
            variants[i].alt.find("N") != std::string::npos)
        {
            continue;
        }
        
        if (variants[i].is_del)
        {
            if (variants[i].pos < target_pos)
            {
                if (variants[i].pos < target_pos)
                {
                    if (target_pos < variants[i].variant_end_pos) tmp_dist = 0; 
                }
                tmp_dist = target_pos - variants[i].variant_end_pos - 1;
            }
            else
            {
                tmp_dist = variants[i].pos - target_pos;
            }
        }
        else
        {    
            tmp_dist = abs(variants[0].pos - target_pos);
        }
    
        if (tmp_dist < min_dist)
        {
            
            if (variants[i].ref != target.ref || variants[i].alt != target.alt)
            { 
                if (indel_only && variants[i].is_substitute) 
                {
                    continue;
                }
                else
                {
                    min_dist = tmp_dist;
                    min_idx = static_cast<int>(i);
                }
            }
        }
    }
    
    return min_dist;
}


void reclassify_candidates(const Variant & target, 
                           std::vector<Read> & targets,
                           std::vector<Read> & candidates)
{
    const int target_aligned_pos = target.pos;
    const std::string target_aligned_ref = target.ref;
    const std::string target_aligned_alt = target.alt;
    
    bool is_target = false;
    std::vector<Read> tmp = {};
    
    for (size_t i = 0; i < candidates.size(); ++i)
    {
        for (auto & v : candidates[i].variants)
        {
            if (v.pos == target_aligned_pos
                && v.ref == target_aligned_ref
                && v.alt == target_aligned_alt)
            {
                candidates[i].target_aligned_pos = target_aligned_pos;
                candidates[i].target_aligned_ref = target_aligned_ref;
                candidates[i].target_aligned_alt = target_aligned_alt;
                
                transfer_elem(targets, candidates, i);
                is_target = true;   
                continue;
            }
        }
        
        if (!is_target)
        {
            transfer_elem(tmp, candidates, i);
        }

        is_target = false;
    }
    std::swap(candidates, tmp);    
}


//TODO decomposiiton of cplx indels
void retarget(Variant & target,
              std::vector<Read> & targets,
              std::vector<Read> & candidates,
              const double low_quality_base_rate_threshold,
              int pos_ambi_thresh,
              int seq_span_thresh,
              bool & is_retargetable)
{
    const int top_kmer_score = candidates[0].kmer_score;
    std::vector<Read> top_candidates;

    for (const auto & candidate : candidates)
    {
        if (candidate.kmer_score == top_kmer_score
            && candidate.dirty_base_rate < low_quality_base_rate_threshold)
        {
            top_candidates.push_back(candidate);
        }
    }
    
    if (top_candidates.empty()) return;
    
    std::sort(top_candidates.begin(), top_candidates.end(),
                [](const Read & a, const Read & b){return a.centrality > b.centrality ? true : false;}); 

    if (!top_candidates[0].may_be_complex 
        || top_candidates[0].variants.empty()
        || top_candidates[0].centrality < 0.2) return;

    int min_idx = -1;
    int dist2gap = dist_to_closest(target, min_idx, top_candidates[0].variants, true);
    
    if (min_idx != -1)
    {
        Variant closest_indel = top_candidates[0].variants[min_idx];
        closest_indel.chrom = target.chrom;
        
        min_idx = -1;
        int dist2near_var = dist_to_closest(closest_indel, min_idx, top_candidates[0].variants, false);
        if (min_idx == -1) return;
        
        if (dist2gap <= pos_ambi_thresh && dist2near_var <= seq_span_thresh)
        {
            target = closest_indel;
            
            reclassify_candidates(target, targets, candidates);           
            is_retargetable = true;
        }     
    }
    else return;
}



void concat_top_ten(std::vector<Read> & lt_tmp,
                    std::vector<Read> & rt_tmp,
                    std::vector<Read> & u_tmp,
                    std::vector<Read> & top_tens)
{
    std::vector<std::vector<Read>> tmps = {lt_tmp, rt_tmp, u_tmp};
    for (const auto & tmp : tmps)
    {
        int  i = 0;
        int tmp_len = std::min(10, static_cast<int>(tmp.size()));
        while (i < tmp_len)
        {
            top_tens.emplace_back(tmp[i]);
            ++i;
        }
    }
    
    sort_by_start(top_tens);
} 


SimplifiedRead merge_to_contig(std::vector<Read> & candidates)
{
    
    std::vector<SimplifiedRead> inputs  = {};
    if (candidates.size() > 29)
    {
        std::vector<Read> lt_tmp = {}, u_tmp = {}, rt_tmp = {};
        for (const auto & read : candidates)
        {
            if (read.clip_ptrn == 'L') 
            {    
                lt_tmp.emplace_back(read);
            }
            else if (read.clip_ptrn == 'R')
            {
                rt_tmp.emplace_back(read);
            }
            else
            {
                u_tmp.emplace_back(read);
            }
        }
    
        sort_by_kmer(lt_tmp);
        sort_by_kmer(rt_tmp);
        sort_by_kmer(u_tmp);
        
        std::vector<Read> top_tens;
        
        concat_top_ten(lt_tmp, rt_tmp, u_tmp, top_tens);
    
        for (const auto & read : top_tens)
        {
            inputs.emplace_back(read.seq, read.base_qualities, -1);
        }
        
        for (auto & read : candidates)
        {
            if (std::find(top_tens.begin(), top_tens.end(), read) != top_tens.end())
            {
                read.used4contig = true;
            }
        }

    }
    else
    {   
        sort_by_start(candidates);
        for (auto & read : candidates)
        {
            read.used4contig = true;
            inputs.emplace_back(read.seq, read.base_qualities, -1);
        }
    }
    
    return merge_reads(inputs);
    
}


void process_unaligned_target(Variant & target, 
                              FastaReference & fr, 
                              const int base_quality_threshold,
                              const double low_quality_base_rate_threshold, 
                              const int match_score,
                              const int mismatch_penalty,
                              const int gap_open_penalty,
                              const int gap_extention_penalty,
                              const size_t kmer_size,
                              const int unspl_loc_ref_start,
                              const std::unordered_map<int, char> & indexed_local_reference,
                              Contig & contig,
                              Reads & targets, 
                              Reads & candidates, 
                              Reads & non_targets)
{
    
    std::string mocked_seq = "";
    std::vector<std::pair<int, int>> mocked_seq_layout;
    std::string ref_seq = "";
    
    kmer_annot_candidates(target, fr, candidates, 
                          mocked_seq, mocked_seq_layout, ref_seq,
                          kmer_size);
    
    std::cout << "here 1" << std::endl;
    sift_by_kmer(target, fr, candidates, non_targets);
    if (candidates.empty()) return;
    
    std::cout << "here 2 " << candidates.size() << std::endl;
    bool is_retargetable = false;                          
    retarget(target, targets, candidates, low_quality_base_rate_threshold, 2, 10, is_retargetable);
    std::cout << "here 3" << std::endl;

    if (is_retargetable)   
    {    
        std::cout << "enter here" << std::endl;
        std::cout << "retarget: " << target.pos << std::endl;
        std::cout << target.ref << std::endl;
        std::cout << target.alt << std::endl;
        //std::cout << candidates[0].kmer_score << std::endl;
        process_aligned_target(target, fr, 
                               base_quality_threshold, low_quality_base_rate_threshold, 
                               match_score, mismatch_penalty, gap_open_penalty, gap_extention_penalty,
                               kmer_size, unspl_loc_ref_start, indexed_local_reference,
                               contig, targets, candidates, non_targets);
        return;
    }

    std::cout << "here 4" << std::endl;
    SimplifiedRead seq_contig = merge_to_contig(candidates);
    std::vector<Read> ref_inputs = {};
    for (auto & candidate : candidates)
    {
        if (candidate.used4contig) ref_inputs.push_back(candidate);
    }
    std::vector<std::pair<int, int>> coord = {};
    std::string ref_contig = construct_ref_contig(ref_inputs, target.chrom, fr, coord);
    std::cout << "here 5" << std::endl;

    Filter filter;
    Alignment aln;
    std::vector<std::pair<int, int>> grid;
    make_grid(gap_open_penalty, gap_extention_penalty, grid); 
    for (const auto & gap_param : grid)
    {
        sw_aln(match_score, mismatch_penalty,
               gap_param.first, gap_param.second, 
               ref_contig, seq_contig.seq, filter, aln);
        std::cout << ref_contig << " " << seq_contig.seq << " " << target.pos << std::endl;
        std::cout << aln.cigar_string << std::endl;
    }
}

