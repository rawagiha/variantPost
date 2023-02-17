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
#include "aligned_target.h"
#include "unaligned_target.h"
#include "read_classifier.h"
#include "fasta/Fasta.h"



/*
struct
{
    bool operator()(const Read & a, 
                    const Read & b, 
                    const std::set<std::string> & target_kmers)
    {
        return (count_kmer_overlap(a.seq, target_kmers) 
                > count_kmer_overlap(b.seq, target_kmers));               
    }

} 
more_similar;
*/

void mock_target_seq(std::string & mocked_seq,
                     std::vector<std::pair<int, int>> & layout,
                     std::string & ref_seq,
                     const Variant & target,
                     FastaReference & fr,
                     const size_t kmer_size)
{
    int lt_frag_end = target.pos;
    int rt_frag_start = lt_frag_end + (target.ref.size() - 1);
    
    size_t k = 3;
    
    std::string lt_frag = fr.getSubSequence(target.chrom, target.pos - kmer_size * k, kmer_size * k);
    layout.emplace_back(0, lt_frag.size());
    std::string mid_frag = target.alt.substr(1, target.alt.size() - 1);
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
    std::sort(candidates.begin(), candidates.end(),
              [](const Read & a, const Read & b){return a.kmer_score > b.kmer_score ? true : false;});
    
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


void retarget(Variant & target,
              std::vector<Read> & targets,
              std::vector<Read> & candidates,
              int pos_ambi_thresh,
              int seq_span_thresh,
              bool & is_retargetable)
{
    const int top_kmer_score = candidates[0].kmer_score;
    std::vector<Read> top_candidates;

    for (auto & candidate : candidates)
    {
        if (candidate.kmer_score == top_kmer_score)
        {
            top_candidates.push_back(candidate);
        }
    }
    
    std::sort(top_candidates.begin(), top_candidates.end(),
                [](const Read & a, const Read & b){return a.centrality > b.centrality ? true : false;}); 

    if (!top_candidates[0].may_be_complex || top_candidates[0].variants.empty()) return;

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
           // Read top = top_candidates[0];
            //top.target_aligned_pos = closest_indel.pos;
            //top.target_aligned_ref = closest_indel.ref;
            //top.target_aligned_alt = closest_indel.alt; 
            //targets.push_back(top);
            
            //closest_indel.chrom = target.chrom;
            target = closest_indel;
            
            reclassify_candidates(target, targets, candidates);           
            is_retargetable = true;
        }     
    }
    else return;
}


void process_unaligned_target(Variant & target, FastaReference & fr, 
                              const int base_quality_threshold,
                              const double low_quality_base_rate_threshold, 
                              std::vector<Read> & targets, 
                              std::vector<Read> & candidates, 
                              std::vector<Read> & non_targets,
                              const size_t kmer_size)
{
    std::string mocked_seq = "";
    std::vector<std::pair<int, int>> mocked_seq_layout;
    std::string ref_seq = "";
    
    kmer_annot_candidates(target, fr, candidates, 
                          mocked_seq, mocked_seq_layout, ref_seq,
                          kmer_size);
                           
    sift_by_kmer(target, fr, candidates, non_targets);

    //exit if kmer_score = 0 for all candidates
    if (!candidates.size()) return; 

    bool is_retargetable = false;
    retarget(target, targets, candidates, 2, 10, is_retargetable);
    
    if (is_retargetable)
    {  
        std::string _c;
        std::cout << target.pos << " " << target.ref << " " << target.alt << " " << target.is_ins << std::endl;
        process_aligned_target(target, fr, base_quality_threshold, low_quality_base_rate_threshold, kmer_size, _c, targets, candidates, non_targets);
    }    
}
