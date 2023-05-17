#include <set>
#include <string>
#include <vector>
#include <cstdlib>
#include <climits>
#include <utility>
#include <cstring>
#include <iterator>
#include <algorithm>

//#include "swlib.h"

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
    std::string mid_frag = target.alt;
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
            target = closest_indel;
            
            reclassify_candidates(target, targets, candidates);           
            is_retargetable = true;
        }     
    }
    else return;
}


inline void add_to_input(
                std::vector<SimplifiedRead> & input, 
                const std::vector<int> & indexes, 
                std::vector<Read> & candidates, 
                int n)
{
    size_t idx_size = indexes.size();
    
    if (idx_size)
    {
    
       
    }  
}

/*
//TODO spliced cases
void prep_reads_to_merge(std::vector<SimplifiedRead> & inputs, std::vector<Read> & candidates)
{
    
    int kmer_score_thresh = static_cast<int>(candidates[0].kmer_score / 2.0);
    
    std::vector<int> lt_clipped = {};
    std::vector<int> rt_clipped = {};
    std::vector<int> unclipped = {};
    for (int i = 0; i < candidates.size(); ++i)
    {
        char clp = candidates[i].clip_ptrn;

        switch (clp)
        {
            case 'L':
                lt_clipped.push_back(i);
                break;
            case 'R':
                rt_clipped.push_back(i);
                break;
            default:
                unclipped.push_back(i);
                break;
        }
    }
}
*/

SimplifiedRead merge_to_fragment(std::vector<Read> & candidates, char clip_ptrn)
{
    std::vector<Read> tmp = {};
    for (auto & read : candidates)
    {
        if (read.kmer_score > 0 && read.clip_ptrn == clip_ptrn)
        {  
            read.used4contig = true;
            tmp.push_back(read);
        }
    }
    
    if (tmp.empty())
    {
        SimplifiedRead empty_read;
        return empty_read;    
    }
    else
    {
        if (clip_ptrn == 'L')
        {
            std::sort(tmp.begin(), tmp.end(),
                [](const Read & a, const Read & b){return a.read_start < b.read_start ? true : false;});    
        }
        else
        {
            std::sort(tmp.begin(), tmp.end(),
                [](const Read & a, const Read & b){return a.read_start < b.read_start ? true : false;});
        }

        std::vector<SimplifiedRead> inputs;
        for (auto & read : tmp)
        {
            inputs.emplace_back(read.seq, read.base_qualities, -1);
        }

        return merge_reads(inputs);
     }
}


SimplifiedRead stitch_fragments(std::vector<SimplifiedRead> & fragments)
{
    bool lt_extension = false;
    size_t frag_size = fragments.size();
     
    if (frag_size == 1)
    {
        return fragments[0];
    }
    else if (frag_size == 2)
    {
        return pairwise_stitch(fragments, lt_extension); 
    }
    else if (frag_size == 3)
    {
        std::vector<SimplifiedRead> tmp = {fragments[0], fragments[1]};
        SimplifiedRead _intermediate = pairwise_stitch(tmp, lt_extension);
        return pairwise_stitch({_intermediate, fragments[2]}, lt_extension);
    }
            
    SimplifiedRead empty_read;
    return empty_read;
}







void process_unaligned_target(Variant & target, 
                              FastaReference & fr, 
                              const int base_quality_threshold,
                              const double low_quality_base_rate_threshold, 
                              const size_t kmer_size,
                              const int unspl_loc_ref_start,
                              const std::unordered_map<int, char> & indexed_local_reference,
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
    
    bool is_retargetable = false;                          
    retarget(target, targets, candidates, 2, 10, is_retargetable);
    
    if (is_retargetable)   
    {    
        std::string _c;
        process_aligned_target(target, fr, 
                               base_quality_threshold, low_quality_base_rate_threshold, 
                               3,2,3,1,
                               kmer_size, unspl_loc_ref_start, indexed_local_reference,
                               _c, targets, candidates, non_targets);
        return;
    }

    sift_by_kmer(target, fr, candidates, non_targets);

    //exit if kmer_score = 0 for all candidates
    if (candidates.empty()) return;
    
    std::vector<char> ptrns = {'R', 'U', 'L'};
    std::vector<SimplifiedRead> fragments = {};
    for (const auto & ptrn : ptrns)
    {
        SimplifiedRead _frag = merge_to_fragment(candidates, ptrn);   
        if (!_frag.empty())
        {
            fragments.push_back(_frag);
            std::cout << ptrn << " " << _frag.seq << std::endl; 
        }
    }
    
    std::vector<Read> refinputs = {};
    for (auto & r : candidates)
    {
        if (r.used4contig) refinputs.push_back(r);
    }
    std::vector<std::pair<int, int>> coord = {};
    std::string aa = construct_ref_contig(refinputs, target.chrom, fr, coord);
    
    std::cout << stitch_fragments(fragments).seq << std::endl;
    std::cout << aa << std::endl;

    const uint8_t match_score = 3, mismatch_penalty = 2;
    const uint8_t gap_open_penalty = 3, gap_extention_penalty = 0;
    Filter filter;
    Alignment aln;
    Aligner aligner(
                match_score, mismatch_penalty,
                gap_open_penalty, gap_extention_penalty
            );

    int32_t mask_len = strlen(stitch_fragments(fragments).seq.c_str()) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;
    aligner.Align(stitch_fragments(fragments).seq.c_str(), 
                  aa.c_str(),
                  aa.size(), 
                  filter, &aln, mask_len); 
    std::cout << aln.cigar_string << std::endl; 
    std::cout << aln.ref_begin << " " << aln.query_begin << std::endl;
}

