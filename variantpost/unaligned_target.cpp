#include <set>
#include <string>
#include <vector> 
#include <utility>
#include <cstring>
#include <iterator>
#include <algorithm>

#include "util.h"
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




void process_unaligned_target(Variant & target, FastaReference & fr, 
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
}

