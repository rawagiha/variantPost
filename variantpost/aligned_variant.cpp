#include <set>
#include <random>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <iterator>

#include "util.h"
#include "swlib.h"
#include "aligned_variant.h"
#include "fasta/Fasta.h"
#include "pileup_parser.h"



// select seed reads for assembly
std::vector<std::pair<std::string, std::string>> get_seed_reads
(
    const std::vector<ParsedRead> & targets,
    std::string & commonest_ptrn,
    ParsedRead & sample_target_read, 
    std::vector<ParsedRead> & targets_used,
    size_t n = 5
);


std::string get_ref_contig(std::vector<ParsedRead> & reads, const std::string & chrom, FastaReference & fr);

std::vector<std::pair<int, int>> decompose_read(const ParsedRead & read);


std::vector<std::pair<int, int>> decompose_contig(const std::vector<std::pair<int, int>> & decomposed_read,
                                                  const ParsedRead & read,
                                                  const std::string & contig);

void expected_number_of_repeats(const std::string & repeat_unit,
                                const std::string & contig,
                                const std::vector<std::pair<int, int>> & decomposed_contig,
                                int & exp_rep_n,
                                std::pair<int, int> & boundary_indexes);


void process_aligned_target(
                            //const std::string & unspliced_local_reference,
                            //const int unspliced_local_reference_start,
                            const std::string & chrom,
                            FastaReference & fr,
                            std::string & contig,
                            // may not need
                            int & target_pos,
                            std::string & target_ref,
                            std::string & target_alt,
                            //
                            std::string & repeat_unit,
                            std::vector<ParsedRead> & targets,
                            std::vector<ParsedRead> & candidates,
                            std::vector<ParsedRead> & non_targets)
{
    ParsedRead sample_target_read;
    std::string commonest_non_ref_ptrn = "";
    std::vector<ParsedRead> targets_with_common_ptrn;
    if (targets.size() == 1) {
        sample_target_read = targets[0];
        targets_with_common_ptrn.push_back(sample_target_read);
        contig = sample_target_read.read_seq;
    }
    else {
        std::vector<std::pair<std::string, std::string>> seeds = get_seed_reads(targets,
                                                                                commonest_non_ref_ptrn,
                                                                                sample_target_read,
                                                                                targets_with_common_ptrn);

        if (seeds.size() > 1) {
            std::vector<std::pair<std::string, std::string>> _others = {seeds.begin() + 1, seeds.end()};
            contig = sw::flatten_reads(seeds[0], _others);
        }
        else {
            contig = seeds[0].first;
        }
    }
    
    // contigs
    std::string ref_contig = get_ref_contig(targets_with_common_ptrn, chrom, fr);

    std::vector<std::pair<int, int>> decomposed_read = decompose_read(sample_target_read);
    std::vector<std::pair<int, int>> decomposed_contig = decompose_contig(decomposed_read, 
                                                                          sample_target_read, 
                                                                          contig);
    
    
    // contig kmers
    size_t k = 10;   
    std::set<std::string> ref_kmers, mut_kmers;
    ref_kmers = make_kmers(ref_contig, k);
    mut_kmers = make_kmers(contig, k);
 
    //informative kmers mut only??
    //std::set<std::string> ref_only, mut_only;
    //std::set_difference(ref_kmers.begin(), ref_kmers.end(), mut_kmers.begin(), mut_kmers.end(), std::inserter(ref_only, ref_only.end()));
    std::set<std::string> mut_only;
    std::set_difference(mut_kmers.begin(), mut_kmers.end(), ref_kmers.begin(), ref_kmers.end(), std::inserter(mut_only, mut_only.end()));    
    
    // repeat condition
    Variant observed_target(sample_target_read.target_aligned_pos, 
                            sample_target_read.target_aligned_ref, 
                            sample_target_read.target_aligned_alt);
    repeat_unit = observed_target.minimal_repeat_unit();
    
    int n = 0;
    std::pair<int, int> boundary_indexes;
    expected_number_of_repeats(repeat_unit, contig, decomposed_contig, n, boundary_indexes);
    
    for (const auto & read : candidates) {
        
        int ref_score = 0; 
        int mut_score = 0;
        
        const char* read_seq = read.read_seq.c_str();
        /*
        for (const auto & ref_mer : ref_only) {
            if(strstr(read_seq, ref_mer.c_str())) --ref_score;
        }
        */ 
        for (const auto & mut_mer : mut_only) {
            if(strstr(read_seq, mut_mer.c_str())) ++mut_score;
        }
        
        //std::cout << read.read_name << " " << ref_score << " " << mut_score << " " << read.dirty_base_rate << std::endl;
        
        if (read.dirty_base_rate < 0.05) { 
            if (mut_score) {
                if (sw::is_compatible(contig, ref_contig, read.read_seq, decomposed_contig, repeat_unit, n, boundary_indexes)) {
                    targets.push_back(read);    
                }
                else if (read.covering_ptrn != 'A') {  //can't remember why != 'A'...
                //non_targets.push_back(read);
                }
             }
             else {
                if (read.covering_ptrn == 'A') non_targets.push_back(read);
             }
        }

    }
    
    std::cout << "target: " << targets.size() << std::endl;
    std::cout << "non_target: " << non_targets.size() << " " << std::endl;
    for (auto & read: non_targets) {
    //for (auto & read: targets) {
        std::cout << read.read_name << " " << read.is_reverse << " " << read.dirty_base_rate <<  std::endl;
    }
    //for (auto & d : spliced_segments) {
    //    std::cout << d.first << "-" << d.second << std::endl;
    //}
}

std::vector<std::pair<std::string, std::string>> get_seed_reads
(
    const std::vector<ParsedRead> & targets, 
    std::string & commonest_ptrn,
    ParsedRead & sample_target_read,
    std::vector<ParsedRead> & tmp,   
    size_t n
)
{   
    std::vector<std::string> non_ref_ptrns;
    for (const auto & read : targets) {
        non_ref_ptrns.push_back(read.non_ref_ptrn_str);
    }
    commonest_ptrn = find_commonest_str(non_ref_ptrns);
    
    std::vector<ParsedRead> targets_with_common_ptrn; 
    for (const auto & read : targets) {
        if (read.non_ref_ptrn_str == commonest_ptrn) {
            targets_with_common_ptrn.push_back(read);
        }
    } 
    
    std::vector<std::pair<std::string, std::string>> seed_candidates;
    if (targets_with_common_ptrn.size() > n) {
        std::shuffle(targets_with_common_ptrn.begin(),
                     targets_with_common_ptrn.end(),
                     std::default_random_engine(123));
        for (auto it = targets_with_common_ptrn.begin(); it != targets_with_common_ptrn.begin() + n; ++it) {
            sample_target_read = (*it);
            tmp.push_back(*it);
            seed_candidates.emplace_back((*it).read_seq, (*it).base_qualities);
        }    
    }
    else {
        for (const auto & read : targets_with_common_ptrn) {
            sample_target_read = read;
            seed_candidates.emplace_back(read.read_seq, read.base_qualities);
        }  
        tmp = targets_with_common_ptrn; 
    }
       
    return seed_candidates;
}       

//simple indel
std::vector<std::pair<int, int>> decompose_read(const ParsedRead & read)
{   
    // with -1 offset
    int i = -1;
    int curr_pos = read.read_start - 1;

    int target_pos = read.target_aligned_pos;

    std::pair<int, int> lt_fragment;
    std::pair<int, int> middle_fragment;
    std::pair<int, int> rt_fragment;

    char op;
    int op_len;
    for (const auto & c : read.cigar_vector) {
        op = c.first;
        op_len = c.second;

       
        switch (op) {
            case 'M':
            case 'S':
            case 'X':
            case '=':
                curr_pos += op_len;
                i += op_len;
                break;
            case 'N':
            case 'D':
                curr_pos += op_len;
                break;
            case 'I':
                i += op_len;
                break;
        }
         
        if (curr_pos == target_pos) {
            lt_fragment.first = 0;
            lt_fragment.second = i + 1;

            middle_fragment.first = i + 1;

            int alt_len = read.target_aligned_alt.size();
            int ref_len = read.target_aligned_ref.size();
            if (alt_len > ref_len) {
                middle_fragment.second = (alt_len - ref_len);
                rt_fragment.first = i + 1 + (alt_len - ref_len);
            }
            else {
                middle_fragment.second = 0;
                rt_fragment.first = i + 1; 
            }

            rt_fragment.second = read.read_seq.size() - (i + alt_len - ref_len);
           
            return {lt_fragment, middle_fragment, rt_fragment};
        }
        
    }

    // target not found (ever happens?)
    return {lt_fragment, middle_fragment, rt_fragment};
}

std::vector<std::pair<int, int>> decompose_contig(const std::vector<std::pair<int, int>> & decomposed_read,
                                                  const ParsedRead & read,
                                                  const std::string & contig)
{   
    size_t expected_mid_len = decomposed_read[1].second;
    size_t lt_contig_start = 0;
    size_t lt_contig_len;
    size_t mid_contig_start;
    size_t mid_contig_len = expected_mid_len;
    size_t rt_contig_start;
    size_t rt_contig_len; 

    std::string seq = read.read_seq;
    std::string lt_seq = seq.substr(decomposed_read[0].first, decomposed_read[0].second);
    std::string mid_seq = seq.substr(decomposed_read[1].first, expected_mid_len); 
    std::string rt_seq = seq.substr(decomposed_read[2].first, decomposed_read[2].second);
    
    size_t lt_seq_start = contig.find(lt_seq);
    if (lt_seq_start == std::string::npos) {
        lt_seq = lt_seq.substr(read.start_offset);
        lt_seq_start = contig.find(lt_seq);
    }
    lt_contig_len = lt_seq_start + lt_seq.size();
    mid_contig_start = lt_contig_len;  

    rt_contig_start = contig.find(rt_seq);
    if (rt_contig_start == std::string::npos) {
        rt_seq = rt_seq.substr(0, rt_seq.size() - read.end_offset);
        rt_contig_start = contig.find(rt_seq);
    }
    rt_contig_len = contig.substr(rt_contig_start).size();  
    
    std::vector<std::pair<int, int>> d {
                                            {lt_contig_start, lt_contig_len},
                                            {mid_contig_start, mid_contig_len},
                                            {rt_contig_start, rt_contig_len}
                                       };
    return d;
}


void expected_number_of_repeats(const std::string & repeat_unit,
                                const std::string & contig,
                                const std::vector<std::pair<int, int>> & decomposed_contig,
                                int & exp_n_rep,
                                std::pair<int, int> & boundary_indexes)
{   
    int lt_bound_idx = decomposed_contig[0].second - 1;
    int rt_bound_idx = decomposed_contig[2].first;

    std::string lt_seq = contig.substr(decomposed_contig[0].first, 
                                       decomposed_contig[0].second);
    std::reverse(lt_seq.begin(), lt_seq.end());
    std::string mid_seq = contig.substr(decomposed_contig[1].first, 
                                        decomposed_contig[1].second);
    std::string rt_seq = contig.substr(decomposed_contig[2].first, 
                                       decomposed_contig[2].second);

    int lt_cnt = count_repeats(repeat_unit, lt_seq);
    int mid_cnt = count_repeats(repeat_unit, mid_seq);
    int rt_cnt = count_repeats(repeat_unit, rt_seq);
    
    int rep_len = repeat_unit.size();
    lt_bound_idx = lt_bound_idx - (rep_len * lt_cnt);
    rt_bound_idx = rt_bound_idx + (rep_len * rt_cnt);
    
    exp_n_rep = (lt_cnt + mid_cnt + rt_cnt);
    boundary_indexes.first = lt_bound_idx;
    boundary_indexes.second = rt_bound_idx;       
} 

std::vector<std::pair<int, int>> contig_coordinates_orig(std::vector<ParsedRead> & reads) 
{
    std::vector<int> starts;
    std::vector<int> ends;
    for (const auto & read : reads) {
        starts.push_back(read.aln_start);
        ends.push_back(read.aln_end);
        //std::cout << read.non_ref_ptrn_str << std::endl;
    }
     
    auto contig_start = std::min_element(starts.begin(), starts.end());
    auto contig_end = std::max_element(ends.begin(), ends.end());
    
    std::vector<std::pair<int, int>> un_spls = reads[0].un_spliced_segments;
    
    std::vector<std::pair<int, int>> ret;
    if (un_spls.size() == 1) { 
        ret.emplace_back(*contig_start, *contig_end);
    } 
    else {             
        for (auto it = un_spls.begin(); it != un_spls.end(); ++it) {
            if (it == un_spls.begin()) {
                ret.emplace_back(*contig_start, (*it).second);
            }
            else if (it == un_spls.end() - 1) {
                ret.emplace_back((*it).first, *contig_end);
            }
            else {
                ret.emplace_back((*it).first, (*it).second);
            } 
        }   
    }
    return ret; 
}

std::string get_ref_contig(std::vector<ParsedRead> & reads, const std::string & chrom, FastaReference & fr)
{
    std::vector<std::pair<int, int>> coordinates = contig_coordinates_orig(reads);
    std::string s = "";
    for (auto & coord : coordinates) {
        s += fr.getSubSequence(chrom, coord.first - 1, (coord.second -  coord.first + 1));
        //std::cout << coord.first << " " << coord.second << std::endl;
    }   

    return s;
}    

                                 
