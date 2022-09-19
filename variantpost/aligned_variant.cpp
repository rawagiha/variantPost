#include <set>
#include <random>
#include <string>
#include <vector>
#include <cstring>
#include <utility>
#include <iterator>
#include <algorithm>

#include "util.h"
#include "swlib.h"
#include "localn.h"
#include "fasta/Fasta.h"
#include "pileup_parser.h"

#include "aligned_variant.h"


std::vector<std::pair<int, int>> decompose_read(const ParsedRead & read);

Reads find_seed_reads(const Reads & targets, const double dirty_thresh, const size_t seed_size)
{   
    Reads clean_targets;
    for (const auto & read : targets) 
    {
        if ((read.dirty_base_rate < dirty_thresh) && (read.read_seq.find("N") == std::string::npos)) 
        {
            clean_targets.push_back(read);
        }
    }
    
    if (clean_targets.empty()) clean_targets = targets;

    std::vector<std::string> non_ref_ptrns;
    for (const auto & read : clean_targets) 
    {
        non_ref_ptrns.push_back(read.non_ref_ptrn_str);
    }
    
    std::string commonest_ptrn = find_commonest_str(non_ref_ptrns);
    
    std::cout << commonest_ptrn << std::endl;
    
    Reads seeds; 
    for (const auto & read : clean_targets) 
    {
        if (read.non_ref_ptrn_str == commonest_ptrn) 
        {
            seeds.push_back(read);
        }
    } 
    
    std::cout << seeds.size() << std::endl;
    
    if (seeds.size() > seed_size) 
    {
        std::cout << "shuff " << std::endl;
        std::shuffle(seeds.begin(), seeds.end(), std::default_random_engine(123));
        
        Reads tmp = {seeds.begin(), seeds.begin() + seed_size};

        
        std::cout << "sort " << std::endl;
        std::sort(tmp.begin(), tmp.end(), 
                  [](const ParsedRead & a, const ParsedRead & b){return a.aln_start < b.aln_start ? true : false;});   
        
        
         
        std::cout << "is this worng" << std::endl;
        return tmp;
    }
    else return seeds;
}


std::vector<std::pair<int, int>> find_coordinates_union(const Reads & reads) 
{
    std::vector<int> starts, ends;
    for (const auto & read : reads)
    {
        starts.push_back(read.aln_start);
        ends.push_back(read.aln_end);
    }
     
    auto contig_start = std::min_element(starts.begin(), starts.end());
    auto contig_end = std::max_element(ends.begin(), ends.end());
    
    std::vector<std::pair<int, int>> un_spls = reads[0].aligned_segments;
    
    std::vector<std::pair<int, int>> coords;
    if (un_spls.size() == 1) 
    { 
        coords.emplace_back(*contig_start, *contig_end);
    } 
    else 
    {             
        for (auto it = un_spls.begin(); it != un_spls.end(); ++it) 
        {
            if (it == un_spls.begin()) 
            {
                coords.emplace_back(*contig_start, (*it).second);
            }
            else if (it == un_spls.end() - 1) 
            {
                coords.emplace_back((*it).first, *contig_end);
            }
            else 
            {
                coords.emplace_back((*it).first, (*it).second);
            } 
        }   
    }
    return coords; 
}


std::string construct_ref_contig(const Reads & reads, const std::string & chrom, FastaReference & fr, std::vector<std::pair<int, int>> & coordinates)
{
    coordinates = find_coordinates_union(reads);
    std::string s = "";
    for (auto & coord : coordinates) 
    {
        s += fr.getSubSequence(chrom, coord.first - 1, (coord.second -  coord.first + 1));
    }   
    return s;
}    


std::vector<std::pair<int, int>> decompose_read(const ParsedRead & read)
{   
    // with -1 offset
    int i = -1;
    int curr_pos = read.read_start - 1;

    int target_pos = read.target_aligned_pos;

    std::pair<int, int> lt_fragment, middle_fragment, rt_fragment;

    char op;
    int op_len;
    for (const auto & c : read.cigar_vector) 
    {
        op = c.first;
        op_len = c.second;
       
        switch (op) 
        {
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
         
        if (curr_pos == target_pos) 
        {
            lt_fragment.first = 0;
            lt_fragment.second = i + 1;

            middle_fragment.first = i + 1;

            int alt_len = read.target_aligned_alt.size();
            int ref_len = read.target_aligned_ref.size();
            if (alt_len > ref_len) 
            {
                middle_fragment.second = (alt_len - ref_len);
                rt_fragment.first = i + 1 + (alt_len - ref_len);
            }
            else 
            {
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



std::vector<std::pair<int, int>> decompose(const MergedRead & mr, const std::string & ref_allele, const std::string & alt_allele)
{
    size_t lt_start = 0;
    size_t lt_len = mr.target_start + 1;
    size_t mid_start = lt_len;
    size_t mid_len = (alt_allele.size() > ref_allele.size() ? alt_allele.size() - 1 : 0);  
    size_t rt_start = mid_start + mid_len;
    size_t rt_len = mr.merged_seq.size() - rt_start;
    
    return {{lt_start, lt_len}, {mid_start, mid_len}, {rt_start, rt_len}};
}


MergedRead preprep(const Reads & targets, const double dirty_thresh, Reads & seeds)
{
    int target_start = 0;
    if (targets.size() == 1)
    {
        seeds = targets;
        target_start = decompose_read(seeds[0])[0].second - 1;   
        
        MergedRead _from_single(seeds[0].read_seq,
                                seeds[0].base_qualities,
                                target_start);
        return _from_single;
    }
    else 
    {
        size_t seed_size = 5;
        seeds = find_seed_reads(targets, dirty_thresh, seed_size);

        std::vector<InputRead> inputs;
        for (const auto & seed : seeds)
        {
            int target_start = decompose_read(seed)[0].second - 1;
            inputs.emplace_back(seed.read_seq, seed.base_qualities, target_start);
              
            // last elem twice
            if (&seed == &seeds.back()) 
            {    
                inputs.emplace_back(seed.read_seq, seed.base_qualities, target_start);
            }
        }

        std::cout << "here" << std::endl;
        MergedRead _from_multiple = merge_reads(inputs);
        std::cout << "here here" << std::endl;
        return _from_multiple;
    }
}


struct Contig
{
private:
    const double dirty_thresh = .1;
public:
    std::string contig_seq;
    std::string contig_qual;
    std::string ref_contig_seq;
    std::string non_ref_ptrn;
     
    double dirty_rate = 0.0;
    bool is_dirty = false;
    size_t seeds_size = 0;
    size_t contig_size = 0;
    double central = 0.0;
    bool lt_clip_seed = false;
    bool rt_clip_seed = false; 
    
    int target_pos;
    std::string target_ref;
    std::string target_alt;
    
    std::vector<std::pair<int, int>> decomposition;
    std::vector<std::pair<int, int>> coordinates;
    
       
    void set_ref_seq(const Reads & seeds, const std::string & chrom, FastaReference & fr)
    {
        ref_contig_seq = construct_ref_contig(seeds, chrom, fr, coordinates);
    }   
    
    void set_dirty_rate(const int qual_thresh)
    {   
        int dirty_bases = 0;
        for (const char & q : contig_qual)
        {
            if ((static_cast<int>(q) - 33) < qual_thresh) ++dirty_bases;   
        }
        
        dirty_rate = static_cast<double>(dirty_bases) / contig_qual.size();
        if (dirty_rate > dirty_thresh) is_dirty = true;
    }
    
    void set_details(const Reads & seeds, const MergedRead & mr)
    {   
        seeds_size = seeds.size();
        ParsedRead sample_seed = seeds[0];
         
        target_pos = sample_seed.target_aligned_pos;
        target_ref = sample_seed.target_aligned_ref;
        target_alt = sample_seed.target_aligned_alt; 
        
        decomposition = decompose(mr, target_ref, target_alt);   
        
        int lt_len = decomposition[0].second;
        int mid_len = decomposition[1].second;
        int rt_len = decomposition[2].second;
        contig_size = lt_len + mid_len + rt_len;

        int closer_to_end = (lt_len < rt_len ? lt_len : rt_len);
        central = static_cast<double>(closer_to_end) / contig_size;
        
        non_ref_ptrn = sample_seed.non_ref_ptrn_str;
        if (non_ref_ptrn.find("lt_clip_end=") != std::string::npos) lt_clip_seed = true;
        if (non_ref_ptrn.find("rt_clip_start=") != std::string::npos) rt_clip_seed = true;
    }

};


Contig set_up_contig(const std::string & chrom, const Reads & targets, const double dirty_rate_thresh, const int qual_thresh, FastaReference & fr)
{
    Reads seeds;
    
    MergedRead mr = preprep(targets, dirty_rate_thresh, seeds);
    
    Contig contig;

    contig.contig_seq = mr.merged_seq;
    contig.contig_qual = mr.merged_qualities;
    contig.set_ref_seq(seeds, chrom, fr);
    contig.set_dirty_rate(qual_thresh);
    contig.set_details(seeds, mr); 
    
    return contig;     
}       

std::set<std::string> diff_kmers(const Contig & contig, const size_t k, const bool is_for_target)
{
    std::set<std::string> mut_kmers = make_kmers(contig.contig_seq, k);
    std::set<std::string> ref_kmers = make_kmers(contig.ref_contig_seq, k);
    
    std::set<std::string> diff;
    if (is_for_target)
    {
        std::set_difference(mut_kmers.begin(), mut_kmers.end(), ref_kmers.begin(), ref_kmers.end(), std::inserter(diff, diff.end()));
    }
    else
    {
        std::set_difference(ref_kmers.begin(), ref_kmers.end(), mut_kmers.begin(), mut_kmers.end(), std::inserter(diff, diff.end()));
    }
    
    return diff;
}

inline void to_left(std::string & unpadded_allele, std::string & lt_seq)
{
    unpadded_allele.pop_back();
    lt_seq.pop_back();
    unpadded_allele = lt_seq.substr(lt_seq.size() - 1, 1) + unpadded_allele;    
}

inline void to_right(std::string & unpadded_allele, std::string & rt_seq)
{
    unpadded_allele.erase(0, 1);
    unpadded_allele += rt_seq.substr(0, 1);
    rt_seq.erase(0, 1);
}

inline bool is_rotatable(const std::string & unpadded_allele)
{
    return (*unpadded_allele.begin() == *(unpadded_allele.end() - 1));
}

void repeat_check(const Variant & trgt, const Contig & contig, const std::string & repeat_unit, int & n_tandem_repeats, bool & is_complete_tandem_repeat, std::pair<int, int> & boundary)
{
    const size_t lt_len = contig.decomposition[0].second;
    const size_t rt_len = contig.decomposition[2].second;
    std::string lt_seq = contig.contig_seq.substr(0, lt_len);
    std::string rt_seq = contig.contig_seq.substr(contig.decomposition[2].first, rt_len);
   
    std::string allele_for_lt = "";
    std::string allele_for_rt = "";
    if (trgt.is_ins)
    {
        allele_for_lt = trgt.alt;
        allele_for_rt = trgt.alt;
    }
    else if (trgt.is_del)
    {
        allele_for_lt = trgt.ref;
        allele_for_rt = trgt.ref;
    }
    
    size_t i = 0;
    while((i < lt_len) & is_rotatable(allele_for_lt)) 
    {
        to_left(allele_for_lt, lt_seq);
        ++i;
    }

    size_t j = 0;
    do
    {
        to_right(allele_for_rt, rt_seq);
        ++j;
    } while((j < rt_len) & is_rotatable(allele_for_rt));

    const size_t lt_bound_end = lt_len - (i + 1);
    const size_t rt_bound_start = contig.decomposition[2].first + j - 1;
    const size_t shiftable_len =  rt_bound_start - lt_bound_end - 1;
    std::string shiftable = contig.contig_seq.substr(lt_bound_end + 1, shiftable_len);
 
    const size_t unit_len = repeat_unit.size();
    if (shiftable_len >= unit_len)
    {
        for (size_t k = 0; k <= (shiftable_len - unit_len); k += unit_len)
        {
            if (shiftable.find(repeat_unit, k) != std::string::npos)
            {
                ++n_tandem_repeats;
                is_complete_tandem_repeat = true;
            }
            else
            {
                is_complete_tandem_repeat = false;
                break;
            }               
        }
        
        size_t remainder = shiftable_len % unit_len;
        if (remainder) is_complete_tandem_repeat = false;
    }

    boundary.first = lt_bound_end;
    boundary.second = rt_bound_start;
}


void process_aligned_target(const std::string & chrom, FastaReference & fr, const int base_quality_threshold, const double low_quality_base_rate_threshold, const int kmer_size, 
                            std::string & _contig, int & target_pos, std::string & target_ref, std::string & target_alt, std::string & _repeat_unit,
                            Reads & targets, Reads & candidates, Reads & non_targets)
{   
    Contig contig = set_up_contig(chrom, targets, low_quality_base_rate_threshold, base_quality_threshold, fr);
    
    const bool is_for_target = true;
    std::set<std::string> informative_kmers = diff_kmers(contig, kmer_size, is_for_target);
    
    Variant obsved_trgt = Variant(contig.target_pos, contig.target_ref, contig.target_alt);
    std::string repeat_unit = obsved_trgt.minimal_repeat_unit();   
    std::string rv_repeat_unit = repeat_unit;
    std::reverse(rv_repeat_unit.begin(), rv_repeat_unit.end());   
     
    std::cout << " before I get here" << std::endl;
    std::cout << contig.contig_seq << " " << contig.decomposition[0].second << std::endl;
    std::cout << contig.ref_contig_seq << std::endl;
    
    int  n_tandem_repeats = 0;
    bool is_complete_tandem_repeat = false;
    std::pair<int, int> repeat_boundary; 
    repeat_check(obsved_trgt, contig, repeat_unit, n_tandem_repeats, is_complete_tandem_repeat, repeat_boundary);
    
    const uint8_t match_score = 3, mismatch_penalty = 0;
    const uint8_t gap_open_penalty = 255, gap_extention_penalty = 255;
    Filter filter;
    Alignment aln;
    Aligner aligner(
                match_score, mismatch_penalty,
                gap_open_penalty, gap_extention_penalty
            );
     
    std::cout << " I get here" << std::endl;
    std::cout << contig.contig_seq << " " << contig.decomposition[0].second << std::endl;
     
    Reads undetermined;
    for (const auto & read : candidates)
    {
        if (is_complete_tandem_repeat & read.has_non_target_in_critical_region) 
        {
            non_targets.push_back(read);
            continue;
        }
       
        int kmer_score = 0;
        const char* read_seq = read.read_seq.c_str();
        for (const auto & kmer : informative_kmers)
        {
            if(strstr(read_seq, kmer.c_str())) ++kmer_score;
        }
        
        const bool is_dirty_query = (read.dirty_base_rate < low_quality_base_rate_threshold);

        if (kmer_score)
        {
            char match_ptrn = match_to_contig(read.read_seq, is_dirty_query,
                                                  contig.contig_seq, contig.ref_contig_seq, contig.decomposition, contig.is_dirty,
                                                  n_tandem_repeats, repeat_unit, rv_repeat_unit, is_complete_tandem_repeat,
                                                  repeat_boundary, filter, aligner, aln);
            switch (match_ptrn)
            {
                case 'T':
                    targets.push_back(read);
                    std::cout << read.read_seq << std::endl;
                    std::cout << read.cigar_string << std::endl;
                    for (auto & v : read.variants)
                    {
                    std::cout << v.pos << " " << v.ref << " " << v.alt;
                    }
                    std::cout<<std::endl;
                    break;
                case 'F':
                    non_targets.push_back(read);
                    break;
                case 'U':
                    undetermined.push_back(read);
                    break;
            }                    
        }
        else
        {
             non_targets.push_back(read);
        }
     
     }   
    
    _contig = contig.contig_seq;
    
    
    Alignment alna;
    
    
    std::cout << "target: " << targets.size() << std::endl;
    std::cout << "non_target: " << non_targets.size() << " " << std::endl;
    std::cout << "undetermined: " << undetermined.size() << std::endl;
     
    for (auto & t : targets)
    {
        //std::cout << t.read_name << " " << t.is_reverse << std::endl;
    }
}




