#include <set>
#include <random>
#include <string>
#include <vector>
#include <climits>
#include <cstring>
#include <utility>
#include <iterator>
#include <algorithm>

#include "util.h"
#include "localn.h"
#include "fasta/Fasta.h"
#include "read_classifier.h"

#include "aligned_target.h"


Reads find_seed_reads(const Reads & targets, const double dirty_thresh, const size_t seed_size)
{   
    Reads clean_targets;
    for (const auto & read : targets) 
    {
        if ((read.dirty_base_rate < dirty_thresh) && (read.seq.find("N") == std::string::npos)) 
        {
            clean_targets.push_back(read);
        }
    }
    
    if (clean_targets.empty()) clean_targets = std::move(targets);

    std::vector<std::string> non_ref_ptrns;
    for (const auto & read : clean_targets) 
    {
        non_ref_ptrns.push_back(read.non_ref_ptrn_str);
    }
    
    std::string commonest_ptrn = find_commonest_str(non_ref_ptrns);
    
    Reads seeds; 
    for (const auto & read : clean_targets) 
    {
        if (read.non_ref_ptrn_str == commonest_ptrn) 
        {
            seeds.push_back(read);
        }
    } 
    
    if (seeds.size() > seed_size) 
    {
        std::shuffle(seeds.begin(), seeds.end(), std::default_random_engine(123));
        
        Reads tmp = {seeds.begin(), seeds.begin() + seed_size};
        
        std::sort(tmp.begin(), tmp.end(), 
                  [](const Read & a, const Read & b){return a.read_start < b.read_start ? true : false;});   
         
        return tmp;
    }
    else return seeds;
}


// input: seed reads 
// input as starts/end vectors -> easier to reuse
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


std::vector<std::pair<int, int>> split_target_read(const Read & read)
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

            rt_fragment.second = read.seq.size() - (i + alt_len - ref_len);
           
            return {lt_fragment, middle_fragment, rt_fragment};
        }
        
    }

    // target not found (ever happens?)
    return {lt_fragment, middle_fragment, rt_fragment};
}


std::vector<std::pair<int, int>> split_merged_target_read(const SimplifiedRead & mr, const std::string & ref_allele, const std::string & alt_allele)
{
    size_t lt_start = 0;
    size_t lt_len = mr.target_start + 1;
    size_t mid_start = lt_len;
    size_t mid_len = (alt_allele.size() > ref_allele.size() ? alt_allele.size() - 1 : 0);  
    size_t rt_start = mid_start + mid_len;
    size_t rt_len = mr.seq.size() - rt_start;
    
    return {{lt_start, lt_len}, {mid_start, mid_len}, {rt_start, rt_len}};
}


SimplifiedRead merge_target_reads(const Reads & targets, const double dirty_thresh, Reads & seeds)
{
    int target_start = 0;
    if (targets.size() == 1)
    {
        seeds = targets;
        target_start = split_target_read(seeds[0])[0].second - 1;   
        
        SimplifiedRead _from_single(seeds[0].seq,
                                    seeds[0].base_qualities,
                                    target_start);
        return _from_single;
    }
    else 
    {
        size_t seed_size = 5;
        seeds = find_seed_reads(targets, dirty_thresh, seed_size);

        std::vector<SimplifiedRead> inputs;
        for (const auto & seed : seeds)
        {
            int target_start = split_target_read(seed)[0].second - 1;
            inputs.emplace_back(seed.seq, seed.base_qualities, target_start);
              
            // last elem twice
            if (&seed == &seeds.back()) 
            {    
                inputs.emplace_back(seed.seq, seed.base_qualities, target_start);
            }
        }

        SimplifiedRead _from_multiple = merge_reads(inputs);
        return _from_multiple;
    }
}


struct RawContig
{
private:
    const double dirty_thresh = .1;
public:
    std::string seq;
    std::string base_qualities;
    std::string ref_seq;
    std::string non_ref_ptrn;
     
    double dirty_rate = 0.0;
    bool is_dirty = false;
    size_t seeds_size = 0;
    size_t contig_size = 0;
    double central = 0.0;
    bool closer_to_lt_end = false;
    bool lt_clip_seed = false;
    bool rt_clip_seed = false;
    bool is_spliced = false;
    
    int target_pos;
    std::string target_ref;
    std::string target_alt;
    
    std::vector<std::pair<int, int>> decomposition;
    std::vector<std::pair<int, int>> coordinates;
          
    void set_ref_seq(const Reads & seeds, const std::string & chrom, FastaReference & fr)
    {
        ref_seq = construct_ref_contig(seeds, chrom, fr, coordinates);
        is_spliced = (coordinates.size() > 1);
    }   
    
    void set_dirty_rate(const int qual_thresh)
    {   
        int dirty_bases = 0;
        for (const char & q : base_qualities)
        {
            if ((static_cast<int>(q) - 33) < qual_thresh) ++dirty_bases;   
        }
        
        dirty_rate = static_cast<double>(dirty_bases) / base_qualities.size();
        if (dirty_rate > dirty_thresh) is_dirty = true;
    }
    
    void set_details(const Reads & seeds, const SimplifiedRead & mr)
    {   
        seeds_size = seeds.size();
        Read sample_seed = seeds[0];
         
        target_pos = sample_seed.target_aligned_pos;
        target_ref = sample_seed.target_aligned_ref;
        target_alt = sample_seed.target_aligned_alt; 
        
        decomposition = split_merged_target_read(mr, target_ref, target_alt);   
        
        int lt_len = decomposition[0].second;
        int mid_len = decomposition[1].second;
        int rt_len = decomposition[2].second;
        contig_size = lt_len + mid_len + rt_len;

        int dist_to_end = 0;
        if (lt_len < rt_len)
        {
            dist_to_end = lt_len;
            closer_to_lt_end = true;
        }
        else
        {
            dist_to_end = rt_len;
        }
        
        central = static_cast<double>(dist_to_end) / contig_size;
        
        non_ref_ptrn = sample_seed.non_ref_ptrn_str;
        if (non_ref_ptrn.find("lt_clip_end=") != std::string::npos) lt_clip_seed = true;
        if (non_ref_ptrn.find("rt_clip_start=") != std::string::npos) rt_clip_seed = true;
    }

    bool is_high_quality()
    {
     return (!lt_clip_seed && !rt_clip_seed && central > .25);   
    }    
};


RawContig set_up_contig(const std::string & chrom, const Reads & targets, const double dirty_rate_thresh, const int qual_thresh, FastaReference & fr)
{
    Reads seeds;
    
    SimplifiedRead mr = merge_target_reads(targets, dirty_rate_thresh, seeds);
    
    RawContig contig;

    contig.seq = mr.seq;
    contig.base_qualities = mr.base_qualities;
    contig.set_ref_seq(seeds, chrom, fr);
    contig.set_dirty_rate(qual_thresh);
    contig.set_details(seeds, mr); 

    return contig;     
}       

/*
//specialized for target
std::set<std::string> diff_kmers(const Contig & contig, const size_t k)
{
    std::set<std::string> mut_kmers = make_kmers(contig.seq, k);
    std::set<std::string> ref_kmers = make_kmers(contig.ref_seq, k);
    
    std::set<std::string> diff;
    std::set_difference(mut_kmers.begin(), mut_kmers.end(), 
                        ref_kmers.begin(), ref_kmers.end(), 
                        std::inserter(diff, diff.end()));
    
    return diff;
}
*/

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

void repeat_check(const Variant & trgt, 
                  const RawContig & contig, 
                  const std::string & repeat_unit, 
                  int & n_tandem_repeats, 
                  bool & is_complete_tandem_repeat, 
                  std::pair<int, int> & boundary)
{
    const size_t lt_len = contig.decomposition[0].second;
    const size_t rt_len = contig.decomposition[2].second;
    std::string lt_seq = contig.seq.substr(0, lt_len);
    std::string rt_seq = contig.seq.substr(contig.decomposition[2].first, rt_len);
    
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
    while((i < lt_len) && is_rotatable(allele_for_lt)) 
    {
        to_left(allele_for_lt, lt_seq);
        ++i;
    }

    size_t j = 0;
    do
    {
        to_right(allele_for_rt, rt_seq);
        ++j;
    } while((j < rt_len) && is_rotatable(allele_for_rt));
    --j; //undo once for rt aln
    
    const size_t lt_bound_end = lt_len - (i + 1);
    const size_t rt_bound_start = contig.decomposition[2].first + j - 1;
    const size_t shiftable_len =  rt_bound_start - lt_bound_end - 1;
    std::string shiftable = contig.seq.substr(lt_bound_end + 1, shiftable_len);
 
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


void classify_candidates(Reads & candidates, const RawContig & contig, 
                         const Filter & filter, const Aligner & aligner, Alignment & alignment, 
                         const std::set<std::string> & target_kmers,
                         const std::set<std::string> & core_kmers,
                         const bool is_complete_tandem_repeat, const int n_tandem_repeats,
                         const std::pair<int, int> & repeat_boundary,
                         const std::string & repeat_unit, const std::string & rv_repeat_unit,
                         const double low_quality_base_rate_threshold,
                         Reads & lt_extenders, Reads & extra_targets, Reads & rt_extenders,
                         Reads & undetermined, Reads & non_targets)
{
    for (size_t i = 0; i < candidates.size(); ++i)
    {
        if (is_complete_tandem_repeat && candidates[i].has_non_target_in_critical_region) 
        {
            //non_targets.push_back(read);
            //transfer_elem(non_targets, candidates, i);
            //continue;
        }
        
        int kmer_score = count_kmer_overlap(candidates[i].seq, target_kmers);

        //covering reads must have core kmers
        if (kmer_score && candidates[i].covering_ptrn == 'A')
        {
            kmer_score = count_kmer_overlap(candidates[i].seq, core_kmers);
        } 
        
        const bool is_dirty_query = (candidates[i].dirty_base_rate < low_quality_base_rate_threshold);

        if (kmer_score)
        {
            char match_ptrn = match_to_contig(candidates[i].seq, is_dirty_query,
                                              contig.seq, contig.ref_seq, contig.decomposition, contig.is_dirty,
                                              n_tandem_repeats, repeat_unit, rv_repeat_unit, is_complete_tandem_repeat,
                                              repeat_boundary, filter, aligner, alignment); 
            
            switch (match_ptrn)
            {
                case 'L':
                    //lt_extenders.push_back(read);
                    transfer_elem(lt_extenders, candidates, i);
                    //std::cout << read.read_name << " " << read.is_reverse << " " << read.cigar_string << std::endl;            
                    break;
                case 'M':
                    transfer_elem(extra_targets, candidates, i);
                    //extra_targets.push_back(read);
                    break;
                case 'R':
                    transfer_elem(rt_extenders, candidates, i);
                    //rt_extenders.push_back(read);
                    break;
                case 'F':
                    //non_targets.push_back(read);
                    transfer_elem(non_targets, candidates, i);
                    break;
                case 'U':
                    //undetermined.push_back(read);
                    transfer_elem(undetermined, candidates, i);
                    break;
            }                    
        }
        else
        {
             //non_targets.push_back(read);
             transfer_elem(non_targets, candidates, i);
        }
     }
     candidates.clear();   
}    



bool is_covered(const std::pair<int, int> & segment, 
                const std::vector<std::pair<int, int>> & coordinates,
                const size_t coord_size)
{   
    for (size_t i = 0; i < coord_size; ++i)
    {
        if (i == 0)
        {
            if (segment.second <= coordinates[i].second) return true;                
        }
        else if (i == (coord_size - 1))
        {
            if (coordinates[i].first <= segment.first) return true;
        }
        else
        {
            if (coordinates[i].first == segment.first 
                && coordinates[i].second == segment.second) return true;
        }
    }
    
    return false;   
}

bool is_compatible_spl_ptrn(const Read & read, const std::vector<std::pair<int, int>> & coordinates)
{
    const size_t coord_size = coordinates.size();
    
    /*
    if (coord_size == 1)
    {
        //unspliced read
    }
    else
    */

    if (coord_size > 1)
    {   
        const int first_skip_start = coordinates[0].second + 1; 
        const int last_skip_start = coordinates[(coord_size - 2)].second + 1;

        // check for 5'and 3' skips
        const size_t n_skips = read.skipped_segments.size();
        for (size_t i = 0; i < n_skips; ++i)
        {
            if (i == 0)
            {
                if (read.skipped_segments[0].first < first_skip_start) return false; 
            }
            
            if (i == (n_skips - 1))
            {
                if (last_skip_start < read.skipped_segments[(n_skips - 1)].first) return false; 
            }
        }
       
        // check if aligned bases are included in coordinates        
        //const size_t n_alns = read.aligned_segments.size(); 
        bool _covered = false;
        for (const auto & segment : read.aligned_segments)
        {               
            _covered = is_covered(segment, coordinates, coord_size);
            if (!_covered) return false;
        } 
        
        /*
        std::vector<int> shared_intron_index; 
        const size_t n_segments = read.skipped_segments.size();

        for (size_t i = 0; i < n_segments; ++i)
        {
            if (std::find(coordinates.begin(), 
                          coordinates.end(), 
                          read.skipped_segments[i]) != coordinates.end())
            { 
                shared_intron_index.push_back(i); 
            }
        }
        */
    }

    return true;

} 

//splice compatibility check
//bool is_splice_compatible (in util?)
//do we need commonest prtn check (as in seed)??
void extend_contig_seq(const RawContig & contig, 
                       const Reads & lt_extenders, 
                       const Reads & rt_extenders, 
                       SimplifiedRead & extended_contig,
                       std::vector<std::pair<int, int>> & extended_coordinates)
{

    //const size_t lt_len = contig.decomposition[0].second;
    //const size_t rt_len = contig.decomposition[2].second;
    
    SimplifiedRead merged_ext;

    if (contig.closer_to_lt_end)
    {
        const size_t n_lt_ext = lt_extenders.size();
        if (n_lt_ext)
        {
            const size_t lt_size = (n_lt_ext > 5 ? 5 : n_lt_ext);
            
            Reads lt_tmp = {lt_extenders.begin(), lt_extenders.begin() + lt_size};

            std::sort(lt_tmp.begin(), lt_tmp.end(),
                  [](const Read & a, const Read & b)
                  {return a.aln_start < b.aln_start ? true : false;});

            std::vector<int> aln_starts;
            std::vector<SimplifiedRead> lt_inputs;
            for (auto & _ext : lt_tmp)
            {
                if (is_compatible_spl_ptrn(_ext, contig.coordinates))
                {
                    lt_inputs.emplace_back(_ext.seq, _ext.base_qualities, -1);
                    aln_starts.push_back(_ext.aln_start);
                }
            }

            if (!lt_inputs.empty())
            {
                merged_ext = merge_reads(lt_inputs);

                std::vector<SimplifiedRead> lt_scaffolds;
                lt_scaffolds.emplace_back(contig.seq, contig.base_qualities, -1);
                lt_scaffolds.emplace_back(merged_ext.seq, merged_ext.base_qualities, -1);
            
                extended_contig = pairwise_stitch(lt_scaffolds, true);

                //update coord
                extended_coordinates[0].first = *(std::min_element(aln_starts.begin(), aln_starts.end()));
            }
        }
    }
    else
    {
        const size_t n_rt_ext = rt_extenders.size();
        
        if (n_rt_ext)
        {
            const size_t rt_size = (n_rt_ext > 5 ? 5 : n_rt_ext); 
             
            Reads rt_tmp = {rt_extenders.end() - rt_size, rt_extenders.end()};

            std::sort(rt_tmp.begin(), rt_tmp.end(), 
                  [](const Read & a, const Read & b)
                  {return a.aln_start < b.aln_start ? false : true;});

            std::vector<int> aln_ends;
            std::vector<SimplifiedRead> rt_inputs;
            for (auto & _ext : rt_tmp)
            {
                if (is_compatible_spl_ptrn(_ext, contig.coordinates))
                {
                    rt_inputs.emplace_back(_ext.seq, _ext.base_qualities, -1);
                    aln_ends.push_back(_ext.aln_end);
                }
            }
            
            
            if (!rt_inputs.empty())
            {
                merged_ext = merge_reads(rt_inputs);
            
                std::vector<SimplifiedRead> rt_scaffolds;
                rt_scaffolds.emplace_back(contig.seq, contig.base_qualities, -1);
                rt_scaffolds.emplace_back(merged_ext.seq, merged_ext.base_qualities, -1);

                extended_contig = pairwise_stitch(rt_scaffolds, false);
                

                //update coord
                extended_coordinates[(extended_coordinates.size() - 1)].second = *(std::max_element(aln_ends.begin(), aln_ends.end()));
            }
        } 
    }
    
    //std::cout << extended_contig.read_seq << " " << extended_contig.base_qualities << std::endl; 
}


void align_to_contig(const std::string & chrom, FastaReference & fr, 
                     SimplifiedRead & extended_contig, 
                     std::vector<std::pair<int, int>> & extended_coordinates)
{
    std::string extended_ref = ""; 
    for (const auto & coord : extended_coordinates)
    {
        extended_ref += fr.getSubSequence(chrom, coord.first - 1, (coord.second -  coord.first + 1));
    }

    std::cout << extended_ref.size() << std::endl;
    std::cout << extended_contig.seq.size() << std::endl;
    std::cout << extended_contig.base_qualities.size() << std::endl;
    
     
    //try gap_ext = 1 first
    //if not as target, try gap_ext=0
    //use whichever larget margin
    //rosenfeld case
    const uint8_t match_score = 3, mismatch_penalty = 2;
    const uint8_t gap_open_penalty = 3, gap_extention_penalty = 1;
    Filter filter;
    Alignment aln;
    Aligner aligner(
                match_score, mismatch_penalty,
                gap_open_penalty, gap_extention_penalty
            );

    int32_t mask_len = strlen(extended_contig.seq.c_str()) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;
    aligner.Align(extended_contig.seq.c_str(), 
                  extended_ref.c_str(), 
                  extended_ref.size(), 
                  filter, &aln, mask_len);

    std::cout << aln.cigar_string << std::endl;
    std::cout << aln.ref_begin << " " << aln.query_begin << std::endl;

    // edit cigar 
    std::vector<int> genomic_pos = expand_coordinates(extended_coordinates);
    std::vector<std::pair<char, int>> cigar_vec = to_cigar_vector(aln.cigar_string);    
    splice_cigar(cigar_vec, aln.ref_begin, genomic_pos, extended_coordinates);
    move_up_insertion(cigar_vec); 
   
    int aln_start = genomic_pos[aln.ref_begin];
    int aln_end = genomic_pos[aln.ref_end];  
    std::cout << aln_start << " " << aln_end << std::endl;
    std::string non;
    std::vector<Variant> aa = find_mapped_variants(aln_start, aln_end, extended_ref.substr(aln.ref_begin), extended_contig.seq.substr(aln.query_begin), extended_contig.base_qualities,  cigar_vec, non);
    for (auto & v : aa)
    {
        std::cout << v.pos << " " << v.ref << " " << v.alt << std::endl;
    }
}

void make_grid(const int gap_open_penalty, const int gap_extension_penalty, std::vector<std::pair<int, int>> & grid)
{
    std::pair<int, int> _default = {gap_open_penalty, gap_extension_penalty};

    //std::vector<std::pair<int, int>> grid = {};
    uint8_t max_gap_open = 5, min_gap_open = 3;
    if (gap_open_penalty > 5) max_gap_open = gap_open_penalty;
    
    for (uint8_t i = min_gap_open; i <= max_gap_open; ++i)
    {      
        grid.emplace_back(i, 1);
        grid.emplace_back(i, 0);
    }
    
    if (std::find(grid.begin(), grid.end(), _default) !=  grid.end())
    {
        // pass
    }
    else
    {
        grid.push_back(_default);
        //return grid;
    }
}


void realn_extended_contig
    (const std::string & chrom, FastaReference & fr,
     SimplifiedRead & extended_contig,
     std::vector<std::pair<int, int>> & extended_coordinates,
     const int match_score,
     const int mismatch_penalty,
     const int gap_open_penalty,
     const int gap_extention_penalty,
     const Variant & target,
     const int unspl_loc_ref_start,
     const std::unordered_map<int, char> & indexed_local_reference,
     std::vector<RealignedGenomicSegment> & realns)
{
    std::string extended_ref = ""; 
    for (const auto & coord : extended_coordinates)
    {
        extended_ref += fr.getSubSequence(chrom, coord.first - 1, (coord.second -  coord.first + 1));
    }
    std::vector<int> genomic_pos = expand_coordinates(extended_coordinates);
    
    
    Filter filter;
    Alignment aln;
    std::vector<Variant> _variants = {};
    std::vector<std::pair<char, int>> _cigar_vec = {};
    std::string variant_quals = "";

    std::vector<std::pair<int, int>> grid;
    make_grid(gap_open_penalty, gap_extention_penalty, grid);  
    for (auto & gap_param : grid)
    {
        sw_aln(match_score, mismatch_penalty, 
               gap_param.first, gap_param.second, 
               extended_ref, extended_contig.seq, filter, aln);
        
        _cigar_vec = to_cigar_vector(aln.cigar_string);
        splice_cigar(_cigar_vec, aln.ref_begin, genomic_pos, extended_coordinates);
        move_up_insertion(_cigar_vec);
        
        int begin_with_offset = aln.query_begin;
        if (_cigar_vec[0].first == 'S') begin_with_offset -= _cigar_vec[0].second;

        _variants = find_mapped_variants(
                            genomic_pos[aln.ref_begin], 
                            genomic_pos[aln.ref_end],
                            extended_ref.substr(aln.ref_begin), 
                            extended_contig.seq.substr(begin_with_offset), 
                            extended_contig.base_qualities.substr(begin_with_offset), 
                            _cigar_vec, variant_quals
                    );
        
        int target_pos = target.pos;
        bool has_target = false;
        for (const auto & v : _variants)
        {
            if (v.is_equivalent(target, unspl_loc_ref_start, indexed_local_reference))
            {
                target_pos = v.pos;
                has_target = true;
            }
        }
        
        RealignedGenomicSegment realn(genomic_pos[aln.ref_begin], 
                                      genomic_pos[aln.ref_end], 
                                      target_pos, 
                                      has_target,
                                      extended_ref.substr(aln.ref_begin), 
                                      extended_contig.seq.substr(begin_with_offset),
                                      extended_contig.base_qualities.substr(begin_with_offset),
                                      _cigar_vec, _variants);
        
        realns.push_back(realn);
    }
}


void find_core_kmers(const RawContig & contig, const int kmer_size, 
                     const std::set<std::string> & differential_kmers, std::set<std::string> & core_kmers)
{
     std::set<std::string> left_diff_kmers = diff_kmers(contig.seq.substr(contig.decomposition[0].first, contig.decomposition[0].second),
                                                        contig.ref_seq,
                                                        kmer_size); 
     
     std::set<std::string> right_diff_kmers = diff_kmers(contig.seq.substr(contig.decomposition[2].first, contig.decomposition[2].second),
                                                         contig.ref_seq,
                                                         kmer_size); 
     
     for (auto & kmer : differential_kmers)
     {
        if (
            (left_diff_kmers.find(kmer) == left_diff_kmers.end()) 
            &&
            (right_diff_kmers.find(kmer) == right_diff_kmers.end())
        ) core_kmers.insert(kmer);
     }

     /* peripheral kmers
     std::set_difference(differential_kmers.begin(), differential_kmers.end(),
                         core_kmers.begin(), core_kmers.end(),
                         std::inserter(peripheral_kmers, peripheral_kmers.end()));      
     */
}


//void process_aligned_target(const std::string & chrom, FastaReference & fr, const int base_quality_threshold, const double low_quality_base_rate_threshold, const int kmer_size, 
//                            std::string & _contig, int & target_pos, std::string & target_ref, std::string & target_alt, std::string & _repeat_unit,
//                            Reads & targets, Reads & candidates, Reads & non_targets)
void process_aligned_target(Variant & target,
                            FastaReference & fr, 
                            const int base_quality_threshold,
                            const double low_quality_base_rate_threshold, 
                            const int match_score,
                            const int mismatch_penalty,
                            const int gap_open_penalty,
                            const int gap_extention_penalty,
                            const int kmer_size,
                            const int unspl_loc_ref_start,
                            const std::unordered_map<int, char> & indexed_local_reference,
                            Contig & contig,
                            Reads & targets, Reads & candidates, Reads & non_targets)

{   
    RawContig raw_contig = set_up_contig(target.chrom, targets, low_quality_base_rate_threshold, base_quality_threshold, fr);
    
    std::set<std::string> informative_kmers = diff_kmers(raw_contig.seq, raw_contig.ref_seq, kmer_size);
    std::set<std::string> core_kmers = {};
    find_core_kmers(raw_contig, kmer_size, informative_kmers, core_kmers);

    // may happen for pure mapping artifacts 
    if (core_kmers.empty())
    {
        transfer_vector(non_targets, targets);
        transfer_vector(non_targets, candidates);
        return;
    }

    
    Variant observed_target = Variant(raw_contig.target_pos, 
                                      raw_contig.target_ref, 
                                      raw_contig.target_alt);
    std::string repeat_unit = observed_target.minimal_repeat_unit();   
    std::string rv_repeat_unit = repeat_unit;
    std::reverse(rv_repeat_unit.begin(), rv_repeat_unit.end());   
    
    int  n_tandem_repeats = 0;
    bool is_complete_tandem_repeat = false;
    std::pair<int, int> repeat_boundary = {}; 
    repeat_check(observed_target, raw_contig, repeat_unit, n_tandem_repeats, is_complete_tandem_repeat, repeat_boundary);
    
    //gap-less alignment
    const uint8_t _gap_open_penalty = 255, _gap_extention_penalty = 255;
    Filter filter;
    Alignment aln;
    Aligner aligner(
                match_score, mismatch_penalty,
                _gap_open_penalty, _gap_extention_penalty
            );
     
    
    Reads lt_extenders, extra_targets, rt_extenders, undetermined;
    
    classify_candidates(candidates, raw_contig, filter, aligner, aln, informative_kmers, core_kmers,
                        is_complete_tandem_repeat, n_tandem_repeats, repeat_boundary,
                        repeat_unit, rv_repeat_unit, low_quality_base_rate_threshold,
                        lt_extenders, extra_targets, rt_extenders, undetermined, non_targets);
    
    
    SimplifiedRead extended(raw_contig.seq, raw_contig.base_qualities);
    std::vector<std::pair<int, int>> ext_coord = raw_contig.coordinates;
    extend_contig_seq(raw_contig, lt_extenders, rt_extenders, extended, ext_coord);
    
    std::vector<RealignedGenomicSegment> realns = {};
    realn_extended_contig(target.chrom, fr, extended, ext_coord, match_score, mismatch_penalty, gap_open_penalty, gap_extention_penalty, target, unspl_loc_ref_start, indexed_local_reference, realns);
    
    make_contig(realns, contig);
    
    transfer_vector(targets, lt_extenders);
    transfer_vector(targets, extra_targets);
    transfer_vector(targets, rt_extenders); 
}




