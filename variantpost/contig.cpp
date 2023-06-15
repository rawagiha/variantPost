#include <random>
#include <utility>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "util.h"
#include "merge.h"
#include "reads.h"
#include "match.h"
#include "contig.h"


void UnalignedContig::furnish(
    const Seq& merged_target_reads,
    const Variant& target
)
{
    seq = merged_target_reads.seq;
    quals = merged_target_reads.base_quals;
    len = seq.size();

    lt_end_idx = merged_target_reads.target_start;
    lt_len = lt_end_idx + 1;
    
    mid_len = (target.alt_len > target.ref_len ? target.alt_len - 1 : 0); 
    
    rt_start_idx = lt_end_idx + 1 + mid_len;
    rt_len = len - rt_start_idx;
}


void unite_coordinates(Coord& coord, const Reads& reads)
{
    if (reads.empty()) return;
    
    std::vector<int> starts, ends;
    for (const auto& read : reads)
    {
        starts.push_back(read.aln_start);
        ends.push_back(read.aln_end);
    }
     
    auto start = std::min_element(starts.begin(), starts.end());
    auto end = std::max_element(ends.begin(), ends.end());
    
    Coord aln_coords = reads[0].aligned_segments;
    
    if (aln_coords.size() == 1) 
    { 
        coord.emplace_back(*start, *end);
    } 
    else 
    {             
        for (auto it = aln_coords.begin(); it != aln_coords.end(); ++it) 
        {
            if (it == aln_coords.begin()) 
            {
                coord.emplace_back(*start, (*it).second);
            }
            else if (it == aln_coords.end() - 1) 
            {
                coord.emplace_back((*it).first, *end);
            }
            else 
            {
                coord.emplace_back((*it).first, (*it).second);
            } 
        }   
    }
}


void find_seed_substitute_reads()
{
//later
}

void find_seed_indel_reads(
    Reads& seeds, 
    Reads& targets,
    const UserParams& user_param, 
    const size_t seed_size
)
{
    Reads clean_targets;
    clean_targets.reserve(targets.size());
    for (const auto& read : targets)
    {
        if (read.nonref_lq_rate < user_param.lq_rate_thresh
            && read.seq.find("N") == std::string::npos)
        {
            clean_targets.push_back(read);
        }
    }
    clean_targets.shrink_to_fit();
    
    //use least dirty one
    if (clean_targets.empty())
    {
        std::sort(
            targets.begin(), targets.end(), [](const Read& a, const Read& b)
            {
                return a.nonref_lq_rate < b.nonref_lq_rate ? true : false;
            }
        );
        clean_targets.push_back(targets[0]);
    }
    
    if (clean_targets.size() == 1)
    {
        std::swap(seeds, clean_targets);
        return;
    }
    
    std::vector<std::string> non_ref_sigs;
    non_ref_sigs.reserve(clean_targets.size());
    for (const auto& read : clean_targets)
    {
        non_ref_sigs.push_back(read.non_ref_signature);
    }
   
    std::string common_sig = find_commonest_str(non_ref_sigs);

    Reads tmp;
    tmp.reserve(clean_targets.size());
    for (const auto& read: clean_targets)
    {
        if (read.non_ref_signature == common_sig)
        {
            tmp.push_back(read);
        }
    }

    sort_by_start(tmp);
    
    if (tmp.size() > seed_size) // seed_size>2
    {     
        std::shuffle(
            tmp.begin() + 1, tmp.end() - 1, std::default_random_engine(123)
        );
        
        seeds.push_back(tmp[0]);
        for (size_t i = 1; i < seed_size - 1; ++i)
        {
            seeds.push_back(tmp[i]);
        }
        seeds.push_back(tmp[tmp.size() - 1]);
        
        sort_by_start(seeds);
    }
    else
    {
        std::swap(seeds, tmp);    
    }
}


int find_target_start(const Read& read)
{
    // with -1 offset
    int i = -1;
    int curr_pos = read.read_start - 1;

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
            default:
                break;
        }
         
        if (curr_pos == read.target_pos) 
        {
            return i;
        }
    }
    return -1;
}


Seq merge_target_seeds(const Reads& seeds) 
{
    std::vector<Seq> inputs;
    for (const auto& seed : seeds)
    {
        inputs.emplace_back(
            seed.seq,
            seed.base_quals,
            find_target_start(seed)
        );
    }

    return merge_reads(inputs);
}


void set_ref_info(
    UnalignedContig& u_contig, 
    const Coord& coord, 
    LocalReference& loc_ref
)
{
    u_contig.coordinates = coord;

    std::string ref_seq;
    for (const auto& c : coord)
    {
        ref_seq += loc_ref.fasta.getSubSequence(
            loc_ref.chrom, c.first - 1, (c.second -  c.first + 1)
        );
    } 
    u_contig.ref_seq = ref_seq;
}

/*
bool is_splice_compatible(const Read& read, const Coord&, mapped_coords)
{
    //prechecked for splicing (e.g. mapped_coords.size() >= 2)
    const size_t coord_size = mapped_coords.size();
    const first_skip_start = mapped_coords[0].second + 1;
    const last_skip_start = mapped_coords[coord_size - 2].second + 1;

    const size_t n_skips = read.skipped_segments.size();
    for (size_t i = 0; i < n_skips; ++i)
    {
    }
}*/

std::string common_splice_ptrn(const Reads& reads)
{
    std::vector<std::string> common_spls;
    common_spls.reserve(reads.size());
    for (const auto& read : reads)
    {
        common_spls.push_back(read.splice_signature);
    }

    return find_commonest_str(common_spls);
}


void make_unaln_contig(
    UnalignedContig& u_contig,
    const Variant& target, 
    Reads& targets, 
    const UserParams& user_params, 
    LocalReference& loc_ref
)
{
    if (target.is_substitute)
    {
    }
    else
    {
        Reads seeds;
        size_t seed_size = 5;
        find_seed_indel_reads(seeds, targets, user_params, seed_size);

        Seq _merged = merge_target_seeds(seeds);
        
        Coord coord;
        unite_coordinates(coord, seeds);
        set_ref_info(u_contig, coord, loc_ref);

        u_contig.furnish(_merged, target);
    }
}


void mock_target_seq(
    std::string& mocked_seq,
    std::string& mocked_ref_seq,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    size_t n = 2;
    int offset = target.is_complex ?  1 : 0;
    
    std::string lt_frag = loc_ref.fasta.getSubSequence(
        loc_ref.chrom,
        target.pos - offset - user_params.kmer_size * n, 
        user_params.kmer_size * n
    );

    std::string mid_frag = "";
    if (target.is_complex)
    {
        mid_frag = target.alt;
    }
    else
    { 
        mid_frag = target.is_ins ? target.alt.substr(1) : "";
    }

    std::string rt_frag = loc_ref.fasta.getSubSequence(
        loc_ref.chrom, 
        target.variant_end_pos - 1, 
        user_params.kmer_size * n
    );
    
    mocked_seq = lt_frag + mid_frag + rt_frag;

    mocked_ref_seq = loc_ref.fasta.getSubSequence(
        loc_ref.chrom,
        target.pos - offset - user_params.kmer_size * n,
        target.ref_len - 1 + offset + user_params.kmer_size * n * 2
    );
}


void prefilter_candidates(
    Reads& candidates,
    Reads& non_targets,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    std::string mocked_seq, mocked_ref_seq;
    mock_target_seq(mocked_seq, mocked_ref_seq, target, user_params, loc_ref);

    Kmers target_kmers = diff_kmers(
        mocked_seq, mocked_ref_seq,  user_params.kmer_size
    );
    
    for (auto& read : candidates)
    {
        read.kmer_score = count_kmer_overlap(read.seq, target_kmers);
    }

    sort_by_kmer(candidates);

    for (Reads::reverse_iterator i = candidates.rbegin(); 
        i != candidates.rend(); ++i)
    {
        if ((*i).kmer_score)
        {
            break; //stop if hits candidate with kmer > 0
        }
        else
        {
            non_targets.insert(
                non_targets.end(), 
                std::make_move_iterator(candidates.rbegin()), 
                std::make_move_iterator(candidates.rbegin() + 1)
            );
            candidates.pop_back();  
        }
    }
}


void prioritize_reads_for_contig_construction(
    Reads& prioritized, 
    const Reads& candidates,
    const UserParams& user_params
)
{
    
    std::string common_spl_ptrn = common_splice_ptrn(candidates);
    
    prioritized.reserve(candidates.size());
    for (const auto& read : candidates)
    {
        if (read.splice_signature == common_spl_ptrn)
        {
            if (read.nonref_lq_rate < user_params.lq_rate_thresh)
            {
                prioritized.push_back(read);
            }
        }
    }
    prioritized.shrink_to_fit();
}


void concat_top_tens(
    Reads& lt_tmp,
    Reads& rt_tmp,
    Reads& u_tmp,
    Reads& top_tens
)
{
    std::vector<Reads> tmps = {lt_tmp, rt_tmp, u_tmp};
    for (auto& tmp : tmps)
    {
        size_t i = 0;
        const size_t size_thresh = 30;
        const size_t tmp_range = std::min(size_thresh, tmp.size());
        while (i < tmp_range)
        {
            transfer_elem(top_tens, tmp, i);
            ++i;
        }
     }

     sort_by_start(top_tens);
}


void flag_contig_member(const Reads& used, Reads& candidates)
{
    for (auto& read : candidates)
    {
        if (std::find(used.begin(), used.end(), read) != used.end())
        {
            read.is_contig_member = true;
        }
    }
}


void suggest_unaln_contig(
    UnalignedContig& u_contig,
    Reads& candidates,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    Reads prioritized;
    
    prioritize_reads_for_contig_construction(
        prioritized, 
        candidates, 
        user_params
    );
    
    //UnalignedContig remains empty
    if (prioritized.empty()) return;
    
    const size_t max_size = 30; // for efficiency
    const size_t search_size = prioritized.size();
    
    Coord coord;
    std::vector<Seq> suggestions;
    suggestions.reserve(max_size);
    
    if (search_size >= max_size)
    {
        Reads lt_tmp, u_tmp, rt_tmp;
        for (size_t i = 0; i < search_size; ++i)
        {
            if (prioritized[i].clip_ptrn == 'L')
            {
                transfer_elem(lt_tmp, prioritized, i);
            }
            else if (prioritized[i].clip_ptrn == 'R')
            {
                transfer_elem(rt_tmp, prioritized, i);
            }
            else
            {
                transfer_elem(u_tmp, prioritized, i);
            }           
        }
        
        sort_by_kmer(lt_tmp);
        sort_by_kmer(rt_tmp);
        sort_by_kmer(u_tmp);
               
        Reads top_tens;       
        concat_top_tens(lt_tmp, rt_tmp, u_tmp, top_tens);

        for (const auto& read : top_tens)
        {
            
            suggestions.emplace_back(read.seq, read.base_quals, -1);
        }
        
        flag_contig_member(top_tens, candidates);  
        unite_coordinates(coord, top_tens);
    }
    else
    {   
        // hack
        Read first = prioritized[0];
        if (prioritized.size() < 4 && first.kmer_score > 3)
        {
            prioritized.push_back(first);
        }
         
        sort_by_start(prioritized);
        
        for (const auto& read : prioritized)
        {
            suggestions.emplace_back(read.seq, read.base_quals, -1);
        }

        flag_contig_member(prioritized, candidates);
        unite_coordinates(coord, prioritized);
    }
    
    suggestions.shrink_to_fit();
    
    Seq merged_suggestions = merge_reads(suggestions);

    u_contig.seq =  merged_suggestions.seq;
    u_contig.quals = merged_suggestions.base_quals;
    
    set_ref_info(u_contig, coord, loc_ref);  
}


