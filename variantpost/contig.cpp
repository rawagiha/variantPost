#include <random>
//#include <iostream>

#include "eval.h"
#include "util.h"
#include "merge.h"
#include "reads.h"
#include "contig.h"


Contig::Contig() {}


Contig::Contig(const Variant& target)
{
    positions = {target.pos};
    ref_bases = {target.ref};
    alt_bases = {target.alt};
}


void Contig::furnish(
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


void extend_ref_coordinate(
    const Variant& target, 
    const UserParams& user_params,
    Coord& coord
)
{  
    // currently only supported for unspliced 
    if (coord.size() != 1) return;

    int ext = 0;
    if (coord[0].first > target.pos)
    {
        ext = coord[0].first -  target.pos + user_params.local_thresh;
        if (target.pos > user_params.local_thresh)
        {
            coord[0].first -= ext;
        }
        else
        {   
            coord[0].first = 0;
        }
    }

    if (coord[0].second < target.pos)
    {
        ext = target.pos - coord[0].second + user_params.local_thresh;
        coord[0].second += ext; // TODO check for chromsom end behavior!!!!!
    }
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
        if (read.overall_lq_rate < user_param.lq_rate_thresh
            && read.seq.find('N') == std::string_view::npos)
        {
            clean_targets.push_back(read);
        }
    }
    clean_targets.shrink_to_fit();
    
    //use most centered one
    if (clean_targets.empty())
    {
        std::sort(
            targets.begin(), targets.end(), [](const Read& a, const Read& b)
            {
                return a.central_score > b.central_score ? true : false;
            }
        );
        clean_targets.push_back(targets[0]);
    }
    
    if (clean_targets.size() == 1)
    {
        std::swap(seeds, clean_targets);
        return;
    }
    
    std::vector<std::string_view> non_ref_sigs;
    non_ref_sigs.reserve(clean_targets.size());
    for (const auto& read : clean_targets)
    {
        non_ref_sigs.push_back(read.non_ref_signature);
       // non_ref_sigs.push_back(read.splice_signature);
    }
   
    std::string_view common_sig = find_commonest_str(non_ref_sigs);
    
    Reads tmp;
    tmp.reserve(clean_targets.size());
    for (const auto& read: clean_targets)
    {
        if (static_cast<std::string_view>(read.non_ref_signature) == common_sig)
        {
            tmp.push_back(read);
        }
    }
    
    sort_by_start(tmp);
    
    std::mt19937 engine(123); 
    if (tmp.size() > seed_size) // seed_size>2
    {     
        std::shuffle(
            tmp.begin(), tmp.end(), engine
        );
        
        for (size_t i = 0; i < seed_size; ++i)
        {
            seeds.push_back(tmp[i]);
        }
        
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
            std::string(seed.seq),
            std::string(seed.base_quals),
            find_target_start(seed)
        );
    }

    return merge_reads(inputs);
}


void set_ref_info(
    Contig& contig, 
    const Coord& coord, 
    LocalReference& loc_ref
)
{
    contig.coordinates = coord;
    
    std::string ref_seq;
    for (const auto& c : coord)
    {
        ref_seq += loc_ref.fasta.getSubSequence(
            loc_ref.chrom, c.first - 1, (c.second -  c.first + 1)
        );
    } 
    contig.ref_seq = ref_seq;
}


std::string_view common_splice_ptrn(const Reads& reads)
{
    std::vector<std::string_view> common_spls;
    common_spls.reserve(reads.size());
    for (const auto& read : reads)
    {
        common_spls.push_back(read.splice_signature);
    }

    return find_commonest_str(common_spls);
}


void make_contig(
    Contig& contig,
    const Variant& target, 
    Reads& targets, 
    const UserParams& user_params, 
    LocalReference& loc_ref
)
{
   Reads seeds;
   size_t seed_size = 5;
   find_seed_indel_reads(seeds, targets, user_params, seed_size);

   Seq _merged = merge_target_seeds(seeds);
        
   Coord coord;
   unite_coordinates(coord, seeds);
   set_ref_info(contig, coord, loc_ref);

   contig.furnish(_merged, target);
}


void mock_target_seq(
    Contig& contig,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    size_t n = 2;
    int offset = target.is_complex ?  1 : 0;
    
    int _start = (target.pos <= int(offset + user_params.kmer_size * n)) 
               ? 0 : target.pos - offset - user_params.kmer_size * n;
    int _len = (target.pos - _start < int(user_params.kmer_size * n)) 
               ? target.pos - _start : user_params.kmer_size * n;
    
    std::string lt_frag = loc_ref.fasta.getSubSequence(
        loc_ref.chrom, _start, _len
    );

    contig.mock_lt_len = lt_frag.size();
    contig.mock_lt_end_idx = contig.mock_lt_len - 1; 
    
    std::string mid_frag = "";
    if (target.is_complex)
    {
        mid_frag = target.alt;
    }
    else
    { 
        mid_frag = target.is_ins ? target.alt.substr(1) : "";
    }

    contig.mock_mid_len = mid_frag.size();
    contig.mock_rt_start_idx = contig.mock_lt_end_idx + contig.mock_mid_len + 1;

    std::string rt_frag = loc_ref.fasta.getSubSequence(
        loc_ref.chrom, 
        target.variant_end_pos - 1, 
        user_params.kmer_size * n
    );
    
    contig.mocked_seq = lt_frag + mid_frag + rt_frag;
    contig.mock_len = contig.mocked_seq.size();
    contig.mock_rt_len = contig.mocked_seq.size() - contig.mock_rt_start_idx;
    
    std::string mocked_q(contig.mock_len, 'F');
    contig.m_quals = mocked_q;
        
    int mock_start = target.pos - offset - user_params.kmer_size * n;
    contig.mocked_ref = loc_ref.fasta.getSubSequence(
        loc_ref.chrom,
        _start,
        contig.mock_len + target.ref_len// for long del
    );
    
    contig.mocked_coord.emplace_back(
        mock_start + 1, mock_start + contig.mock_len + target.ref_len
    );
}


// used for non complex inputs
void prefilter_candidates(
    Contig& contig,
    Reads& candidates,
    Reads& non_targets,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    mock_target_seq(contig, target, user_params, loc_ref);

    Kmers target_kmers;
    diff_kmers(
        contig.mocked_seq, contig.mocked_ref,  user_params.kmer_size, target_kmers
    );

    for (auto& read : candidates)
    {
        read.kmer_score = count_kmer_overlap(read.seq, target_kmers);
        
        // kmer becomes zero for immediate variant clusters
        /* 
        if (read.kmer_score == 0)
        {
            if (!read.dist_to_non_target || !read.dist_to_clip)
            {
                std::cout << "get pseudo" << std::endl;
                read.kmer_score = 1; //pseudo count   
            }    
        }*/
              
    }
    
    sort_by_kmer(candidates);

    for (Reads::reverse_iterator i = candidates.rbegin(); 
        i != candidates.rend(); ++i)
    {
        if ((*i).kmer_score > 0)
        {
            break; //stop if candidate with kmer > 0
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
    
    candidates.shrink_to_fit();
}


void extract_indels(
    std::vector<Variant>& indels,
    const std::vector<Variant>& decomposed
)
{                     
    for (const auto& v : decomposed)
    {
        if (!v.is_substitute) indels.push_back(v);
    }
}


void prefilter_cplx_candidates(
    Contig& contig,
    Reads& candidates,
    Reads& non_targets,
    bool& has_pass_candidates,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref,
    std::vector<Variant>& decomposed
)
{   
    mock_target_seq(contig, target, user_params, loc_ref);
    
    Kmers target_kmers;
    diff_kmers(
        contig.mocked_seq, contig.mocked_ref,  user_params.kmer_size, target_kmers
    );
   
    std::vector<Variant> indels;
    to_simple_variants(user_params, contig, loc_ref, decomposed);   
    
    //extract_indels(indels, decomposed);    

    std::vector<Variant> shared;
    for (auto& read : candidates)
    {
        if (read.may_be_complex && !read.variants.empty())
        {
            find_shared_variants(shared, read.variants, decomposed);
            if (shared == decomposed)
            {
                read.kmer_score = 255; //pseudo count
            }           
        }
        shared.clear(); 
       
        if (read.kmer_score <= 0)
        { 
            read.kmer_score = count_kmer_overlap(read.seq, target_kmers);           
        }
    } 
    
    sort_by_kmer(candidates);
    
    for (Reads::reverse_iterator i = candidates.rbegin(); 
        i != candidates.rend(); ++i)
    {
        if ((*i).kmer_score > 0)
        {
            has_pass_candidates = true;
            break; //stop if candidate with kmer > 0
        }
        else
        {
            (*i).is_deprioritized = true;
            
            /* judge on kmer = 0 may be too harsh
            non_targets.insert(
                non_targets.end(), 
                std::make_move_iterator(candidates.rbegin()), 
                std::make_move_iterator(candidates.rbegin() + 1)
            );
            candidates.pop_back();*/  
        }
    }
    
    //candidates.shrink_to_fit();
}


void prioritize_reads_for_contig_construction(
    Reads& prioritized, 
    const Reads& _candidates,
    const UserParams& user_params
)
{
    Reads candidates;
    for (auto & read : _candidates)
    {
        if (!read.is_deprioritized) candidates.push_back(read);     
    }
    
    std::string_view common_spl_ptrn = common_splice_ptrn(candidates);    
    
    bool is_cplx_gap_matched = (candidates[0].kmer_score == 255);
     
    prioritized.reserve(candidates.size());
    for (const auto& read : candidates)
    {
        if (static_cast<std::string_view>(read.splice_signature)== common_spl_ptrn)
        {
            if (is_cplx_gap_matched)
            {
                if (read.kmer_score == 255) prioritized.push_back(read);
            }
            else if (read.overall_lq_rate < user_params.lq_rate_thresh)
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
        const size_t size_thresh = 5;
        const size_t tmp_range = std::min(size_thresh, tmp.size());
        while (i < tmp_range)
        {
            transfer_elem(top_tens, tmp, i);
            ++i; 
        }
     }

     sort_by_start(top_tens);
}


void suggest_contig(
    const Variant& target,
    Contig& contig,
    Reads& candidates,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    contig.by_kmer_suggestion = true;   
    
    Reads prioritized;
    prioritize_reads_for_contig_construction(
        prioritized, 
        candidates, 
        user_params
    );
    
    //Contig remains empty
    if (prioritized.empty()) return;
    
    const size_t max_size = 10; // for efficiency
    const size_t search_size = prioritized.size();
    
    Coord coord;
    std::vector<Seq> suggestions;
    suggestions.reserve(max_size);
    
    if (search_size >= max_size)
    {
        /*
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
        */          
        Reads top_tens;       
        //concat_top_tens(lt_tmp, rt_tmp, u_tmp, top_tens);
        
        sort_by_kmer(prioritized);
        for (size_t i= 0; i < max_size; ++i)
        {
            transfer_elem(top_tens, prioritized, i);
        }
        
        //sort_by_start(top_tens);
        for (const auto& read : top_tens)
        {
            
            suggestions.emplace_back(
                std::string(read.seq), 
                std::string(read.base_quals), 
                -1
            );
        }
        
        unite_coordinates(coord, top_tens);
    }
    else
    {   
        sort_by_kmer(prioritized);   
        
        if (prioritized.size() < 3)
        {
            prioritized.push_back(prioritized[prioritized.size() - 1]);
        }

        for (const auto& read : prioritized)
        {
            suggestions.emplace_back(
                std::string(read.seq), 
                std::string(read.base_quals), 
                -1
            );
        }

        unite_coordinates(coord, prioritized);
    }
    
    suggestions.shrink_to_fit();
    
    Seq merged_suggestions = merge_reads(suggestions);
    
    contig.n_seeds = suggestions.size();
    contig.seq = merged_suggestions.seq;
    contig.quals = merged_suggestions.base_quals;
   
    extend_ref_coordinate(target, user_params, coord);
    
    set_ref_info(contig, coord, loc_ref); 
}


bool is_successful_extension(
    const int orig_start,
    const int ext_start,
    const int orig_end,
    const int ext_end,
    std::string_view orig_seq,
    std::string_view ext_seq
)
{
    if (orig_seq.size() >= ext_seq.size()) return false;
    //if (ext_seq.find(orig_seq) == std::string_view::npos) return false;
    if (orig_end - orig_start >= ext_end - ext_start) return false;
    
    return true; 
}


void prep_lt_input(
    Reads& lt_matches,
    int& ext_coord_start,
    std::vector<Seq>& inputs,
    const int ext_sz = 3
)
{
    sort_by_start(lt_matches);
    
    int i = ext_sz;
    const int lt_last = static_cast<int>(lt_matches.size()) - 1;
    while (i >= 0)
    {
        if (lt_last - i >= 0)
        {
            inputs.emplace_back(
                std::string(lt_matches[lt_last - i].seq),
                lt_matches[lt_last - i].base_quals, -1
            );
                    
            if (lt_matches[lt_last - i].aln_start < ext_coord_start)
            {
                ext_coord_start = lt_matches[lt_last - i].aln_start;
            }
        }
        --i; 
    }
}


void prep_rt_input(
    Reads& rt_matches,
    int& ext_coord_end,
    std::vector<Seq>& inputs,
    const size_t ext_sz = 3
)
{
    sort_by_start(rt_matches);
    
    size_t rt_last = rt_matches.size() - 1;
    size_t j = 0;
    while (j <= rt_last && j <= ext_sz)
    {
        inputs.emplace_back(
            std::string(rt_matches[j].seq),
            rt_matches[j].base_quals, -1
        );

        if (ext_coord_end < rt_matches[j].aln_end)
        {
            ext_coord_end =  rt_matches[j].aln_end;
        }
        ++j;
    }
}


void extend_contig(
    const char eval,
    Contig& contig,
    Reads& lt_matches,
    Reads& rt_matches,
    LocalReference& loc_ref
)
{   
    int ext_coord_start = contig.coordinates.front().first;
    int ext_coord_end = contig.coordinates.back().second;
    
    Coord ext_coord;
    std::vector<Seq> inputs;
    switch (eval)
    { 
        case 'L':
        {    
            if (lt_matches.empty()) return;
            
            prep_lt_input(lt_matches, ext_coord_start, inputs);            
            inputs.emplace_back(contig.seq, contig.quals, -1);
            ext_coord.emplace_back(ext_coord_start, ext_coord_end);
            break; 
        }
        case 'R':
        {
            if (rt_matches.empty()) return;
           
            inputs.emplace_back(contig.seq, contig.quals, -1);
            prep_rt_input(rt_matches, ext_coord_end, inputs);
            ext_coord.emplace_back(ext_coord_start, ext_coord_end);
            break;
       }
       case 'E': 
       {
            if (lt_matches.empty())
            {
                inputs.emplace_back(contig.seq, contig.quals, -1);
            }
            else
            {    
                prep_lt_input(lt_matches, ext_coord_start, inputs);        
                inputs.emplace_back(contig.seq, contig.quals, -1);
            }
            
            if (!rt_matches.empty())
            {
                prep_rt_input(rt_matches, ext_coord_end, inputs);
            }
            ext_coord.emplace_back(ext_coord_start, ext_coord_end);
            break;
        }      
    }

    Seq extended = merge_reads(inputs);
    
    bool is_passed = is_successful_extension(
        contig.coordinates.front().first, ext_coord_start,
        contig.coordinates.back().second, ext_coord_end, contig.seq, extended.seq
    );
     
    if (is_passed)
    {    
        contig.seq = extended.seq;
        contig.quals = extended.base_quals;
        contig.coordinates = ext_coord;

        set_ref_info(contig, ext_coord, loc_ref);
    }
}
