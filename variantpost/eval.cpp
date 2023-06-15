#include <cmath>
#include <string>
#include <vector> 
#include <climits>

#include "eval.h"
#include "util.h"
#include "merge.h"
#include "contig.h"


void local_alignment(
    const uint8_t match_score,
    const uint8_t mismatch_penal,
    const uint8_t gap_open_penal,
    const uint8_t gap_ext_penal,
    const std::string& seq,
    const std::string& ref_seq,
    const Filter& filter,
    Alignment& alignment
)
{
    Aligner aligner(
        match_score, 
        mismatch_penal, 
        gap_open_penal, 
        gap_ext_penal
    ); 
       
    int32_t mask_len = strlen(seq.c_str()) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;

    aligner.Align(
        seq.c_str(), 
        ref_seq.c_str(), 
        ref_seq.size(),
        filter, 
        &alignment, 
        mask_len
    ); 
}


void gap_penal_grid(
    const int gap_open_penal, 
    const int gap_extension_penal, 
    std::vector<std::pair<int, int>>& grid
)
{
    std::pair<int, int> _default = {gap_open_penal, gap_extension_penal};

    uint8_t max_gap_open = 5, min_gap_open = 3;
    if (gap_open_penal > max_gap_open) max_gap_open = gap_open_penal;
    
    for (uint8_t i = min_gap_open; i <= max_gap_open; ++i)
    {      
        grid.emplace_back(i, 1);
        grid.emplace_back(i, 0);
    }
    
    if (std::find(grid.begin(), grid.end(), _default) ==  grid.end())
    {
        grid.push_back(_default);
    }
}


bool has_mapped_ends(
    UserParams& user_params,
    std::vector<std::pair<char, int>>& cigar_vec
)
{
    size_t last_idx = cigar_vec.size() - 1;
    if (cigar_vec[0].first == '='
        && cigar_vec[0].second >= user_params.local_thresh
        && cigar_vec[last_idx].first == '='
        && cigar_vec[last_idx].second >= user_params.local_thresh)
    {
        return true;
    }
    else return false;
}


int dist_from_target(
    const int lt_lim, 
    const int target_pos,
    const int rt_lim, 
    const LocalReference& loc_ref,
    Variant& variant
)
{ 
    int d;
    if (variant.pos < target_pos)
    {
        variant.set_rightmost_pos(loc_ref);
        int lt_d =  lt_lim - (variant.rpos + variant.ref_len - 1)   
        d = (lt_d < 0) ? 0 : lt_d;
    }
    else
    {
        variant.set_leftmost_pos(loc_ref);
        int rt_d = variant.lpos - rt_lim
        d = (rt_d < 0) ? 0 : rt_d;
    } 

    return d;
}

int nearest_variant_idx(
    const int lt_lim,
    const int target_pos,
    const int rt_lim, 
    const LocalReference& loc_ref,
    std::vector<Variant>& variants
)
{
    int idx = 0, tmp_dist = INT_MAX, dist = INT_MAX;
    for (int i = 0; i < variants.size(); ++i)
    {
        tmp_dist = dist_from_target(
            lt_lim, 
            target_pos, 
            rt_lim, 
            loc_ref, 
            variant
        );

        if (tmp_dist < dist)  
        {
            dist = tmp_dist;
            idx = i;
        } 
    }
    
    return idx; 
}

bool isok(
    const Variant& target,
    const std::vector<std::pair<char, int>>& cigar_vec,
    const std::vector<Variant>& variants,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    if (variants.empty()) return false;
    
    bool is_well_ended = has_mapped_ends(user_params, cigar_vec);    

    int _idx = nearest_variant_idx(
        target.lpos,
        target.pos,
        target.rpos + target.ref_len - 1,
        loc_ref,
        variants
    );
    
    Variant maybe_target = variants[_idx];

    if (target.is_equivalent(maybe_target, loc_ref) 
        && is_well_ended) return true;
          
    if (std::abs(maybe_target.pos - target.pos) > user_params.local_thresh) 
    {
        return false;
    }                 
}

//to be rename...
void eval_by_aln(
    const UnalignedContig& u_contig,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    std::vector<std::pair<int, int>> gap_penals;
    gap_penal_grid(
        user_params.gap_open_penal,
        user_params.gap_ext_penal,
        gap_penals
    );

    std::vector<int> genomic_pos = expand_coordinates(u_contig.coordinates);
    
    Filter filter;
    Alignment aln;
    std::vector<Variant> _variants = {};
    std::vector<std::pair<char, int>> _cigar_vec = {};
    std::string variant_quals = "";

    for (const auto& gap_penal : gap_penals)
    { 
        local_alignment(
            user_params.match_score,
            user_params.mismatch_penal,
            gap_penal.first,
            gap_penal.second,
            u_contig.seq,
            u_contig.ref_seq, 
            filter,
            aln
        );

        _cigar_vec = to_cigar_vector(aln.cigar_string);
        splice_cigar(_cigar_vec, aln.ref_begin, genomic_pos, u_contig.coordinates);
        move_up_insertion(_cigar_vec);
   
        int begin_with_offset = aln.query_begin;
        if (_cigar_vec[0].first == 'S') begin_with_offset -= _cigar_vec[0].second;
        
        
        parse_variants(
            genomic_pos[aln.ref_begin], 
            u_contig.ref_seq.substr(aln.ref_begin), 
            u_contig.seq.substr(begin_with_offset), 
            u_contig.quals.substr(begin_with_offset), 
            _cigar_vec, loc_ref.dict, _variants, variant_quals, true
        );
        
        
        std::cout << aln.cigar_string << std::endl;;
        for (auto& v : _variants)
        {
            std::cout << v.pos << " " << v.ref << " " << v.alt << " -- ";
        }
        std::cout << std::endl;
        
        _variants.clear();    
    }

} 
