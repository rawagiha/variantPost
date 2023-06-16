#include <cmath>
#include <string>
#include <vector> 
#include <climits>
#include <utility>

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


void to_simple_variants(
    const UserParams& user_params,
    const Contig& contig,
    LocalReference& loc_ref,
    std::vector<Variant>& simples
)
{
    Filter filter;
    Alignment alignment;
    local_alignment(
        user_params.match_score,
        user_params.mismatch_penal,
        user_params.gap_open_penal,
        user_params.gap_ext_penal,
        contig.mocked_seq,
        contig.mocked_ref,
        filter,
        alignment
    );

    std::vector<std::pair<char, int>> cigar_vec = to_cigar_vector(
        alignment.cigar_string
    );
    
    move_up_insertion(cigar_vec);
    
    std::vector<int> pos_vec = expand_coordinates(contig.mocked_coord); 
    
    int query_begin = alignment.query_begin;
    if (cigar_vec[0].first == 'S')
    {
        query_begin -= cigar_vec[0].second;
    }
    
    bool include_clips = false;   
    std::string _tmp = ""; //won't need
    parse_variants(
        pos_vec[alignment.ref_begin],
        contig.mocked_ref.substr(alignment.ref_begin),
        contig.mocked_seq.substr(query_begin),
        contig.mocked_seq.substr(query_begin), //for mocking
        cigar_vec,
        loc_ref.dict,
        simples,
        _tmp,
        include_clips
    );
}


void postprocess_alignment(
    AlnResult& rslt, 
    const std::vector<int>& pos_vec,
    const Contig& contig,
    const LocalReference& loc_ref,
    const Alignment& alignment)
{
    rslt.cigar_str = alignment.cigar_string;
    rslt.cigar_vec = to_cigar_vector(rslt.cigar_str);
    splice_cigar(
        rslt.cigar_vec, 
        alignment.ref_begin, 
        pos_vec, 
        contig.coordinates
    );
    move_up_insertion(rslt.cigar_vec);

    int query_begin = alignment.query_begin;
    if (rslt.cigar_vec[0].first == 'S')
    {      
        query_begin -= rslt.cigar_vec[0].second;
    }

    rslt.genomic_start_pos = pos_vec[alignment.ref_begin];
    rslt.seq = contig.seq.substr(query_begin);
    rslt.ref_seq = contig.ref_seq.substr(alignment.ref_begin);
    rslt.quals = contig.quals.substr(query_begin);
    
    bool include_clips = true;
    std::string _tmp = ""; //won't need
    parse_variants(
        rslt.genomic_start_pos,
        rslt.ref_seq,
        rslt.seq,
        rslt.quals,
        rslt.cigar_vec,
        loc_ref.dict,
        rslt.variants,
        _tmp,
        include_clips
    );
}


void gap_penal_grid(
    const int gap_open_penal, 
    const int gap_extension_penal, 
    std::vector<std::pair<int, int>>& grid
)
{
    std::pair<int, int> _default = {gap_open_penal, gap_extension_penal};

    uint8_t max_gap_open = 10, min_gap_open = 2; //configureab;e
    if (gap_open_penal > max_gap_open) max_gap_open = gap_open_penal;
    
    for (uint8_t i = min_gap_open; i <= max_gap_open; i += 2)
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
    const UserParams& user_params,
    const std::vector<std::pair<char, int>>& cigar_vec
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


void annot_end_mappping(AlnResult& rslt, const UserParams& user_params)
{
    size_t last_idx = rslt.cigar_vec.size() - 1;
    if (rslt.cigar_vec[0].first == '=')
    {
        rslt.lt_ref_mapped_len = rslt.cigar_vec[0].second;
        rslt.lt_well_ref_mapped = (
            rslt.lt_ref_mapped_len >= user_params.local_thresh
        ) ? true : false;
    }

    if (rslt.cigar_vec[last_idx].first == '=')
    {
        rslt.rt_ref_mapped_len = rslt.cigar_vec[last_idx].second; 
        rslt.rt_well_ref_mapped = (
            rslt.rt_ref_mapped_len >= user_params.local_thresh
        ) ? true : false;   
    }

    rslt.is_well_ref_mapped = (
        rslt.lt_well_ref_mapped && rslt.rt_well_ref_mapped
    ) ? true : false;
}


bool has_query(
    const Variant& query, 
    const std::vector<Variant>& variants,
    const LocalReference& loc_ref)
{ 
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
        int lt_d =  lt_lim - (variant.rpos + variant.ref_len - 1);   
        d = (lt_d < 0) ? 0 : lt_d;
    }
    else
    {
        variant.set_leftmost_pos(loc_ref);
        int rt_d = variant.lpos - rt_lim;
        d = (rt_d < 0) ? 0 : rt_d;
    } 

    return d;
}


void annot_nearest_variant(
    const int lt_lim,
    const int target_pos,
    const int rt_lim, 
    const LocalReference& loc_ref,
    std::vector<Variant>& variants,
    int& idx,
    int& dist
)
{
    int tmp_dist = INT_MAX;
    for (size_t i = 0; i < variants.size(); ++i)
    {
        tmp_dist = dist_from_target(
            lt_lim, 
            target_pos, 
            rt_lim, 
            loc_ref, 
            variants[i]
        );

        if (tmp_dist < dist)  
        {
            dist = tmp_dist;
            idx = i;
        } 
    }
}


void eval_by_variant(
    AlnResult& rslt,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    if (rslt.variants.empty()) return;

    annot_end_mappping(rslt, user_params);

    annot_nearest_variant(
        target.lpos,
        target.pos,
        target.rpos + target.ref_len - 1,
        loc_ref,
        rslt.variants,
        rslt.idx_to_closest,
        rslt.dist_to_closest
    );

    if (rslt.dist_to_closest > user_params.local_thresh) return;
    
    Variant maybe_target = rslt.variants[rslt.idx_to_closest];
    
    rslt.is_passed = true;
    
    if (target.is_equivalent(maybe_target, loc_ref)) 
    {    
        rslt.has_target = true;
        if (rslt.is_well_ref_mapped) rslt.terminate_search = true;
    }

    //possible cases reach here
    // 1) simple long indels with partial contig
    //    -> do partial contig?
    // 2) simple indels in readends -> realigned differently
    //    -> retarget
} 


void eval_by_variant_lst(
    AlnResult& rslt,
    const std::vector<Variant>& simples,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    if (rslt.variants.empty()) return;

    annot_end_mappping(rslt, user_params);

    size_t total_match = 0, indel_match = 0;
    for (const auto& simple : simples)
    {
        for (const auto& variant : rslt.variants)
        {
            if (simple.is_equivalent(variant, loc_ref))
            {
                ++total_match;
                if (!simple.is_substitute) ++indel_match;
                break; //exit from the inner loop
            } 
        }
    }
    
    if (!total_match) return;
    
    if (total_match == simples.size())
    {
        rslt.has_target = true;
        rslt.is_passed = true;
        rslt.terminate_search = true;
    } 
    else if (!indel_match)
    {   
        rslt.is_passed = true;
    }
}


//to be rename...
void eval_by_aln(
    const Contig& contig,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    
    std::vector<Variant> simples;
    if (target.is_complex)
    {
        to_simple_variants(user_params, contig, loc_ref, simples);
    }
    
    
    // evaluation by grid search
    std::vector<std::pair<int, int>> gap_penals;
    gap_penal_grid(
        user_params.gap_open_penal,
        user_params.gap_ext_penal,
        gap_penals
    );

    std::vector<int> pos_vec = expand_coordinates(contig.coordinates);
    
    Filter filter;
    Alignment aln;
    std::vector<AlnResult> rslts;
    rslts.reserve(gap_penals.size());
    for (const auto& gap_penal : gap_penals)
    { 
        AlnResult rslt;
        
        local_alignment(
            user_params.match_score,
            user_params.mismatch_penal,
            gap_penal.first,
            gap_penal.second,
            contig.seq,
            contig.ref_seq, 
            filter,
            aln
        );

        postprocess_alignment(rslt, pos_vec, contig, loc_ref, aln);
        
        if (target.is_complex)
        {
            eval_by_variant_lst(rslt, simples, user_params, loc_ref);
        }
        else
        {
            eval_by_variant(rslt, target, user_params, loc_ref); 
        }

        if (rslt.terminate_search)
        {    
            std::cout << "annot contig for alingment " << std::endl;
            return;   
        }
        else if (rslt.is_passed)
        {
            rslts.push_back(rslt);
        }

    }
    // if made (target guaranteed)  -> do most stable -> check for extension -> ext -> annot align
    // if suggested -> use most stable with check has or not
    // if made => 
    
} 
