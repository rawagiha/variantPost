#include <cmath>
#include <string>
#include <vector> 
#include <climits>
#include <utility>
#include <string_view>

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
    
    //bool include_clips = false;   
    std::string _tmp = ""; //won't need
    parse_variants(
        pos_vec[alignment.ref_begin],
        contig.mocked_ref.substr(alignment.ref_begin),
        contig.mocked_seq.substr(query_begin),
        contig.mocked_seq.substr(query_begin), //for mocking
        cigar_vec,
        loc_ref.dict,
        simples,
        _tmp
    );

    std::sort(simples.begin(), simples.end());
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
    
    //bool include_clips = true;
    std::string _tmp = ""; //won't need
    parse_variants(
        rslt.genomic_start_pos,
        rslt.ref_seq,
        rslt.seq,
        rslt.quals,
        rslt.cigar_vec,
        loc_ref.dict,
        rslt.variants,
        _tmp
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


inline bool pass_gap_check(AlnResult& rslt, const UserParams& user_params)
{
    if (rslt.variants.empty()) return false;
    
    annot_end_mappping(rslt, user_params);
    if (rslt.is_well_ref_mapped && !has_gaps(rslt.cigar_str))
    {
        return false;
    }
    
    return true;
}


void eval_by_variant(
    AlnResult& rslt,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
    if (!pass_gap_check(rslt, user_params)) return;
    
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

} 


void eval_by_variant_lst(
    AlnResult& rslt,
    const std::vector<Variant>* p_decomposed,
    const UserParams& user_params,
    LocalReference& loc_ref,
    int& alternative_pos 
)
{
    if (!pass_gap_check(rslt, user_params)) return;

    std::vector<Variant> shared;
    find_shared_variants(shared, rslt.variants, *p_decomposed);
    
    size_t total_match = shared.size(), indel_match = 0;
    
    //find indel with biggest change
    int alternative_indel_len = 0;
    for (const auto& v : shared)
    {
        if (!v.is_substitute)
        {
            ++indel_match;
            if (alternative_indel_len <= v.ref_len + v.alt_len)    
            {    
                alternative_pos = v.pos;
                alternative_indel_len = (v.ref_len + v.alt_len);
            }
         }
    }
    /*
    for (const auto& simple : *p_decomposed)
    {
        
        
        for (const auto& variant : rslt.variants)
        {
            if (simple.is_equivalent(variant, loc_ref))
            {
                ++total_match;
                if (!simple.is_substitute) 
                {    
                    ++indel_match;
                    if (alternative_indel_len <= (simple.ref_len + simple.alt_len))
                    {
                        alternative_pos = simple.pos;
                        alternative_indel_len = (simple.ref_len + simple.alt_len);
                    }
                
                }
                break; //exit from the inner loop
            } 
        }
    }*/
    
    if (!total_match) return;
    
    if (total_match == p_decomposed->size())
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


void annot_alignment(Contig& contig, const AlnResult& rslt)
{
    int genomic_pos = rslt.genomic_start_pos; 

    contig.seq = rslt.seq;
    contig.len = contig.seq.size();
    contig.ref_seq = rslt.ref_seq;
    contig.quals = rslt.quals;
    
    int op_len;
    char op = '\0';
    size_t seq_idx = 0, ref_idx = 0, v_idx = 0;    
    for (const auto& cigar : rslt.cigar_vec)
    {  
        op = cigar.first;
        op_len = cigar.second;
       
        switch (op)
        {
            case '=':
            case 'X':
                for (int i = 0; i < op_len; ++i)
                {                   
                    contig.positions.push_back(genomic_pos);
                    contig.ref_bases.push_back(contig.ref_seq.substr(ref_idx, 1));
                    contig.alt_bases.push_back(contig.seq.substr(seq_idx, 1));
                    contig.base_quals.push_back(contig.quals.substr(seq_idx, 1));
                    
                    ++genomic_pos;
                    ++ref_idx;
                    ++seq_idx;

                    if (op == 'X') 
                    {    
                        ++v_idx;
                    }
                }
                break;
            case 'I':
                contig.alt_bases.pop_back();
                contig.alt_bases.push_back(contig.seq.substr(seq_idx - 1, op_len + 1));
                contig.base_quals.pop_back();
                contig.base_quals.push_back(contig.quals.substr(seq_idx - 1, op_len + 1));
                
                seq_idx += op_len;
                ++v_idx;
                break;
            case 'D':
                contig.ref_bases.pop_back();
                contig.ref_bases.push_back(rslt.variants[v_idx].ref);

                ref_idx += op_len;
                genomic_pos += op_len;
                ++v_idx;
                break;
            case 'N':
                contig.skip_starts.push_back(genomic_pos);
                genomic_pos += op_len;
                contig.skip_ends.push_back(genomic_pos - 1);
                break;
            case 'S':
                seq_idx += op_len;
                break;
        }
    }
}


void update_contig_layout(Contig& contig, const int pos)
{
    auto it = std::find(
        contig.positions.begin(), 
        contig.positions.end(), 
        pos
    );
    
    size_t idx;
    //keep this block for furtehr cheking
    if (it == contig.positions.end()) 
    {
        int diff = INT_MAX;
        size_t tmp_idx = 0;
        for (size_t i = 0; i < contig.positions.size(); ++i)
        {
            if (contig.ref_bases[i].size() != contig.alt_bases[i].size()
                && diff > std::abs(contig.positions[i] - pos))
            { 
                diff = std::abs(contig.positions[i] - pos);
                tmp_idx = i;
            } 
        }
        
        idx = tmp_idx;    
    }    
    
    idx = std::distance(contig.positions.begin(), it);

    contig.lt_len = 0;
    contig.mid_len =0;
    for (size_t i = 0; i <= idx; ++i)
    {
        if (i == idx) 
        {
            size_t _len = contig.alt_bases[i].size();
            if (_len > 1) contig.mid_len = (_len - 1);
            contig.lt_len += 1;
        }
        else
        {
            contig.lt_len += contig.alt_bases[i].size();
        }
    }
    contig.lt_end_idx = contig.lt_len - 1;
    contig.rt_start_idx = contig.lt_len + contig.mid_len;
    contig.rt_len = contig.len - contig.rt_start_idx; 
}


char eval_by_aln(
    Contig& contig,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref,
    const std::vector<Variant>* p_decomposed
)
{
    
    //for (auto& h : *p_decomposed) std::cout << h.pos << " " << h.ref << " " << h.alt << std::endl;
    
    //std::vector<Variant> simples;
    //if (target.is_complex)
    //{
    //    to_simple_variants(user_params, contig, loc_ref, simples);
    //}
    
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
        
        int target_pos = target.pos;
        //if (target.is_complex)
        if (p_decomposed != NULL)
        {
            //target_pos may change -> one of decomposed simple indels
            eval_by_variant_lst(rslt, p_decomposed, user_params, loc_ref, target_pos);
        }
        else
        { 
            eval_by_variant(rslt, target, user_params, loc_ref);
        }

        if (rslt.terminate_search)
        {    
            annot_alignment(contig, rslt);
            update_contig_layout(contig, target_pos);
            //std::cout << "terminated " << std::endl;
            return 'A';   
        }
        else if (rslt.is_passed)
        {
            rslts.push_back(rslt);
        }
    }
    
    if (rslts.empty()) return 'C';

    std::vector<std::string_view> aln_ptrns;
    for (const auto& rslt : rslts)
    {
        aln_ptrns.push_back(rslt.cigar_str);
    }
    
    //most stable alignment
    std::string_view common_ptrn = find_commonest_str(aln_ptrns);
    auto it = std::find(aln_ptrns.begin(), aln_ptrns.end(), common_ptrn);
    auto rslt = rslts[std::distance(aln_ptrns.begin(), it)];
    
    if (contig.by_kmer_suggestion)
    {
        if (rslt.has_target) 
        {
            annot_alignment(contig, rslt);
            update_contig_layout(contig, target.pos);
            //std::cout << "assembled " << std::endl;
            return 'A';
        }
        //count check by mocked seq
        //std::cout << "mock used" << std::endl;
        return 'B';   
    }
     
    //spliced -> no extension
    for (const auto& c : rslt.cigar_vec)
    {
        if (c.first == 'N') return 'A';
    }
    
    //do extension
    if (!rslt.lt_well_ref_mapped && rslt.rt_well_ref_mapped)
    {
        return 'L';
    }
    else if (rslt.lt_well_ref_mapped && !rslt.rt_well_ref_mapped)
    {
        return 'R';
    }
    else return 'E';
}


void aln_extended_contig( 
    Contig& contig,
    const Variant& target,
    const UserParams& user_params,
    LocalReference& loc_ref
)
{
   std::vector<int> pos_vec = expand_coordinates(contig.coordinates);
    
   Filter filter;
   Alignment aln;
   AlnResult rslt; 
    
   local_alignment(
       user_params.match_score,
       user_params.mismatch_penal,
       user_params.gap_open_penal,
       user_params.gap_ext_penal,
       contig.seq,
       contig.ref_seq, 
       filter,
       aln
   );

   postprocess_alignment(rslt, pos_vec, contig, loc_ref, aln);
   annot_alignment(contig, rslt);
   update_contig_layout(contig, target.pos);
}
