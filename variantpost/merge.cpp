#include <deque>
#include <string>
#include <vector>
#include <numeric>
#include <string.h>
#include <algorithm>

#include "util.h"
#include "merge.h"
#include "ssw/ssw_cpp.h"


static char BASES[5] = {'A', 'C', 'G', 'T', 'N'};


Seq::Seq() {}


Seq::Seq(
    const std::string& seq,
    const std::string& base_quals,
    int target_start
) : seq(seq),
    base_quals(base_quals),
    target_start(target_start) {}


bool Seq::empty() const
{
    return seq.empty();
}


template <typename T>
std::vector<int> get_max_indices(const T& arr)
{
    std::vector<int> indices;

    auto it = std::max_element(std::begin(arr), std::end(arr));
    while (it != std::end(arr)) 
    {
        indices.push_back(std::distance(std::begin(arr), it));
        it = std::find(std::next(it), std::end(arr), *it);
    }
    return indices;
} 


template<typename T>
double average(const std::vector<T>& v)
{
    if (v.empty()) return 0.0;
    
    const double cnt = static_cast<double>(v.size());
    
    return std::accumulate(v.begin(), v.end(), 0.0) / cnt;
}


/*
inline bool has_gaps(std::string& cigar_string)
{   
    return (cigar_string.find("I") != std::string::npos 
         || cigar_string.find("D") != std::string::npos);
}
*/

struct Overlap
{   
    int index;
    int ref_start;
    int ref_end;
    int query_start;
    int query_end;
    int target_start;

    Overlap
    (   const int index, 
        const int ref_start, 
        const int ref_end, 
        const int query_start, 
        const int query_end,
        int target_start
    ) : index(index), 
        ref_start(ref_start), ref_end(ref_end), 
        query_start(query_start), query_end(query_end), 
        target_start(target_start) {}
};


struct BaseCount
{
    size_t a = 0, c = 0, g = 0, t = 0, n = 0;
    std::vector<int> aq, cq, gq, tq, nq;
    std::vector<int> start_checks;
    
    void add(char base, char _qual, int is_target_start=0)
    {
        start_checks.push_back(is_target_start);
          
        int qual = static_cast<int>(_qual); 
         
        switch(base)
        {
            case 'A':
                ++a;
                aq.push_back(qual);   
                break;
            case 'C':
                ++c;
                cq.push_back(qual);
                break;
            case 'G':
                ++g;
                gq.push_back(qual);
                break;
            case 'T':
                ++t;
                tq.push_back(qual);
                break;
            default:
                ++n;
                nq.push_back(qual);
                break; 
        }
    }

    std::pair<char, char> get_consensus() const 
    {
        size_t counts[] = {a, c, g, t, n};
        double quals[] = {
            average(aq), 
            average(cq), 
            average(gq), 
            average(tq), 
            average(nq)
        };

        std::vector<int> max_cnt_indices = get_max_indices(counts);

        int max_idx = -1;
        double max_q = 0.0, tmp_q = 0.0;
        for (const auto idx : max_cnt_indices)
        {
            tmp_q = quals[idx];
            if (max_q <= tmp_q) 
            {
                max_q = tmp_q;
                max_idx = idx;
            }     
        } 
        
        return {BASES[max_idx], static_cast<char>(quals[max_idx])};
    }
};


void find_overlaps(std::vector<Overlap>& overlaps, const std::vector<Seq>& inputs)
{
    const uint8_t match_score = 2, mismatch_penalty = 2;
    const uint8_t gap_open_penalty = 255, gap_extention_penalty = 255;

    Filter filter;
    Alignment aln;   
    Aligner aligner(
        match_score, 
        mismatch_penalty, 
        gap_open_penalty, 
        gap_extention_penalty
    );
    
    int32_t mask_len = 0;
    int32_t prev_start = 0, prev_end = inputs[0].seq.size();
    size_t j = 0, n = inputs.size();
    for (size_t i = 0; i < n - 1;)
    {
        
        j = i + 1;
   
        mask_len = strlen(inputs[j].seq.c_str()) / 2;
        mask_len = mask_len < 15 ? 15 : mask_len;
            
        // align next seq (as query) vs current seq (as ref)
        aligner.Align(
            inputs[j].seq.c_str(),                 
            inputs[i].seq.c_str(), 
            inputs[i].seq.size(), 
            filter, 
            &aln, 
            mask_len
        );
            
        
        if (has_gaps(aln.cigar_string)) return;
        if (aln.ref_begin >=  prev_end || prev_start >= aln.ref_end) return; 
        
        overlaps.emplace_back(
            i, 
            aln.ref_begin, 
            aln.ref_end, 
            aln.query_begin, 
            aln.query_end, 
            inputs[i].target_start
        );
        
        i = j;
        prev_start = aln.query_begin;
        prev_end = aln.query_end;
    }
}


void merge_to_right(
    const int prev_start, 
    const int curr_start,
    const int prev_end, 
    const int curr_end,
    const int target_start,
    const std::string& seq,
    const std::string& seq_qual,
    int& merge_start, 
    int& merge_end,
    std::deque<BaseCount>& base_cnts,
    bool& is_target_start
)
{
    merge_start += (curr_start - prev_start);
    merge_end = base_cnts.size(); 
    const int _merge_end = merge_end;
    int merge_idx = 0;

    for (int j = 0; j <= (prev_end - curr_start); ++j)
    {
        if (merge_start + j < merge_end) 
        {    
            if (target_start == curr_start + j) is_target_start = 1;
            base_cnts[merge_start + j].add(
                seq[curr_start + j], 
                seq_qual[curr_start + j], 
                is_target_start
            );
            is_target_start = 0;
        }
        
        merge_idx = merge_start + j;
                    
    }
                                        
    for (int j = 1; j <= (curr_end - prev_end); ++j)
    {
        if (target_start == prev_end + j) is_target_start = 1;
                            
        if (merge_idx + j < _merge_end)
        {
            base_cnts[merge_idx + j].add(
                seq[prev_end + j], 
                seq_qual[prev_end + j], 
                is_target_start
            );
        }
        else
        {
            BaseCount bc;
            bc.add(
                seq[prev_end + j], 
                seq_qual[prev_end + j], 
                is_target_start
            );
            base_cnts.push_back(bc);
        }                           
        
        is_target_start = 0;
    }
}


Seq merge_reads(const std::vector<Seq>& inputs)
{
    if (inputs.size() == 1)
    {
        Seq _merged(
            inputs[0].seq, 
            inputs[0].base_quals, 
            inputs[0].target_start
        );
        
        return _merged;
    }
    
    std::deque<BaseCount> base_cnts;
    std::vector<Overlap> overlaps; 
    find_overlaps(overlaps, inputs);
    
    if (overlaps.empty())
    {
        Seq _merged(
            inputs[0].seq, 
            inputs[0].base_quals, 
            inputs[0].target_start
        );
        
        return _merged;
    }

    int merge_start = 0, merge_end = 0;
    int curr_start = 0, curr_end = 0;
    int prev_start = 0, prev_end = 0;
    int target_start = -1;
    int n_overlap_reads = overlaps.size();
    bool is_target_start = false;

    for (int i = 0; i < n_overlap_reads; ++i)
    {
        curr_start = overlaps[i].ref_start;
        curr_end = overlaps[i].ref_end;
        
        target_start = overlaps[i].target_start;
        
        std::string seq = inputs[overlaps[i].index].seq;
        std::string seq_qual = inputs[overlaps[i].index].base_quals;

        if (i == 0)
        {
            //aligned segment only// 
            for (int j = 0; j <= (curr_end - curr_start); ++j)
            {
                if (target_start == curr_start + j) is_target_start = 1;
                
                BaseCount bc;
                bc.add(
                    seq[curr_start + j], 
                    seq_qual[curr_start + j], 
                    is_target_start
                );
                base_cnts.push_back(bc);
                
                is_target_start = 0;
            }
        }
        else 
        {
            if (prev_start <= curr_start)
            {
                merge_to_right(
                    prev_start, 
                    curr_start, 
                    prev_end, 
                    curr_end, 
                    target_start, 
                    seq, seq_qual, 
                    merge_start, 
                    merge_end, 
                    base_cnts, 
                    is_target_start
                );
            }
            else
            {
                int _ps = prev_start;                
                merge_to_right(
                    prev_start, 
                    _ps, 
                    prev_end, 
                    curr_end, 
                    target_start, 
                    seq, 
                    seq_qual, 
                    merge_start, 
                    merge_end, 
                    base_cnts, 
                    is_target_start
                );
                
                int lt_move = prev_start - curr_start;
                const int _merge_start = merge_start;
                for (int j = 0; j <= lt_move; ++j)
                {
                    if (target_start == prev_start - j) is_target_start = 1;
                        
                    if (_merge_start - j >= 0)
                    {
                        base_cnts[_merge_start - j].add(
                            seq[prev_start - j], 
                            seq_qual[prev_start - j], 
                            is_target_start
                        );
                        
                        if (j) merge_start -= 1;
                    }
                    else
                    {
                        if (target_start == prev_start - j) is_target_start = 1;

                        BaseCount bc;
                        bc.add(
                            seq[prev_start - j], 
                            seq_qual[prev_start - j], 
                            is_target_start
                        );
                        base_cnts.push_front(bc);
                        is_target_start = 0;
                        merge_start = 0;
                    }    
                } 
            }
        }
        
        prev_start = overlaps[i].query_start;
        prev_end = overlaps[i].query_end;   
    }

    std::string merged_read = "";
    std::string merged_qualities = "";
    std::vector<int> n_start_checks;
    for (const auto& b : base_cnts)
    {
        merged_read += b.get_consensus().first;
        merged_qualities += b.get_consensus().second;
        n_start_checks.push_back(
            std::accumulate(b.start_checks.begin(), b.start_checks.end(), 0)
        );
    }

    int target_start_idx = get_max_indices(n_start_checks)[0];   
    
    Seq __merged(merged_read, merged_qualities, target_start_idx);   
    
    return __merged; 
} 
