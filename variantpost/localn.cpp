#include <deque>
#include <random>
#include <string>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <string.h>

#include "util.h"
#include "localn.h"
#include "ssw/ssw_cpp.h"

static char BASES[5] = {'A', 'C', 'G', 'T', 'N'};


template <typename T>
std::vector<int> get_max_indices(const T & arr)
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
double average(const std::vector<T>  & v)
{
    if (v.empty()) return 0.0;
    const double cnt = static_cast<double>(v.size());
    return std::accumulate(v.begin(), v.end(), 0.0) / cnt;
}


inline bool has_gaps(std::string & cigar_string)
{   
    return (cigar_string.find("I") != std::string::npos 
         || cigar_string.find("D") != std::string::npos);
}


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
        target_start(target_start){}
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
        double quals[] = {average(aq), average(cq), average(gq), average(tq), average(nq)};

        std::vector<int> max_cnt_indices = get_max_indices(counts);

        int max_idx;
        if (max_cnt_indices.size() > 1)
        {
            std::vector<int> max_qual_indices = get_max_indices(quals);
            max_idx = max_qual_indices[0];
        }
        else max_idx = max_cnt_indices[0];
        
        return {BASES[max_idx], static_cast<char>(quals[max_idx])};
    }

};


std::vector<Overlap> find_overlaps(const std::vector<InputRead> & inputs)
{
    std::vector<Overlap> overlaps;   
    
    const uint8_t match_score = 3, mismatch_penalty = 10;
    const uint8_t gap_open_penalty = 255, gap_extention_penalty = 255;

    Filter filter;
    Alignment aln;   
    Aligner aligner(
                match_score, mismatch_penalty, 
                gap_open_penalty, gap_extention_penalty
            );
    
    size_t n = inputs.size(); 
    for (size_t i = 0; i < n - 1;)
    {
        
        size_t j = i;
        bool is_gapped = true;
        do
        {
            ++j;
            int32_t mask_len = strlen(inputs[j].read_seq.c_str()) / 2;
            mask_len = mask_len < 15 ? 15 : mask_len;
            
            // align next seq (as query) vs current seq (as ref)
            aligner.Align(
                        inputs[j].read_seq.c_str(),                 
                        inputs[i].read_seq.c_str(), inputs[i].read_seq.size(), 
                        filter, &aln, mask_len
                    );
            
            std::string cigar_string = aln.cigar_string;
            if (!has_gaps(cigar_string)) is_gapped = false;
        }   
        while (j < n - 1 && is_gapped);
        
        overlaps.emplace_back(
                        i, 
                        aln.ref_begin, aln.ref_end, aln.query_begin, aln.query_end, 
                        inputs[i].target_start
                    );
        
        i = j;
            
    }

    return overlaps;
}


// 
InputRead pairwise_stitch(const std::vector<InputRead> & inputs, bool lt_extention)
{
    std::vector<Overlap> overlaps = find_overlaps(inputs);   

    if (lt_extention)
    {
        std::string lt_ext = inputs[1].read_seq.substr(0, overlaps[0].query_start);
        std::string lt_ext_q = inputs[1].base_qualities.substr(0, overlaps[0].query_start);

        InputRead lt_extended{lt_ext + inputs[0].read_seq, lt_ext_q + inputs[0].base_qualities, -1};

        return lt_extended;
    }
    else
    {
        std::string rt_ext = inputs[1].read_seq.substr(overlaps[0].query_end + 1);
        std::string rt_ext_q = inputs[1].base_qualities.substr(overlaps[0].query_end + 1);

        InputRead rt_extended{inputs[0].read_seq + rt_ext, inputs[0].base_qualities + rt_ext_q, -1};

        return rt_extended;
    }
}


MergedRead merge_reads(const std::vector<InputRead> & inputs)
{
    if (inputs.size() == 1)
    {
        MergedRead _merged(inputs[0].read_seq, inputs[0].base_qualities, -1);
        return _merged;
    }
    
    std::deque<BaseCount> base_cnts;
    std::vector<Overlap> overlaps = find_overlaps(inputs);
    
    int merge_start = 0, merge_end = 0;
    int curr_start = 0, curr_end = 0;
    int prev_start = 0, prev_end = 0;
    int target_start = -1;

    bool is_last = false, is_target_start = false;
    int n_overlap_reads = overlaps.size();
    for (int i = 0; i < n_overlap_reads; ++i)
    {
        is_last = (i == (n_overlap_reads - 1));
        
        curr_start = overlaps[i].ref_start;
        curr_end = overlaps[i].ref_end;
        target_start = overlaps[i].target_start;
        
        std::string seq = inputs[overlaps[i].index].read_seq;
        std::string seq_qual = inputs[overlaps[i].index].base_qualities;
             
        if (!i)
        {
            for (int j = 0; j <= (curr_end - curr_start); ++j)
            {
                BaseCount bc;
                if (target_start == curr_start + j) is_target_start = 1;
                bc.add(seq[curr_start + j], seq_qual[curr_start + j], is_target_start);
                base_cnts.push_back(bc);
                is_target_start = 0;
            }
        }
        else  
        {
            if (prev_start <= curr_start)
            {
                merge_start += (curr_start - prev_start);
                merge_end = base_cnts.size(); 

                if (curr_start <= prev_end) 
                {
                    for (int j = 0; j <= (prev_end - curr_start); ++j)
                    {
                        if (merge_start + j < merge_end) 
                        {    
                            if (target_start == curr_start + j) is_target_start = 1;
                            base_cnts[merge_start + j].add(seq[curr_start + j], seq_qual[curr_start + j], is_target_start);
                            is_target_start = 0;
                        }
                    }
                    
                    // right extension (meaningful if prev_end < curr_end)
                    if (!is_last)
                    {
                        for (int j = 1; j <= (curr_end - prev_end); ++j)
                        {
                            BaseCount bc;
                            if (target_start == prev_end + j) is_target_start = 1;
                            bc.add(seq[prev_end + j], seq_qual[prev_end + j], is_target_start);
                            base_cnts.push_back(bc);
                            is_target_start = 0;
                        }
                    }
                  
                }
                else
                {
                    //construction fail 
                    //should be checked earlier? 
                }
            }
            else
            {
                if(!is_last)
                {
                    for (int j = -1; (curr_start - prev_start) <= j; --j)
                    {
                        BaseCount bc;
                        if (target_start == prev_start + j) is_target_start = 1;
                        bc.add(seq[prev_start + j], seq_qual[prev_start + j], is_target_start);
                        base_cnts.push_front(bc);
                        is_target_start = 0;
                    }
                    merge_start = prev_start - curr_start;
                }
                else merge_start = 0;
                
                merge_end = base_cnts.size();
                
                if (prev_start <= curr_end)
                {
                    for (int j = 0; j <= (curr_end - prev_start); ++j)
                    {
                        if (merge_start + j < merge_end) 
                        {    
                            if (target_start == prev_start + j) is_target_start = 1;
                            base_cnts[merge_start + j].add(seq[prev_start + j], seq_qual[prev_start + j], is_target_start);
                            is_target_start = 0;
                        }
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
    for (const auto & b : base_cnts)
    {
        merged_read += b.get_consensus().first;
        merged_qualities += b.get_consensus().second;
        n_start_checks.push_back(std::accumulate(b.start_checks.begin(), b.start_checks.end(), 0));
    }
    
    /*
    int p = 0;
    for (const auto & b : base_cnts)
    {   
        for (const auto q : b.start_checks)
        {
            std::cout << q << ",";
        }
        std::cout << "  " << p << " " <<std::endl;
        ++p;
    }
    */
     
    int target_start_idx = get_max_indices(n_start_checks)[0];   
    
    MergedRead mr(merged_read, merged_qualities, target_start_idx);   
    
    return mr; 
}


//complete contig
char match_to_contig
(
    const std::string & query, const bool is_dirty_query,
    const std::string & contig_seq, const std::string & ref_contig_seq, 
    const std::vector<std::pair<int, int>> & decomposed_contig, const bool is_dirty_contig,
    const int n_tandem_repeats, const std::string & repeat_unit, 
    const std::string & rv_repeat_unit, const bool is_complete_tandem_repeat,
    const std::pair<int, int> & repeat_boundary,
    const Filter & filter,
    const Aligner & aligner, 
    Alignment & alignment
)          
{
    //contig layout
    const int lt_end_idx = decomposed_contig[0].second - 1;
    const int rt_start_idx = decomposed_contig[2].first;
    const int contig_len = contig_seq.size();
    const bool is_del = (rt_start_idx == (lt_end_idx + 1));
    
    int32_t mask_len = strlen(query.c_str()) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;
            
    aligner.Align(query.c_str(), contig_seq.c_str(), contig_len, filter, &alignment, mask_len);
   
    const int ref_start = alignment.ref_begin;
    const int ref_end = alignment.ref_end;
    
    
    // still has gap
    std::string cigar_string = alignment.cigar_string;
    size_t found_ins = cigar_string.find("I");
    size_t found_del = cigar_string.find("D");
    if (found_ins != std::string::npos || found_del != std::string::npos) return 'F';
       
    // doesn't overlap the target region
    if (ref_end <= lt_end_idx || rt_start_idx <= ref_start) return 'F';
    
    // del needs to be cov
    if (is_del)
    {
        if ((ref_start <= lt_end_idx) && (rt_start_idx <= ref_end)) { /*pass*/ }
        else return 'F'; 
    }
    
    // checks for complete tandem repeat  
    int lt_rep = 0, rt_rep = 0;
    if (is_complete_tandem_repeat) 
    { 
        const int boundary_start = repeat_boundary.first;
        const int read_start = alignment.query_begin;

        // must be bounded
        if (boundary_start < ref_start || ref_end < repeat_boundary.second) return 'F';
        
        std::string lt_fragment = query.substr(0, (boundary_start + 1 - ref_start) + read_start);
        std::reverse(lt_fragment.begin(), lt_fragment.end());
        lt_rep = count_repeats(rv_repeat_unit, lt_fragment);

        std::string rt_fragment = query.substr((boundary_start + 1 - ref_start) + read_start);
        rt_rep = count_repeats(repeat_unit, rt_fragment);

        // repeat as expected
        if (n_tandem_repeats != (lt_rep + rt_rep)) return 'F';  
    } 
    
    //critical region
    // Extra one base (N) + boundary (B) + repeat (R): NBRRRRRRRBN  
    int lt_crit_bound = lt_end_idx - (lt_rep + 1)* repeat_unit.size();
    int rt_crit_bound = rt_start_idx + (rt_rep + 1)* repeat_unit.size();

    //inserted region
    int lt_ins_bound = -1;
    int rt_ins_bound = -1;
    if (!is_del)
    {
        lt_ins_bound = lt_end_idx; 
        rt_ins_bound = rt_start_idx;
    }
    
    char op;
    int op_len;
    int i = ref_start;
    int n_crit_mismatched_bases = 0;  
    int n_ins_matched_bases = 0;
    std::vector<std::pair<char, int>> cigar_vec = to_cigar_vector(cigar_string);
    for (const auto & c : cigar_vec) 
    {
         op = c.first;
         op_len = c.second;
         
         switch (op) 
         {
            case '=':
                //i += op_len;
                for (int j = 0; j <  op_len; ++j)
                {
                    if (lt_ins_bound < i && i <rt_ins_bound)
                    {
                        ++n_ins_matched_bases;
                    }
                    ++i;
                }

                break;
            case 'X':
                for (int j = 0; j < op_len; ++j)
                {
                    if (lt_crit_bound <= i && i <= rt_crit_bound)
                    {
                        ++n_crit_mismatched_bases;
                    }
                    //if (lt_ins_bound < i && i <rt_ins_bound)
                    //{
                    //    ++n_ins_mismatched_bases;
                    //}
                    ++i;
                }
                break;
            case 'S':
                if (lt_crit_bound <= i && i <= rt_crit_bound)
                {
                    ++n_crit_mismatched_bases;
                }
                //if (lt_ins_bound < i && i <rt_ins_bound)
                //{
                //    ++n_ins_mismatched_bases;           
                //}
                break;
            default:
                break;
         }        
            
    }
    
    if (is_complete_tandem_repeat && (n_crit_mismatched_bases)) return 'F';
        
    if (!is_del)
    {
        int ins_len = (rt_start_idx - lt_end_idx);
        double frac_covered_ins = n_ins_matched_bases / static_cast<double>(ins_len);
        
        //smale -> adjustable
        if (ins_len < 4)
        {
            if (frac_covered_ins < 1.0) return 'F';
        }
        else 
        {    //simscore
            if (frac_covered_ins < 0.66) return 'F'; 
        }  
    }
   
    // check for extenders
    if ((ref_start < lt_end_idx && rt_start_idx < ref_end) 
        &&
        (cigar_string.find('X') == std::string::npos)) 
    {
        if (ref_start == 0) 
        {
            //left extender
            return 'L';
        }
       
        if (ref_end == (contig_len - 1))
        {
            //right extender
            return 'R';
            //std::cout << query << " " << ref_contig_seq << std::endl;
        }
    }          
    
    return 'M';
}


/*
void sw::merge_reads(std::vector<std::string> & seqs, std::vector<std::string> & seq_quals)
{
    std::deque<BaseCnt> base_cnts;
    std::vector<Overlap> overlaps = find_overlaps(seqs);

    int merge_start = 0, merge_end = 0;

*/
/*
sw::Alignment::Alignment( uint16_t __alignment_score, int32_t  __ref_begin,
                          int32_t  __ref_end, int32_t  __query_begin, int32_t  __query_end,
                          const std::string & __cigar_string )
{
    alignment_score = __alignment_score;
    ref_begin = __ref_begin;
    ref_end = __ref_end;
    query_begin = __query_begin;
    query_end = __query_end;
    cigar_string = __cigar_string;
}


sw::Alignment sw::align( const std::string & ref,
                         const std::string & query,
                         const uint8_t & match_score,
                         const uint8_t & mismatch_penalty,
                         const uint8_t & gap_open_penalty,
                         const uint8_t & gap_extending_penalty )
{

    int32_t mask_len = strlen( query.c_str() ) / 2;
    mask_len = mask_len < 15 ? 15 : mask_len;

    StripedSmithWaterman::Alignment __alignment;

    StripedSmithWaterman::Aligner aligner( match_score,
                                           mismatch_penalty,
                                           gap_open_penalty,
                                           gap_extending_penalty );

    StripedSmithWaterman::Filter filter;

    aligner.Align( query.c_str(), ref.c_str(), ref.size(), filter, &__alignment,
                   mask_len );


    sw::Alignment alignment( __alignment.sw_score, __alignment.ref_begin,
                             __alignment.ref_end, __alignment.query_begin, __alignment.query_end,
                             __alignment.cigar_string );

    return alignment;
}


//@Function:
//      Parse CIGAR to a list. 60=2D8I32=1S -> [60=, 2D, 8I, 32=, 1S]
std::vector<std::string> decompose_cigar_string( const std::string &
        cigar_string )
{
    std::regex cigar_pattern( R"(\d+[MIDNSHPX=])" );

    std::regex_iterator<std::string::const_iterator> cigar_itr(
        cigar_string.begin(),
        cigar_string.end(), cigar_pattern );
    std::regex_iterator<std::string::const_iterator> end_of_itr;

    std::vector<std::string> cigarette;
    while ( cigar_itr != end_of_itr ) {
        cigarette.push_back( ( *cigar_itr ).str() );
        ++cigar_itr;
    }

    return cigarette;
}


inline std::vector<std::string> _decompose_cigar_string(const std::string & cigar_string)
{
    std::vector<std::string> cigarette;
    
    size_t pos = 0;
    size_t newpos = 0;
    const size_t len = cigar_string.size();
    
    while (pos < len) {
        newpos = cigar_string.find_first_of("MIDNSHPX=", pos) + 1;
        cigarette.push_back(cigar_string.substr(pos, newpos - pos));
        pos = newpos;
    }
    
    return cigarette; 
}


inline bool is_gap( std::string cigar_itr )
{
    return (
               ( cigar_itr.find( "I" ) != std::string::npos )
               ||
               ( cigar_itr.find( "D" ) != std::string::npos )
           );
}

//@Function:
//      Check if there exist consecutive gaps (I or D)
//      [60=, 2D, 8I, 32=, 1S] -> true
//      [60=, 6I, 2X, 32=, 1S] -> false
bool has_consecutive_gap(const std::vector<std::string> & cigarette)
{
    bool prev_is_gap = false;
    for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {
        if ( ( prev_is_gap ) && is_gap( ( *itr ) ) ) {
            return true;
        }
        prev_is_gap = is_gap( ( *itr ) );
    }
    return false;
}


//@Function:
//      Merge gaps into a signle gap of specified type
//      [4I, 2I] -> "6I"
std::string concat_gaps( const std::vector<std::string> & cigarette,
                         std::string  gap_type )
{
    std::string gaps;

    if ( !cigarette.empty() ) {
        uint16_t total_gap_len = 0;

        for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
                itr != cigarette.end(); ++itr ) {
            total_gap_len += std::stoi( ( *itr ).substr( 0, ( *itr ).length() - 1 ) );
        }

        gaps = std::to_string( total_gap_len ) + gap_type;
    }

    return gaps;
}


//@Function:
//      Merge consecutive gaps and make insertion come first
//      [4=, 2I, 2D, 1I, 3=, 3D, 1I, 2D, 4I, 4=]
//      -> [4=, 3I, 2D, 3=, 5I, 5D, 4=]
void edit_cigar( std::vector<std::string> & cigarette )
{
    std::vector<std::string> tmp, ins, del;

    bool prev_is_gap = false;
    for ( std::vector<std::string>::const_iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {
        if ( is_gap( ( *itr ) ) ) {
            if ( ( *itr ).find( "I" ) != std::string::npos ) {
                ins.push_back( *itr );
            }
            else {
                del.push_back( *itr );
            }
            prev_is_gap = true;
        }
        else {
            if ( prev_is_gap ) {
                std::string merged_ins = concat_gaps( ins, "I" );
                if ( !merged_ins.empty() ) {
                    tmp.push_back( merged_ins );
                }
                std::string merged_del = concat_gaps( del, "D" );
                if ( !merged_del.empty() ) {
                    tmp.push_back( merged_del );
                }
                tmp.push_back( *itr );
                ins.clear();
                del.clear();
            }
            else {
                tmp.push_back( *itr );
            }

            prev_is_gap = false;
        }
    }
    std::swap( cigarette, tmp );
}

char sw::match_to_contig(const std::string & query, const bool is_dirty_query,
                         const std::string & contig_seq, const std::string & ref_contig_seq, const std::vector<std::pair<int, int>> & decomposed_contig, const bool is_dirty_contig,
                         const int n_tandem_repeats, const std::string & repeat_unit, const std::string & rv_repeat_unit, const bool is_complete_tandem_repeat,
                         const std::pair<int, int> & repeat_boundary)

char sw::is_compatible(const std::string & contig,
                       const std::string & ref_contig,
                       const std::string & query,
                       const std::vector<std::pair<int, int>> & decomposed_contig,
                       const std::string & repeat_unit,
                       const std::string & reversed_repeat_unit,
                       const int expected_num_repeats,
                       const bool is_complete_tandem_repeat,
                       const std::pair<int, int> & boundary_indexes,
                       const bool is_dirty) 


{
    //contig layout
    const int lt_end_idx = decomposed_contig[0].second - 1;
    const int rt_start_idx = decomposed_contig[2].first;
      
    //gap-less aln
    const int match_score = 3;
    const int mismatch_penalty = 10;
    const int gap_open_penalty = 255;
    const int gap_extention_penalty = 255;
    
    //TODO -> do not creat alingment object each time. pass by referece.
    sw::Alignment alignment = sw::align(contig_seq, query, match_score, mismatch_penalty,
                                        gap_open_penalty, gap_extention_penalty);
    const int ref_start = alignment.ref_begin;
    const int ref_end = alignment.ref_end;
    
    // still has gap
    std::string cigar_string = alignment.cigar_string;
    size_t found_ins = cigar_string.find("I");
    size_t found_del = cigar_string.find("D");
    if ((found_ins != std::string::npos) || (found_del != std::string::npos)) return 'F';
       
    // doesn't overlap the target region
    if ((ref_end <= lt_end_idx) || (rt_start_idx <= ref_start)) return 'F';
   
    //TODO -> require covering both sides for deletions.
    //     for ins, check lt-covering, rt-coveirng or both
    
    
    // checks for complete tandem repeat  
    if (is_complete_tandem_repeat) 
    { 
        const int boundary_start = repeat_boundary.first;
        const int read_start = alignment.query_begin;

        // must be bounded
        if ((boundary_start < ref_start) || (ref_end < repeat_boundary.second)) return 'F';
        
        std::string lt_fragment = query.substr(0, (boundary_start + 1 - ref_start) + read_start);
        std::reverse(lt_fragment.begin(), lt_fragment.end());
        int lt_rep = count_repeats(rv_repeat_unit, lt_fragment);

        std::string rt_fragment = query.substr((boundary_start + 1 - ref_start) + read_start);
        int rt_rep = count_repeats(repeat_unit, rt_fragment);

        // must repeat exactly n times 
        if (n_tandem_repeats != (lt_rep + rt_rep)) return 'F';         
                
    }  
    
    // check with alignment against referece 
    sw::Alignment ref_alignment = sw::align(ref_contig_seq, query, match_score, mismatch_penalty, gap_open_penalty, gap_extention_penalty);
    
    if (alignment.alignment_score > ref_alignment.alignment_score) 
    {
        return 'T';
    }
    else
    {
        if (is_dirty_query || is_dirty_contig) 
        {   
            return 'U';
        } 
        else
        {
            return 'F';
        }
    }
    
    char op;
    int op_len;
    int i = ref_start;
    int j = 0;
    int n_mismatched_bases = 0;  
    std::vector<std::pair<char, int>> cigar_vec = to_cigar_vector(cigar_string);
    for (const auto & c : cigar_vec) 
    {
         op = c.first;
         op_len = c.second;
         
         switch (op) 
         {
            case '=':
                i += op_len;
                j = i;
                break;
            case 'X':
                while ((lt_end_idx <= j) && (j <= rt_start_idx) && (j <= (i + op_len))) 
                {
                    ++n_mismatched_bases;
                    ++j;
                }
                 
                i += op_len;     
                j = i;
                break;
            default:
                break;
         }        
            
    }
    
    int critical_aln_start = std::max(ref_start, lt_end_idx);
    int critical_aln_end = std::min(rt_start_idx, ref_end);
    int n_total_bases = critical_aln_end - critical_aln_start + 1;
    
    //mismatche rate in critical region
    double mitmatch_rate = n_mismatched_bases / static_cast<double>(n_total_bases); 
    if (mitmatch_rate > 0.333) {
        return 'F';
    }
    else {
        sw::Alignment ref_alignment = sw::align(ref_contig, query, match_score, mismatch_penalty,
                                                gap_open_penalty, gap_extention_penalty);
        
        if (alignment.alignment_score > ref_alignment.alignment_score) {
            return 'T';
        }
        else {
            if (is_dirty) return 'T';
            else return 'U';
        }
        
    }
    
}

std::vector<sw::ParsedVariant> sw::find_variants( const sw::Alignment &
        alignment,
        const std::string & ref,
        const std::string & query,
        const uint32_t & genomic_ref_start )
{

    std::vector<std::string> cigarette = _decompose_cigar_string(
            alignment.cigar_string );
    std::vector<sw::ParsedVariant> variants;

    if ( has_consecutive_gap( cigarette ) ) {
        edit_cigar( cigarette );
    }

    uint32_t genomic_pos = ( genomic_ref_start > 0 ) ? genomic_ref_start - 1 : 0;

    uint32_t ref_idx = alignment.ref_begin, query_idx = alignment.query_begin;

    for ( std::vector<std::string>::iterator itr = cigarette.begin();
            itr != cigarette.end(); ++itr ) {

        char operation = ( *itr ).back(); // cigar operation
        uint16_t op_len = std::stoi( ( *itr ).substr( 0,
                                     ( *itr ).length() - 1 ) ); // operation length

        if ( operation == 'I' ) {
            sw::ParsedVariant ins;

            ins.is_indel = true;
            ins.is_ins = true;
            ins.is_del = false;

            ins.lt_ref = ref.substr( 0, ref_idx );
            ins.lt_query = query.substr( 0, query_idx );
            ins.ins_seq = query.substr( query_idx, op_len );
            ins.variant_len = op_len;
            ins.rt_ref = ref.substr( ref_idx );
            ins.rt_query = query.substr( query_idx + op_len );

            ins.lt_clipped_segment = query.substr( 0, alignment.query_begin );
            ins.rt_clipped_segment = query.substr( alignment.query_end + 1,
                                                   query.length() - alignment.query_end );

            ins.genomic_pos = genomic_pos;

            variants.push_back( ins );

            query_idx += op_len;

        }
        else if ( operation == 'D' ) {
            sw::ParsedVariant del;

            del.is_indel = true;
            del.is_ins = false;
            del.is_del = true;

            del.lt_ref = ref.substr( 0, ref_idx );
            del.lt_query = query.substr( 0, query_idx );
            del.del_seq = ref.substr( ref_idx, op_len );
            del.variant_len = op_len;
            del.rt_ref = ref.substr( ref_idx + op_len );
            del.rt_query = query.substr( query_idx );

            del.lt_clipped_segment = query.substr( 0, alignment.query_begin );
            del.rt_clipped_segment = query.substr( alignment.query_end + 1,
                                                   query.length() - alignment.query_end );

            del.genomic_pos = genomic_pos;

            variants.push_back( del );

            ref_idx += op_len;
            genomic_pos += op_len;
        }
        else if ( operation == 'X' ) {
            sw::ParsedVariant smv;      // single or multi-nucleotide variant

            smv.is_indel = false;
            smv.is_ins = false;
            smv.is_del = false;

            smv.lt_ref = ref.substr( 0, ref_idx );
            smv.lt_query = query.substr( 0, query_idx );
            smv.ref_base = ref.substr( ref_idx, op_len );
            smv.alt_base = query.substr( query_idx, op_len );
            smv.variant_len = op_len;
            smv.rt_ref = ref.substr( ref_idx + op_len );
            smv.rt_query = query.substr( query_idx + op_len );

            smv.lt_clipped_segment = query.substr( 0, alignment.query_begin );
            smv.rt_clipped_segment = query.substr( alignment.query_end + 1,
                                                   query.length() - alignment.query_end );

            smv.genomic_pos = genomic_pos;

            variants.push_back( smv );

            ref_idx += op_len;
            query_idx += op_len;
            genomic_pos += op_len;
        }
        else if ( operation == 'S' ) {
            //pass for softclips
        }
        else {
            ref_idx += op_len;
            query_idx += op_len;
            genomic_pos += op_len;
        }

    }

    return variants;
}



namespace bases {
char bases[5] = {'A', 'C', 'G', 'T', 'N'};
}

static char bases[5] = {'A', 'C', 'G', 'T', 'N'};


//@function
//  return indices of max elems
//@return
//  {max_elem_indices}
template <typename T>
std::vector<int> get_max_indices(const T & arr)
{
    std::vector<int> indices;

    auto it = std::max_element(std::begin(arr), std::end(arr));
    while (it != std::end(arr)) {
        indices.push_back(std::distance(std::begin(arr), it));
        it = std::find(std::next(it), std::end(arr), *it);
    }
    return indices;
}

template<typename T>
double average(const std::vector<T>  & v)
{
    if (v.empty()) return 0.0;
    const double cnt = static_cast<double>(v.size());
    return std::accumulate(v.begin(), v.end(), 0.0) / cnt;
}


struct BaseCount {

    //double a = 0.0, c = 0.0, g = 0.0, t = 0.0, n = 0.0;
    std::vector<double> a, c, g, t, n;
    
    double a_cnt = 0.0001, c_cnt = 0.0001, g_cnt = 0.0001, t_cnt = 0.0001,
           n_cnt = 0.0001;

    BaseCount() {}

    BaseCount(char base, char qual)
    {
        double _qual = static_cast<double>(qual);
        
        switch (base) {
        case 'A':
            //a = _qual;
            a.push_back(_qual);
            ++a_cnt;
            break;
        case 'C':
            //c = _qual;
            c.push_back(_qual);
            ++c_cnt;
            break;
        case 'G':
            //g = _qual;
            g.push_back(_qual);
            ++g_cnt;
            break;
        case 'T':
            //t = _qual;
            t.push_back(_qual);
            ++t_cnt;
            break;
        default:
            _qual = 33;
            n.push_back(_qual);
            ++n_cnt;
            break;
        }
    }

    void add(char base, char qual)
    {
        int _qual = static_cast<int>(qual);
        
        switch (base) {
        case 'A':
            //a += _qual;
            a.push_back(_qual);
            ++a_cnt;
            break;
        case 'C':
            //c += _qual;
            c.push_back(_qual);
            ++c_cnt;
            break;
        case 'G':
            //g += _qual;
            g.push_back(_qual);
            ++g_cnt;
            break;
        case 'T':
            //t += _qual;
            t.push_back(_qual);
            ++t_cnt;
            break;
        default:
            _qual = 33;
            n.push_back(_qual);
            ++n_cnt;
            break;
        }
    }

    std::pair<char, char> get_consensus()
    {
        double counts[] = {a_cnt, c_cnt, g_cnt, t_cnt, n_cnt};
        //double quals[] = {a, c, g, t, n};
        double quals[] = {average(a), average(c), average(g), average(t), average(n)};


        std::vector<int> max_cnt_indices = get_max_indices(counts);

        int max_idx;
        if (max_cnt_indices.size() > 1) {
            std::vector<int> max_qual_indices = get_max_indices(quals);
            max_idx = max_qual_indices[0];
        }
        else {
            max_idx = max_cnt_indices[0];
        }

        //std::vector<char> ret = {bases[max_idx], static_cast<char>(quals[max_idx] / counts[max_idx])};
        //std::pair<char, char> ret {bases[max_idx], static_cast<char>(quals[max_idx] / counts[max_idx])};
        std::pair<char, char> ret {bases[max_idx], static_cast<char>(quals[max_idx])};  

        return ret;
    }
};


void update(std::deque<BaseCount> & consensus,  const int & read2_begin,
            const std::string & lt_ext, const std::string & lt_qual,
            const std::string & mread2, const std::string & mqual2,
            const std::string & rt_ext, const std::string & rt_qual)
{
    if (read2_begin == 0) {

        // update middle part
        const auto lt_len = lt_ext.size();
        for (size_t i = lt_len; (i < consensus.size() && (i - lt_len) < mread2.size() - 1); ++i) {
            consensus[i].add(mread2[i - lt_len], mqual2[i - lt_len]);
        }

        // rt extention
        for (size_t i = 0; i < rt_ext.size(); ++i) {
            consensus.emplace_back(BaseCount(rt_ext[i], rt_qual[i]));
        }
    }
    else {
        // update middle part
        for (size_t i = 0; i < mread2.size(); ++i) {
            consensus[i].add(mread2[i], mqual2[i]);
        }

        //lt extention
        for (int i = lt_ext.size() - 1; i >= 0; --i) {
            consensus.emplace_front(BaseCount(lt_ext[i], lt_qual[i]));
        }
    }
}


std::pair<std::string, std::string> get_consensus_contig(std::deque<BaseCount> & consensus)
{
    std::string consensus_seq;
    std::string consensus_qual;
    for (std::deque<BaseCount>::iterator itr = consensus.begin();
            itr != consensus.end(); ++itr) {
        std::pair<char, char> c = (*itr).get_consensus();
        consensus_seq += c.first;
        consensus_qual += c.second;
    }

    std::pair<std::string, std::string> ret {consensus_seq, consensus_qual};

    return ret;
}


void pairwise_merge(std::vector<BaseCount> & consensus,
                    const std::pair<std::string, std::string> & read)
{ 
    std::vector<std::string> c = get_consensus_contig(consensus);

    std::string seq1 = c[0];
    std::string qual1 = c[1];

    std::string seq2 = read.first;
    std::string qual2 = read.second;
    
    size_t seq2_len = seq2.size()
    
    // gap-less alignment
    const int match_score = 3;
    const int mismatch_penalty = 2;
    const int gap_open_penalty = std::max(seq1.size(), seq2_len);
    const int gap_extention_penalty = 1;

    sw::Alignment aln = sw::align(seq1, seq2, 
                                  match_score, mismatch_penalty, 
                                  gap_open_penalty, gap_extention_penalty);
    
    uint16_t aln_score = aln.alignment_score;
    uint16_t matched_segment_len = aln_score / match_score 
    double matched_proportion = ((double)matched_segment_len) / ((double)seq2_len);
    if (matched_proportion > 0.95) { 
        
       const int seq1_begin = aln.ref_begin;
       const int seq1_end = aln.ref_end;
       const int seq2_begin = aln.query_begin;
       const int seq2_end = aln.query_end;

       


    }
    else {
    }
}




void pairwise_stitch2(...)
{
    
    //thresh
    int n_exact_match_thresh = 10;

    
    std::string reference;
    std::string query;
    
    // gap-less alignment
    const int match_score = 3;
    const int mismatch_penalty = 2;
    const int gap_open_penalty = 255;
    const int gap_extention_penalty = 1;

    sw::Alignment aln = sw::align(reference, query, match_score, mismatch_penalty, gap_open_penalty, gap_extention_penalty);
    
    std::string cigar_string = aln.cigar_string;
    if ((cigar_string.find("I") != std::string::npos) || (cigar_string.find("D") != std::string::npos)) 
    {    
        return; //gapped -> not stitchable
    }

    int longest_exact_match_len = 0;
    std::vector<std::pair<char, int>> cigar_vec = to_cigar_vector(cigar_string);   
    for (const auto & c : cigar_vec) 
    {
        if (c.first == '=') 
        {    
            if (longest_exact_match_len < c.second) 
            {    
                longest_exact_match_len = c.second;
            }
        }      
    }    
    if (longest_exact_match_len < n_exact_match_thresh)
    {   
        return; // too short exact match -> not stitchable
    }
    
    const int r_begin = aln.ref_begin;
    const int r_end = aln.ref_end
    const int q_begin = aln.query_begin;
    const int q_end = aln.query_end;
    
    const bool is_r_start_aligned = (!r_begin);
    const bool is_q_start_aligned = (!q_begin);
    const bool is_r_end_aligned = (r_end == (reference.size() - 1));
    const bool is_q_end_aligned = (q_end == (query.size() - 1));

    size_t mid_len = q_end - q_begin + 1;
    
    std::string lt_ext, mid, rt_ext;
    if (is_r_start_aligned) 
    {
        lt_ext = query.substr(0, q_begin);
        //do    
    }
    else if (is_q_end_aligned)
    { 
        
    } 


    
    if (is_r_end_aligned)
    {
        rt_ext = query.substr(q_end + 1); 
        is_rt_extendable = true;
    }
}   


struct BaseCnt
{
    size_t a = 0, c = 0, g = 0, t = 0, n = 0;
    std::vector<int> aq, cq, gq, tq, nq;
    
    void add(char base, char _qual)
    {
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
        double quals[] = {average(aq), average(cq), average(gq), average(tq), average(nq)};

        std::vector<int> max_cnt_indices = get_max_indices(counts);

        int max_idx;
        if (max_cnt_indices.size() > 1)
        {
            std::vector<int> max_qual_indices = get_max_indices(quals);
            max_idx = max_qual_indices[0];
        }
        else max_idx = max_cnt_indices[0];
        
        return {bases[max_idx], static_cast<char>(quals[max_idx])};
    }
};

struct Overlap
{   
    int index;
    int ref_start;
    int ref_end;
    int query_start;
    int query_end;

    Overlap(const int index, 
            const int ref_start, 
            const int ref_end, 
            const int query_start, 
            const int query_end) : index(index), ref_start(ref_start), ref_end(ref_end), query_start(query_start), query_end(query_end) {}
};


std::vector<Overlap> find_overlaps(std::vector<std::string> & seqs)
{
     
    std::vector<std::string> ttt = {"GTGTATAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGAGAAGCAACATCTCCGAAAGCAAACAAGCCAACAAGGNNNN",
                                    "TTTTTTTTTTTTTTTTTTTTAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGAGAAGCAACATCTCCGAAAGCCAACAAGCCAACAAGGAAATCCTCGN",
                                    "CTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCGCGTCGCTATCAAGAGAAGCAACATCTCCGAAAGCCAACAAGCCAACAAGGAAATCCTCGATGAAGC",
                                    "AAAGTTAAAATTCCCGTCGCTATCAAGAGAAGCAACATCTCCGAAAGCCACCAAGCCAACAAGGAAATCCTCGATGAAGCCTACGTGATGGCCAGCGTGG",
                                    "AAAGTTAAAATTCCCGTCGCTATCAAGAGAAGCAACATCTCCGAAAGCCACCAAGCCAACAAGGAAATCCTCGATGAAGCCTACGTGATGGCCAGCGTGG"};

    
    
    std::vector<std::string> ttt = {"TACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCCTCCTTCCTCTGAAAATATGTCTTGTCACTTGTGACTTGAATGTAGACAGAAAACCTCCTAAAAACTATCTTTGCTTCTCTTCA",
                                    "GAGTGTGAGGCGTATTATACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCCTCCTTCCTCTGAACATATGTCTTGTCACTTGTGACTTGAATGTAGACAGAAAACCTCCTAAAAACT",
                                    "AGAATGAGTGTGAGGCGTATTATACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCCTCCTTCCTCTGAAAATATGTCTTGTCACTTGTGACTTGAATGTAGACAGAAAACCTCCTAA",
                                    "GATTGATCGAGGGTAAATGTGTCTTCAAGATTCTACAACAGATTCTCTCGTCAGGGGGTTGAGAATGAGTGTGAGGCGTATTATACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCC",
                                    "GATTGATCGAGGGTAAATGTGTCTTCAAGATTCTACAACAGATTCTCTCGTCAGGGGGTTGAGAATGAGTGTGAGGCGTATTATACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCC"};

    

     
    std::vector<std::string> ttt = { "NAANANNANANANANNCCNANAGAATGAGTGTGAGGCGTATTATACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCCTCCTTCCTCTGAAAATATGTCTTGTCACTTGTGACTTGAATGTAGACAGAAAACCTCCTAA",
                                    "NTTNTACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCCTCCTTCCTCTGAAAATATGTCTTGTCACTTGTGACTTGAATGTAGACAGAAAACCTCCTAAAAACTATCTTTGCTTCTCTTCA",
                                     "NGGNGNNTGATTGATCGAGGGTAAATGTGTCTTCAAGATTCTACAACAGATTCTCTCGTCAGGGGGTTGAGAATGAGTGTGAGGCGTATTATACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCC",
                                     "NCCGAGTGTGAGGCGTATTATACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCCTCCTTCCTCTGAAAATATGTCTTGTCACTTGTGACTTGAATGTAGACAGAAAACCTCCTAAAAACT", 
                                     "NCCGAGTGTGAGGCGTATTATACCATAGCCGCCTAGCTTATCAAGAGTTTATCTTTCTCTTTTCTTGCTACAAACCCAGTGCATTTCCTCCTTCCTCTGAAAATATGTCTTGTCACTTGTGACTTGAATGTAGACAGAAAACCTCCTAAAAACT"};

    
    seqs = ttt;
    
    
    std::vector<Overlap> overlaps;   
    
    const int match_score = 3;
    const int mismatch_penalty = 10;
    const int gap_open_penalty = 255;
    const int gap_extention_penalty = 255;

    StripedSmithWaterman::Alignment alignment;   
    StripedSmithWaterman::Aligner aligner(match_score, mismatch_penalty, gap_open_penalty, gap_extention_penalty);
    StripedSmithWaterman::Filter filter;
    
    int n = seqs.size();
    for (size_t i = 0; i < n - 1;)
    {
        int32_t mask_len = strlen(seqs[i+1].c_str() ) / 2;
        mask_len = mask_len < 15 ? 15 : mask_len;
        int j = i;
        int k = i;
        bool is_gapped = true;
        //flag for supporting alignment prev_start <= curr_start  prev_start <= curr_end
        do
        {
            ++j;
            aligner.Align(seqs[j].c_str(), seqs[i].c_str(), seqs[i].size(), filter, &alignment, mask_len);
            std::string cigar_string = alignment.cigar_string;
            if (cigar_string.find("I") == std::string::npos || cigar_string.find("D") == std::string::npos) is_gapped = false;
            std::cout << i << " " << j << " " << cigar_string << std::endl; 
        }   
        while ((j < n -1) && is_gapped);
        i = j;

        overlaps.emplace_back(k, alignment.ref_begin, alignment.ref_end, alignment.query_begin, alignment.query_end);
            
    }

    return overlaps;
}

void sw::merge_reads(std::vector<std::string> & seqs, std::vector<std::string> & seq_quals)
{
    std::deque<BaseCnt> base_cnts;
    std::vector<Overlap> overlaps = find_overlaps(seqs);

    int merge_start = 0, merge_end = 0;
    int curr_start = 0, curr_end = 0;
    int prev_start = 0, prev_end = 0;
    
    bool is_last = false;
    int n_overlap_reads = overlaps.size();
    for (int i = 0; i < n_overlap_reads; ++i)
    {
        is_last = (i == (n_overlap_reads - 1));
        
        curr_start = overlaps[i].ref_start;
        curr_end = overlaps[i].ref_end;
        
        std::string seq = seqs[overlaps[i].index];
        std::string seq_qual = seq_quals[overlaps[i].index];
             
        if (!i)
        {
            for (int j = 0; j <= (curr_end - curr_start); ++j)
            {
                BaseCnt bc;
                bc.add(seq[curr_start + j], seq_qual[curr_start + j]);
                base_cnts.push_back(bc);
            }
        }
        else  
        {
            if (prev_start <= curr_start)
            {
                merge_start += (curr_start - prev_start);
                merge_end = base_cnts.size(); 

                if (curr_start <= prev_end) 
                {
                    for (int j = 0; j <= (prev_end - curr_start); ++j)
                    {
                        if (merge_start + j < merge_end) 
                        {    
                            base_cnts[merge_start + j].add(seq[curr_start + j], seq_qual[curr_start + j]);
                        }
                    }
                    
                    // right extension (meaningful if prev_end < curr_end)
                    if (!is_last)
                    {
                        for (int j = 1; j <= (curr_end - prev_end); ++j)
                        {
                            BaseCnt bc;
                            bc.add(seq[prev_end + j], seq_qual[prev_end + j]);
                            base_cnts.push_back(bc);
                        }
                    }
                  
                }
                else
                {
                    //construction fail 
                    //should be checked earlier? 
                }
            }
            else
            {
                if(!is_last)
                {
                    for (int j = -1; (curr_start - prev_start) <= j; --j)
                    {
                        BaseCnt bc;
                        bc.add(seq[prev_start + j], seq_qual[prev_start + j]);
                        base_cnts.push_front(bc);
                    }
                    merge_start = prev_start - curr_start;
                }
                else merge_start = 0;
                
                merge_end = base_cnts.size();
                
                if (prev_start <= curr_end)
                {
                    for (int j = 0; j <= (curr_end - prev_start); ++j)
                    {
                        if (merge_start + j < merge_end) 
                        {    
                            base_cnts[merge_start + j].add(seq[prev_start + j], seq_qual[prev_start + j]);
                        }
                    }
                }    

            }
            
            std::cout << prev_start << " " << curr_start << " " << merge_start << std::endl;    
        }
        
        // query serves as ref in the next iteration
        // memorize the aln start/end for the next
        prev_start = overlaps[i].query_start;
        prev_end = overlaps[i].query_end;       
    }
    for (auto & i : base_cnts)
    {
        std::cout << i.get_consensus().first;
    }
    std::cout << std::endl;     
}

void pairwise_stitch(std::deque<BaseCount> & consensus,
//                     const std::vector<std::string> & v)
                     const std::pair<std::string, std::string> & v)   
{

    std::pair<std::string, std::string> c = get_consensus_contig(consensus);
    
    std::string read1 = c.first;
    std::string qual1 = c.second;

    std::string read2 = v.first;
    std::string qual2 = v.second;

    // gap-less alignment
    const int match_score = 3;
    const int mismatch_penalty = 2;
    const int gap_open_penalty = 255;
    const int gap_extention_penalty = 255;

    sw::Alignment aln = sw::align(read1, read2, match_score, mismatch_penalty,
                                  gap_open_penalty, gap_extention_penalty);

    
    std::vector<std::string> cigar_vec = _decompose_cigar_string(aln.cigar_string);
    
    const int read1_begin = aln.ref_begin;
    const int read1_end = aln.ref_end;
    const int read2_begin = aln.query_begin;
    const int read2_end = aln.query_end;
    
    const bool is_read1_end_covered = (read1_end == (read1.size() - 1));
     
    
    bool is_stitchable = true;
    std::string lt_ext, lt_qual, rt_ext, rt_qual, mread1, mqual1, mread2, mqual2;
    if (read2_begin == 0 && is_read1_end_covered) {
        lt_ext = read1.substr( 0, read1_begin );
        lt_qual = qual1.substr( 0, read1_begin );
        rt_ext = read2.substr( read2_end + 1 );
        rt_qual = qual2.substr( read2_end + 1 );

        mread2 = read2.substr( 0, read2_end + 1 );
        mqual2 = qual2.substr( 0, read2_end + 1 );

    }
    else if ( read1_begin == 0 ) {
        lt_ext = read2.substr( 0, read2_begin );
        lt_qual = qual2.substr( 0, read2_begin );
        rt_ext = read1.substr( read1_end + 1 );
        rt_qual = qual1.substr( read1_end + 1 );
        
        size_t m_len = read2_end - read2_begin + 1;
        mread2 = read2.substr( read2_begin, m_len );
        mqual2 = qual2.substr( read2_begin, m_len );
    }
    else {
        is_stitchable = false;
        // not stitchable -> reconsider
    }


    // check for stichability needed
    if (is_stitchable)
    {    
        update( consensus, read2_begin, lt_ext, lt_qual, mread2, mqual2, rt_ext,
                rt_qual );
    }
   
}

std::pair<std::string, std::string> sw::flatten_reads(const std::pair<std::string, std::string> seed_read,
                                                      const std::vector<std::pair<std::string, std::string>> & reads)
{
    std::deque<BaseCount>  consensus;
    std::string seed_bases = seed_read.first;
    std::string seed_quals = seed_read.second;

    for ( size_t i = 0; i < seed_bases.size(); ++i ) {
        consensus.emplace_back(BaseCount( seed_bases[i], seed_quals[i]));
    }


    for (std::vector<std::pair<std::string, std::string>>::const_iterator itr = reads.begin();
            itr != reads.end(); ++itr ) {
        pairwise_stitch( consensus, *itr );
    }


    std::string seq = "";
    std::string qual = "";
    for ( size_t i = 0; i < consensus.size(); ++i ) {
        BaseCount c = consensus[i];
        seq += c.get_consensus().first;
        qual += c.get_consensus().second;
    }

    return {seq, qual};
}

*/
