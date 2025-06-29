#include "pileup.h"

Pileup::Pileup(const std::vector<std::string>& read_names,
               const std::vector<bool>& are_reverse, 
               const std::vector<std::string>& cigar_strs,
               const std::vector<int>& aln_starts,
               const std::vector<int>& aln_ends,
               const std::vector<std::string>& read_seqs,
               const std::vector<std::vector<int>>& quals,
               const std::vector<int>& mapqs,
               const std::vector<bool>& are_from_first_bam,
               const UserParams& params,
               LocalReference& loc_ref,
               Variant& target)
{
    const int qc_start = target.lpos - params.local_thresh;
    const int qc_end = target.end_pos + params.local_thresh;
    
    size_t n_reads = read_names.size(); reads.reserve(n_reads);
    for (size_t i = 0; i < n_reads; ++i)
    {
        if (mapqs[i] < params.mapq_thresh) continue; //can do this at Python?
        
        reads.emplace_back(read_names[i], are_reverse[i], cigar_strs[i],
                           aln_starts[i], aln_ends[i], read_seqs[i],
                           quals[i], mapqs[i], are_from_first_bam[i]);

        auto& read = reads[i];
        read.setReference(loc_ref); 
        
        if (read.is_na_ref) continue;
        read.setVariants(loc_ref); 

        read.parseCoveringPattern(loc_ref, target);
        if (read.covering_ptrn == 'C') continue;
        
        read.parseLocalPattern(loc_ref, target);
        read.qualityCheck(qc_start, qc_end, params); 
        
        if (loc_ref.has_flankings)
            read.isStableNonReferenceAlignment(loc_ref); 
        
        // search complex target by string match
        if (target.is_complex         && 
            read.covering_ptrn == 'A' && read.is_stable_non_ref)
                read.hasTargetComplexVariant(loc_ref, target);
        // ranking
        if (read.has_target)
        {
            read.setSignatureStrings(); read.rank = 's';
        }
        else if (read.has_local_events)
        {
            read.setSignatureStrings(); read.rank = 'u';
        }
        else
            read.rank = 'n'; // here local ref only?      
    }
}

void Pileup::setLocalHaploTypes(LocalReference& loc_ref)
{
    
} 
