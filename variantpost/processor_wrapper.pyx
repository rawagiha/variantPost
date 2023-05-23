from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

#cdef extern from "util.h":
cdef extern from "pileup_processor.h" namespace "pp":
    
    cdef cppclass ProcessedPileup:
        
        ProcessedPileup() except +
        ProcessedPileup(vector[int] &,
                        vector[string] &,
                        vector[string] &,
                        vector[string] &, 
                        int,
                        string &,
                        string &,
                        vector[string] &,
                        vector[bool_t] &,
                        vector[bool_t] &,
                        vector[bool_t] &) except +
        
        vector[int] positions
        vector[string] ref_bases
        vector[string] alt_bases
        vector[string] base_quals
        int target_pos
        string ref, alt
        vector[string] read_names
        vector[bool_t] are_reverse                             
        vector[bool_t] are_target
        vector[bool_t] are_from_first_bam
    
    
    
    ProcessedPileup  _process_pileup(
            string &,
            string &,
            int,
            string &,
            string &,
            int,
            int,
            float,
            int,
            int,
            int,
            int,
            int, 
            int,
            int,
            #string &, 
            vector[string] &,
            vector[bool_t] &, 
            vector[string] &,
            vector[int] &, 
            vector[int] &,
            vector[string] &,
            #vector[string] &,           
            vector[vector[int]] &, 
            vector[int] &,
            vector[bool_t] &
    )


cdef object process_pileup(
    string & fastafile,
    string & chrom,
    int pos, 
    string & ref, 
    string & alt,
    int mapping_quality_threshold,
    int base_quality_threshold,
    float low_quality_base_rate_threshold,
    int match_score,
    int mismatch_penalty,
    int gap_open_penalty,
    int gap_extention_penalty,
    int kmer_size,
    int unspliced_local_reference_start,
    int unspliced_local_reference_end,
    #string & unspliced_local_reference,
    vector[string] & read_names,
    vector[bool_t] & are_reverse,
    vector[string] & cigar_strings,
    vector[int] & aln_starts,
    vector[int] & aln_ends,
    vector[string] & read_seqs,
    #vector[string] & ref_seqs,
    vector[vector[int]] & quals,
    vector[int] & mapqs,
    vector[bool_t] & are_first_bam,
):

    res = _process_pileup(
        fastafile,
        chrom,
        pos, 
        ref, 
        alt,
        mapping_quality_threshold,
        base_quality_threshold,
        low_quality_base_rate_threshold,
        match_score,
        mismatch_penalty,
        gap_open_penalty,
        gap_extention_penalty,
        kmer_size,
        unspliced_local_reference_start,
        unspliced_local_reference_end,
        #unspliced_local_reference,
        read_names,
        are_reverse,
        cigar_strings,
        aln_starts,
        aln_ends,
        read_seqs,
        #ref_seqs,
        quals,
        mapqs,
        are_first_bam
    )
    
    #res is not a python obj..
    
    #print(res.positions, res.ref_bases, res.target_pos, res.ref, res.alt) 
    return res.positions, res.ref_bases, res.target_pos, res.ref, res.alt
