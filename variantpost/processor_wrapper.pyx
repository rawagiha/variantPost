from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

cdef extern from "processor.h" namespace "pp":
    string  process_pileup(
            string &,
            string &,
            int,
            string &,
            string &,
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


cdef string test_it(
    string & fastafile,
    string & chrom,
    int pos, 
    string & ref, 
    string & alt,
    int base_quality_threshold,
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

    res = process_pileup(
        fastafile,
        chrom,
        pos, 
        ref, 
        alt,
        base_quality_threshold,
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
    
    return res
