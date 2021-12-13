from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

cdef extern from "pileup_parser.h" namespace "pileup":
    void  parse_pileup(
            int, 
            string &, 
            vector[string] &,
            vector[bool_t] &, 
            vector[string] &,
            vector[int] &, 
            vector[int] &,
            vector[string] &,
            vector[string] &,           
            vector[vector[int]] &, 
            vector[int] & 
    )


def test_it(
    unspliced_local_reference_start,
    unspliced_local_reference,
    read_names,
    are_reverse,
    cigar_strings,
    aln_starts,
    aln_ends,
    read_seqs,
    ref_seqs,
    quals,
    mapqs,
):

    parse_pileup(
        unspliced_local_reference_start,
        unspliced_local_reference,
        read_names,
        are_reverse,
        cigar_strings,
        aln_starts,
        aln_ends,
        read_seqs,
        ref_seqs,
        quals,
        mapqs,
    )
