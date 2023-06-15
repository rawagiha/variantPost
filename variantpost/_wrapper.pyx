from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

from collections import OrderedDict

cdef extern from "search.h":

    cdef cppclass SearchResult:

        SearchResult() except +
        SearchResult(vector[int] &,
                        vector[string] &,
                        vector[string] &,
                        vector[string] &,
                        vector[int] &,
                        vector[int] &,
                        int,
                        string &,
                        string &,
                        vector[string] &,
                        vector[bool_t] &,
                        vector[int] &,
                        vector[bool_t] &) except +

        vector[int] positions
        vector[string] ref_bases
        vector[string] alt_bases
        vector[string] base_quals
        vector[int] skip_starts
        vector[int] skip_ends
        int target_pos
        string ref, alt
        vector[string] read_names
        vector[bool_t] are_reverse
        vector[int] target_statuses
        vector[bool_t] are_from_first_bam



    SearchResult _search_target(
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
            int,
            vector[string] &,
            vector[bool_t] &,
            vector[string] &,
            vector[int] &,
            vector[int] &,
            vector[string] &,
            vector[vector[int]] &,
            vector[int] &,
            vector[bool_t] &
    )


class AnnotatedRead(object):
    def __init__(self, read_name, is_reverse, is_first_bam, target_status):
        self.read_name = read_name
        self.is_reverse = is_reverse
        self.is_first_bam = is_first_bam
        self.target_status = target_status


cdef object search_target(
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
     int local_threshold,
     int unspliced_local_reference_start,
     int unspliced_local_reference_end,
     vector[string] & read_names,
     vector[bool_t] & are_reverse,
     vector[string] & cigar_strings,
     vector[int] & aln_starts,
     vector[int] & aln_ends,
     vector[string] & read_seqs,
     vector[vector[int]] & quals,
     vector[int] & mapqs,
     vector[bool_t] & are_first_bam,
):

    res = _search_target(
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
        local_threshold,
        unspliced_local_reference_start,
        unspliced_local_reference_end,
        read_names,
        are_reverse,
        cigar_strings,
        aln_starts,
        aln_ends,
        read_seqs,
        quals,
        mapqs,
        are_first_bam,
    )
    
    contig_dict = OrderedDict()
    for pos, ref_base, alt_base, base_qual in zip(
        res.positions, res.ref_bases, res.alt_bases, res.base_quals
    ):
        contig_dict[pos] = (ref_base.decode("utf-8"), alt_base.decode("utf-8"), base_qual)

    annot_reads = []
    for read_name, is_reverse, target_status, is_first_bam in zip(
        res.read_names, res.are_reverse, res.target_statuses, res.are_from_first_bam
    ):
        annot_reads.append(AnnotatedRead(read_name, is_reverse, is_first_bam, target_status))
    
    skips = [(start, end) for start, end in zip(res.skip_starts, res.skip_ends)]
    
    return contig_dict, skips, annot_reads
