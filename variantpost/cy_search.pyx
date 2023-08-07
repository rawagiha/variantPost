import random
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t
from collections import OrderedDict

cdef extern from "search.h":

    cdef cppclass SearchResult:

        SearchResult() except +
        SearchResult(
            vector[int]&,
            vector[string]&,
            vector[string]&,
            vector[string]&,
            vector[int]&,
            vector[int]&,
            int,
            string&,
            string&,
            vector[string]&,
            vector[bool_t]&,
            vector[int]&,
            vector[bool_t]&
        ) except +

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
        bool_t is_retargeted


    void _search_target(
        SearchResult&,
        string&,
        string&,
        int,
        string&,
        string&,
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
        int,
        vector[string]&,
        vector[bool_t]&,
        vector[string]&,
        vector[int]&,
        vector[int]&,
        vector[string]&,
        vector[vector[int]]&,
        vector[int]&,
        vector[bool_t]&
    )


cdef inline bint is_qualified_read(read, bint exclude_duplicates):
    
    if exclude_duplicates:
        if read.cigarstring and (not read.is_duplicate) and  (not read.is_secondary) and (not read.is_supplementary) and read.reference_end:
            return True
    else:
        if read.cigarstring and (not read.is_secondary) and (not read.is_supplementary) and read.reference_end:
            return True
    
    return False


def fetch_reads(bam, chrom, pos, chrom_len, window, exclude_duplicates, fetched_reads, est_cov, is_secondary):
    reads = bam.fetch(
        chrom, max(0, pos - window), min(pos + window, chrom_len), until_eof=False
    )
    
    for read in reads:
        if is_qualified_read(read, exclude_duplicates):
            fetched_reads.append((read, is_secondary))
            
            if read.reference_start <= pos <= read.reference_end:
                est_cov += 1 
    
    return est_cov


cdef inline void pack_to_lists(
    tuple read_tuple, 
    vector[string]& read_names, 
    vector[bool_t]& are_reverse, 
    vector[string]& cigar_strings, 
    vector[int]& aln_starts, 
    vector[int]& aln_ends, 
    vector[string]& read_seqs, 
    vector[string]& ref_seqs, 
    vector[vector[int]]& qual_seqs, 
    vector[int]& mapqs, 
    vector[bool_t]& is_primary, 
):
    cdef object read = read_tuple[0]
    read_names.push_back(read.query_name.encode())
    are_reverse.push_back(read.is_reverse)
    cigar_strings.push_back(read.cigarstring.encode())
    aln_starts.push_back(read.reference_start + 1)
    aln_ends.push_back(read.reference_end)
    read_seqs.push_back(read.query_sequence.encode())
    qual_seqs.push_back(read.query_qualities)
    mapqs.push_back(read.mapping_quality)
    
    if read_tuple[1]:
        is_primary.push_back(False)
    else:
        is_primary.push_back(True)


class AnnotatedRead(object):
    def __init__(self, read_name, is_reverse, is_first_bam, target_status):
        self.read_name = read_name
        self.is_reverse = is_reverse
        self.is_first_bam = is_first_bam
        self.target_status = target_status


cdef object search_target(
     object bam,
     object second_bam,
     int chrom_len,
     bint exclude_duplicates,
     int window,
     string  fastafile,
     str  chrom,
     int pos,
     string  ref,
     string  alt,
     int mapping_quality_threshold,
     int base_quality_threshold,
     float low_quality_base_rate_threshold,
     int downsample_thresh,
     int match_score,
     int mismatch_penalty,
     int gap_open_penalty,
     int gap_extention_penalty,
     int kmer_size,
     int local_threshold,
     int retarget_threshold,
     int unspliced_local_reference_start,
     int unspliced_local_reference_end,
):
    cdef SearchResult rslt     
    cdef int est_cov = 0
    cdef list fetched_reads = []
    
    # read_fetching from first bam
    est_cov = fetch_reads(
        bam, 
        chrom, 
        pos, 
        chrom_len, 
        window, 
        exclude_duplicates, 
        fetched_reads, 
        est_cov, 
        False
    )  
    
    if second_bam:
        est_cov = fetch_reads(
            second_bam, 
            chrom, 
            pos, 
            chrom_len, 
            window, 
            exclude_duplicates, 
            fetched_reads, 
            est_cov, 
            True
        )  
    
    if est_cov > downsample_thresh:
        n_sample = int(len(fetched_reads) * (downsample_thresh / est_cov))
        random.seed(123)
        fetched_reads = random.sample(fetched_reads,  n_sample)
    
    cdef int buff_size = len(fetched_reads);

    cdef vector[string] read_names 
    read_names.reserve(buff_size)
    cdef vector[bool_t] are_reverse 
    are_reverse.reserve(buff_size)
    cdef vector[string] cigar_strings 
    cigar_strings.reserve(buff_size)
    cdef vector[int] aln_starts 
    aln_starts.reserve(buff_size)
    cdef vector[int] aln_ends 
    aln_ends.reserve(buff_size)
    cdef vector[string] read_seqs
    read_seqs.reserve(buff_size) 
    cdef vector[string] ref_seqs
    ref_seqs.reserve(buff_size) 
    cdef vector[vector[int]] qual_seqs
    qual_seqs.reserve(buff_size) 
    cdef vector[int] mapqs 
    mapqs.reserve(buff_size)
    cdef vector[bool_t] are_first_bam
    are_first_bam.reserve(buff_size)
    
    for read in fetched_reads:
        pack_to_lists(
            read, 
            read_names, 
            are_reverse, 
            cigar_strings,
            aln_starts, 
            aln_ends, 
            read_seqs, 
            ref_seqs, 
            qual_seqs, 
            mapqs, 
            are_first_bam
        )
    
    _search_target(
        rslt,
        fastafile,
        chrom.encode(),
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
        retarget_threshold,
        unspliced_local_reference_start,
        unspliced_local_reference_end,
        read_names,
        are_reverse,
        cigar_strings,
        aln_starts,
        aln_ends,
        read_seqs,
        qual_seqs,
        mapqs,
        are_first_bam,
    )
    
    contig_dict = OrderedDict()
    for pos, ref_base, alt_base, base_qual in zip(
        rslt.positions, rslt.ref_bases, rslt.alt_bases, rslt.base_quals
    ):
        contig_dict[pos] = (ref_base.decode("utf-8"), alt_base.decode("utf-8"), base_qual.decode("utf-8"))
    
    annot_reads = []
    for read_name, is_reverse, target_status, is_first_bam in zip(
        rslt.read_names, rslt.are_reverse, rslt.target_statuses, rslt.are_from_first_bam
    ):
        annot_reads.append(AnnotatedRead(read_name.decode("utf-8"), is_reverse, is_first_bam, target_status))
    
    skips = [(start, end) for start, end in zip(rslt.skip_starts, rslt.skip_ends)]
    
    return contig_dict, skips, rslt.read_names,  rslt.are_reverse, rslt.target_statuses, rslt.are_from_first_bam, rslt.is_retargeted
