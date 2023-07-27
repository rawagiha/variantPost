#import time
import random

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

from collections import OrderedDict
#from pysam.libcalignmentfile cimport AlignmentFile
#from pysam.libcalignedsegment cimport AlignedSegment

RAW_DEPTH = 0

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
        bool_t is_retargeted


    void _search_target(
            SearchResult &,
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



cdef inline bint is_qualified_read(read, bint exclude_duplicates):
    
    if exclude_duplicates:
        if read.cigarstring and (not read.is_duplicate) and  (not read.is_secondary) and (not read.is_supplementary) and read.reference_end:
            return True
    else:
        if read.cigarstring and (not read.is_secondary) and (not read.is_supplementary) and read.reference_end:
            return True
    
    return False


cdef object fetch_reads(bam, str chrom, int pos, int chrom_len, int window, bint exclude_duplicates):
    out = []
    reads = bam.fetch(
        chrom, max(0, pos - window), min(pos + window, chrom_len), until_eof=False
    )
    
    return [read for read in reads if is_qualified_read(read, exclude_duplicates)]


def downsampler(chrom, pos, bam, downsample_thresh, reads):
    """Downsample reads if depth exceeds downsample_thresh
    """
    depth = bam.count(chrom, pos - 1, pos)

    if depth > downsample_thresh:
        pileup_size = len(reads)
        random.seed(123)
        reads = random.sample(reads, int(pileup_size * (downsample_thresh / depth)))
        sample_factor = pileup_size / len(reads)
    else:
        sample_factor = 1.0

    return reads, sample_factor



cdef inline void pack_to_lists(
    object read, 
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
    bint is_secondary
):
    read_names.push_back(read.query_name.encode())
    are_reverse.push_back(read.is_reverse)
    cigar_strings.push_back(read.cigarstring.encode())
    aln_starts.push_back(read.reference_start + 1)
    aln_ends.push_back(read.reference_end)
    read_seqs.push_back(read.query_sequence.encode())
    qual_seqs.push_back(read.query_qualities)
    mapqs.push_back(read.mapping_quality)
    
    if is_secondary:
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
     int downsample_threshold,
     string  fastafile,
     str  chrom,
     int pos,
     string  ref,
     string  alt,
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
):

    cdef SearchResult rslt
    #cdeobject AlignedSegment read
    #tt = time.time()
        
    cdef int buff_size = 0;
    first_reads = fetch_reads(bam, chrom, pos, chrom_len, window, exclude_duplicates)
    buff_size += len(first_reads)
    
    if second_bam:
        second_reads = fetch_reads(second_bam, chrom, pos, chrom_len, window, exclude_duplicates)
        buff_size += len(second_reads)

    #print("I/O by pysam", time.time() - tt)

    #tt = time.time()
    #if downsample_threshold < 0:
    #    sample_factor = 1.0
    #else:
    #    reads, sample_factor = downsampler(chrom, pos, bam, downsample_threshold, reads)

    #cdef int n = len(reads)
    
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

    #tt = time.time()
    
    for read in first_reads:
        pack_to_lists(read, read_names, are_reverse, cigar_strings,
                     aln_starts, aln_ends, read_seqs, ref_seqs, qual_seqs, mapqs, are_first_bam, False)
    
    if second_bam:
        for _read in second_reads:
            pack_to_lists(_read, read_names, are_reverse, cigar_strings,
                          aln_starts, aln_ends, read_seqs, ref_seqs, qual_seqs, mapqs, are_first_bam, True)
    
    #print("prep--", time.time() - tt)
    
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
    #contig_dict = list()
    for pos, ref_base, alt_base, base_qual in zip(
        rslt.positions, rslt.ref_bases, rslt.alt_bases, rslt.base_quals
    ):
        contig_dict[pos] = (ref_base.decode("utf-8"), alt_base.decode("utf-8"), base_qual.decode("utf-8"))
        #contig_dict.append([pos, ref_base.decode("utf-8"), alt_base.decode("utf-8"), base_qual.decode("utf-8")])       
    
    annot_reads = []
    for read_name, is_reverse, target_status, is_first_bam in zip(
        rslt.read_names, rslt.are_reverse, rslt.target_statuses, rslt.are_from_first_bam
    ):
        annot_reads.append(AnnotatedRead(read_name.decode("utf-8"), is_reverse, is_first_bam, target_status))
    
    skips = [(start, end) for start, end in zip(rslt.skip_starts, rslt.skip_ends)]
    
    return contig_dict, skips, rslt.read_names,  rslt.are_reverse, rslt.target_statuses, rslt.are_from_first_bam, rslt.is_retargeted
    #return contig_dict, skips, annot_reads
