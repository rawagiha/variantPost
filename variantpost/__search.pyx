# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: cdivision=True
# cython: nonecheck=False

import random
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t
from collections import OrderedDict

from .long_reads import shorten_read

cdef extern from "search.h":

    cdef cppclass SearchResult:

        SearchResult() except +

        vector[int] positions
        vector[string] ref_bases
        vector[string] alt_bases
        #vector[string] base_quals
        #vector[int] skip_starts
        #vector[int] skip_ends
        int retarget_pos
        string ref, alt
        #vector[string] read_names
        #vector[bool_t] are_reverse
        vector[int] target_statuses
        #vector[bool_t] are_from_first_bam
        #vector[string] trans_vars
        #bool_t is_retargeted

    void _search_target(
        SearchResult&,
        string&,
        string&,
        int,
        string&,
        string&,
        #int,
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
        vector[string]&,
        #vector[int]&,
        vector[bool_t]&,
        bool_t
    )

# utilities
cdef inline int int_max(int a, int b) nogil: return a if a > b else b
cdef inline int int_min(int a, int b) nogil: return a if a < b else b

cdef inline bint is_qualified_read(read, bint exclude_duplicates):
    cdef int flag = read.flag
    
    # 0x100: secondary, 0x800: supplementary, 0x4: unmapped
    # 0x4: unmapped -> this tests if read has a cigerstring
    # 0x100 | 0x800 | 0x4 -> 0x904
    if (flag & 0x904):
        return False
    # 0x400: duplicate
    if exclude_duplicates and (flag & 0x400):
        return False
    if read.reference_start == -1:
        return False
    return True


def fetch_reads(
    object bam, 
    str chrom, 
    int pos, 
    int chrom_len, 
    int window, 
    bint exclude_duplicates, 
    list fetched_reads, 
    int est_cov, 
    bint is_secondary
):
    cdef int start = int_max(0, pos - window)
    cdef int end = int_min(pos + window, chrom_len)
    
    reads = bam.fetch(chrom, start, end, until_eof=False)
    
    cdef object read
    cdef int ref_start, ref_end

    for read in reads:
        if is_qualified_read(read, exclude_duplicates):        
            fetched_reads.append((read, is_secondary))
            
            ref_start = read.reference_start
            ref_end = read.reference_end
            if ref_start <= pos <= ref_end:
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
    vector[string]& qual_seqs, 
    #vector[int]& mapqs, 
    vector[bool_t]& is_primary,
    #vector[string]& cb,             #for isoanalysis
    int unspliced_local_reference_start,
    int unspliced_local_reference_end,
    int window_len,
    int k,
    list tags,
    bint get_tag_flag
):
    cdef object read = read_tuple[0]
    cdef bint secondary_flag = read_tuple[1]
    cdef str q_seq = read.query_sequence
    cdef str q_qual = read.query_qualities_str

    if len(q_seq) < window_len:
        read_names.push_back(<string>read.query_name.encode('ascii'))
        are_reverse.push_back(<bool_t>read.is_reverse)
        cigar_strings.push_back(<string>read.cigarstring.encode('ascii'))
        aln_starts.push_back(<int>read.reference_start + 1)
        aln_ends.push_back(<int>read.reference_end)
        read_seqs.push_back(<string>q_seq.encode('ascii'))

        if not q_qual:
            qual_seqs.push_back((b'F' * len(q_seq)))
        else:
            qual_seqs.push_back(<string>q_qual.encode('ascii'))
    
        is_primary.push_back(secondary_flag)
        if get_tag_flag:
            tags.append(read.get_tags())
    else:
        #TODO This is basically long-read support. Work on later
        res = shorten_read(
            read, 
            unspliced_local_reference_start, 
            unspliced_local_reference_end,
            k
        )
       
        if res:
            read_names.push_back(read.query_name.encode('ascii'))
            are_reverse.push_back(read.is_reverse)
            cigar_strings.push_back(res[2].encode('ascii'))
            aln_starts.push_back(res[3])
            aln_ends.push_back(res[4])
            read_seqs.push_back(res[0].encode('ascii'))
            qual_seqs.push_back(res[1].encode('ascii'))
            is_primary.push_back(secondary_flag)
            if get_tag_flag:
                tags.append(read.get_tags())


cpdef object search_target(
     object bam,
     object second_bam,
     int chrom_len,
     bint exclude_duplicates,
     int window,
     string  fastafile,
     str chrom,
     str bam_chrom, 
     int pos,
     string ref,
     string alt,
     int mapping_quality_threshold,
     int base_quality_threshold,
     float low_quality_base_rate_threshold,
     int downsample_thresh,
     int match_score,
     int mismatch_penalty,
     int gap_open_penalty,
     int gap_extention_penalty,
     int kmer_size,
     int dimer_window,
     int local_threshold,
     int unspliced_local_reference_start,
     int unspliced_local_reference_end,
     int k,
     bint get_tags_flag=True # user option add later
):
    cdef SearchResult rslt     
    cdef int est_cov = 0
    cdef list fetched_reads = []
    cdef bool_t has_second = False

    # read_fetching from first bam
    est_cov = fetch_reads(bam, bam_chrom, pos, chrom_len, window, 
                          exclude_duplicates, fetched_reads, est_cov, 
                          False if second_bam else True)  
    
    if second_bam:
        has_second = True
        est_cov = fetch_reads(second_bam, bam_chrom, pos, chrom_len, window, 
                              exclude_duplicates, fetched_reads, est_cov, True)  
    
    if downsample_thresh > 0 and est_cov > downsample_thresh:
        n_sample = int(len(fetched_reads) * (downsample_thresh / est_cov))
        random.seed(123)
        fetched_reads = random.sample(fetched_reads, n_sample)
    
    cdef int buff_size = len(fetched_reads);
    cdef vector[string] read_names, cigar_strings, read_seqs, qual_seqs 
    cdef vector[bool_t] are_reverse, are_first_bam
    cdef vector[int] aln_starts, aln_ends

    read_names.reserve(buff_size)
    are_reverse.reserve(buff_size)
    cigar_strings.reserve(buff_size)
    aln_starts.reserve(buff_size)
    aln_ends.reserve(buff_size)
    read_seqs.reserve(buff_size)
    qual_seqs.reserve(buff_size)
    are_first_bam.reserve(buff_size)
   
    cdef list tags = []
    cdef int widow_len = unspliced_local_reference_end - unspliced_local_reference_start
    cdef tuple r_tup
    cdef object r_obj

    for r_tup in fetched_reads:
        r_obj = r_tup[0]
        if r_obj.mapping_quality < mapping_quality_threshold:
            continue
        
        pack_to_lists(r_tup, read_names, are_reverse, cigar_strings,
                      aln_starts, aln_ends, read_seqs, qual_seqs, 
                      are_first_bam, unspliced_local_reference_start,
                      unspliced_local_reference_end, widow_len, k, tags, get_tags_flag)
         
    _search_target(
        rslt,
        fastafile,
        chrom.encode('ascii'),
        pos,
        ref,
        alt,
        #mapping_quality_threshold,  #check if this is still meaningful
        base_quality_threshold,
        low_quality_base_rate_threshold,
        match_score,
        mismatch_penalty,
        gap_open_penalty,
        gap_extention_penalty,
        kmer_size,
        dimer_window,
        local_threshold,
        #retarget_threshold,
        unspliced_local_reference_start,
        unspliced_local_reference_end,
        read_names,
        are_reverse,
        cigar_strings,
        aln_starts,
        aln_ends,
        read_seqs,
        qual_seqs,
        #mapqs,
        are_first_bam,
        has_second,
    )
    
    contig_dict = OrderedDict()
    cdef int p
    cdef string refb, altb
    for i in range(rslt.positions.size()):
        p = rslt.positions[i]
        refb = rslt.ref_bases[i]
        altb = rslt.alt_bases[i]
        contig_dict[p] = (refb.decode('ascii'), altb.decode('ascii'), "F")
    
    #TODO hard code "F", is this okay?? What if this base in actual data is low qual
    #here is prev version: keep for now
        #qual_chars = base_qual.decode("utf-8")
        #qual = -1
        #if len(qual_chars) == 1:
        #    qual = ord(qual_chars) - 33
        #else:
        #    quals = [ord(c) - 33 for c in qual_chars]
        #    qual = statistics.median(quals)
        
    
    # KEEP this for later dev
    #annot_reads = []
    #for read_name, is_reverse, target_status, is_first_bam in zip(
    #    rslt.read_names, rslt.are_reverse, rslt.target_statuses, rslt.are_from_first_bam
    #):
    #    annot_reads.append(
    #        AnnotatedRead(read_name.decode("utf-8"), is_reverse, is_first_bam, target_status)
    #    )
    
    #skips = [(start, end) for start, end in zip(rslt.skip_starts, rslt.skip_ends)]
    
    return (
        contig_dict,
        rslt.target_statuses,
        are_reverse,
        are_first_bam,
        tags,
        rslt.retarget_pos,
        rslt.ref.decode("ascii"),
        rslt.alt.decode("ascii")
     #   cb
    )
    #return (
    #    contig_dict, 
    #    skips, 
    #    rslt.read_names, 
    #    rslt.are_reverse, 
    #    rslt.target_statuses, 
    #    rslt.are_from_first_bam, 
    #    rslt.is_retargeted, 
    #    rslt.retarget_pos,
    #    rslt.trans_vars
    #)
