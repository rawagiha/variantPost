import random
#import statistics
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t
from collections import OrderedDict

from .long_reads import shorten_read

cdef extern from "search.h":

    cdef cppclass SearchResult:

        SearchResult() except +

        #vector[int] positions
        #vector[string] ref_bases
        #vector[string] alt_bases
        #vector[string] base_quals
        #vector[int] skip_starts
        #vector[int] skip_ends
        #int retarget_pos
        #string ref, alt
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


cdef inline bint is_qualified_read(read, bint exclude_duplicates):
    
    if exclude_duplicates:
        if (
            read.cigarstring 
            and not read.is_duplicate 
            and not read.is_secondary 
            and not read.is_supplementary 
            and read.reference_end
        ):
            return True
    else:
        if (
            read.cigarstring 
            and not read.is_secondary 
            and not read.is_supplementary 
            and read.reference_end
        ):
            return True
    
    return False


def fetch_reads(
    bam, 
    chrom, 
    pos, 
    chrom_len, 
    window, 
    exclude_duplicates, 
    fetched_reads, 
    est_cov, 
    is_secondary
):
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
    vector[string]& qual_seqs, 
    #vector[int]& mapqs, 
    vector[bool_t]& is_primary,
    #vector[string]& cb,             #for isoanalysis
    int unspliced_local_reference_start,
    unspliced_local_reference_end,
    int window_len,
    int k,
    list tags
):
    cdef object read = read_tuple[0]
    
    if len(read.query_sequence) < window_len:
        read_names.push_back(read.query_name.encode())
        are_reverse.push_back(read.is_reverse)
        cigar_strings.push_back(read.cigarstring.encode())
        aln_starts.push_back(read.reference_start + 1)
        aln_ends.push_back(read.reference_end)
        read_seqs.push_back(read.query_sequence.encode())
        if not read.query_qualities_str:
            qual_seqs.push_back(('F'*len(read.query_sequence)).encode())
        else:
            qual_seqs.push_back(read.query_qualities_str.encode())
    
        if read_tuple[1]:
            is_primary.push_back(True) # to be renamed
        else:
            is_primary.push_back(False)
        tags.append(read.get_tags())
    else:
        res = shorten_read(
            read, 
            unspliced_local_reference_start, 
            unspliced_local_reference_end,
            k
        )
       
        if res:
            read_names.push_back(read.query_name.encode())
            are_reverse.push_back(read.is_reverse)
            cigar_strings.push_back(res[2].encode())
            aln_starts.push_back(res[3])
            aln_ends.push_back(res[4])
            read_seqs.push_back(res[0].encode())
            qual_seqs.push_back(res[1].encode())
            #mapqs.push_back(read.mapping_quality)

            if read_tuple[1]:
                is_primary.push_back(True)
            else:
                is_primary.push_back(False)
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
):
    cdef SearchResult rslt     
    cdef int est_cov = 0
    cdef list fetched_reads = []
    cdef bool_t has_second = False

    # read_fetching from first bam
    est_cov = fetch_reads(
        bam, 
        bam_chrom, 
        pos, 
        chrom_len, 
        window, 
        exclude_duplicates, 
        fetched_reads, 
        est_cov, 
        False if second_bam else True
    )  
    
    if second_bam:
        has_second = True
        est_cov = fetch_reads(
            second_bam, 
            bam_chrom, 
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
    cdef vector[string] qual_seqs
    qual_seqs.reserve(buff_size) 
    #cdef vector[int] mapqs 
    #mapqs.reserve(buff_size)
    cdef vector[bool_t] are_first_bam

    #isoseq#
    #cdef vector[string] cb
    #cb.reserve(buff_size)

    #opposite!!!
    are_first_bam.reserve(buff_size)
    
    tags = []
    cdef int widow_len = unspliced_local_reference_end - unspliced_local_reference_start
    for read in fetched_reads:
        
        if read[0].mapping_quality < mapping_quality_threshold:
            continue
        
        pack_to_lists(
            read, 
            read_names, 
            are_reverse, 
            cigar_strings,
            aln_starts, 
            aln_ends, 
            read_seqs, 
            qual_seqs, 
            #mapqs, 
            are_first_bam,
            #cb,
            unspliced_local_reference_start,
            unspliced_local_reference_end,
            widow_len,
            k,
            tags
        )
         
    _search_target(
        rslt,
        fastafile,
        chrom.encode(),
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
    
   # contig_dict = OrderedDict()
   # for pos, ref_base, alt_base, base_qual in zip(
   #     rslt.positions, rslt.ref_bases, rslt.alt_bases, rslt.base_quals
   # ):
        #qual_chars = base_qual.decode("utf-8")
        #qual = -1
        #if len(qual_chars) == 1:
        #    qual = ord(qual_chars) - 33
        #else:
        #    quals = [ord(c) - 33 for c in qual_chars]
        #    qual = statistics.median(quals)
        
    #    contig_dict[pos] = (
    #        ref_base.decode("utf-8"), 
    #        alt_base.decode("utf-8"), 
    #        base_qual.decode("utf-8")
    #    )
    
    #annot_reads = []
    #for read_name, is_reverse, target_status, is_first_bam in zip(
    #    rslt.read_names, rslt.are_reverse, rslt.target_statuses, rslt.are_from_first_bam
    #):
    #    annot_reads.append(
    #        AnnotatedRead(read_name.decode("utf-8"), is_reverse, is_first_bam, target_status)
    #    )
    
    #skips = [(start, end) for start, end in zip(rslt.skip_starts, rslt.skip_ends)]
    
    return (
        rslt.target_statuses,
        are_reverse,
        are_first_bam,
        tags
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
