import re
import random
import numpy as np
import time
import cython

#from libcpp.string cimport string
#from libcpp.vector cimport vector
#from libcpp cimport bool as bool_t

RAW_DEPTH=0
cigar_ptrn = re.compile(b"[0-9]+[MIDNSHPX=]")

def edit_chrom_prefix(chrom, bam):
    """Add/delete "chr" prefix for incomaptibe nomenclature
    """
    chrom_names_in_bam = bam.references

    if chrom in chrom_names_in_bam:
        return chrom
    else:
        if chrom.startswith("chr"):
            return chrom.replace("chr", "")
        else:
            return "chr" + chrom


cdef bint is_qualified_read(object read, bint exclude_duplicates):
    
    if exclude_duplicates:
        if read.cigarstring and (not read.is_duplicate) and  (not read.is_secondary) and (not read.is_supplementary) and read.reference_end:
            return True
    else:
        if read.cigarstring and (not read.is_secondary) and (not read.is_supplementary) and read.reference_end:
            return True
    
    return False

cdef object fetch_reads(object bam, str chrom, int pos, int chrom_len, int window, bint exclude_duplicates):
    reads = bam.fetch(
        chrom, max(0, pos - window), min(pos + window, chrom_len), until_eof=False
    )
    #return (read for read in reads if is_qualified_read(read, exclude_duplicates, i))
    return (read for read in reads if is_qualified_read(read, exclude_duplicates))


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


def get_spliced_reference_seq(
    chrom,
    aln_start,
    cigar_list,
    fasta,
):
    __pos = aln_start - 1  # 0-based

    consuming_operations = (b"M", b"X", b"D")
    ref_seq = ""
    for c in cigar_list:
        op, op_len = c[-1:], int(c[:-1])
        if op in consuming_operations:
            ref_seq += fasta.fetch(chrom, __pos, __pos + op_len)
            __pos += op_len
        elif op == b"N":
            __pos += op_len
        else:
            # non consuming operations (b"I", b"S", b"H", b"P")
            pass
   
    return ref_seq.encode()

def make_qual_seq(qual_arr):
    a = np.frombuffer(qual_arr, dtype=np.int8)
    a += 33
    return a.tobytes().decode("utf-8")




cdef void  process(object read, 
                         str chrom, 
                         object fasta, 
                         list read_names, 
                         list are_reverse, 
                         list cigar_strings, 
                         list aln_starts, 
                         list aln_ends, 
                         list read_seqs, 
                         list ref_seqs, 
                         list qual_seqs, 
                         list mapqs, 
                         list is_primary, bint is_secondary):
    cdef bytes cigar_string = read.cigarstring.encode()
    
    read_names.append(read.query_name.encode())
    are_reverse.append(read.is_reverse)
    
    cdef int aln_start = read.reference_start + 1
    
    cigar_strings.append(cigar_string)
    aln_starts.append(aln_start)
    aln_ends.append(read.reference_end)
    read_seqs.append(read.query_sequence.encode())
    if b"N" in cigar_string:
        cigar_list = cigar_ptrn.findall(cigar_string)
        ref_seqs.append(get_spliced_reference_seq(chrom, aln_start, cigar_list, fasta))
    else:
        ref_seqs.append(b"")
    qual_seqs.append(read.query_qualities)
    mapqs.append(read.mapping_quality)
    
    if is_secondary:
        is_primary.append(False)
    else:
        is_primary.append(True)


def preprocess(
    chrom,
    pos,
    chrom_len,
    bam,
    second_bam,
    unspliced_local_reference,
    unspliced_local_reference_start,
    fasta,
    exclude_duplicates,
    window,
    downsample_thresh,
):
    chrom = edit_chrom_prefix(chrom, bam)
    
    tt = time.time()
    reads = fetch_reads(bam, chrom, pos, chrom_len, window, exclude_duplicates)
    print("I/O by pysam", time.time() - tt)

    tt = time.time()
    if downsample_thresh < 0:
        sample_factor = 1.0
    else:
        reads, sample_factor = downsampler(chrom, pos, bam, downsample_thresh, reads)
    
    read_names = []  #
    are_reverse = [] #
    cigar_strings = [] #
    aln_starts = [] #
    aln_ends = []     #
    read_seqs = []  #
    ref_seqs = []
    qual_seqs = [] #
    mapqs = [] #
    are_first_bam = [] # from primary bam
    
    #cdef object read
    cdef bytes cigar_string
    cdef int aln_start, aln_end
    
    for read in reads: 
        
        process(read, chrom, fasta, read_names, are_reverse, cigar_strings, 
                aln_starts, aln_ends, read_seqs, ref_seqs, qual_seqs, mapqs, are_first_bam, False)
    
    if second_bam:
        reads = fetch_reads(bam, chrom, pos, chrom_len, window, exclude_duplicates)
        
        for read in reads:
            process(read, chrom, fasta, read_names, are_reverse, cigar_strings,
                    aln_starts, aln_ends, read_seqs, ref_seqs, qual_seqs, mapqs, are_first_bam, True)
    
    
    return (
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
