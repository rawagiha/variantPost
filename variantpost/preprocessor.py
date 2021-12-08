import re
import random
import numpy as np

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


def is_qualified_read(read, exclude_duplicates):
    if exclude_duplicates:
        return all(
            (
                (not read.is_duplicate),
                (not read.is_secondary),
                read.cigarstring,
                read.reference_start,
            )
        )
    else:
        return all(((not read.is_secondary), read.cigarstring))


def fetch_reads(bam, chrom, pos, chrom_len, window, exclude_duplicates):
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


def get_read_wise_reference_seq(
    chrom,
    aln_start,
    aln_end,
    unspliced_local_reference,
    unspliced_local_reference_start,
    cigar_string,
    cigar_list,
    fasta,
):
    """
    unspliced_local_reference_start: 0-based
    aln_start and aln_end: 1-based
    """
    __pos = aln_start - 1  # 0-based

    if not b"N" in cigar_string:
        start_idx = __pos - unspliced_local_reference_start
        ref_seq = unspliced_local_reference[start_idx : start_idx + (aln_end - __pos)]
    else:
        consuming_operations = (b"M", b"X", b"D")
        ref_seq = ""
        for c in cigar_list:
            op, op_len = c[-1], int(c[:-1])
            if op in consuming_operations:
                ref_seq += fasta.fetch(chrom, __pos, __pos + op_len)
                __pos += op_len
            elif op == b"N":
                __pos += op_len
            else:
                # non consuming operations (b"I", b"S", b"H", b"P")
                pass

    return ref_seq.encode("utf-8")


def make_qual_seq(qual_arr):
    a = np.frombuffer(qual_arr, dtype=np.int8)
    a += 33
    return a.tobytes().decode("utf-8")


def preprocess(
    chrom,
    pos,
    chrom_len,
    bam,
    unspliced_local_reference,
    unspliced_local_reference_start,
    fasta,
    exclude_duplicates,
    window,
    downsample_thresh,
):
    chrom = edit_chrom_prefix(chrom, bam)
    reads = fetch_reads(bam, chrom, pos, chrom_len, window, exclude_duplicates)

    if downsample_thresh < 0:
        sample_factor = 1.0
    else:
        reads, sample_factor = downsampler(chrom, pos, bam, downsample_thresh, reads)

    read_names = []
    are_reverse = []
    cigar_strings = []
    cigar_lists = []
    aln_starts = []
    clipped_starts = []
    aln_ends = []
    clipped_ends = []
    read_seqs = []
    ref_seqs = []
    qual_seqs = []
    mapqs = []
    for read in reads:
        read_names.append(read.query_name.encode("utf-8"))
        are_reverse.append(read.is_reverse)

        cigar_string = read.cigarstring.encode("utf-8")
        cigar_list = cigar_ptrn.findall(cigar_string)

        aln_start = read.reference_start + 1
        start_offset = int(cigar_list[0][:-1]) if cigar_list[0].endswith(b"S") else 0
        clipped_start = aln_start - start_offset

        aln_end = read.reference_end
        if aln_end is None:
            aln_end = aln_start + sum(
                int(c[:-1])
                for c in cigar_list
                if c[-1] in (b"M", b"N", b"D", b"=", b"X")
            )

        end_offset = int(cigar_list[-1][:-1]) if cigar_list[-1].endswith(b"S") else 0
        clipped_end = aln_end + end_offset

        cigar_strings.append(cigar_string)
        cigar_lists.append(cigar_list)
        aln_starts.append(aln_start)
        clipped_starts.append(clipped_start)
        aln_ends.append(aln_end)
        clipped_ends.append(clipped_end)
        read_seqs.append(read.query_sequence.encode("utf-8"))
        ref_seqs.append(
            get_read_wise_reference_seq(
                chrom,
                aln_start,
                aln_end,
                unspliced_local_reference,
                unspliced_local_reference_start,
                cigar_string,
                cigar_list,
                fasta,
            )
        )
        qual_seqs.append(make_qual_seq(read.query_qualities))
        mapqs.append(read.mapping_quality)
    return (
        read_names,
        are_reverse,
        cigar_strings,
        cigar_lists,
        aln_starts,
        clipped_starts,
        aln_ends,
        clipped_ends,
        read_seqs,
        ref_seqs,
        qual_seqs,
        mapqs,
    )
