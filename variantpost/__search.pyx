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

# 既存の long_reads はそのまま利用（ここがボトルネックになる可能性はあるが、ロジック維持のため）
from .long_reads import shorten_read

cdef extern from "search.h":
    cdef cppclass SearchResult:
        SearchResult() except +
        vector[int] positions
        vector[string] ref_bases
        vector[string] alt_bases
        int ppos
        string ref, alt, pref, palt, pltseq, prtseq
        vector[int] target_statuses

    void _search_target(
        SearchResult&, string&, string&, int, string&, string&,
        int, float, int, int, int, int, int, int, int, int, int,
        vector[string]&, vector[bool_t]&, vector[string]&,
        vector[int]&, vector[int]&, vector[string]&,
        vector[string]&, vector[bool_t]&, bool_t
    )

# --- Utilities ---
cdef inline int int_max(int a, int b) nogil: return a if a > b else b
cdef inline int int_min(int a, int b) nogil: return a if a < b else b

cdef inline bint is_qualified_read_fast(int flag, int ref_start, bint exclude_duplicates) nogil:
    # 0x904: secondary (0x100), supplementary (0x800), unmapped (0x4)
    if (flag & 0x904): return False
    if exclude_duplicates and (flag & 0x400): return False
    if ref_start == -1: return False
    return True

# --- Main Function ---
cpdef object search_target(
     object bam,
     object second_bam,
     int chrom_len,
     bint exclude_duplicates,
     int window,
     string fastafile,
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
     bint get_tags_flag=True
):
    cdef SearchResult rslt     
    cdef int start = int_max(0, pos - window)
    cdef int end = int_min(pos + window, chrom_len)
    cdef int window_len = unspliced_local_reference_end - unspliced_local_reference_start

    # C++ Vectors (ストリームライン用)
    cdef vector[string] read_names, cigar_strings, read_seqs, qual_seqs 
    cdef vector[bool_t] are_reverse, are_first_bam
    cdef vector[int] aln_starts, aln_ends
    
    # Python 側の管理
    cdef list tags = []
    cdef list temp_reads = [] # ダウンサンプリング用の一時保持
    cdef int est_cov = 0
    cdef object read
    cdef bint is_sec_bam

    # 1. Fetch & Filter (First Pass)
    # 複数の BAM を効率的に回す
    bams = [(bam, False)]
    if second_bam:
        bams.append((second_bam, True))

    for target_bam, is_sec_bam in bams:
        # fetch は Python オブジェクトを返すが、このイテレータを直接回す
        for read in target_bam.fetch(bam_chrom, start, end, until_eof=False):
            # 属性を一度だけ取得 (Python 呼び出しを最小化)
            r_flag = read.flag
            r_ref_start = read.reference_start
            
            if is_qualified_read_fast(r_flag, r_ref_start, exclude_duplicates):
                r_ref_end = read.reference_end
                if r_ref_start <= pos <= r_ref_end:
                    est_cov += 1
                
                # mapping_quality チェックもここで行う
                if read.mapping_quality >= mapping_quality_threshold:
                    temp_reads.append((read, is_sec_bam))

    # 2. Downsampling
    if downsample_thresh > 0 and est_cov > downsample_thresh:
        n_sample = int(len(temp_reads) * (downsample_thresh / est_cov))
        random.seed(123)
        temp_reads = random.sample(temp_reads, n_sample)

    # 3. Packing (Second Pass) - ここを極限までストリームライン化
    cdef int buff_size = len(temp_reads)
    read_names.reserve(buff_size)
    are_reverse.reserve(buff_size)
    cigar_strings.reserve(buff_size)
    aln_starts.reserve(buff_size)
    aln_ends.reserve(buff_size)
    read_seqs.reserve(buff_size)
    qual_seqs.reserve(buff_size)
    are_first_bam.reserve(buff_size)

    cdef object r_obj
    cdef bint secondary_flag
    cdef str q_seq, q_qual

    for r_obj, secondary_flag in temp_reads:
        q_seq = r_obj.query_sequence
        
        if q_seq is not None and len(q_seq) < window_len:
            # 短いリードの高速パッキング
            read_names.push_back(<string>r_obj.query_name.encode('ascii'))
            are_reverse.push_back(<bool_t>r_obj.is_reverse)
            cigar_strings.push_back(<string>r_obj.cigarstring.encode('ascii'))
            aln_starts.push_back(<int>r_obj.reference_start + 1)
            aln_ends.push_back(<int>r_obj.reference_end)
            read_seqs.push_back(<string>q_seq.encode('ascii'))
            
            q_qual = r_obj.query_qualities_str
            if not q_qual:
                qual_seqs.push_back(<string>(b'F' * len(q_seq)))
            else:
                qual_seqs.push_back(<string>q_qual.encode('ascii'))
            
            are_first_bam.push_back(secondary_flag)
            if get_tags_flag:
                tags.append(r_obj.get_tags())
        else:
            # Long-read 処理 (shorten_read は Python 関数のまま呼び出し)
            res = shorten_read(r_obj, unspliced_local_reference_start, unspliced_local_reference_end, k)
            if res:
                read_names.push_back(<string>r_obj.query_name.encode('ascii'))
                are_reverse.push_back(<bool_t>r_obj.is_reverse)
                cigar_strings.push_back(<string>res[2].encode('ascii'))
                aln_starts.push_back(<int>res[3])
                aln_ends.push_back(<int>res[4])
                read_seqs.push_back(<string>res[0].encode('ascii'))
                qual_seqs.push_back(<string>res[1].encode('ascii'))
                are_first_bam.push_back(secondary_flag)
                if get_tags_flag:
                    tags.append(r_obj.get_tags())

    # 4. C++ への受け渡し
    _search_target(
        rslt, fastafile, chrom.encode('ascii'), pos, ref, alt,
        base_quality_threshold, low_quality_base_rate_threshold,
        match_score, mismatch_penalty, gap_open_penalty, gap_extention_penalty,
        kmer_size, dimer_window, local_threshold,
        unspliced_local_reference_start, unspliced_local_reference_end,
        read_names, are_reverse, cigar_strings, aln_starts, aln_ends,
        read_seqs, qual_seqs, are_first_bam, (second_bam is not None)
    )
    
    # 5. 結果の構築
    contig_dict = OrderedDict()
    cdef size_t i
    for i in range(rslt.positions.size()):
        contig_dict[rslt.positions[i]] = (
            rslt.ref_bases[i].decode('ascii'), 
            rslt.alt_bases[i].decode('ascii'), 
            "F"
        )
    
    return (
        contig_dict, rslt.target_statuses, are_reverse, are_first_bam, tags,
        rslt.ppos, rslt.pref.decode("ascii"), rslt.palt.decode("ascii"),
        rslt.pltseq.decode("ascii"), rslt.prtseq.decode("ascii"),
        rslt.ref.decode("ascii"), rslt.alt.decode("ascii")
    )
