import pysam

IS_REF_MAPPED = [0, 7, 8]
IS_REF_SKIPPED = [2, 3]
IS_REF_CONSUMING = (0, 2, 3, 7, 8)
IS_QRY_CONSUMING = (0, 1, 4, 7, 8)

CIGAR = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"}


def find_cut_start_idx(cigar_tuple, aln_start, aln_end, cut_start, cut_end, k):

    if cut_end <= aln_start or aln_end <= cut_start:
        return (-1,)

    start_idx, end_idx = -1, -1
    start_span_start, start_span_end = -1, -1
    end_span_start, end_span_end = -1, -1

    starts_with_sclip = True if cigar_tuple[0][0] == 4 else False
    ends_with_sclip = True if cigar_tuple[len(cigar_tuple) - 1][0] == 4 else False

    sclip_starter, sclip_ender = False, False
    curr_pos, prev_pos = aln_start, aln_start
    for i, c in enumerate(cigar_tuple):
        op, op_len = c[0], c[1]

        prev_pos = curr_pos
        move_len = op_len if op in IS_REF_CONSUMING else 0
        curr_pos += move_len

        # first cut-start mapped pass
        if start_idx == -1 and cut_start < curr_pos:
            if op in IS_REF_MAPPED:
                if i == 1 and starts_with_sclip and cut_start < aln_start and k == 1:
                    start_idx = 0
                    sclip_starter = True
                else:
                    start_idx = i
                start_span_start = prev_pos
                start_span_end = curr_pos

        # first cut-start passed and still overlapped
        if start_idx != -1 and prev_pos <= cut_end:
            # mapped and further
            if op in IS_REF_MAPPED and end_idx <= i:
                if i == len(cigar_tuple) - 2 and ends_with_sclip and aln_end < cut_end and k == 1:
                    end_idx = len(cigar_tuple) - 1
                    sclip_ender = True
                else:
                    end_idx = i
                end_span_start = prev_pos
                end_span_end = curr_pos - 1

    return (
        start_idx,
        start_span_start,
        start_span_end,
        end_idx,
        end_span_start,
        end_span_end,
        sclip_starter,
        sclip_ender,
    )


def shorten_cigar(
    cigar_tuple,
    cut_start,
    cut_end,
    start_idx,
    start_span_start,
    start_span_end,
    end_idx,
    end_span_start,
    end_span_end,
):
    new_aln_start = max(start_span_start, cut_start)
    new_aln_end = min(end_span_end, cut_end)

    if start_idx == end_idx:
        return (
            "{}{}".format(
                (new_aln_end + 1 - new_aln_start), CIGAR[cigar_tuple[start_idx][0]]
            ),
            new_aln_start,
            new_aln_end,
        )
    else:
        shortend = []

        if cigar_tuple[start_idx][0] == 4:
            shortend.append(
                "{}{}".format(
                    cigar_tuple[start_idx][1], CIGAR[cigar_tuple[start_idx][0]]
                )
            )
        else:
            shortend.append(
                "{}{}".format(
                    (start_span_end - new_aln_start), CIGAR[cigar_tuple[start_idx][0]]
                )
            )

        for i in range(start_idx + 1, end_idx):
            shortend.append("{}{}".format(cigar_tuple[i][1], CIGAR[cigar_tuple[i][0]]))

        if cigar_tuple[end_idx][0] == 4:
            shortend.append(
                "{}{}".format(cigar_tuple[end_idx][1], CIGAR[cigar_tuple[end_idx][0]])
            )
        else:
            shortend.append(
                "{}{}".format(
                    (new_aln_end + 1 - end_span_start), CIGAR[cigar_tuple[end_idx][0]]
                )
            )

    return "".join(shortend), new_aln_start, new_aln_end


def read_cut(
    seq,
    quals,
    cigar_tuple,
    orig_aln_start,
    new_aln_start,
    orig_aln_end,
    new_aln_end,
    include_lt_clip,
    include_rt_clip,
):
    seq_start_idx, seq_end_idx = -1, -1
    curr_qry_idx = 0
    curr_pos, prev_pos = orig_aln_start, orig_aln_start
    for i, c in enumerate(cigar_tuple):
        op, op_len = c[0], c[1]

        prev_pos = curr_pos
        move_len = op_len if op in IS_REF_CONSUMING else 0
        curr_pos += move_len

        if prev_pos <= new_aln_start <= curr_pos:
            seq_start_idx = curr_qry_idx + new_aln_start - prev_pos

        if prev_pos <= new_aln_end <= curr_pos:
            seq_end_idx = curr_qry_idx + new_aln_end - prev_pos

        # if seq_start_idx > -1 and seq_end_idx > -1:
        #    break

        qry_move = op_len if op in IS_QRY_CONSUMING else 0
        curr_qry_idx += qry_move

    # starts_with_sclip = True if cigar_tuple[0][0] == 4 else False
    # ends_with_sclip = True if cigar_tuple[len(cigar_tuple) - 1][0] == 4 else False

    _start, _end = seq_start_idx, seq_end_idx + 1

    if include_lt_clip:
        _start = 0

    if include_rt_clip:
        _end = len(seq)

    # return seq[seq_start_idx : seq_end_idx+1], quals[seq_start_idx : seq_end_idx+1]
    return seq[_start:_end], quals[_start:_end]


def shorten_read(read, cut_start, cut_end, k):

    res = find_cut_start_idx(
        read.cigartuples,
        read.reference_start + 1,
        read.reference_end,
        cut_start,
        cut_end,
        k
    )

    check = (i == -1 for i in res)
    if any(check):
        return

    (
        start_idx,
        start_span_start,
        start_span_end,
        end_idx,
        end_span_start,
        end_span_end,
        include_lt_clip,
        include_rt_clip,
    ) = res

    shortend_cigar, new_aln_start, new_aln_end = shorten_cigar(
        read.cigartuples,
        cut_start,
        cut_end,
        start_idx,
        start_span_start,
        start_span_end,
        end_idx,
        end_span_start,
        end_span_end,
    )

    shortend_read_seq, shortend_qual_vec = read_cut(
        read.query_sequence,
        read.query_qualities,
        read.cigartuples,
        read.reference_start + 1,
        new_aln_start,
        read.reference_end,
        new_aln_end,
        include_lt_clip,
        include_rt_clip,
    )

    return (
        shortend_read_seq,
        shortend_qual_vec,
        shortend_cigar,
        new_aln_start,
        new_aln_end,
    )
