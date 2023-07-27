import numpy as np
from difflib import SequenceMatcher
from collections import OrderedDict, Counter

# from .variant import Variant


def _phase(
    contig_dict,
    skips,
    target_pos,
    is_indel,
    local_thresh,
    match_penal,
    max_common_substr_len,
):
    if len(contig_dict) == 0:
        return None

    # contig_dict = to_tight_layout(contig_dict, target_pos)
    snvs, indels, actual_target = profile_non_refs(contig_dict, target_pos, is_indel)

    if actual_target.is_null():
        return None

    contig_dict, snvs, indels = crop_contig(
        contig_dict, skips, snvs, indels, actual_target.pos
    )

    if not snvs and not indels:
        trim_contig(contig_dict, actual_target.pos, actual_target.end_pos)
    elif not indels:
        lt_end = find_peak(
            contig_dict, actual_target, snvs, match_penal, local_thresh, True
        )
        rt_end = find_peak(
            contig_dict, actual_target, snvs, match_penal, local_thresh, False
        )
        trim_contig(contig_dict, lt_end, rt_end)
    else:
        remove_common_substrs(contig_dict, actual_target.pos, max_common_substr_len)
        remove_unclustered_snvs(
            contig_dict, actual_target, snvs, match_penal, local_thresh
        )

    return greedy_phasing(contig_dict)


class NonReferenceEvent(object):
    def __init__(self, pos, ref, alt, qual):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.end_pos = pos + len(ref) - 1
        self.qual = qual

    def __eq__(self, other):
        pos_eq = self.pos == other.pos
        ref_eq = self.ref == other.ref
        alt_eq = self.alt == other.alt

        return all([pos_eq, ref_eq, alt_eq])

    def is_null(self):
        if self.ref == self.alt == "N":
            return True
        else:
            return False


# def _rt_aln(prev_ref, prev_alt, prev_qual, curr_ref, curr_alt, curr_qual):
#    if len(prev_ref) > len(prev_alt):
#        if prev_ref[1] == curr_ref:
#            return (prev_ref[1 : ] + curr_ref, curr_ref, curr_qual)
#    else:
#        if prev_alt[1] == curr_ref:
#            return (curr_ref, prev_alt[1 : ] + curr_ref, curr_qual +  prev_qual[1 : ])

#    return None

# def _update(contig_dicti_lst, i, alned):
#    contig_dicti_lst[i - 1][1][0] = contig_dicti_lst[i - 1][1][0][0]
#    contig_dicti_lst[i - 1][1][1] = contig_dicti_lst[i - 1][1][1][0]
#    contig_dicti_lst[i - 1][1][2] = contig_dicti_lst[i - 1][1][2][0]
#    contig_dicti_lst[i][1][0] = alned[0]
#    contig_dicti_lst[i][1][1] = alned[1]
#    contig_dicti_lst[i][1][2] = alned[2]


# def _to_tight_layout(contig_dict, target_pos):
#    data = [[k, list(v)] for k,v in contig_dict.items()]
#
#    prev_pos, prev_ref, prev_alt, prev_qual = data[0], data[0][1][0], data[0][1][1], data[0][1][2]
#    for i in range(1, len(data)- 1):
#        curr_pos, curr_ref, curr_alt, curr_qual = data[i][0], data[i][1][0], data[i][1][1], data[i][1][2]
#        if curr_ref == curr_alt:
#            if len(prev_ref) != len(prev_alt):
#                alned = rt_aln(prev_ref, prev_alt, prev_qual, curr_ref, curr_alt, curr_qual)
#
#                if alned:
#                    update(data, i, alned)
#        prev_pos, prev_ref, prevr_alt, prev_qual = curr_pos, curr_ref, curr_alt, curr_qual
#
#    return OrderedDict((i[0], (i[1][0], i[1][1], i[1][2])) for i in data)


# def rt_aln(newlst, prev_ref, prev_alt, prev_qual, curr_ref, curr_alt, curr_qual):
#    if len(prev_ref) > len(prev_alt):
#        if prev_ref[1] == curr_ref:

#
# def to_tight_layout(cntg_lst, target_pos):
#
#    _tight = list()
#
#    prev_pos, prep_ref, prev_alt, prev_qual = cntg_lst[0]
#    for i in range(1, len(cntg_lst) -1):
#        curr_pos, curr_ref, curr_alt, curr_qual = cntg_lst[i]
#        if curr_ref == curr_alt:
#            if len(prev_ref) != len(prev_alt):


def parse_contig_for_exons(contig_dict, skips):
    contig_start = list(contig_dict.items())[0][0]
    contig_end = list(contig_dict.items())[-1][0]

    exon_start = contig_start
    exon_end = contig_end
    exons = []
    for skip in skips:
        exon_end = skip[0] - 1
        exons.append((exon_start, exon_end))
        exon_start = skip[1] + 1
        exon_end = contig_end

    exons.append((exon_start, exon_end))

    return contig_start, contig_end, exons


def profile_non_refs(contig_dict, target_pos, is_indel):
    actual_pos = -1
    actual_event = NonReferenceEvent(-1, "N", "N", 0)
    pos_diff = np.inf

    snvs, indels = [], []
    for pos, v in contig_dict.items():
        ref, alt, qual = v[0], v[1], v[2]

        if ref != alt:

            non_ref = NonReferenceEvent(pos, ref, alt, qual)

            if abs(pos - target_pos) < pos_diff:
                if is_indel:
                    if len(ref) != len(alt):
                        actual_pos = pos
                        pos_diff = abs(pos - target_pos)
                        actual_event = non_ref
                else:
                    if len(ref) == len(alt) and ref != alt:
                        actual_pos = pos
                        pos_diff = abs(pos - target_pos)
                        actual_event = non_ref

            if len(ref) == len(alt):
                snvs.append(non_ref)
            else:
                indels.append(non_ref)
    
    if is_indel and (actual_event in indels):
        indels.remove(actual_event)

    if not is_indel and (actual_event in snvs):
        snvs.remove(actual_event)

    return snvs, indels, actual_event


def crop_contig(contig_dict, skips, snvs, indels, target_pos):
    """
    crop to target exon and trim outmost "N" bases
    """

    lt_lim, rt_lim, exons = parse_contig_for_exons(contig_dict, skips)
    for exon in exons:
        if exon[0] <= target_pos <= exon[1]:
            lt_lim = exon[0] - 1
            rt_lim = exon[1] + 1

    lt_disqualified, rt_disqualified = [], []
    for pos, v in contig_dict.items():
        if "N" in v[0] or "N" in v[1]:
            if pos < target_pos:
                lt_disqualified.append(pos)
            else:
                rt_disqualified.append(pos)

    if lt_disqualified:
        lt_lim = max(lt_lim, max(lt_disqualified))

    if rt_disqualified:
        rt_lim = min(rt_lim, min(rt_disqualified))

    contig_dict = {pos: v for pos, v in contig_dict.items() if lt_lim < pos < rt_lim}
    snvs = [snv for snv in snvs if lt_lim < snv.pos < rt_lim]
    indels = [indel for indel in indels if lt_lim < indel.pos < rt_lim]

    return contig_dict, snvs, indels


def loss(i, match_penal, local_thresh):
    return -1 * match_penal * min(i / local_thresh, 1)


def find_peak(contig_dict, target, snvs, match_penal, local_thresh, is_left):

    if is_left:
        loci = [pos for pos, _data in contig_dict.items() if pos <= target.pos][::-1]
        snv_loci = [snv.pos for snv in snvs if snv.pos < target.pos]
    else:
        loci = [pos for pos, _data in contig_dict.items() if pos > target.end_pos]
        snv_loci = [snv.pos for snv in snvs if snv.pos > target.pos]

    peak_locus = -np.inf if is_left else np.inf
    if not loci:
        return peak_locus

    score, gain = 0.0, 1.0
    scores = []
    for i, locus in enumerate(loci):
        if locus in snv_loci:
            score += gain
        else:
            score += loss(i, match_penal, local_thresh)

        scores.append(score)

    peak_score = max(scores)
    if peak_score > 0.0:
        peak_idx = [idx for idx, _score in enumerate(scores) if _score == peak_score][
            -1
        ]
        peak_locus = loci[peak_idx]
        score = peak_score

    if peak_locus == -np.inf:
        peak_locus = target.pos
    elif peak_locus == np.inf:
        peak_locus = target.end_pos

    return peak_locus


def trim_contig(contig_dict, lt_end, rt_end):
    for pos in contig_dict.copy():
        if pos < lt_end:
            del contig_dict[pos]
        elif rt_end < pos:
            del contig_dict[pos]


def greedy_phasing(contig_dict):
    phased_pos = list(contig_dict.keys())[0]
    phased_ref, phased_alt = "", ""
    for pos, _data in contig_dict.items():
        phased_ref += _data[0]
        phased_alt += _data[1]

    return phased_pos, phased_ref, phased_alt


def ext_common_substr(start, contig_dict):
    common_start, common_end = start, start
    common_substr = [common_start]
    for pos, _data in contig_dict.items():
        if pos > start and _data[0] == _data[1]:
            common_start = pos
            common_substr.append(pos)
        elif pos > common_start and _data[0] != _data[1]:
            common_end = pos
            common_substr.append(pos)
            break

    if not common_substr:
        common_substr = [common_start, common_end]

    return common_substr


def near_rt_pos(pos, contig_dict, contig_end):
    found = False
    while not found and pos < contig_end:
        pos += 1
        found = contig_dict.get(pos, False)

    return pos


def near_lt_pos(pos, contig_dict, is_left):
    if is_left:
        not_found = True
    else:
        not_found = False if contig_dict.get(pos, None) else True

    contig_start = list(contig_dict.keys())[0]
    while not_found and contig_start < pos:
        pos -= 1
        _data = contig_dict.get(pos, False)
        if _data:
            not_found = False

            # del involved
            if len(_data[0]) > 1:
                pos += len(_data[0])
    return pos


def common_substrs(contig_dict):
    commons = []
    _keys = list(contig_dict.keys())

    contig_pos, contig_end = _keys[0], _keys[-1]
    while contig_pos < contig_end:
        common_substr = ext_common_substr(contig_pos, contig_dict)
        commons.append(common_substr)

        contig_pos = near_rt_pos(common_substr[-1], contig_dict, contig_end)

    return commons


def del_commons(contig_dict, commons, max_common_str_len, is_left):
    if not is_left:
        commons[::-1]

    deletables = []
    for substr in commons:
        if substr[0] == substr[-1]:
            start = substr[0]
        else:
            start = near_lt_pos(substr[0], contig_dict, is_left)

        end = substr[-1]

        if end - start > max_common_str_len:
            if is_left:
                deletables.append(end)
            else:
                deletables.append(start)

    if deletables:
        _keys = list(contig_dict.keys())

        if is_left:
            lim = max(deletables)
            for locus in _keys:
                if locus < lim:
                    del contig_dict[locus]
        else:
            lim = min(deletables)
            for locus in _keys:
                if locus > lim:
                    del contig_dict[locus]


def remove_common_substrs(contig_dict, pos, max_common_str_len):
    commons = common_substrs(contig_dict)
    lt_commons = [substr for substr in commons if substr[1] < pos]
    rt_commons = [substr for substr in commons if pos <= substr[0]]
    del_commons(contig_dict, lt_commons, max_common_str_len, True)
    del_commons(contig_dict, rt_commons, max_common_str_len, False)


def get_far_indel(contig_dict, is_left):
    if is_left:
        for pos, _data in contig_dict.items():
            if len(_data[0]) != len(_data[1]):
                return NonReferenceEvent(pos, _data[0], _data[1], _data[2])
    else:
        tmp_pos, tmp_ref, tmp_aln, tmp_qual = -1, "N", "N", "N"
        for pos, _data in contig_dict.items():
            if len(_data[0]) != len(_data[1]) and pos > tmp_pos:
                tmp_pos = pos
                tmp_ref = _data[0]
                tmp_aln = _data[1]
                tmp_qual = _data[2]

        if tmp_ref != "N" and tmp_aln != "N":
            return NonReferenceEvent(tmp_pos, tmp_ref, tmp_aln, tmp_qual)


def remove_unclustered_snvs(contig_dict, target, snvs, match_penal, local_thresh):

    # [(pos, (ref, aln, qual))]
    alns = list(contig_dict.items())

    # closed with indels
    if len(alns[0][1][0]) != len(alns[0][1][1]) and len(alns[-1][1][0]) != len(
        alns[-1][1][1]
    ):
        return None

    # lt process
    lt_far_indel = get_far_indel(contig_dict, True)
    if not lt_far_indel:
        lt_far_indel = target

    lt_peak_locus = find_peak(
        contig_dict, lt_far_indel, snvs, match_penal, local_thresh, True
    )

    # rt process
    rt_far_indel = get_far_indel(contig_dict, False)
    if not rt_far_indel:
        rt_far_indel = target

    rt_peak_locus = find_peak(
        contig_dict, rt_far_indel, snvs, match_penal, local_thresh, False
    )

    trim_contig(contig_dict, lt_peak_locus, rt_peak_locus)
