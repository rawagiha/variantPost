import math
import statistics
from difflib import SequenceMatcher
from collections import OrderedDict, Counter


def _phase(
    contig_dict,
    skips,
    target_pos,
    is_indel,
    local_thresh,
    base_qual_thresh,
    match_penal,
    max_common_substr_len,
):
    if len(contig_dict) == 0:
        return None

    snvs, indels, actual_target = profile_non_refs(contig_dict, target_pos, is_indel)
    if actual_target.is_null():
        return None

    contig_dict, snvs, indels = crop_contig(
        contig_dict, skips, snvs, indels, actual_target.pos, base_qual_thresh
    )

    if not snvs and not indels:
        trim_contig(contig_dict, actual_target.pos, actual_target.end_pos)
        return greedy_phasing(contig_dict)

    # to capture cases where phasable only after right alignment
    lt_max = lt_max_lim(actual_target, indels, max_common_substr_len, contig_dict)
    phasable_dist = phasable_dist_to_mismatch(match_penal, local_thresh, contig_dict)
    rt_min = rt_min_lim(
        actual_target, snvs, indels, phasable_dist, max_common_substr_len, contig_dict
    )

    if not indels:
        lt_end = find_peak(
            contig_dict, actual_target, snvs, match_penal, local_thresh, True
        )

        rt_end = find_peak(
            contig_dict, actual_target, snvs, match_penal, local_thresh, False
        )
        rt_end = max(rt_end, rt_min)

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


def phasable_dist_to_mismatch(match_penal, local_thresh, contig_dict):
    score = 0
    for i in range(len(contig_dict)):
        score += loss(i, match_penal, local_thresh)
        if score < -1.0:
            return i

    return len(contig_dict)


def is_rt_shiftable(non_ref, contig_dict):
    if len(non_ref.ref) == len(non_ref.alt):
        return False

    longer = non_ref.alt if len(non_ref.alt) > len(non_ref.ref) else non_ref.ref
    shorter = non_ref.alt if len(non_ref.alt) < len(non_ref.ref) else non_ref.ref

    # complex
    if longer[: len(shorter)] != shorter:
        return False

    rt_data = contig_dict.get(non_ref.end_pos + 1, None)
    if rt_data:
        if rt_data[0] != rt_data[1]:
            return False

        return longer[1] == rt_data[0]
    else:
        return False


def has_rt_variants(subject_pos, snvs, indels):
    for snv in snvs:
        if subject_pos < snv.pos:
            return True

    for indel in indels:
        if subject_pos < indel.pos:
            return True

    return False


def rt_aln_pos(non_ref, contig_dict, contig_end):
    curr_pos = non_ref.pos
    curr_allele = non_ref.alt if len(non_ref.alt) > len(non_ref.ref) else non_ref.ref
    variant_len = len(non_ref.ref)
    is_alignable = True
    while is_alignable and curr_pos < contig_end:
        variant_end = curr_pos + variant_len
        data = contig_dict.get(variant_end, None)
        if data:
            if curr_allele[1] == data[0]:
                curr_pos += 1
                curr_allele = curr_allele[1:] + data[0]
            else:
                return curr_pos
        else:
            return curr_pos

    return curr_pos


def lt_max_lim(target, indels, max_common_substr_len, contig_dict):
    min_lim = math.inf
    lt_indels = [indel for indel in indels if indel.pos < target.pos]

    if lt_indels:
        nearest_indel = lt_indels[-1]
        if target.pos - nearest_indel.pos >= max_common_substr_len:
            rt_pos = rt_aln_pos(nearest_indel, contig_dict, target.pos)
            if rt_pos - target.pos < max_common_substr_len:
                return rt_pos
    return min_lim


def rt_min_lim(target, snvs, indels, phasable_dist, max_common_substr_len, contig_dict):
    min_lim = -1

    nearest_snv_pos, nearest_indel_pos = math.inf, math.inf
    if is_rt_shiftable(target, contig_dict):

        contig_end = next(reversed(contig_dict))

        # snvs, indels sorte by construction
        for snv in snvs:
            if target.pos < snv.pos:
                nearest_snv_pos = snv.pos
                break

        for indel in indels:
            if target.pos < indel.pos:
                nearest_indel_pos = indel.pos
                break

        # note: math.inf < math.inf -> False
        if nearest_snv_pos < nearest_indel_pos:
            if nearest_snv_pos - target.pos >= phasable_dist:
                rt_pos = rt_aln_pos(target, contig_dict, contig_end)
                if nearest_snv_pos - rt_pos < phasable_dist:
                    return nearest_snv_pos

        elif nearest_indel_pos < nearest_snv_pos:
            if nearest_indel_pos - target.pos >= max_common_substr_len:
                rt_pos = rt_aln_pos(target, contig_dict, contig_end)
                if nearest_indel_pos - rt_pos < max_common_substr_len:
                    return nearest_indel_pos

    return min_lim


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
    pos_diff = math.inf

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


def to_numeric_qual(char_qual):
    return statistics.median([ord(c) - 33 for c in char_qual])


def crop_contig(contig_dict, skips, snvs, indels, target_pos, base_qual_thresh):
    """
    crop to target exon and trim outmost "N" & low-qualbases
    """

    lt_lim, rt_lim, exons = parse_contig_for_exons(contig_dict, skips)
    for exon in exons:
        if exon[0] <= target_pos <= exon[1]:
            lt_lim = exon[0] - 1
            rt_lim = exon[1] + 1

    lt_disqualified, rt_disqualified = [], []
    for pos, v in contig_dict.items():
        if "N" in v[0] or "N" in v[1] or to_numeric_qual(v[2]) < base_qual_thresh:
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

    peak_locus = -math.inf if is_left else math.inf
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

    if peak_locus == -math.inf:
        peak_locus = target.pos
    elif peak_locus == math.inf:
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

    if lt_far_indel.pos <= target.pos:
        lt_peak_locus = find_peak(
            contig_dict, lt_far_indel, snvs, match_penal, local_thresh, True
        )
    else:
        lt_peak_locus = -math.inf

    # rt process
    rt_far_indel = get_far_indel(contig_dict, False)
    if not rt_far_indel:
        rt_far_indel = target

    if rt_far_indel.pos >= target.pos:
        rt_peak_locus = find_peak(
            contig_dict, rt_far_indel, snvs, match_penal, local_thresh, False
        )
    else:
        rt_peak_locus = math.inf

    trim_contig(contig_dict, lt_peak_locus, rt_peak_locus)
