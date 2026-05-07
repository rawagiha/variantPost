import math

def loss(i, match_penal, local_thresh):
    return -1 * match_penal * min(i / local_thresh, 1)

class NonReferenceEvent:
    # __slots__ を用いることで辞書割り当てを回避し、属性アクセスを劇的に高速化・省メモリ化
    __slots__ = ['pos', 'ref', 'alt', 'end_pos', 'qual']

    def __init__(self, pos, ref, alt, qual):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.end_pos = pos + len(ref) - 1
        self.qual = qual

    def __eq__(self, other):
        if not isinstance(other, NonReferenceEvent):
            return False
        return self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def is_null(self):
        return self.ref == "N" and self.alt == "N"


def _phase(
    contig_dict,
    skips,
    target_pos,
    is_indel,
    local_thresh,
    base_qual_thresh,
    trans_vars,
    match_penal,
    max_common_substr_len,
):
    if not contig_dict:
        return None

    snvs, indels, actual_target = profile_non_refs(contig_dict, target_pos, is_indel)
    if actual_target.is_null():
        return None

    contig_dict, snvs, indels = crop_contig(
        contig_dict,
        skips,
        snvs,
        indels,
        actual_target.pos,
        base_qual_thresh,
        trans_vars,
    )

    if not snvs and not indels:
        trim_contig(contig_dict, actual_target.pos, actual_target.end_pos)
        return greedy_phasing(contig_dict)

    lt_max = lt_max_lim(actual_target, indels, max_common_substr_len, contig_dict)
    phasable_dist = phasable_dist_to_mismatch(match_penal, local_thresh, contig_dict)
    rt_min = rt_min_lim(
        actual_target, snvs, indels, phasable_dist, max_common_substr_len, contig_dict
    )

    if not indels:
        lt_end = find_peak(contig_dict, actual_target, snvs, match_penal, local_thresh, True)
        rt_end = find_peak(contig_dict, actual_target, snvs, match_penal, local_thresh, False)
        rt_end = max(rt_end, rt_min)

        trim_contig(contig_dict, lt_end, rt_end)
    else:
        remove_common_substrs(contig_dict, actual_target.pos, max_common_substr_len)
        remove_unclustered_snvs(contig_dict, actual_target, snvs, match_penal, local_thresh)
    
    return greedy_phasing(contig_dict)


def phasable_dist_to_mismatch(match_penal, local_thresh, contig_dict):
    score = 0.0
    for i in range(len(contig_dict)):
        # loss() のインライン化による関数呼び出しオーバーヘッドの削減
        score -= match_penal * min(i / local_thresh, 1.0)
        if score < -1.0:
            return i
    return len(contig_dict)


def is_rt_shiftable(non_ref, contig_dict):
    len_ref, len_alt = len(non_ref.ref), len(non_ref.alt)
    if len_ref == len_alt:
        return False

    if len_alt > len_ref:
        longer, shorter = non_ref.alt, non_ref.ref
    else:
        longer, shorter = non_ref.ref, non_ref.alt

    if not longer.startswith(shorter):
        return False

    rt_data = contig_dict.get(non_ref.end_pos + 1)
    if rt_data and rt_data[0] == rt_data[1] and longer[1] == rt_data[0]:
        return True
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
    
    while curr_pos < contig_end:
        variant_end = curr_pos + variant_len
        data = contig_dict.get(variant_end)
        if data and curr_allele[1] == data[0]:
            curr_pos += 1
            curr_allele = curr_allele[1:] + data[0]
        else:
            break
    return curr_pos


def lt_max_lim(target, indels, max_common_substr_len, contig_dict):
    # indels は位置順と仮定し、リストの末尾から探索することで速度向上
    for indel in reversed(indels):
        if indel.pos < target.pos:
            if target.pos - indel.pos >= max_common_substr_len:
                rt_pos = rt_aln_pos(indel, contig_dict, target.pos)
                if rt_pos - target.pos < max_common_substr_len:
                    return rt_pos
            break
    return math.inf


def rt_min_lim(target, snvs, indels, phasable_dist, max_common_substr_len, contig_dict):
    if not is_rt_shiftable(target, contig_dict):
        return -1

    contig_end = next(reversed(contig_dict.keys()))

    nearest_snv_pos = next((s.pos for s in snvs if target.pos < s.pos), math.inf)
    nearest_indel_pos = next((i.pos for i in indels if target.pos < i.pos), math.inf)

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

    return -1


def profile_non_refs(contig_dict, target_pos, is_indel):
    actual_event = NonReferenceEvent(-1, "N", "N", 0)
    pos_diff = math.inf

    snvs, indels = [], []
    for pos, (ref, alt, qual) in contig_dict.items():
        if ref == alt:
            continue

        non_ref = NonReferenceEvent(pos, ref, alt, qual)
        is_cur_indel = len(ref) != len(alt)
        diff = abs(pos - target_pos)

        if diff < pos_diff and is_indel == is_cur_indel:
            pos_diff = diff
            actual_event = non_ref

        if is_cur_indel:
            indels.append(non_ref)
        else:
            snvs.append(non_ref)

    # list.remove の線形探索を避け、リスト内包表記で再構築
    if is_indel:
        indels = [x for x in indels if x != actual_event]
    else:
        snvs = [x for x in snvs if x != actual_event]

    return snvs, indels, actual_event


def to_numeric_qual(char_qual):
    # 重い statistics.median を排除し、自前でソート・中央値算出を行う
    if not char_qual:
        return 0
    quals = sorted(ord(c) - 33 for c in char_qual)
    n = len(quals)
    mid = n // 2
    if n % 2 == 0:
        return (quals[mid - 1] + quals[mid]) / 2.0
    return quals[mid]


def crop_contig(contig_dict, skips, snvs, indels, target_pos, base_qual_thresh, trans_vars):
    keys = list(contig_dict.keys())
    contig_start, contig_end = keys[0], keys[-1]

    lt_lim, rt_lim = contig_start, contig_end

    # エクソンの検索を最適化
    exon_start = contig_start
    for skip in skips:
        exon_end = skip[0] - 1
        if exon_start <= target_pos <= exon_end:
            lt_lim, rt_lim = exon_start - 1, exon_end + 1
            break
        exon_start = skip[1] + 1
    else:
        if exon_start <= target_pos <= contig_end:
            lt_lim, rt_lim = exon_start - 1, contig_end + 1

    lt_disqualified_max = -math.inf
    rt_disqualified_min = math.inf

    if trans_vars:
        for v_str in trans_vars:
            v_ = v_str.decode("utf-8").split("_")
            pos_ = int(v_[0])
            res = contig_dict.get(pos_)
            if res and v_[2] == res[1]:
                if pos_ < target_pos:
                    if pos_ > lt_disqualified_max: lt_disqualified_max = pos_
                elif pos_ < rt_disqualified_min:
                    rt_disqualified_min = pos_

    for pos, (ref, alt, qual) in contig_dict.items():
        if "N" in ref or "N" in alt or to_numeric_qual(qual) < base_qual_thresh:
            if pos < target_pos:
                if pos > lt_disqualified_max: lt_disqualified_max = pos
            elif target_pos < pos:
                if pos < rt_disqualified_min: rt_disqualified_min = pos

    if lt_disqualified_max != -math.inf:
        lt_lim = max(lt_lim, lt_disqualified_max)
    if rt_disqualified_min != math.inf:
        rt_lim = min(rt_lim, rt_disqualified_min)

    # 高速な辞書/リストの再構築
    new_contig_dict = {pos: v for pos, v in contig_dict.items() if lt_lim < pos < rt_lim}
    new_snvs = [snv for snv in snvs if lt_lim < snv.pos < rt_lim]
    new_indels = [indel for indel in indels if lt_lim < indel.pos < rt_lim]

    return new_contig_dict, new_snvs, new_indels


def find_peak(contig_dict, target, snvs, match_penal, local_thresh, is_left):
    if is_left:
        loci = [pos for pos in reversed(contig_dict.keys()) if pos <= target.pos]
        snv_loci = {snv.pos for snv in snvs if snv.pos < target.pos}
    else:
        loci = [pos for pos in contig_dict if pos > target.end_pos]
        snv_loci = {snv.pos for snv in snvs if snv.pos > target.pos}

    if not loci:
        return -math.inf if is_left else math.inf

    score = 0.0
    peak_score = -math.inf
    peak_locus = loci[0]

    # リストを複数回舐めるのをやめ、1パスでピークを探索
    for i, locus in enumerate(loci):
        if locus in snv_loci:
            score += 1.0
        else:
            score -= match_penal * min(i / local_thresh, 1.0)

        # >= にすることで、元のコードの `[-1]` (最後に出現した最大値) と同じ挙動を保証
        if score >= peak_score:
            peak_score = score
            peak_locus = locus

    if peak_score <= 0.0:
        peak_locus = -math.inf if is_left else math.inf

    if peak_locus == -math.inf:
        return target.pos
    elif peak_locus == math.inf:
        return target.end_pos

    return peak_locus


def trim_contig(contig_dict, lt_end, rt_end):
    # .copy() を避けて削除すべきキーだけをリスト化してから削除
    to_delete = [pos for pos in contig_dict if pos < lt_end or rt_end < pos]
    for pos in to_delete:
        del contig_dict[pos]


def greedy_phasing(contig_dict):
    if not contig_dict:
        return None, "", ""
    
    phased_pos = next(iter(contig_dict))
    # += による重い文字列結合を .join() へ変更（最速）
    phased_ref = "".join(v[0] for v in contig_dict.values())
    phased_alt = "".join(v[1] for v in contig_dict.values())

    return phased_pos, phased_ref, phased_alt


def common_substrs(contig_dict):
    commons = []
    keys = list(contig_dict.keys())
    if not keys:
        return commons

    n_keys = len(keys)
    contig_end = keys[-1]
    idx = 0

    while idx < n_keys:
        start_pos = keys[idx]
        common_start = start_pos
        common_end = start_pos
        common_substr = [start_pos]

        j = idx + 1
        while j < n_keys:
            pos = keys[j]
            _data = contig_dict[pos]
            if _data[0] == _data[1]:
                common_start = pos
                common_substr.append(pos)
                j += 1
            elif pos > common_start and _data[0] != _data[1]:
                common_end = pos
                common_substr.append(pos)
                break
            else:
                break
        
        if len(common_substr) == 1:
            common_substr = [common_start, common_end]
        commons.append(common_substr)

        next_pos = common_substr[-1]
        idx = j
        while idx < n_keys and keys[idx] <= next_pos:
            idx += 1

    return commons


def near_lt_pos(pos, contig_dict, is_left):
    not_found = True if is_left else (pos not in contig_dict)
    
    keys = list(contig_dict.keys())
    if not keys:
        return pos
        
    contig_start = keys[0]
    while not_found and contig_start < pos:
        pos -= 1
        _data = contig_dict.get(pos)
        if _data:
            not_found = False
            if len(_data[0]) > 1:
                pos += len(_data[0])
    return pos


def del_commons(contig_dict, commons, max_common_str_len, is_left):
    if not commons:
        return
        
    deletables = []
    for substr in commons:
        start = substr[0] if substr[0] == substr[-1] else near_lt_pos(substr[0], contig_dict, is_left)
        end = substr[-1]

        if end - start > max_common_str_len:
            deletables.append(end if is_left else start)

    if deletables:
        if is_left:
            lim = max(deletables)
            to_del = [k for k in contig_dict if k < lim]
        else:
            lim = min(deletables)
            to_del = [k for k in contig_dict if k > lim]
            
        for k in to_del:
            del contig_dict[k]


def remove_common_substrs(contig_dict, pos, max_common_str_len):
    commons = common_substrs(contig_dict)
    lt_commons = [substr for substr in commons if substr[1] < pos]
    rt_commons = [substr for substr in commons if pos <= substr[0]]
    del_commons(contig_dict, lt_commons, max_common_str_len, True)
    del_commons(contig_dict, rt_commons, max_common_str_len, False)


def get_far_indel(contig_dict, is_left):
    # 走査方向に応じて dict.items() か reversed() を使い分け、最初に見つけた時点で早期リターン
    iterator = contig_dict.items() if is_left else reversed(contig_dict.items())
    for pos, (ref, alt, qual) in iterator:
        if len(ref) != len(alt):
            return NonReferenceEvent(pos, ref, alt, qual)
    return None


def remove_unclustered_snvs(contig_dict, target, snvs, match_penal, local_thresh):
    if not contig_dict:
        return None

    first_v = next(iter(contig_dict.values()))
    last_v = next(reversed(contig_dict.values()))
    
    if len(first_v[0]) != len(first_v[1]) and len(last_v[0]) != len(last_v[1]):
        return None

    lt_far_indel = get_far_indel(contig_dict, True) or target
    lt_peak_locus = find_peak(contig_dict, lt_far_indel, snvs, match_penal, local_thresh, True) if lt_far_indel.pos <= target.pos else -math.inf

    rt_far_indel = get_far_indel(contig_dict, False) or target
    rt_peak_locus = find_peak(contig_dict, rt_far_indel, snvs, match_penal, local_thresh, False) if rt_far_indel.pos >= target.pos else math.inf

    trim_contig(contig_dict, lt_peak_locus, rt_peak_locus)
