BASES = {"A", "C", "G", "T", "N"}

from pysam import FastaFile


class Variant(object):
    def __init__(self, chrom, pos, ref, alt, reference, window=50):

        self.chrom = chrom
        self.pos = pos  # 1-based

        self.ref = ref
        self.alt = alt

        if not set(list(self.ref)) <= BASES or not set(list(self.alt)) <= BASES:
            raise ValueError("Found characters other than 'A', 'C', 'G', 'T', 'N'.")

        self.reference = reference

        if len(ref) < window:
            self.window = window
        else:
            self.window = len(ref) + window

        self.reference_len = reference.get_reference_length(self.chrom)

        # 1-based
        self.unspliced_local_reference_start = max(0, pos - self.window * 4) + 1
        self.unspliced_local_reference_end = min(
            pos + self.window * 4, self.reference_len
        )

    @property
    def is_indel(self):
        return len(self.ref) != len(self.alt)

    @property
    def is_del(self):
        return len(self.ref) > len(self.alt)

    @property
    def is_ins(self):
        return len(self.ref) < len(self.alt)

    @property
    def is_simple_indel(self):
        if self.is_indel:
            shorter_allele = self.ref if self.is_ins else self.alt
            longer_allele = self.ref if self.is_del else self.alt
            return shorter_allele == longer_allele[: len(shorter_allele)]

        return False

    @property
    def is_normalized(self):
        if self.ref[-1] != self.alt[-1]:
            if self.ref[0] != self.alt[0]:
                return True
            elif len(self.ref) == 1 or len(self.alt) == 1:
                return True

        return False

    def normalize(self, inplace=False):

        if inplace:
            i = self
        else:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)

        condition_1 = i.ref[-1] == i.alt[-1] != "N"
        lhs = i.reference.fetch(i.chrom, i.unspliced_local_reference_start, i.pos - 1)[
            ::-1
        ]

        n = 0
        lhs_len = len(lhs)
        while condition_1 and n < lhs_len:

            left_base = lhs[n]

            i.ref = left_base + i.ref[:-1]
            i.alt = left_base + i.alt[:-1]
            i.pos -= 1

            condition_1 = i.ref[-1] == i.alt[-1] != "N"

            n += 1

        condition_2 = i.ref[0] == i.alt[0]
        condition_3 = len(i.ref) > 1 and len(i.alt) > 1
        while condition_2 and condition_3:
            i.ref = i.ref[1:]
            i.alt = i.alt[1:]
            i.pos += 1
            condition_2 = i.ref[0] == i.alt[0]
            condition_3 = len(i.ref) > 1 and len(i.alt) > 1

        if inplace:
            return None
        else:
            return i

    def __eq__(self, other):
        if (
            self.chrom == other.chrom
            and self.pos == other.pos
            and self.ref == other.ref
            and self.alt == other.alt
        ):
            return True
        else:
            if self.chrom != other.chrom:
                return False

            i, j = self.normalize(), other.normalize()

            return (
                i.chrom == j.chrom
                and i.pos == j.pos
                and i.ref == j.ref
                and i.alt == j.alt
            )

    def __hash__(self):
        i = self.normalize()
        return hash((i.chrom, i.pos, i.ref, i.alt))

    def __getstate__(self):
        return (self.chrom, self.pos, self.ref, self.alt, self.reference.filename)

    def __setstate__(self, state):
        self.chrom = state[0]
        self.pos = state[1]
        self.ref = state[2]
        self.alt = state[3]

        self.reference = FastaFile(state[4])

    def __bool__(self):
        if self.ref == self.alt:
            return False
        else:
            return True

    def lt_pos(self):
        if not self.is_indel:
            return self.pos
        else:
            if self.is_simple_indel:
                curr_allele = self.alt if self.is_ins else self.ref
                curr_pos = self.pos
                if curr_allele[0] != curr_allele[-1]:
                    return curr_pos

                lt_flank = self.reference.fetch(
                    self.chrom, self.unspliced_local_reference_start - 1, curr_pos
                )[::-1]
                i, n = 0, len(lt_flank)
                while i < n:
                    tmp_allele = lt_flank[i] + curr_allele[:-1]
                    i += 1

                    if tmp_alelle[0] != tmp_alelle[-1]:
                        break
                    curr_allele = tmp_allele

                return curr_pos - i
            else:
                return self.pos

    def rt_pos(self):
        if not self.is_indel:
            return self.pos + len(self.ref) - 1
        else:
            if self.is_simple_indel:
                if self.is_ins:
                    curr_pos = self.pos
                    curr_allele = self.alt
                else:
                    curr_pos = self.pos + len(self.ref) - 1
                    curr_allele = self.ref

                rt_flank = self.reference.fetch(
                    self.chrom, curr_pos, self.unspliced_local_reference_end
                )
                i, n = 0, len(rt_flank)
                while i < n:
                    tmp_alelle = curr_allele[1:] + rt_flank[i]
                    i += 1

                    if tmp_alelle[0] != tmp_alelle[-1]:
                        break
                    curr_allele = tmp_alelle

                return curr_pos + i - 1
            else:
                return self.pos + len(self.ref) - 1

    def query_vcf(self, vcf, contig_name="", match_by="normalization"):
        vcf_chrom = contig_name if contig_name else self.chrom

        lt_pos, rt_pos = self.lt_pos, self.rt_pos
        if match_by == "normalization":
            search_start, search_end = lt_pos - 1, rt_pos
        else:
            search_start, search_end = (
                self.unspliced_local_reference_start - 1,
                self.unspliced_local_reference_end,
            )

        vcf_entries = vcf.fetch(vcf_chrom, search_start, search_end)

        hits = []
        if match_by == "normalization":
            for _entry in vcf_entries:
                for _alt in _entry.alts:
                    if self == Variant(
                        self.chrom, _entry.pos, _entry.ref, _alt, self.reference
                    ):
                        hits.append(_entry)
                        break
        else:
            for _entry in vcf_entries:
                ref_len = len(_entry.ref)
                for _alt in _entry.alts:
                    if ref_len == len(_aln):
                        if lt_pos <= _entry.pos <= rt_pos:
                            hits.append(_entry)
                            break
                    else:
                        v = Variant(
                            self.chrom, _entry.pos, _entry.ref, _alt, self.reference
                        )
                        _lt_pos, _rt_pos = v.lt_pos, v.rt_pos

                        if lt_pos <= _lt_pos <= rt_pos or lt_pos <= _rt_pos <= rt_pos:
                            hits.append(_entry)
                            break

        return hits


#
#
#
#
#        matchbys = ["normalization", "??"]
