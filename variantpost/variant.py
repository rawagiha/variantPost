from typing import Optional, List, Any, Union
from collections import namedtuple
from pysam import FastaFile, VariantFile

BASESET = frozenset({"A", "C", "G", "T", "N"})


MatchedRecord = namedtuple(
    "MatchedRecord",
    [
        "chrom",
        "pos",
        "id",
        "ref",
        "alts",
        "qual",
        "filter",
        "info",
        "format",
        "samples",
    ],
)

def to_tuple(_entry: Any) -> MatchedRecord:
    """pysam's VariantRecard -> MatchedRecord"""
    return MatchedRecord(
        _entry.chrom,
        _entry.pos,
        _entry.id,
        _entry.ref,
        _entry.alts,
        _entry.qual,
        _entry.filter,
        _entry.info,
        _entry.format,
        _entry.samples,
    )


class Variant(object):
    """This class abstracts a variant relative to a linear reference genome.
    Equality holds between :class:`~variantpost.Variant` objects
    if they are identical in the normalized form (equivalent alignments)

    Parameters
    ----------
    chrom : string
        chromosome name.

    pos : integer
        1-based genomic position.

    ref : string
        VCF-style reference allele.

    alt : string
        VCF-style alternative allele.

    reference : pysam.Fasta
        reference FASTA file supplied as
        `pysam.FastaFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.FastaFile>`__ object.

    Raises
    ------
    ValueError
        if the input alleles contain letters other than A, C, G, T and N.

    """
    
    __slots__ = (
        'chrom', 'pos', 'ref', 'alt', 'cpos', 'cref', 'calt', 'phased_as_complex', 'reference', 'window',
        'reference_len', 'k', 'unspliced_local_reference_start', 'unspliced_local_reference_end'
    )

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, reference: FastaFile, window: int = 50):
        
        # input validation
        if not set(ref).issubset(BASESET) or not set(alt).issubset(BASESET):
            raise ValueError("Found characters other than 'A', 'C', 'G', 'T', 'N'.")
        
        self.chrom = chrom
        self.pos = pos  # 1-based
        self.ref = ref
        self.alt = alt

        self.cpos = pos
        self.cref = ref
        self.calt = alt
        self.phased_as_complex = False

        self.reference = reference

        ref_len = len(ref)
        self.window = window if ref_len < window else ref_len + window
        self.reference_len = reference.get_reference_length(self.chrom)

        self.k = 1 if ref_len > window * 4  else 4 # coefficient for ref window 

        # 1-based
        self.unspliced_local_reference_start = max(
            0, pos - max(150, window * self.k)
        ) + 1
        
        self.unspliced_local_reference_end = min(
            pos + self.window * self.k, self.reference_len
        )

    @property
    def is_indel(self) -> bool:
        """True for insertions or deletions. False otherwise."""
        return len(self.ref) != len(self.alt)

    @property
    def is_del(self) -> bool:
        """True for deletions. False for insertions or substitutions.   
        """
        return len(self.ref) > len(self.alt)

    @property
    def is_ins(self) -> bool:
        """True for insetions. False for deletions or substitutions."""
        return len(self.ref) < len(self.alt)

    @property
    def is_simple_indel(self) -> bool:
        """True for indel that is not complex (co-occurrence of insertion and deletion).
        False for complex indels or substitutions.
        """
        ref_len = len(self.ref)
        alt_len = len(self.alt)

        if ref_len == alt_len:
            return False

        shorter_allele = self.ref if ref_len < alt_len else self.alt
        longer_allele = self.ref if ref_len > alt_len else self.alt

        # test prefix by startswith
        return longer_allele.startswith(shorter_allele)

    
    @property
    def indel_seq(self):
        """returns the inserted/deleted sequence for non-complex indels. None for substitutions.
        """
        if self.is_ins:
            return self.alt[len(self.ref) :]
        elif self.is_del:
            return self.ref[len(self.alt) :]
        else:
            return ""

    @property
    def is_normalized(self) -> bool:
        """True if left-aligened and the allele representations are minimal."""
        ref, alt = self.ref, self.alt
        if ref[-1] != alt[-1]:
            if ref[0] != alt[0]:
                return True
            if len(ref) == 1 or len(alt) == 1:
                return True
        return False

    def left_flank(self, window=50, normalize=False):
        """extracts the left-flanking reference sequence. See also :meth:`~indelpost.Variant.right_flank`.

        Parameters
        ----------
        window : integer
            extract the reference sequence [variant_pos - window, variant_pos].
        normalize : bool
            if True, the normalized indel position is used as the end of the flanking sequence.
        """
        if normalize:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference, skip_validation=True)
        else:
            i = self

        if i.is_non_complex_indel():
            pos = i.pos
        else:
            pos = i.pos - 1

        lt_flank = i.reference.fetch(i.chrom, max(0, pos - window), pos)

        return lt_flank


    def right_flank(self, window=50, normalize=False):
        """extracts the right-flanking reference sequence. See also :meth:`~indelpost.Variant.left_flank`.

        Parameters
        ----------
        window : integer
            extract the reference sequence [variant_end_pos, variant_end_pos + window].
        normalize : bool
            if True, the normalized indel position is used as the start of the flanking sequence.
        """
        if normalize:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference, skip_validation=True)
        else:
            i = self

        ref_lim = i.reference.get_reference_length(i.chrom)
        if i.is_non_complex_indel() and i.is_ins:
            rt_flank = i.reference.fetch(i.chrom, i.pos, min(i.pos + window, ref_lim))
        else:
            if i.is_non_complex_indel() and i.is_del:
                event_len = len(i.indel_seq)
            else:
                event_len = len(i.ref) - 1
            rt_flank = i.reference.fetch(i.chrom, i.pos + event_len, min(i.pos + event_len + window, ref_lim))

        return rt_flank


    def normalize(self, inplace: bool = False) -> Optional['Variant']:
        """left-aligns in genomic position and minimize the allele represention.

        Parameters
        ---------
            inplace : bool
                normalizes this (self) :class:`~variantpost.Variant` object if True.
                Otherwise, returns a normalized copy of the object. Default to False.
        """
        i = self if inplace else Variant(self.chrom, self.pos, self.ref, self.alt, self.reference)
        
        # binding to local variables (avoid frequent property access)
        ref, alt, pos = i.ref, i.alt, i.pos

        condition_1 = ref[-1] == alt[-1] and ref[-1] != "N"
        if condition_1:
            lhs = i.reference.fetch(i.chrom, i.unspliced_local_reference_start, pos - 1)
            n = len(lhs) - 1

            while condition_1 and n >= 0:
                left_base = lhs[n]
                ref = left_base + ref[:-1]
                alt = left_base + alt[:-1]
                pos -= 1
                condition_1 = ref[-1] == alt[-1] and ref[-1] != "N"
                n -= 1

        condition_2 = ref[0] == alt[0]
        condition_3 = len(ref) > 1 and len(alt) > 1
        
        condition_2 = ref[0] == alt[0]
        condition_3 = len(ref) > 1 and len(alt) > 1

        while condition_2 and condition_3:
            ref = ref[1:]
            alt = alt[1:]
            pos += 1
            condition_2 = ref[0] == alt[0]
            condition_3 = len(ref) > 1 and len(alt) > 1

        i.ref, i.alt, i.pos = ref, alt, pos

        if inplace:
            return None
        return i

    def __eq__(self, other : Any) -> bool:
        if not isinstance(other, Variant):
            return False
        
        if self.chrom != other.chrom:
            return False

        if self.pos == other.pos and self.ref == other.ref and self.alt == other.alt:
            return True

        i = self.normalize()
        j = other.normalize()

        return i.pos == j.pos and i.ref == j.ref and i.alt == j.alt

    def __hash__(self) -> int:
        i = self.normalize()
        return hash((i.chrom, i.pos, i.ref, i.alt))

    def __getstate__(self) -> tuple:
        return (self.chrom, self.pos, self.ref, self.alt, self.reference.filename)

    def __setstate__(self, state) -> None:
        self.chrom = state[0]
        self.pos = state[1]
        self.ref = state[2]
        self.alt = state[3]
        self.reference = FastaFile(state[4])

    def __bool__(self) -> bool:
        return self.ref != self.alt

    def left_pos(self) -> int:
        """returns a left-aligned genomic position."""
        if not self.is_indel or not self.is_simple_indel:
            return self.pos

        curr_allele = self.alt if self.is_ins else self.ref
        curr_pos = self.pos
        
        if curr_allele[0] != curr_allele[-1]:
            return curr_pos

        lt_flank = self.reference.fetch(
            self.chrom, self.unspliced_local_reference_start - 1, curr_pos
        )
        
        idx = len(lt_flank) - 1
        n_shifted = 0
        
        while idx >= 0:
            tmp_allele = lt_flank[idx] + curr_allele[:-1]
            if tmp_allele[0] != tmp_allele[-1]:
                break
            curr_allele = tmp_allele
            idx -= 1
            n_shifted += 1

        return curr_pos - n_shifted

    def right_pos(self) -> int:
        """returns the variant-end position after right-aligned."""
        if not self.is_indel or not self.is_simple_indel:
            return self.pos + len(self.ref) - 1

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
            tmp_allele = curr_allele[1:] + rt_flank[i]
            i += 1
            if tmp_allele[0] != tmp_allele[-1]:
                break
            curr_allele = tmp_allele

        return curr_pos + i - 1
    
    
    def is_non_complex_indel(self):
        """returns True only if non-complex indel (False if complex indel or substitution).
        """
        i = self.normalize()
        ref, alt = i.ref, i.alt
        if len(ref) == len(alt):
            return False
        
        if ref[0] != alt[0]:
            return False
        
        the_shorter = ref if i.is_ins else alt
        if len(the_shorter) > 1:
            return False
        
        return True    

    def count_repeats(self, by_repeat_unit=True):
        """counts indel repeats in the flanking reference sequences. The search window is
        defined by :meth:`~indelpost.Variant.left_flank` and :meth:`~indelpost.Variant.right_flank`.
        
        Parameters
        ----------
        by_repeat_unit : bool
            count by the smallest tandem repeat unit. For example, the indel sequence "ATATATAT" has
            tandem units "ATAT" and "AT". The occurrence of "AT" will be counted if True (default).
        """ 

        if self.is_non_complex_indel():
            seq = self.indel_seq
        else:
            seq = self.alt
        
        if by_repeat_unit:
            seq = to_minimal_repeat_unit(seq)
        
        lt_flank = self.left_flank()
        lt_repeat = repeat_counter(seq, lt_flank[::-1])
        
        rt_flank = self.indel_seq[len(seq):] 
        rt_flank += self.right_flank()
        rt_repeat = repeat_counter(seq, rt_flank)
        
        return lt_repeat + rt_repeat

    def query_vcf(self, vcf: VariantFile, chrom_name: Optional[str] = None, match_by_equivalence: bool = True) -> List[MatchedRecord]:
        """returns a `list <https://docs.python.org/3/library/stdtypes.html#list>`__ of :class:`MatchedRecord`.
        

        Parameters
        ----------
            vcf : pysam.VariantFile
                VCF file to be queried.
                Supply as `pysam.VariantFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantFile>`__ object. The VCF file muet be indexed.

            chrom_name : string
                specify an alias chromosome name if the VCF file uses a chromosome nomenclature different from the reference genome in :class:`~variantpost.Variant`.
                If not specified (default), the nomenclature in the FASTA file in :class:`~variantpost.Variant` will be used.

            match_by_equivalence : bool
                queries the VCF records by normalization if True (default). Otherwise, positionally overlapping records will be returned.
        
        
        :class:`MatchedRecord` is a `namedtuple <https://docs.python.org/3/library/collections.html#collections.namedtuple>`__ 
        with the following fields

        - **chrom** - VCF CHROM field. 
        - **pos** - VCF POS field.
        - **id** - VCF ID field.
        - **ref** - VCF REF field.
        - **alts** - VCF ALT field as `tuple <https://docs.python.org/3/library/stdtypes.html#tuples>`__. May contain multiple alleles.
        - **qual** - VCF QUAL field.
        - **filter** - VCF FILTER field.
        - **info** - VCF INFO field. Values are accessible by keys defined in the header.
        - **format** - VCF FORMAT field.
        - **samples** - VCF genotype field. Genotypes are accessible by using sample names as key.
        
        """

        vcf_chrom = chrom_name or self.chrom
        lt_pos, rt_pos = self.left_pos(), self.right_pos()

        if match_by_equivalence:
            search_start, search_end = lt_pos - 1, rt_pos
            self_norm = self.normalize() 
        else:
            search_start, search_end = (
                self.unspliced_local_reference_start - 1,
                self.unspliced_local_reference_end,
            )

        vcf_entries = vcf.fetch(vcf_chrom, search_start, search_end)
        hits = []

        for _entry in vcf_entries:
            if not _entry.alts:
                continue

            ref_len = len(_entry.ref)

            for _alt in _entry.alts:
                if not set(_alt).issubset(BASESET):
                    continue

                if match_by_equivalence:
                    if self.pos == _entry.pos and self.ref == _entry.ref and self.alt == _alt:
                        hits.append(to_tuple(_entry))
                        break

                    tmp_variant = Variant(self.chrom, _entry.pos, _entry.ref, _alt, self.reference)
                    tmp_norm = tmp_variant.normalize()
                    
                    if (self_norm.pos == tmp_norm.pos and 
                        self_norm.ref == tmp_norm.ref and 
                        self_norm.alt == tmp_norm.alt):
                        hits.append(to_tuple(_entry))
                        break
                else:
                    if ref_len == len(_alt):
                        if lt_pos <= _entry.pos <= rt_pos:
                            hits.append(to_tuple(_entry))
                            break
                    else:
                        v = Variant(self.chrom, _entry.pos, _entry.ref, _alt, self.reference)
                        _lt_pos, _rt_pos = v.left_pos(), v.right_pos()

                        if lt_pos <= _lt_pos <= rt_pos or lt_pos <= _rt_pos <= rt_pos:
                            hits.append(to_tuple(_entry))
                            break

        return hits

def to_minimal_repeat_unit(seq):
    """Find repeat unit in indel sequence
    """
    mid = int(len(seq) / 2)
    min_unit = seq
    
    j = 1
    found = False
    
    while j <= mid and not found:
        tandems = [seq[i : i + j] for i in range(0, len(seq), j)]
        if len(set(tandems)) == 1:
            found = True
            min_unit = list(set(tandems))[0]
        j += 1
    
    return min_unit

def repeat_counter(query_seq, flank_seq):
    """
    """
    qlen, flen = len(query_seq), len(flank_seq)
    count = 0
    
    if flen < qlen:
        return count
    try:
        for i in range(0, flen, qlen):
            if flank_seq[i  : i + qlen] == query_seq:
                count += 1
            else:
                break
    except:
        return -1
    
    return count


def to_tuple(_entry):
    MatchedRecord = namedtuple(
        "MatchedRecord",
        [
            "chrom",
            "pos",
            "id",
            "ref",
            "alts",
            "qual",
            "filter",
            "info",
            "format",
            "samples",
        ],
    )

    return MatchedRecord(
        _entry.chrom,
        _entry.pos,
        _entry.id,
        _entry.ref,
        _entry.alts,
        _entry.qual,
        _entry.filter,
        _entry.info,
        _entry.format,
        _entry.samples,
    )
