from collections import namedtuple
from pysam import FastaFile

BASES = {"A", "C", "G", "T", "N"}
BASESET = set(BASES)

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

    def __init__(self, chrom, pos, ref, alt, reference, window=50):

        self.chrom = chrom
        self.pos = pos  # 1-based

        self.ref = ref
        self.alt = alt

        if not set(list(self.ref)) <= BASES or not set(list(self.alt)) <= BASES:
            raise ValueError("Found characters other than 'A', 'C', 'G', 'T', 'N'.")

        self.reference = reference

        ref_len = len(ref)
        if ref_len < window:
            self.window = window
        else:
            self.window = ref_len + window

        self.reference_len = reference.get_reference_length(self.chrom)

        self.k = 1 if ref_len > window * 4  else 4

        # 1-based
        self.unspliced_local_reference_start = max(
            0, pos - max(150, window * self.k)
        ) + 1
        
        self.unspliced_local_reference_end = min(
            pos + self.window * self.k, self.reference_len
        )

    @property
    def is_indel(self):
        """True for insertions or deletions. False otherwise."""
        return len(self.ref) != len(self.alt)

    @property
    def is_del(self):
        """True for deletions. False for insertions or substitutions.   
        """
        return len(self.ref) > len(self.alt)

    @property
    def is_ins(self):
        """True for insetions. False for deletions or substitutions."""
        return len(self.ref) < len(self.alt)

    @property
    def is_simple_indel(self):
        """True for indel that is not complex (co-occurrence of insertion and deletion).
        False for complex indels or substitutions.
        """
        if self.is_indel:
            shorter_allele = self.ref if self.is_ins else self.alt
            longer_allele = self.ref if self.is_del else self.alt
            return shorter_allele == longer_allele[: len(shorter_allele)]

        return False

    @property
    def is_normalized(self):
        """True if left-aligened and the allele representations are minimal."""
        if self.ref[-1] != self.alt[-1]:
            if self.ref[0] != self.alt[0]:
                return True
            elif len(self.ref) == 1 or len(self.alt) == 1:
                return True

        return False

    def normalize(self, inplace=False):
        """left-aligns in genomic position and minimize the allele represention.

        Parameters
        ---------
            inplace : bool
                normalizes this (self) :class:`~variantpost.Variant` object if True.
                Otherwise, returns a normalized copy of the object. Default to False.
        """
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

    def left_pos(self):
        """returns a left-aligned genomic position."""

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

                    if tmp_allele[0] != tmp_allele[-1]:
                        break
                    curr_allele = tmp_allele

                return curr_pos - i + 1
            else:
                return self.pos

    def right_pos(self):
        """returns the variant-end position after right-aligned."""
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

    def query_vcf(self, vcf, chrom_name=None, match_by_equivalence=True):
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
        vcf_chrom = chrom_name if chrom_name else self.chrom

        lt_pos, rt_pos = self.left_pos(), self.right_pos()

        if match_by_equivalence:
            search_start, search_end = lt_pos - 1, rt_pos
        else:
            search_start, search_end = (
                self.unspliced_local_reference_start - 1,
                self.unspliced_local_reference_end,
            )

        vcf_entries = vcf.fetch(vcf_chrom, search_start, search_end)

        hits = []
        if match_by_equivalence:
            for _entry in vcf_entries:
                
                # may be empty (ALT = .)
                if not _entry.alts:
                    continue

                for _alt in _entry.alts:
                    
                    _bases = set(list(_alt))
                    
                    # may contain other than ACTGN
                    if not _bases <= BASESET:
                        continue;

                    if self == Variant(
                        self.chrom, _entry.pos, _entry.ref, _alt, self.reference
                    ):
                        hits.append(to_tuple(_entry))
                        break
        else:
            for _entry in vcf_entries:
                ref_len = len(_entry.ref)
                
                # may be empty (ALT = .)
                if not _entry.alts:
                    continue                
                
                for _alt in _entry.alts:
                    
                    _bases = set(list(_alt))
                    
                    # may contain other than ACTGN
                    if not _bases <= BASESET:
                        continue;
                     
                    if ref_len == len(_alt):
                        if lt_pos <= _entry.pos <= rt_pos:
                            hits.append(to_tuple(_entry))
                            break
                    else:
                        v = Variant(
                            self.chrom, _entry.pos, _entry.ref, _alt, self.reference
                        )
                        _lt_pos, _rt_pos = v.left_pos(), v.right_pos()

                        if lt_pos <= _lt_pos <= rt_pos or lt_pos <= _rt_pos <= rt_pos:
                            hits.append(to_tuple(_entry))
                            break

        return hits


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
