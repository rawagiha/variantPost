from collections import namedtuple


from .variant import Variant
from .phaser import _phase, loss
from variantpost.__search import search_target


class VariantAlignment(object):
    """This class accepts the target variant as :class:`~variantpost.Variant` and the BAM file
        as `pysam.AlignmentFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile>`__
        to process the variant alignment.


    Parameters
    ----------
    variant : Variant
         :class:`~variantpost.Variant` object representing the target variant.

    bam : pysam.AlignmentFile
        BAM file supplied as
        `pysam.AlignmentFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile>`__ object.

    second_bam : pysam.AlignmentFile
        A second BAM file for paired analysis. Default: None.

    chrom_name : string
        Specify an alias chromosome name if the BAM file uses a chromosome nomenclature different from
        the reference used in :class:`~variantpost.Variant`. If not specified (default), the nomenclature
        in the reference will be used.

    mapping_quality_threshold : integer
        A mininum mapping quality to be analized. Default 1.

    base_quality_threshold : integer
        Non-reference base-calls with a Phred-scale quality score below the threshold are labeled low quality.

    low_quality_base_rate_threshold : float

    downsample_threshold : integer

    match_score : integer

    mismatch_penalty : integer

    gap_open_penalty : integer

    gap_extension_penalty : integer

    kmer_size : integer

    local_threshold : integer

    match_penalty_for_phasing : float

    """

    def __init__(
        self,
        variant,
        bam,
        second_bam=None,
        chrom_name=None,
        exclude_duplicates=True,
        mapping_quality_threshold=1,
        base_quality_threshold=30,
        low_quality_base_rate_threshold=0.1,
        downsample_threshold=5000,
        match_score=3,
        mismatch_penalty=2,
        gap_open_penalty=3,
        gap_extension_penalty=1,
        kmer_size=32,
        local_threshold=20,
        match_penalty_for_phasing=0.5,
    ):
        if not variant.is_normalized:
            variant.normalize(inplace=True)

        self.variant = variant
        self.chrom = variant.chrom
        self.bam_chrom = chrom_name if chrom_name else variant.chrom
        self.target_pos = variant.pos
        self.target_is_indel = variant.is_indel
        self.window = variant.window
        self.reference = variant.reference
        self.local_thresh = local_threshold
        self.has_second = second_bam

        retarget_thresh = _retarget_thresh(
            match_penalty_for_phasing, local_threshold, self.window
        )

        # interact with c++ code
        (
            self.contig_dict,
            self.skips,
            self.read_names,
            self.are_reverse,
            self.target_status,
            self.are_first_bam,
            self.is_retargeted,
        ) = search_target(
            bam,
            second_bam,
            variant.reference_len,
            exclude_duplicates,
            self.window,
            variant.reference.filename,
            self.chrom,
            self.bam_chrom,
            variant.pos,
            variant.ref.encode(),
            variant.alt.encode(),
            mapping_quality_threshold,
            base_quality_threshold,
            low_quality_base_rate_threshold,
            downsample_threshold,
            match_score,
            mismatch_penalty,
            gap_open_penalty,
            gap_extension_penalty,
            kmer_size,
            local_threshold,
            retarget_thresh,
            variant.unspliced_local_reference_start,
            variant.unspliced_local_reference_end,
        )

        self.is_with_target = any([status == 1 for status in self.target_status])

    def count_alleles(self):
        """returns :class:`AlleleCount` as `namedtuple <https://docs.python.org/3/library/collections.html#collections.namedtuple>`__ of read counts.
       :class:`AlleleCount` has the following fields accessible by attribute:
       
        
        - "s": count of read names supporting the variant.
        - "n": count of read names not supporting the variant. 
        - "u": count of read names undetermined to be supporting/non-supporting

        Strand breakdowns are also available by:
        
        - "s_fw": count of forward reads supportingt the variant.
        - "s_rv": count of reverse reads supportingt the variant. 
        - ...

        For paired analysis, :class:`PairedAlleleCount` is returnd and has the following fields: 

        - "first": :class:`AlleleCount` for the first BAM file.
        - "second": :class:`AlleleCount` for the second BAM file.

        """
        if self.has_second:
            return self._paired_count()
        else:
            return self._unpaired_count()

    def _unpaired_count(self):
        sf, sr, nf, nr, uf, ur = ([] for i in range(6))

        for read_name, status, is_rv in zip(
            self.read_names, self.target_status, self.are_reverse
        ):

            flags = (status, is_rv)

            if flags == (1, True):
                sr.append(read_name)
            elif flags == (1, False):
                sf.append(read_name)
            elif flags == (0, True):
                nr.append(read_name)
            elif flags == (0, False):
                nf.append(read_name)
            elif flags == (-1, True):
                ur.append(read_name)
            elif flags == (-1, False):
                uf.append(read_name)

        return fill_cnt_data(sf, sr, nf, nr, uf, ur)

    def _paired_count(self):

        sf1, sf2, sr1, sr2, nf1, nf2, nr1, nr2, uf1, uf2, ur1, ur2 = (
            [] for i in range(12)
        )
        for read_name, status, is_rv, is_first in zip(
            self.read_names, self.target_status, self.are_reverse, self.are_first_bam
        ):

            flags = (status, is_rv, is_first)

            if flags == (1, True, True):
                sr1.append(read_name)
            elif flags == (1, True, False):
                sr2.append(read_name)
            elif flags == (1, False, True):
                sf1.append(read_name)
            elif flags == (1, False, False):
                sf2.append(read_name)
            elif flags == (0, True, True):
                nr1.append(read_name)
            elif flags == (0, True, False):
                nr2.append(read_name)
            elif flags == (0, False, True):
                nf1.append(read_name)
            elif flags == (0, False, False):
                nf2.append(read_name)
            elif flags == (-1, True, True):
                ur1.append(read_name)
            elif flags == (-1, True, False):
                ur2.append(read_name)
            elif flags == (-1, False, True):
                uf1.append(read_name)
            elif flags == (-1, False, False):
                uf2.append(read_name)

        PairedAlleleCount = namedtuple("PairedAlleleCount", ["first", "second"])

        pac = PairedAlleleCount(
            fill_cnt_data(sf1, sr1, nf1, nr1, uf1, ur1),
            fill_cnt_data(sf2, sr2, nf2, nr2, uf2, ur2),
        )

        return pac

    def phase(self, match_penalty_for_phasing=0.5, max_common_substr_len=15):
        if self.is_retargeted:
            self.target_is_indel = True

        try:
            phased = _phase(
                self.contig_dict,
                self.skips,
                self.target_pos,
                self.target_is_indel,
                self.local_thresh,
                match_penalty_for_phasing,
                max_common_substr_len,
            )

            if phased and self.is_with_target:
                return Variant(
                    self.chrom, phased[0], phased[1], phased[2], self.reference
                ).normalize()
            elif self.is_with_target:
                return self.variant
            else:
                ref_base = self.reference.fetch(
                    self.chrom, self.target_pos - 1, self.target_pos
                )
                return Variant(
                    self.chrom, self.target_pos, ref_base, ref_base, self.reference
                )
        except:
            if self.is_with_target:
                return self.variant
            else:
                ref_base = self.reference.fetch(
                    self.chrom, self.target_pos - 1, self.target_pos
                )
                return Variant(
                    self.chrom, self.target_pos, ref_base, ref_base, self.reference
                )


def fill_cnt_data(sf, sr, nf, nr, uf, ur):
    AlleleCount = namedtuple(
        "AlleleCount", ["s", "s_fw", "s_rv", "n", "n_fw", "n_rv", "u", "u_fw", "u_rv"]
    )
    ac = AlleleCount(
        len(set(sf + sr)),
        len(sf),
        len(sr),
        len(set(nf + nr)),
        len(nf),
        len(nr),
        len(set(uf + ur)),
        len(uf),
        len(ur),
    )

    return ac


def _retarget_thresh(match_penal, local_thresh, window):
    score = 0
    for i in range(window):
        score += loss(i, match_penal, local_thresh)
        if score < -1.0:
            return i
