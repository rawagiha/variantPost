from .preprocessor import preprocess
from variantpost.pileup_parser_wrapper import test_it 

class VariantAlignment(object):
    def __init__(
        self, variant, bam, window=200, exclude_duplicates=True, downsample_thresh=-1
    ):

        (
            chrom,
            pos,
            chrom_len,
            unspliced_local_reference,
            unspliced_local_reference_start,
            reference,
        ) = (
            variant.chrom,
            variant.pos,
            variant.reference_len,
            variant.unspliced_local_reference,
            variant.unspliced_local_reference_start,
            variant.reference,
        )
        
        #processed reads for c++ wrapper
        reads = preprocess(
            chrom,
            pos,
            chrom_len,
            bam,
            unspliced_local_reference,
            unspliced_local_reference_start,
            reference,
            exclude_duplicates,
            window,
            downsample_thresh,
        )

        test_it(unspliced_local_reference_start, unspliced_local_reference.encode("utf-8"), * reads)
