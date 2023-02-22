from .preprocessor import preprocess
from variantpost.processor_wrapper cimport test_it

import time

class VariantAlignment(object):
    def __init__(
        self, variant, bam, 
        second_bam=None, 
        window=200, 
        exclude_duplicates=True, 
        downsample_threshold=-1, 
        mapping_quality_threshold=1, 
        base_quality_threshold=20, 
        low_quality_base_rate_threshold=0.05, 
        kmer_size=32, 
    ):

        (
            chrom,
            pos,
            ref,
            alt,
            chrom_len,
            unspliced_local_reference,
            unspliced_local_reference_start,
            unspliced_local_reference_end,
            reference,
        ) = (
            variant.chrom,
            variant.pos,
            variant.ref,
            variant.alt,
            variant.reference_len,
            variant.unspliced_local_reference,
            variant.unspliced_local_reference_start,
            variant.unspliced_local_reference_end,            
            variant.reference,
        )
        
        t = time.time()

        ##########
        fastafile = variant.reference.filename
        #########
        
        #processed reads for c++ wrapper
        reads = preprocess(
            chrom,
            pos,
            chrom_len,
            bam,
            second_bam,
            unspliced_local_reference,
            unspliced_local_reference_start,
            reference,
            exclude_duplicates,
            window,
            downsample_threshold,
        )
       
        tt = time.time()
        print(tt -t, "preprosess")

        res = test_it(fastafile,
                      chrom.encode(), 
                      pos, ref.encode(), alt.encode(),
                      mapping_quality_threshold,
                      base_quality_threshold, 
                      low_quality_base_rate_threshold,
                      kmer_size,
                      unspliced_local_reference_start, 
                      unspliced_local_reference_end, 
                      #unspliced_local_reference.encode("utf-8"), 
                      reads[0], reads[1], reads[2], reads[3], reads[4], reads[5], reads[6], reads[7], reads[8]) #refseq removed

        print(time.time() - tt, "c++ time")
        print(time.time() - t, "total varaln --- {}".format(res))
