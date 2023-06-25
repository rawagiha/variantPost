from .phaser import phase
from variantpost.cy_search cimport search_target
from .variant import Variant

#import time

class VariantAlignment(object):
    def __init__(
        self, variant, bam, 
        second_bam=None, 
        window=200, 
        exclude_duplicates=True, 
        downsample_threshold=-1, 
        mapping_quality_threshold=-1, 
        base_quality_threshold=20, 
        low_quality_base_rate_threshold=0.1, 
        match_score=3,
        mismatch_penalty=2,
        gap_open_penalty=3,
        gap_extention_penalty=1, 
        kmer_size=24, 
        local_threshold=20
    ):
        if not variant.is_normalized:
            variant.normalize(inplace=True)
            variant = Variant(
                        variant.chrom, 
                        variant.pos, 
                        variant.ref, 
                        variant.alt, 
                        variant.reference
                    )

        (
            chrom,
            pos,
            ref,
            alt,
            chrom_len,
            unspliced_local_reference_start,
            unspliced_local_reference_end,
            reference,
        ) = (
            variant.chrom,
            variant.pos,
            variant.ref,
            variant.alt,
            variant.reference_len,
            variant.unspliced_local_reference_start,
            variant.unspliced_local_reference_end,            
            variant.reference,
        )
        
        #t = time.time()

        ##########
        fastafile = variant.reference.filename
        #########
        
        ##tt = time.time()
        #print(tt -t, "preprosess")

        # interact with c++ code
        self.contig_dict, self.skips, self.annotated_reads = search_target(
                bam,
                second_bam,
                chrom_len,
                exclude_duplicates,
                window,
                downsample_threshold,
                fastafile,
                chrom, 
                pos, ref.encode(), alt.encode(),
                mapping_quality_threshold,
                base_quality_threshold, 
                low_quality_base_rate_threshold,
                match_score,
                mismatch_penalty,
                gap_open_penalty,
                gap_extention_penalty,
                kmer_size,
                local_threshold,
                unspliced_local_reference_start, 
                unspliced_local_reference_end, 
        )

        
    def count_alleles(self):
            sup, non_sup, fwrv, fs = [], [], [], []
            for i in self.annotated_reads:
                if i.target_status == 1:
                    sup.append(i.read_name)
                elif i.target_status == 0:
                    non_sup.append(i.read_name)
                    
            return (len(sup), len(non_sup))    
        
        #print(time.time() - tt, "c++ time")
        
        #print(annotated_reads)
        #phase(contig_dict, skips,  pos)
        
        #print(time.time() - t, "total varaln --- {}".format("aho"))


