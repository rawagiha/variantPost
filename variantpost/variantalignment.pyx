from .phaser import phase
from .preprocessor import preprocess
from variantpost.processor_wrapper cimport process_pileup

import time

class VariantAlignment(object):
    def __init__(
        self, variant, bam, 
        second_bam=None, 
        window=200, 
        exclude_duplicates=True, 
        downsample_threshold=-1, 
        mapping_quality_threshold=-1, 
        base_quality_threshold=20, 
        low_quality_base_rate_threshold=0.05, 
        match_score=3,
        mismatch_penalty=2,
        gap_open_penalty=3,
        gap_extention_penalty=1, 
        kmer_size=16, 
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
        preprocessed_pileup = preprocess(
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

        # interact with c++ code
        (
            contig_dict,
            skips,
            annotated_reads,
        ) = process_pileup(
                fastafile,
                chrom.encode(), 
                pos, ref.encode(), alt.encode(),
                mapping_quality_threshold,
                base_quality_threshold, 
                low_quality_base_rate_threshold,
                match_score,
                mismatch_penalty,
                gap_open_penalty,
                gap_extention_penalty,
                kmer_size,
                unspliced_local_reference_start, 
                unspliced_local_reference_end, 
                preprocessed_pileup.read_names,
                preprocessed_pileup.are_reverse,
                preprocessed_pileup.cigar_strings,
                preprocessed_pileup.aln_starts,
                preprocessed_pileup.aln_ends,
                preprocessed_pileup.read_seqs,
                preprocessed_pileup.qual_seqs,
                preprocessed_pileup.mapqs,
                preprocessed_pileup.are_first_bam
        )

        #for a in annotated_reads:
        #    if not a.target_status:
        #        print(a.read_name, "non-tar")
        #    if a.target_status < 0:
        #        print(a.read_name, "undeter")
        
        print(time.time() - tt, "c++ time")
        
        phase(contig_dict, skips,  pos)
        
        print(time.time() - t, "total varaln --- {}".format("aho"))


