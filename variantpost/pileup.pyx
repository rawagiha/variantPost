def fetch_reads(bam, chrom, pos):
    return bam.fetch(chro, pos - 10, pos + 10)
