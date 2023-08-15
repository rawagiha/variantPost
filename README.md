variantPost supports short-read analysis for SNVs, MNVs, and indels.  


### Installation
variantPost depends on [cython](https://cython.org/) and [c++17](https://en.cppreference.com/w/cpp/17).
Please build gcc 7 or higher (tested on gcc 9.1.0).


### Usage
The API works with [pysam] (https://github.com/pysam-developers)
usage is similar to indelPost

```
import pysam
from variantpost import Variant, VariantAlignment

reference = pysam.FastaFile("/path/to/ref.fa")

bam_1 = pysam.AlignmentFile("/path/to/bam_1.bam")
bam_2 = pysam.AlignmentFile("/path/to/bam_1.bam")

v = Variant("chrN", 123, "A", "ATCG", reference)

# unpaired
valn = VariantAlignment(v, bam_1)

# paired 
valn = VariantAlignment(v, bam_1, bam_2)

# allele count

cnt = valn.count_alleles()

# if unpaired, cnt is an AlleleCount obj
print(cnt.s) # supporting reads fw and rv (unique count) 
print(cnt.s_fw, cnt.s_rv) # fw/rv breakdown
print(cnt.n) # non-supporting reads fw and rv 
print(cnt.n_fw, cnt.n_rv) # fw/rv breakdown
print(cnt.u) # undetermined reads (low qual etc...)
print(cnt.u_fw, cnt.u_rv) # fw/rv breakdown

# if paired, cnt is a PairedAlleleCount obj, which is like (AlleleCount obj for 1, AlleleCount obj for 2)

cnt1 = cnt.first # AlleleCount obj for bam_1
cnt2 = cnt.second # AlleleCount obj for bam_2


# phasing

phased = valn.phase() 
```

Please report anything of concern. 
