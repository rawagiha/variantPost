# indelinside
<p align="left">
  <img src="./indelinside_logo.svg" alt="indel inside" width="150">
</p>

indelinside is a command-line tool to reanalyze somatic indels on locally personalized genome for indel signature analysis.
The algorithm will:
* harmonize indel representations that may be different across multiple variant callers.
* construct local diploid haplotypes from normal BAM file along with a haplotype carrying the target somatic indel.  
* infer from which germline haplotyes the target indel is derived.    
* realign the target indel haplotype to the inferred germline haplotype for personalization.




variantPost supports tumor/normal-paired analyis for cancer genomics applications.

Visit [documentation](https://variantpost.readthedocs.io/en/latest/) for detail.

### Installation
variantPost requires a Linux machine with a gcc compiler supporting for [c++17](https://en.cppreference.com/w/cpp/17).

To install
```
pip install git+https://github.com/rawagiha/variantPost
```

Upon installation, [cython](https://cython.org/) and [pysam](https://github.com/pysam-developers)
will also be installed if not pre-installed. 

### Usage
[documentation](https://variantpost.readthedocs.io/en/latest/)

### Acknowledgements
variantPost internally uses the following packages. I thank the developers for making them freely available. 
- [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
- [fastahack](https://github.com/ekg/fastahack)
- [Ratcliff-Obershelp algorithm](https://github.com/wernsey/miscsrc)
