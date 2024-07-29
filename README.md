# variantPost
(under development)

variantPost is a Python library for short variant processing via realignment and read-based phasing to resolve alignment ambiguities.
By importing the library, users write their own scripts to solve alignment-sentive problems such as:
* compare SNVs, MNVs, and indels that may differently be called by multiple variant callers (e.g., complex indels).
* compare short variant alignments in multiple mappings (e.g., match DNA variants to RNA-Seq to check expression, where the DNA/RNA alignments may be different).  
* construct a complex indel or MNV from a simple short variant by read-based phasing.    
* count reads supporting the target variant from BAM file by realignment.
* pull variant records matching the target varitnt from VCF file.

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
