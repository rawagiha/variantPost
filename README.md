variantPost supports short-read analysis for SNVs, MNVs, and indels.  
(under development)

### Installation
variantPost depends on [cython](https://cython.org/) and [c++17](https://en.cppreference.com/w/cpp/17).
Please build with gcc 7 or higher.


### Usage
The API works with [pysam](https://github.com/pysam-developers)
usage is similar to indelPost

[documentation](https://variantpost.readthedocs.io/en/latest/)

### Acknowledgements
variantPost internally uses the following packages.
- [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
- [fastahack](https://github.com/ekg/fastahack)
