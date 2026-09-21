<p align="center">
  <img src="./indelinside_logo.svg" alt="indelinside logo" width="180">
</p>

<h1 align="center">indelinside</h1>

<p align="center">
  <b>Reanalyzing Somatic Indels on Locally Personalized Genomes for Indel Signature Analysis</b>
</p>

<br>

**`indelinside`** is a command-line tool that eliminates mapping ambiguities and alignment noise when profiling somatic indel signatures by realigning variants against a **locally personalized germline background**.

---

### 💡 Core Workflow

1. **Harmonization**
   Standardizes inconsistent indel representations across different variant callers into a unified canonical form.

2. **Local Haplotype Assembly**
   Constructs local diploid germline haplotypes (**`hap1`** & **`hap2`**) from the normal BAM file, alongside the somatic indel haplotype (**`hap0`**).

3. **Origin Inference**
   Determines precisely which germline haplotype (**`hap1`** or **`hap2`**) the somatic indel was derived from.

4. **Personalized Realignment**
   Realigns the somatic haplotype (**`hap0`**) to its inferred parent germline background—delivering ultra-clean inputs for downstream signature matrices.

<br>

<p align="center">
  <img src="./fig_repo.svg" alt="indelinside algorithm workflow" width="720">
</p>

<br>

# indelinside
<p align="left">
  <img src="./indelinside_logo.svg" alt="indel inside" width="150">
</p>

indelinside is a command-line tool to reanalyze somatic indels on locally personalized genome for indel signature analysis.
The algorithm will:
* harmonize indel representations that may be different across multiple variant callers.
* construct local diploid haplotypes (hap1 & hap2) from normal BAM file along with a haplotype carrying the target somatic indel (hap0).  
* infer from which germline haplotyes the target indel is derived.    
* realign the target indel haplotype to the inferred germline haplotype for personalization.

<p align="center">
    <img src="./fig_repo.svg" alt="algorithm" width="700">
</p>


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
