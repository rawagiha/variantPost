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

### Core Workflow

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

### Installation
Available as a part of variantPost library:
```
pip install git+https://github.com/rawagiha/variantPost
```

---

### Basic Command

Multiple VCF files from different variant callers can be passed directly to the `-v` / `--vcf` option:

```bash
indelinside personalize \
  -t tumor.wgs.bam \
  -n normal.wgs.bam \
  -r reference.fa \
  -v caller_a.vcf caller_b.vcf caller_c.vcf \
  -o output.txt

### Acknowledgements
variantPost internally uses the following packages. I thank the developers for making them freely available. 
- [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
- [fastahack](https://github.com/ekg/fastahack)
