<p align="center">
  <img src="./indelinside_logo.svg" alt="indelinside logo" width="180">
</p>

<h1 align="center">indelinside</h1>

<p align="center">
  <b>Reanalyzing Somatic Indels on Locally Personalized Genomes for Indel Signature Analysis</b>
</p>

<br>

**`indelinside`** is a command-line tool that resolves mapping ambiguities for indel signature analysis by realigning variants against a **locally personalized germline background**.

---

### Core Workflow

1. **Harmonization**
   Unifies inconsistent indel representations across different variant callers.

2. **Local Haplotype Assembly**
   Constructs local diploid germline haplotypes (**`hap1`** & **`hap2`**) from the normal BAM file, alongside the somatic indel haplotype (**`hap0`**) from the tumor BAM file.

3. **Origin Inference**
   Determines which germline haplotype (**`hap1`** or **`hap2`**) the somatic indel was derived from.

4. **Personalized Realignment**
   Realigns the somatic haplotype (**`hap0`**) to its inferred parent germline background. Depending on the context, the indel may be personalized to a substitute, a different indel class, or a complex indel—delivering a polished input for downstream signature analysis.

<br>

<p align="center">
  <img src="./fig_repo.svg" alt="indelinside algorithm workflow" width="720">
</p>

<br>

---

### Installation

`indelinside` is distributed as a command-line utility within the `variantPost` package:

```bash
pip install git+[https://github.com/rawagiha/variantPost](https://github.com/rawagiha/variantPost)
```

---

### Usage

The pipeline consists of two steps: local personalization of somatic indels followed by feature extraction into an indel signature matrix.

#### Step 1: Personalize Alignments

Input BAM files and arbitrary numbers of VCFs from different variant callers into `personalize`. Variants are automatically consolidated, harmonized, and realigned against inferred germline haplotypes.

```bash
indelinside personalize \
  -t tumor.wgs.bam \
  -n normal.wgs.bam \
  -r reference.fa \
  -v caller_a.vcf caller_b.vcf caller_c.vcf \
  -o indelinside.out.txt
```

#### Step 2: Generate Indel Signature Matrix

Generate a COSMIC-compatible indel signature matrix from the personalized output. Filter consensus calls across callers using `-c` / `--consensus_level`. `-c N` selects indels called by `N` or more callers (`-c 1` is the union of all callers). 

```bash
indelinside matrix \
  -i indelinside.out.txt \
  -c 2 \
  --sample_name my_sample
```

The output matrix can be used as input for external signature analysis tools such as [SigProfilerAssignment](https://github.com/SigProfilerSuite/SigProfilerAssignment).  

```python
from SigProfilerAssignment import Analyzer

Analyzer.cosmic_fit(
    samples="my_sample.indel_83_matrix.txt",  # COSMIC-compatible 83-indel channel matrix  
    output=/path/to/output_dir,
    input_type="matrix", 
    context_type="ID",  
    collapse_to_SBS96=False,
    signature_database=/path/to/COSMIC_signature_database # restrict IDs as needed
)
```

<br>

### Acknowledgements
variantPost internally uses the following packages. I thank the developers for making them freely available. 
- [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
- [fastahack](https://github.com/ekg/fastahack)
