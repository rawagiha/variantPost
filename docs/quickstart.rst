Quick Start
===========

Overview
--------

The indelinside workflow consists of two main steps:

1. Personalize somatic indels using local germline haplotypes.
2. Generate a COSMIC-compatible indel signature matrix.

Step 1: Personalize Somatic Indels
----------------------------------

Run the ``personalize`` command:

.. code-block:: bash

   indelinside personalize \
       -t tumor.wgs.bam \
       -n normal.wgs.bam \
       -r reference.fa \
       -v caller_a.vcf caller_b.vcf caller_c.vcf \
       -o indelinside.out.txt \
       -p 12

Arguments:

* ``-t``: Tumor BAM file.
* ``-n``: Normal BAM file.
* ``-r``: Reference genome FASTA file.
* ``-v``: One or more somatic VCF files.
* ``-o``: Output file.
* ``-p``: Number of processing threads.

The pipeline consolidates variants from multiple callers,
harmonizes their representations, and performs
local haplotype-based personalization.

Step 2: Generate an Indel Signature Matrix
------------------------------------------

Use the output from the personalization step:

.. code-block:: bash

   indelinside matrix \
       -i indelinside.out.txt \
       -c 2 \
       --sample_name my_sample

The ``-c`` option specifies the consensus level.

For example, ``-c 2`` selects indels supported by at
least two variant callers.

The output is a COSMIC-compatible 83-channel indel
signature matrix.

Downstream Signature Analysis
-----------------------------

The resulting matrix can be used with external
mutational signature analysis tools, such as
SigProfilerAssignment.

For advanced control over individual processing steps,
see :doc:`python_api/index`.Quick Start
===========

Overview
--------

The indelinside workflow consists of two main steps:

1. Personalize somatic indels using local germline haplotypes.
2. Generate a COSMIC-compatible indel signature matrix.

Step 1: Personalize Somatic Indels
----------------------------------

Run the ``personalize`` command:

.. code-block:: bash

   indelinside personalize \
       -t tumor.wgs.bam \
       -n normal.wgs.bam \
       -r reference.fa \
       -v caller_a.vcf caller_b.vcf caller_c.vcf \
       -o indelinside.out.txt \
       -p 12

Arguments:

* ``-t``: Tumor BAM file.
* ``-n``: Normal BAM file.
* ``-r``: Reference genome FASTA file.
* ``-v``: One or more somatic VCF files.
* ``-o``: Output file.
* ``-p``: Number of processing threads.

The pipeline consolidates variants from multiple callers,
harmonizes their representations, and performs
local haplotype-based personalization.

Step 2: Generate an Indel Signature Matrix
------------------------------------------

Use the output from the personalization step:

.. code-block:: bash

   indelinside matrix \
       -i indelinside.out.txt \
       -c 2 \
       --sample_name my_sample

The ``-c`` option specifies the consensus level.

For example, ``-c 2`` selects indels supported by at
least two variant callers.

The output is a COSMIC-compatible 83-channel indel
signature matrix.

Downstream Signature Analysis
-----------------------------

The resulting matrix can be used with external
mutational signature analysis tools, such as
SigProfilerAssignment.

For advanced control over individual processing steps,
see :doc:`python_api/index`.
