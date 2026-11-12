Installation
============

Requirements
------------

indelinside requires Python 3.11 or later.

Install from PyPI
-----------------

The recommended installation method is pip:

.. code-block:: bash

   pip install variantpost

This package installs both the ``variantpost`` Python library
and the ``indelinside`` command-line interface.

Verify Installation
-------------------

Check that the CLI is available:

.. code-block:: bash

   indelinside --help

You should see the available commands, including
``personalize`` and ``matrix``.

Dependencies
------------

Please ensure that the following input files are available
before running the pipeline:

* Tumor BAM file
* Normal BAM file
* Reference genome FASTA file
* One or more somatic VCF files

For BAM and FASTA files, make sure the corresponding
index files are available when required.
