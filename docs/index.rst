indelinside
===========

Somatic Indel Personalization and Signature Analysis

.. image:: https://img.shields.io/badge/docs-Read%20the%20Docs-blue
   :target: https://indelinside.readthedocs.io/

Overview
--------

**indelinside** is a command-line tool for somatic
insertion and deletion (indel) personalization
and mutational signature analysis.

It is built on top of the ``variantpost`` Python library.

The recommended way to use this project is through
the ``indelinside`` command-line interface.

For users who need more granular control over
the underlying algorithms, ``variantpost`` provides
a Python API for building custom workflows.

Key Features
------------

* Somatic indel representation harmonization
* Local haplotype-based indel personalization
* COSMIC-compatible indel mutational signature matrices
* Python API for advanced and custom workflows

Getting Started
---------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Command-Line Interface

   cli/index

.. toctree::
   :maxdepth: 2
   :caption: Advanced Usage

   python_api/index

.. toctree::
   :maxdepth: 1
   :caption: Reference

   references
