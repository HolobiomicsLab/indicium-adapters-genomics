"""Bundled data resources for the genomics adapter.

``edam_index.json`` ships with the wheel. ``gene_index.tsv`` does not: it is
regenerated from the HGNC complete set by ``build_gene_index.py``, a
maintainer-only script that is never run at install or import time.
"""
