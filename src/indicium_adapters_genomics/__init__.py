"""Genomics DomainAdapter for the indicium ecosystem.

Exposes :class:`GenomicsAdapter` via the ``indicium.adapters`` entry point, so
``indicium_adapters.discover_adapters()`` finds it on install and ASB's
classification vocabulary picks up ``primary_domain="genomics"``.

Two file projections sit alongside the adapter, mirroring how the sibling
metabolomics package projects mzTab-M rows:
:func:`~indicium_adapters_genomics.vcf.parse_vcf_to_claims` and
:func:`~indicium_adapters_genomics.gff3.parse_gff3_to_claims`.
"""

from ._genes import GENE_INDEX, GeneIndex, GeneMatch
from .genomics import GenomicsAdapter
from .gff3 import parse_gff3_to_claims
from .vcf import parse_vcf_to_claims

__all__ = [
    "GENE_INDEX",
    "GeneIndex",
    "GeneMatch",
    "GenomicsAdapter",
    "__version__",
    "parse_gff3_to_claims",
    "parse_vcf_to_claims",
]

__version__ = "0.3.0"
