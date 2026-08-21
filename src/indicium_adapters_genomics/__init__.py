"""Genomics DomainAdapter for the indicium ecosystem.

Exposes :class:`GenomicsAdapter` via the ``indicium.adapters`` entry point, so
``indicium_adapters.discover_adapters()`` finds it on install and ASB's
classification vocabulary picks up ``primary_domain="genomics"``.
"""

from .genomics import GenomicsAdapter

__all__ = ["GenomicsAdapter", "__version__"]

__version__ = "0.1.0"
