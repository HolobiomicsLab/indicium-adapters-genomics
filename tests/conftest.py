"""Shared pytest fixtures."""
from __future__ import annotations

from pathlib import Path

import pytest

from indicium_adapters_genomics import GeneIndex, GenomicsAdapter

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def adapter() -> GenomicsAdapter:
    """The adapter with whatever gene index the installation ships."""
    return GenomicsAdapter()


@pytest.fixture
def gene_index() -> GeneIndex:
    """A two-symbol index read from the fixture TSV (real HGNC records)."""
    return GeneIndex.from_tsv(FIXTURES / "gene_index.tsv")


@pytest.fixture
def gene_adapter(gene_index: GeneIndex) -> GenomicsAdapter:
    """An adapter wired to the fixture nomenclature index."""
    return GenomicsAdapter(gene_index)


@pytest.fixture
def base_claim() -> dict:
    """A minimal 5-slot claim dict for adapter input."""
    return {
        "subject": "HBB",
        "qualifier": "located_on",
        "object": "11",
        "ontology_terms": {},
        "evidence_id": "fixture-001",
    }
