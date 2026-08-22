"""Contract for the CURIE prefix table.

The table is the only place a genomics CURIE learns what it expands to, so the
tests here guard the two ways it can silently go wrong: a CURIE minted under a
prefix the table does not define, and a prefix whose IRI cannot be concatenated
with a local part.
"""
from __future__ import annotations

from urllib.parse import urlparse

import pytest

from indicium_adapters_genomics import _edam, _identifiers
from indicium_adapters_genomics._namespaces import GENOMICS_PREFIXES

PREFIX_TERMINATORS = (":", "/", "_", "#")

# The namespaces a genomics claim is expected to reach for. Named explicitly so
# dropping one is a test failure rather than a quiet loss of resolution.
REQUIRED_PREFIXES = (
    "ENA", "SRA", "GEO", "ArrayExpress", "BioProject",
    "dbSNP", "ClinVar", "ClinVar.record", "ClinVar.submission",
    "ENSEMBL", "RefSeq", "HGNC", "NCBIGene",
    "GO", "SO", "NCBITaxon", "edam",
)

# One expanded IRI per identifier family, checked whole rather than by prefix.
# The way this table rots is a prefix that stays declared while its expansion
# stops resolving the local parts the grammar mints under it.
EXPANSIONS = (
    ("rs334", "https://identifiers.org/dbsnp:rs334"),
    ("VCV000012345", "https://www.ncbi.nlm.nih.gov/clinvar/variation/VCV000012345"),
    ("RCV000012345", "https://identifiers.org/clinvar.record:RCV000012345"),
    ("SCV000012345", "https://identifiers.org/clinvar.submission:SCV000012345"),
    ("ENSG00000139618", "https://identifiers.org/ensembl:ENSG00000139618"),
    ("NM_000546.6", "https://identifiers.org/refseq:NM_000546.6"),
    ("HGNC:4827", "https://identifiers.org/hgnc:4827"),
    ("SO:0000704", "http://purl.obolibrary.org/obo/SO_0000704"),
)


def test_every_required_namespace_is_declared():
    missing = [prefix for prefix in REQUIRED_PREFIXES if prefix not in GENOMICS_PREFIXES]
    assert not missing, f"undeclared namespaces: {missing}"


def test_every_base_iri_is_http_and_ends_curie_safely():
    for prefix, base in GENOMICS_PREFIXES.items():
        assert urlparse(base).scheme in ("http", "https"), f"{prefix} is not http(s)"
        assert base.endswith(PREFIX_TERMINATORS), f"{prefix} does not end CURIE-safely: {base}"


def test_every_identifier_authority_prefix_can_be_expanded():
    """A grammar that mints a CURIE the table cannot expand is a dead identifier."""
    undeclared = [p for p in _identifiers.authorities() if p not in GENOMICS_PREFIXES]
    assert not undeclared, f"minted but not expandable: {undeclared}"


def test_the_edam_prefix_matches_what_the_edam_index_mints():
    """`edam` is lower-case here because that is the prefix `_edam` emits."""
    curie = _edam.lookup("comparative genomics")
    assert curie is not None
    assert curie.split(":", 1)[0] in GENOMICS_PREFIXES


@pytest.mark.parametrize(("text", "expected"), EXPANSIONS)
def test_a_minted_curie_expands_to_the_authoritys_own_iri(text, expected):
    prefix, local = _identifiers.resolve(text).split(":", 1)
    assert GENOMICS_PREFIXES[prefix] + local == expected


def test_base_iris_are_distinct():
    assert len(set(GENOMICS_PREFIXES.values())) == len(GENOMICS_PREFIXES)
