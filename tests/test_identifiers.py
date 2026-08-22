"""Contract for structural identifier recognition.

Every grammar here is published by the authority that issues the identifier, so
the tests come in pairs: a real identifier of that family resolves, and a string
that merely looks close does not. The second half matters more -- a wrong CURIE
is worse than no CURIE, and the way this module could rot is by loosening a
pattern until ordinary prose starts matching.
"""
from __future__ import annotations

import pytest

from indicium_adapters_genomics import _identifiers
from indicium_adapters_genomics._namespaces import GENOMICS_PREFIXES

RESOLVING = [
    ("rs334", "dbSNP:rs334"),
    ("VCV000012345", "ClinVar:VCV000012345"),
    ("RCV000012345.7", "ClinVar.record:RCV000012345.7"),
    ("SCV000012345", "ClinVar.submission:SCV000012345"),
    ("ENSG00000139618", "ENSEMBL:ENSG00000139618"),
    ("ENST00000380152.8", "ENSEMBL:ENST00000380152.8"),
    ("ENSMUSG00000017146", "ENSEMBL:ENSMUSG00000017146"),
    ("NM_000546.6", "RefSeq:NM_000546.6"),
    ("NC_000011.10", "RefSeq:NC_000011.10"),
    ("NP_000537.3", "RefSeq:NP_000537.3"),
    ("NZ_CP012345.1", "RefSeq:NZ_CP012345.1"),
    ("HGNC:4827", "HGNC:4827"),
    ("SO:0000704", "SO:0000704"),
    ("GO:0006355", "GO:0006355"),
    ("NCBITaxon:9606", "NCBITaxon:9606"),
]

NOT_RESOLVING = [
    "the reference genome",
    "variant calling",
    "rs",
    "rsID",
    "ENSG0000013961",           # ten digits, not eleven
    "ENSX00000139618",          # X is not a feature-type letter
    "NM_1234",                  # too few digits for a RefSeq accession
    "QQ_000546.6",              # not a RefSeq accession prefix
    "NZ_000546",                # NZ_ takes a WGS project code, not bare digits
    "VCV12345",                 # too few digits for a ClinVar accession
    "SO:704",                   # OBO local parts are seven digits
    "",
    "   ",
]


@pytest.mark.parametrize(("text", "expected"), RESOLVING)
def test_published_identifiers_resolve(text, expected):
    assert _identifiers.resolve(text) == expected


@pytest.mark.parametrize("text", NOT_RESOLVING)
def test_lookalikes_and_prose_do_not_resolve(text):
    assert _identifiers.resolve(text) is None


def test_an_identifier_inside_a_sentence_does_not_resolve():
    """Matching is full-string: scanning prose belongs to asb-domain-core."""
    assert _identifiers.resolve("carriers of rs334 and rs335") is None


def test_surrounding_whitespace_is_tolerated():
    assert _identifiers.resolve("  rs334\n") == "dbSNP:rs334"


def test_case_is_not_relaxed():
    """The authorities issue a fixed case; relaxing it would match prose."""
    assert _identifiers.resolve("RS334") is None
    assert _identifiers.resolve("ensg00000139618") is None


def test_the_three_clinvar_families_do_not_share_one_prefix():
    """identifiers.org registers VCV, RCV and SCV separately; one prefix would
    expand two of the three to an IRI their registry pattern rejects."""
    prefixes = {
        _identifiers.resolve(text).split(":", 1)[0]
        for text in ("VCV000012345", "RCV000012345", "SCV000012345")
    }
    assert prefixes == {"ClinVar", "ClinVar.record", "ClinVar.submission"}


def test_the_sequence_ontology_resolver_accepts_only_so():
    """GFF3's type column admits SO alone, so the general resolver is too wide."""
    assert _identifiers.resolve_sequence_ontology("SO:0000704") == "SO:0000704"
    assert _identifiers.resolve("GO:0006355") == "GO:0006355"
    assert _identifiers.resolve_sequence_ontology("GO:0006355") is None
    assert _identifiers.resolve_sequence_ontology("mRNA") is None


def test_every_grammar_names_an_authority():
    authorities = _identifiers.authorities()
    assert authorities
    assert all(name.strip() for name in authorities.values())


def test_every_resolvable_family_mints_a_declared_prefix():
    for text, expected in RESOLVING:
        prefix = expected.split(":", 1)[0]
        assert prefix in GENOMICS_PREFIXES, f"{text} mints undeclared prefix {prefix}"
