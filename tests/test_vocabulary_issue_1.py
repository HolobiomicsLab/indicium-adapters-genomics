"""Issue #1: the in-domain topics the harvest carries must be declared.

The corpus harvest tagged genomics rows with EDAM topics this adapter never
uttered, so ``domain_admission`` could name no field for them and half the
harvest went unjudged. These tests pin the terms that closed the gap, and pin
the one candidate deliberately left out.
"""

from __future__ import annotations

import re

SLUG = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")

ADDED_SUBDOMAINS = (
    "genetics",
    "phylogenetics",
    "population-genetics",
    "phylogeny",
    "gwas-study",
    "phylogenomics",
    "copy-number-variation",
    "molecular-genetics",
)

ADDED_TECHNIQUES = (
    "genotype-and-phenotype",
    "sequence-assembly",
    "oxford-nanopore",
    "pacbio",
)


def test_the_measured_subdomains_are_declared(adapter):
    missing = [term for term in ADDED_SUBDOMAINS if term not in adapter.subdomains]
    assert not missing, f"issue #1 subdomains still undeclared: {missing}"


def test_the_measured_techniques_are_declared(adapter):
    missing = [term for term in ADDED_TECHNIQUES if term not in adapter.techniques]
    assert not missing, f"issue #1 techniques still undeclared: {missing}"


def test_the_bare_term_mapping_stays_out(adapter):
    """A single common word matches far more topics than it identifies.

    ``domain_admission.names`` admits a row when every word of a declared term
    appears in one of its topics, so "mapping" alone would admit "Read mapping",
    "Genetic mapping" and "Mapping" from any field at all.
    """
    assert "mapping" not in adapter.techniques
    assert "mapping" not in adapter.subdomains


def test_every_declared_term_is_a_lowercase_hyphen_slug(adapter):
    bad = [t for t in (*adapter.subdomains, *adapter.techniques) if not SLUG.match(t)]
    assert not bad, f"not lowercase hyphen-separated slugs: {bad}"


def test_no_term_is_declared_twice_anywhere_in_the_adapter(adapter):
    terms = [*adapter.subdomains, *adapter.techniques]
    duplicates = sorted({t for t in terms if terms.count(t) > 1})
    assert not duplicates, f"declared more than once: {duplicates}"
