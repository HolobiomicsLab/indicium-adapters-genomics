"""Protocol-conformance contract for the genomics adapter.

The invariants every indicium DomainAdapter must hold, plus the two this
adapter exists to satisfy for ASB: it declares ``primary_domain`` (without it a
genomics paper classifies as "other"), and it resolves concepts only by exact
match (a wrong CURIE is worse than no CURIE).
"""

from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import urlparse

from indicium_adapters import CURRENT_API_VERSION, DomainAdapter

from indicium_adapters_genomics import parse_gff3_to_claims, parse_vcf_to_claims

SNAKE_CASE = re.compile(r"^[a-z][a-z0-9_]*$")
SLUG = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
PREFIX_TERMINATORS = (":", "/", "_", "#")


def test_satisfies_the_domain_adapter_protocol(adapter):
    assert isinstance(adapter, DomainAdapter)
    assert adapter.api_version == CURRENT_API_VERSION


def test_declares_primary_domain_so_asb_does_not_classify_as_other(adapter):
    assert adapter.primary_domain == "genomics"
    assert adapter.domain_id == "genomics"


def test_qualifiers_are_snake_case(adapter):
    """Perspicacité's asb:domainQualifier SHACL shape enforces this pattern."""
    bad = [q for q in adapter.qualifiers if not SNAKE_CASE.match(q)]
    assert not bad, f"non-snake_case qualifiers: {bad}"


def test_vocabulary_terms_are_slugs(adapter):
    for term in (*adapter.subdomains, *adapter.techniques):
        assert SLUG.match(term), f"not a kebab-case slug: {term}"


def test_vocabulary_has_no_duplicates(adapter):
    assert len(set(adapter.subdomains)) == len(adapter.subdomains)
    assert len(set(adapter.techniques)) == len(adapter.techniques)


def test_ontology_prefixes_expand_to_curie_safe_iris(adapter):
    for prefix, base in adapter.ontology_prefixes.items():
        parsed = urlparse(base)
        assert parsed.scheme in ("http", "https"), f"{prefix} is not http(s)"
        assert base.endswith(PREFIX_TERMINATORS), f"{prefix} does not end CURIE-safely: {base}"


def test_extraction_context_names_the_domain_and_its_qualifiers(adapter):
    context = adapter.extraction_context()
    assert "genomics" in context.lower()
    for qualifier in ("sequenced_on", "assembled_with", "variant_called_by"):
        assert qualifier in context


# --- enrichment ---------------------------------------------------------------


def test_enrich_claim_is_pure(adapter):
    claim = {"subject": "Genomics", "predicate": "uses", "object": "GATK"}
    before = dict(claim)
    adapter.enrich_claim(claim)
    assert claim == before, "enrich_claim mutated its input"


def test_known_label_resolves_to_an_edam_curie(adapter):
    out = adapter.enrich_claim({"subject": "Genomics", "object": "GATK"})
    assert out["ontology_terms"]["subject"].startswith("edam:")
    assert out["ontology_term_justifications"]["subject"] == "semapv:LexicalMatching"


def test_lookup_is_case_and_whitespace_insensitive(adapter):
    spaced = adapter.enrich_claim({"subject": "  gEnOmIcS  "})
    exact = adapter.enrich_claim({"subject": "Genomics"})
    assert spaced["ontology_terms"]["subject"] == exact["ontology_terms"]["subject"]


def test_unresolved_text_is_recorded_not_silently_dropped(adapter):
    out = adapter.enrich_claim({"subject": "a phrase that is not an EDAM concept"})
    assert out["ontology_terms_unmatched"]["subject"] == "a phrase that is not an EDAM concept"
    assert "ontology_terms" not in out


def test_substring_of_a_concept_does_not_resolve(adapter):
    """A prefix of a real index entry must not pick that entry up — a wrong
    CURIE is worse than no CURIE. Both strings below are proper substrings of
    concepts the shipped index does carry, so the test fails if matching ever
    loosens; a string absent from the index altogether would prove nothing."""
    for prefix in ("comparative genom", "genomic"):
        assert adapter.enrich_claim({"subject": prefix}).get("ontology_terms") is None


def test_an_already_resolved_slot_is_left_to_its_owner(adapter):
    claim = {"subject": "Genomics", "ontology_terms": {"subject": "SO:0000704"}}
    out = adapter.enrich_claim(claim)
    assert out["ontology_terms"]["subject"] == "SO:0000704"


def test_empty_and_non_string_slots_are_skipped(adapter):
    out = adapter.enrich_claim({"subject": "", "object": None, "qualifier": 42})
    assert "ontology_terms" not in out
    assert "ontology_terms_unmatched" not in out


def test_the_claim_shapes_this_domain_needs_have_a_qualifier(adapter):
    """Variant calls, expression contrasts, assembly metrics and annotation."""
    for qualifier in (
        "variant_called_by",
        "genotyped_as",
        "has_consequence",
        "differentially_expressed_in",
        "has_assembly_metric",
        "annotated_by",
        "located_on",
        "has_feature",
    ):
        assert qualifier in adapter.qualifiers


def test_the_projections_only_emit_declared_qualifiers(adapter):
    """A parser inventing a qualifier would put claims outside the vocabulary."""
    fixtures = Path(__file__).parent / "fixtures"
    claims = parse_vcf_to_claims(fixtures / "sample.vcf")
    claims += parse_gff3_to_claims(fixtures / "sample.gff3")
    emitted = {claim["qualifier"] for claim in claims}
    assert emitted <= adapter.qualifiers
