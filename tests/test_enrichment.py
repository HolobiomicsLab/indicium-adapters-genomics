"""Contract for ``GenomicsAdapter.enrich_claim``.

Three resolvers feed one slot, so the tests fix the order they run in, the
SEMAPV justification each earns, and the purity guarantee the whole method rests
on. Enrichment runs over claim dicts that are already extracted and validated;
nothing here reads article text.
"""
from __future__ import annotations

from indicium_adapters_genomics import GenomicsAdapter


def test_a_variant_identifier_resolves_as_a_composite_match(adapter):
    """The text already *is* the authority's identifier, so it is not a lexical guess."""
    out = adapter.enrich_claim({"subject": "rs334", "qualifier": "located_on"})
    assert out["ontology_terms"]["subject"] == "dbSNP:rs334"
    assert out["ontology_term_justifications"]["subject"] == "semapv:CompositeMatching"


def test_a_clinvar_accession_resolves(adapter):
    out = adapter.enrich_claim({"object": "VCV000012345"})
    assert out["ontology_terms"]["object"] == "ClinVar:VCV000012345"


def test_a_gene_symbol_resolves_to_hgnc_with_ensembl_as_an_xref(gene_adapter):
    out = gene_adapter.enrich_claim({"subject": "HBB", "object": "11"})
    assert out["ontology_terms"]["subject"] == "HGNC:4827"
    assert out["xrefs"]["subject"]["ENSEMBL"] == "ENSEMBL:ENSG00000244734"
    assert out["ontology_term_justifications"]["subject"] == "semapv:LexicalMatching"


def test_an_unresolved_gene_symbol_is_recorded_not_dropped(gene_adapter):
    out = gene_adapter.enrich_claim({"subject": "NOT_A_GENE"})
    assert out["ontology_terms_unmatched"]["subject"] == "NOT_A_GENE"


def test_identifier_grammar_wins_over_the_edam_labels(gene_adapter):
    """An rsID must never fall through to a lexical index."""
    out = gene_adapter.enrich_claim({"subject": "rs334"})
    assert out["ontology_terms"]["subject"].startswith("dbSNP:")


def test_an_edam_label_still_resolves_after_the_new_resolvers(gene_adapter):
    out = gene_adapter.enrich_claim({"subject": "comparative genomics"})
    assert out["ontology_terms"]["subject"].startswith("edam:")


def test_enrich_claim_does_not_mutate_nested_maps(gene_adapter):
    """Regression: copying only the outer dict leaks writes back into the caller."""
    terms: dict[str, str] = {}
    claim = {"subject": "HBB", "ontology_terms": terms, "xrefs": {}}
    gene_adapter.enrich_claim(claim)
    assert terms == {}
    assert claim["xrefs"] == {}


def test_pre_existing_xrefs_on_another_slot_survive(gene_adapter):
    claim = {"subject": "HBB", "xrefs": {"object": {"HGNC": "HGNC:11998"}}}
    out = gene_adapter.enrich_claim(claim)
    assert out["xrefs"]["object"] == {"HGNC": "HGNC:11998"}
    assert "HGNC" in out["xrefs"]["subject"]


def test_a_declared_qualifier_is_not_reported_as_unresolved_text(adapter):
    """Regression: it put an unmatched entry on every claim the adapter saw."""
    out = adapter.enrich_claim({"subject": "rs334", "qualifier": "located_on"})
    assert "qualifier" not in out.get("ontology_terms_unmatched", {})


def test_a_qualifier_outside_the_declared_set_is_still_resolved(adapter):
    """Skipping is for this adapter's own vocabulary, not for the whole slot."""
    out = adapter.enrich_claim({"qualifier": "comparative genomics"})
    assert out["ontology_terms"]["qualifier"].startswith("edam:")


def test_the_shipped_adapter_reports_whether_nomenclature_is_installed():
    """The index may be absent; callers must be able to see that it is."""
    assert isinstance(GenomicsAdapter().gene_index.is_available, bool)
