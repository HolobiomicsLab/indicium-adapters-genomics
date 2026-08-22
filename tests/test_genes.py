"""Contract for gene-symbol resolution.

The bundled index is deliberately absent from the repository (see ``_genes``),
so these tests drive the loader from the fixture TSV. The behaviour that must
never regress is the last test: an index with no data reports that it has none
instead of returning the same "not found" a populated index would.
"""
from __future__ import annotations

from indicium_adapters_genomics import GeneIndex

HBB_ROW = {
    "symbol": "HBB",
    "hgnc": "HGNC:4827",
    "ensembl": "ENSG00000244734",
    "ncbigene": "3043",
}


def test_fixture_index_loads_only_rows_with_a_symbol(gene_index):
    assert len(gene_index) == 2
    assert gene_index.is_available


def test_symbol_resolves_to_the_hgnc_id_as_principal(gene_index):
    match = gene_index.resolve("HBB")
    assert match.principal == "HGNC:4827"


def test_symbol_carries_every_cross_reference(gene_index):
    match = gene_index.resolve("TP53")
    assert match.xrefs == {
        "HGNC": "HGNC:11998",
        "ENSEMBL": "ENSEMBL:ENSG00000141510",
        "NCBIGene": "NCBIGene:7157",
    }


def test_lookup_retries_upper_cased(gene_index):
    assert gene_index.resolve("hbb").principal == "HGNC:4827"


def test_unlisted_symbol_returns_none(gene_index):
    assert gene_index.resolve("NOT_A_GENE") is None


def test_blank_symbol_returns_none(gene_index):
    assert gene_index.resolve("   ") is None


def test_rows_can_be_supplied_directly():
    index = GeneIndex(rows={"HBB": HBB_ROW})
    assert index.resolve("HBB").xrefs["ENSEMBL"] == "ENSEMBL:ENSG00000244734"


def test_a_row_with_no_identifier_resolves_to_none():
    index = GeneIndex(rows={"HBB": {"symbol": "HBB", "name": "haemoglobin subunit beta"}})
    assert index.resolve("HBB") is None


def test_an_empty_index_reports_that_it_is_unavailable():
    """An absent index must stay distinguishable from a gene the authority omits."""
    index = GeneIndex(rows={})
    assert not index.is_available
    assert len(index) == 0
    assert index.resolve("HBB") is None
