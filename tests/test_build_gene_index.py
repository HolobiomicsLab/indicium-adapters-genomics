"""Tests for the pure part of the maintainer index builder.

The download path is not exercised — it hits HGNC. What is tested is the row
filter, because shipping withdrawn records would resolve symbols the authority
no longer recognises.
"""
from __future__ import annotations

from indicium_adapters_genomics import GeneIndex
from indicium_adapters_genomics.data.build_gene_index import (
    OUT_COLUMNS,
    _to_row,
    _write_index,
)

APPROVED = {
    "hgnc_id": "HGNC:4827",
    "symbol": "HBB",
    "name": "hemoglobin subunit beta",
    "status": "Approved",
    "ensembl_gene_id": "ENSG00000244734",
    "entrez_id": "3043",
}


def test_an_approved_record_projects_onto_the_output_schema():
    assert _to_row(APPROVED) == {
        "symbol": "HBB",
        "hgnc": "HGNC:4827",
        "ensembl": "ENSG00000244734",
        "ncbigene": "3043",
        "name": "hemoglobin subunit beta",
    }


def test_a_withdrawn_record_is_dropped():
    assert _to_row({**APPROVED, "status": "Entry Withdrawn"}) is None


def test_a_record_without_a_symbol_is_dropped():
    assert _to_row({**APPROVED, "symbol": "  "}) is None


def test_a_record_missing_cross_references_still_projects():
    row = _to_row({"symbol": "HBB", "hgnc_id": "HGNC:4827", "status": "Approved"})
    assert row["ensembl"] == ""


def test_write_index_emits_a_header_and_only_approved_rows(tmp_path):
    source = tmp_path / "hgnc.txt"
    header = "\t".join(("hgnc_id", "symbol", "name", "status", "ensembl_gene_id", "entrez_id"))
    source.write_text(
        f"{header}\n"
        "HGNC:4827\tHBB\themoglobin subunit beta\tApproved\tENSG00000244734\t3043\n"
        "HGNC:31\tABC\told name\tEntry Withdrawn\t\t\n",
        encoding="utf-8",
    )
    target = tmp_path / "gene_index.tsv"
    assert _write_index(source, target) == 1
    lines = target.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(OUT_COLUMNS)
    assert len(lines) == 2


def test_rows_from_the_builder_load_back_into_the_index(tmp_path):
    """The builder's output schema and the loader's expectations must agree."""
    source = tmp_path / "hgnc.txt"
    header = "\t".join(("hgnc_id", "symbol", "name", "status", "ensembl_gene_id", "entrez_id"))
    source.write_text(
        f"{header}\nHGNC:4827\tHBB\thb beta\tApproved\tENSG00000244734\t3043\n",
        encoding="utf-8",
    )
    target = tmp_path / "gene_index.tsv"
    _write_index(source, target)
    assert GeneIndex.from_tsv(target).resolve("HBB").principal == "HGNC:4827"


def test_output_schema_is_the_documented_one():
    assert OUT_COLUMNS == ("symbol", "hgnc", "ensembl", "ncbigene", "name")
