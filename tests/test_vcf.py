"""Tests for the VCF -> asb:Claim projection."""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from indicium_adapters_genomics import parse_vcf_to_claims

FIXTURE = Path(__file__).parent / "fixtures" / "sample.vcf"


def _by_qualifier(claims: list[dict], qualifier: str) -> list[dict]:
    return [claim for claim in claims if claim.get("qualifier") == qualifier]


def test_parses_fixture_without_error():
    claims = parse_vcf_to_claims(FIXTURE)
    assert isinstance(claims, list)
    assert len(claims) >= 1


def test_every_record_yields_a_located_on_claim():
    located = _by_qualifier(parse_vcf_to_claims(FIXTURE), "located_on")
    assert [claim["object"] for claim in located] == ["11", "17", "7"]


def test_a_record_without_an_id_gets_a_positional_label():
    located = _by_qualifier(parse_vcf_to_claims(FIXTURE), "located_on")
    assert located[2]["subject"] == "7:117559590ATCT>A"


def test_an_rsid_yields_an_identifies_claim_with_a_composite_justification():
    """The ID column is written by the calling pipeline, so the CURIE is curated."""
    identifies = _by_qualifier(parse_vcf_to_claims(FIXTURE), "identifies")
    first = identifies[0]
    assert first["object"] == "dbSNP:rs334"
    assert first["ontology_terms"]["object"] == "dbSNP:rs334"
    assert first["ontology_term_justifications"]["object"] == "semapv:CompositeMatching"


def test_a_clinvar_accession_in_the_id_column_also_identifies():
    identifies = _by_qualifier(parse_vcf_to_claims(FIXTURE), "identifies")
    objects = [claim["object"] for claim in identifies]
    assert objects == ["dbSNP:rs334", "ClinVar:VCV000012345"]


def test_every_entry_of_a_semicolon_separated_id_column_identifies(tmp_path):
    """Regression: only the first of VCF's ID list reached an identifies claim."""
    vcf = tmp_path / "multi-id.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "11\t5227002\trs334;VCV000012345\tT\tA\t60\tPASS\t.\n",
        encoding="utf-8",
    )
    identifies = _by_qualifier(parse_vcf_to_claims(vcf), "identifies")
    assert [claim["object"] for claim in identifies] == [
        "dbSNP:rs334", "ClinVar:VCV000012345",
    ]
    assert {claim["subject"] for claim in identifies} == {"rs334"}


def test_an_unrecognised_id_entry_is_skipped_not_guessed_at(tmp_path):
    vcf = tmp_path / "local-id.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "11\t5227002\tinhouse_0001;rs334\tT\tA\t60\tPASS\t.\n",
        encoding="utf-8",
    )
    identifies = _by_qualifier(parse_vcf_to_claims(vcf), "identifies")
    assert [claim["object"] for claim in identifies] == ["dbSNP:rs334"]


def test_a_missing_id_yields_no_identifies_claim():
    claims = parse_vcf_to_claims(FIXTURE)
    assert len(_by_qualifier(claims, "identifies")) == 2


def test_uncalled_genotypes_are_skipped():
    """`./.` is an absent call, not a genotype worth asserting."""
    genotypes = _by_qualifier(parse_vcf_to_claims(FIXTURE), "genotyped_as")
    assert [claim["object"] for claim in genotypes] == ["0/1", "1/1", "0/1", "0/0"]


def test_a_genotype_claim_names_its_sample_in_the_evidence_id():
    genotypes = _by_qualifier(parse_vcf_to_claims(FIXTURE), "genotyped_as")
    assert genotypes[0]["evidence_id"] == "vcf:sample:11:5227002:NA00001"


def test_the_gt_subfield_is_read_by_position_not_assumed_first():
    """The fixture's first record uses `GT:DP`, so a naive split would read `35`."""
    genotypes = _by_qualifier(parse_vcf_to_claims(FIXTURE), "genotyped_as")
    assert genotypes[0]["object"] == "0/1"


def test_a_row_shorter_than_the_header_does_not_raise():
    """The second fixture record omits the NA00002 column entirely."""
    genotypes = _by_qualifier(parse_vcf_to_claims(FIXTURE), "genotyped_as")
    assert all(claim["object"] for claim in genotypes)


def test_gzip_input_produces_the_same_claims(tmp_path):
    compressed = tmp_path / "sample.vcf.gz"
    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(FIXTURE.read_text(encoding="utf-8"))
    assert parse_vcf_to_claims(compressed) == parse_vcf_to_claims(FIXTURE)


def test_claims_are_stable_across_calls():
    assert parse_vcf_to_claims(FIXTURE) == parse_vcf_to_claims(FIXTURE)


def test_a_file_with_no_records_yields_no_claims(tmp_path):
    empty = tmp_path / "empty.vcf"
    empty.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\n", encoding="utf-8")
    assert parse_vcf_to_claims(empty) == []


def test_missing_file_raises_filenotfound():
    with pytest.raises(FileNotFoundError):
        parse_vcf_to_claims(Path("/nonexistent/sample.vcf"))
