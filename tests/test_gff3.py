"""Tests for the GFF3 -> asb:Claim projection."""
from __future__ import annotations

from pathlib import Path

import pytest

from indicium_adapters_genomics import parse_gff3_to_claims

FIXTURE = Path(__file__).parent / "fixtures" / "sample.gff3"


def _by_qualifier(claims: list[dict], qualifier: str) -> list[dict]:
    return [claim for claim in claims if claim.get("qualifier") == qualifier]


def test_parses_fixture_without_error():
    claims = parse_gff3_to_claims(FIXTURE)
    assert isinstance(claims, list)
    assert len(claims) >= 1


def test_every_feature_yields_a_located_on_claim():
    located = _by_qualifier(parse_gff3_to_claims(FIXTURE), "located_on")
    assert [claim["subject"] for claim in located] == [
        "gene00001", "mRNA00001", "exon00001",
    ]


def test_an_so_accession_type_is_carried_onto_the_subject():
    located = _by_qualifier(parse_gff3_to_claims(FIXTURE), "located_on")
    gene = located[0]
    assert gene["ontology_terms"]["subject"] == "SO:0000704"
    assert gene["ontology_term_justifications"]["subject"] == "semapv:CompositeMatching"


def test_a_type_written_as_a_term_name_is_not_guessed_at():
    """`mRNA` needs an SO label index to resolve; inventing a CURIE is worse."""
    located = _by_qualifier(parse_gff3_to_claims(FIXTURE), "located_on")
    assert located[1]["ontology_terms"] == {}


def test_the_source_column_yields_an_annotated_by_claim():
    annotated = _by_qualifier(parse_gff3_to_claims(FIXTURE), "annotated_by")
    assert [claim["object"] for claim in annotated] == ["Ensembl", "Ensembl"]


def test_an_absent_source_yields_no_annotated_by_claim():
    """The exon row carries `.`, which the spec defines as "no value"."""
    annotated = _by_qualifier(parse_gff3_to_claims(FIXTURE), "annotated_by")
    assert "exon00001" not in [claim["subject"] for claim in annotated]


def test_parent_attributes_yield_has_feature_claims():
    nested = _by_qualifier(parse_gff3_to_claims(FIXTURE), "has_feature")
    assert [(c["subject"], c["object"]) for c in nested] == [
        ("gene00001", "mRNA00001"),
        ("mRNA00001", "exon00001"),
    ]


def test_the_fasta_section_is_not_parsed_as_features():
    assert len(parse_gff3_to_claims(FIXTURE)) == 7


def test_evidence_ids_locate_the_feature_in_the_file():
    located = _by_qualifier(parse_gff3_to_claims(FIXTURE), "located_on")
    assert located[0]["evidence_id"] == "gff3:sample:ctg123:1000-9000"


def test_multiple_parents_each_yield_a_claim(tmp_path):
    gff3 = tmp_path / "multi.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "ctg1\t.\texon\t1\t10\t.\t+\t.\tID=e1;Parent=t1,t2\n",
        encoding="utf-8",
    )
    nested = _by_qualifier(parse_gff3_to_claims(gff3), "has_feature")
    assert [claim["subject"] for claim in nested] == ["t1", "t2"]


def test_attribute_values_are_percent_decoded(tmp_path):
    """The spec escapes reserved characters; a raw `%2C` would be a wrong label."""
    gff3 = tmp_path / "escaped.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "ctg1\tmaker\tgene\t1\t10\t.\t+\t.\tID=gene%20one\n",
        encoding="utf-8",
    )
    located = _by_qualifier(parse_gff3_to_claims(gff3), "located_on")
    assert located[0]["subject"] == "gene one"


def test_an_escaped_comma_inside_a_parent_id_is_not_a_separator(tmp_path):
    """Regression: decoding before the comma split turned one Parent into two."""
    gff3 = tmp_path / "escaped-parent.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "ctg1\t.\texon\t1\t10\t.\t+\t.\tID=e1;Parent=t%2C1,t2\n",
        encoding="utf-8",
    )
    nested = _by_qualifier(parse_gff3_to_claims(gff3), "has_feature")
    assert [claim["subject"] for claim in nested] == ["t,1", "t2"]


def test_a_non_so_curie_in_the_type_column_is_not_carried(tmp_path):
    """Column three admits SO alone; a GO CURIE there is not a feature type."""
    gff3 = tmp_path / "wrong-authority.gff3"
    gff3.write_text(
        "##gff-version 3\nctg1\t.\tGO:0006355\t1\t10\t.\t+\t.\tID=f1\n",
        encoding="utf-8",
    )
    located = _by_qualifier(parse_gff3_to_claims(gff3), "located_on")
    assert located[0]["ontology_terms"] == {}


def test_a_feature_without_id_or_name_gets_a_positional_label(tmp_path):
    gff3 = tmp_path / "bare.gff3"
    gff3.write_text("##gff-version 3\nctg1\t.\texon\t5\t9\t.\t+\t.\t.\n", encoding="utf-8")
    located = _by_qualifier(parse_gff3_to_claims(gff3), "located_on")
    assert located[0]["subject"] == "ctg1:5-9"


def test_claims_are_stable_across_calls():
    assert parse_gff3_to_claims(FIXTURE) == parse_gff3_to_claims(FIXTURE)


def test_missing_file_raises_filenotfound():
    with pytest.raises(FileNotFoundError):
        parse_gff3_to_claims(Path("/nonexistent/sample.gff3"))
