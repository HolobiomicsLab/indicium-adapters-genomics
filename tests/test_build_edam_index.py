"""Tests for the pure part of the maintainer EDAM index builder.

The OLS4 fetch is not exercised — it hits the network. What is tested is the
projection from an OLS term record onto the index and the subtopic list, because
that projection is what makes the vocabulary EDAM's rather than an editor's.
"""

from __future__ import annotations

from indicium_adapters_genomics.data.build_edam_index import (
    DOMAIN_ROOT,
    _entries,
    _index_of,
    _slug,
    _subtopics_of,
)

GENOMICS = {"iri": DOMAIN_ROOT, "label": "Genomics", "synonyms": ["Whole genomes"]}
METAGENOMICS = {
    "iri": "http://edamontology.org/topic_3174",
    "label": "Metagenomics",
    "synonyms": ["Ecogenomics", "Shotgun metagenomics"],
}
POPULATION = {
    "iri": "http://edamontology.org/topic_3796",
    "label": "Population genomics",
    "synonyms": None,
}


def test_a_term_contributes_its_label_and_every_synonym():
    assert _entries(METAGENOMICS) == {
        "metagenomics": "edam:topic_3174",
        "ecogenomics": "edam:topic_3174",
        "shotgun metagenomics": "edam:topic_3174",
    }


def test_a_term_without_synonyms_still_contributes_its_label():
    assert _entries(POPULATION) == {"population genomics": "edam:topic_3796"}


def test_entries_are_keyed_the_way_lookup_normalises():
    entry = _entries({"iri": DOMAIN_ROOT, "label": "  Whole   Genomes "})
    assert entry == {"whole genomes": "edam:topic_0622"}


def test_the_index_is_sorted_so_a_rebuild_diffs_cleanly():
    index = _index_of([POPULATION, METAGENOMICS])
    assert list(index) == sorted(index)


def test_a_multi_word_label_becomes_a_kebab_case_slug():
    assert _slug("Population genomics") == "population-genomics"


def test_the_domain_root_is_excluded_from_its_own_subtopics():
    assert _subtopics_of([GENOMICS, METAGENOMICS]) == ["metagenomics"]


def test_subtopics_are_deduplicated_and_sorted():
    assert _subtopics_of([POPULATION, METAGENOMICS, METAGENOMICS]) == [
        "metagenomics",
        "population-genomics",
    ]
