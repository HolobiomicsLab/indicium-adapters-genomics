"""Integrity contract for the shipped EDAM concept index.

The index cannot be checked against EDAM offline, so what is pinned here is
everything that can go wrong without a network call: a CURIE that is not an
EDAM accession, a key the lookup normaliser could never produce, and — the one
that bit this package — a subdomain in the classification vocabulary that no
longer names a concept the index can resolve.
"""

from __future__ import annotations

import re

from indicium_adapters_genomics import _edam

EDAM_CURIE = re.compile(r"^edam:(topic|operation|data|format)_\d{4}$")

INDEX = _edam._payload()["index"]


def test_every_curie_is_an_edam_accession():
    bad = sorted(curie for curie in INDEX.values() if not EDAM_CURIE.match(curie))
    assert not bad, f"not EDAM accessions: {bad}"


def test_every_key_is_already_normalised():
    """A key the normaliser cannot produce is a lookup that can never fire."""
    unreachable = sorted(key for key in INDEX if _edam._normalise(key) != key)
    assert not unreachable, f"unreachable keys: {unreachable}"


def test_the_recorded_roots_are_edam_iris():
    assert _edam.roots()
    for iri in _edam.roots():
        assert iri.startswith(_edam.EDAM_PREFIX), iri


def test_every_subtopic_resolves_through_the_index():
    """Binds the classification vocabulary to the concept index, so a subdomain
    cannot survive an EDAM release that renamed or dropped its topic."""
    unresolved = [
        slug for slug in _edam.subtopics() if _edam.lookup(slug.replace("-", " ")) is None
    ]
    assert not unresolved, f"subdomains with no EDAM concept: {unresolved}"


def test_the_domain_root_is_not_listed_as_its_own_subtopic():
    assert "genomics" not in _edam.subtopics()


def test_size_reports_the_whole_index():
    assert _edam.size() == len(INDEX)
