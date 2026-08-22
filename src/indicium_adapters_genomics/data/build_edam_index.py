"""Maintainer script: regenerate ``edam_index.json`` from EDAM via EBI OLS4.

The index is the label spine ``_edam.py`` matches claim text against, and the
same walk emits the ``subtopics`` list that becomes ``GenomicsAdapter.subdomains``.
It is cut from a recorded set of EDAM roots: each root contributes itself and its
direct subclasses, and each of those contributes its preferred label plus every
exact synonym EDAM records. Nothing is hand-added -- a label that EDAM does not
carry resolves to a CURIE EDAM never issued, which is the wrong-CURIE failure the
whole package is built to avoid.

Usage (never run at install or import time)::

    python -m indicium_adapters_genomics.data.build_edam_index

Re-run on an EDAM release. The roots below are the only vocabulary decision in
the file; everything else is read from EDAM.
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Final

OLS_TERMS_URL: Final[str] = "https://www.ebi.ac.uk/ols4/api/ontologies/edam/terms"

EDAM_IRI_PREFIX: Final[str] = "http://edamontology.org/"

EDAM_CURIE_PREFIX: Final[str] = "edam"

# Genomics. It names the domain itself, so it is the one collected term that is
# not also a subdomain of it.
DOMAIN_ROOT: Final[str] = f"{EDAM_IRI_PREFIX}topic_0622"

# Metatranscriptomics and Epistasis: genomics papers claim both, but EDAM files
# them under Transcriptomics and Genotype-and-phenotype, so neither is reachable
# from the domain root and each has to be named.
ADDITIONAL_ROOTS: Final[tuple[str, ...]] = (
    f"{EDAM_IRI_PREFIX}topic_3941",
    f"{EDAM_IRI_PREFIX}topic_3974",
)

ROOTS: Final[tuple[str, ...]] = (DOMAIN_ROOT, *ADDITIONAL_ROOTS)

SOURCE: Final[str] = "EDAM (edamontology.org/EDAM.csv); entries verified against EBI OLS4"

CHILD_PAGE_SIZE: Final[int] = 200

REQUEST_TIMEOUT_SECONDS: Final[int] = 60


def _normalise(text: str) -> str:
    """Case-fold and collapse whitespace, matching ``_edam.lookup``'s key form."""
    return " ".join(str(text).lower().split())


def _fetch(url: str) -> dict:
    """GET ``url`` and return the decoded JSON body."""
    with urllib.request.urlopen(url, timeout=REQUEST_TIMEOUT_SECONDS) as response:
        return json.load(response)


def _term(iri: str) -> dict:
    """Return the OLS4 term record for one EDAM ``iri``."""
    query = urllib.parse.urlencode({"iri": iri})
    terms = _fetch(f"{OLS_TERMS_URL}?{query}").get("_embedded", {}).get("terms", [])
    if not terms:
        raise RuntimeError(f"EDAM has no term {iri}")
    return terms[0]


def _children(term: dict) -> list[dict]:
    """Return the direct subclasses of ``term``, or [] when it is a leaf."""
    link = term["_links"].get("hierarchicalChildren", {}).get("href")
    if not link:
        return []
    page = _fetch(f"{link}?size={CHILD_PAGE_SIZE}")
    return page.get("_embedded", {}).get("terms", [])


def _entries(term: dict) -> dict[str, str]:
    """Map every label and exact synonym of ``term`` to its ``edam:`` CURIE."""
    local = term["iri"].removeprefix(EDAM_IRI_PREFIX)
    curie = f"{EDAM_CURIE_PREFIX}:{local}"
    names = [term["label"], *(term.get("synonyms") or [])]
    return {_normalise(name): curie for name in names if str(name).strip()}


def _slug(label: str) -> str:
    """Return ``label`` as the kebab-case slug the classification vocabulary uses."""
    return "-".join(_normalise(label).split())


def _collect(roots: tuple[str, ...]) -> list[dict]:
    """Return every EDAM term the index covers: each root and its direct subclasses."""
    terms: dict[str, dict] = {}
    for iri in roots:
        root = _term(iri)
        for term in (root, *_children(root)):
            terms.setdefault(term["iri"], term)
    return list(terms.values())


def _index_of(terms: list[dict]) -> dict[str, str]:
    """Map every label and exact synonym across ``terms`` to its ``edam:`` CURIE."""
    index: dict[str, str] = {}
    for term in terms:
        index.update(_entries(term))
    return dict(sorted(index.items()))


def _subtopics_of(terms: list[dict]) -> list[str]:
    """Return the preferred-label slug of every term except the domain root."""
    return sorted({_slug(term["label"]) for term in terms if term["iri"] != DOMAIN_ROOT})


def main() -> None:
    """Rebuild ``edam_index.json`` next to this script from :data:`ROOTS`."""
    terms = _collect(ROOTS)
    payload = {
        "index": _index_of(terms),
        "subtopics": _subtopics_of(terms),
        "roots": list(ROOTS),
        "source": SOURCE,
    }
    target = Path(__file__).parent / "edam_index.json"
    target.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
    print(f"[edam] wrote {target} with {len(payload['index']):,} entries", file=sys.stderr)


if __name__ == "__main__":
    main()
