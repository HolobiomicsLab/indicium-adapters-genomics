"""Structural recognition of genomics entity identifiers.

Genomics claim slots routinely carry an identifier rather than a phrase --
``rs334``, ``ENSG00000139618``, ``NM_000546.6``, ``VCV000012345``. Each of those
is issued under a published grammar, so recognising one is a syntactic fact
about the string, not a guess from a vocabulary. Every spec below names the
authority that defines its shape:

===============  =========================================================
Prefix           Authority for the grammar
===============  =========================================================
``dbSNP``        NCBI dbSNP reference-SNP identifiers (``rs`` + digits)
``ClinVar``      NCBI ClinVar variation accessions (``VCV`` + digits)
``ClinVar.record``      NCBI ClinVar records (``RCV`` + digits)
``ClinVar.submission``  NCBI ClinVar submissions (``SCV`` + digits)
``ENSEMBL``      Ensembl stable-ID scheme (ENS + species code + type letter)
``RefSeq``       NCBI RefSeq accession-prefix table
``HGNC``         HUGO Gene Nomenclature Committee identifiers
``SO`` ``GO``    OBO Foundry CURIEs, passed through unchanged
===============  =========================================================

Two rules keep this from becoming a second accession scanner (ADR 0001):

* Every pattern is **full-match against one claim-slot value**. Nothing here
  scans prose, so it cannot compete with ``asb-domain-core``'s canonical
  ``accession_scanner``, which owns *deposition accessions mined from text*.
* Nothing here is a **repository** accession. GEO / SRA / ENA / ArrayExpress /
  BioProject appear in ``_namespaces.py`` as expansions only.

No gene symbol is recognised here, because a symbol has no distinguishing
shape -- ``SET``, ``MAX`` and ``REST`` are ordinary English words. Symbols need
an authority index, which is ``_genes.py``'s job.
"""

from __future__ import annotations

import re
from typing import Final, NamedTuple


class _IdentifierSpec(NamedTuple):
    """One identifier grammar: which prefix it mints, and who defines it."""

    prefix: str
    authority: str
    pattern: re.Pattern[str]


def _curie_spec(prefix: str, local: str, authority: str) -> _IdentifierSpec:
    """Build a spec for a CURIE that is passed through as written.

    ``prefix`` is the CURIE prefix, ``local`` a regex for the local part, and
    ``authority`` the body that issues the identifier. Returns the spec.
    """
    return _IdentifierSpec(prefix, authority, re.compile(rf"{prefix}:(?P<id>{local})"))


_OBO_LOCAL: Final[str] = r"\d{7}"

# NCBI RefSeq accession prefixes taking an all-digit local part: genomic
# (AC/NC/NG/NT/NW), transcript (NM/NR/XM/XR) and protein (AP/NP/XP/YP/WP).
_REFSEQ_PREFIXES: Final[str] = "AC|NC|NG|NT|NW|NM|NR|XM|XR|AP|NP|XP|YP|WP"

# NZ_ is the one RefSeq family whose local part is not all digits: a WGS
# accession carries the 2-4 letter project code first, as NZ_CP012345 does.
_REFSEQ_ID: Final[str] = (
    rf"(?:(?:{_REFSEQ_PREFIXES})_\d{{6,}}|NZ_[A-Z]{{2,4}}\d{{6,}})(?:\.\d+)?"
)

_SO_PREFIX: Final[str] = "SO"

# Ensembl stable id: ENS, an optional species code, a feature-type letter
# (E xon, F amily, G ene, P rotein, R egulatory, T ranscript), 11 digits.
_ENSEMBL_ID: Final[str] = r"ENS[A-Z]{0,6}[EFGPRT]\d{11}(?:\.\d+)?"

_SPECS: Final[tuple[_IdentifierSpec, ...]] = (
    _IdentifierSpec("dbSNP", "NCBI dbSNP", re.compile(r"(?P<id>rs\d+)")),
    _IdentifierSpec("ClinVar", "NCBI ClinVar", re.compile(r"(?P<id>VCV\d{9,}(?:\.\d+)?)")),
    _IdentifierSpec(
        "ClinVar.record", "NCBI ClinVar", re.compile(r"(?P<id>RCV\d{9,}(?:\.\d+)?)")
    ),
    _IdentifierSpec(
        "ClinVar.submission", "NCBI ClinVar", re.compile(r"(?P<id>SCV\d{9,}(?:\.\d+)?)")
    ),
    _IdentifierSpec("ENSEMBL", "Ensembl", re.compile(rf"(?P<id>{_ENSEMBL_ID})")),
    _IdentifierSpec("RefSeq", "NCBI RefSeq", re.compile(rf"(?P<id>{_REFSEQ_ID})")),
    _curie_spec("HGNC", r"\d+", "HUGO Gene Nomenclature Committee"),
    _curie_spec(_SO_PREFIX, _OBO_LOCAL, "Sequence Ontology"),
    _curie_spec("GO", _OBO_LOCAL, "Gene Ontology"),
    _curie_spec("OBI", _OBO_LOCAL, "Ontology for Biomedical Investigations"),
    _curie_spec("GENO", _OBO_LOCAL, "Genotype Ontology"),
    _curie_spec("UBERON", _OBO_LOCAL, "Uberon anatomy ontology"),
    _curie_spec("EFO", _OBO_LOCAL, "EBI Experimental Factor Ontology"),
    _curie_spec("NCBITaxon", r"\d+", "NCBI Taxonomy"),
)


def resolve(text: str) -> str | None:
    """Return the CURIE for ``text`` when it *is* an identifier, else None.

    ``text`` is one claim-slot value. Matching is case-sensitive and anchored to
    the whole stripped string: the authorities all issue their identifiers in a
    fixed case, and a looser match would mint a confident-looking wrong CURIE
    out of ordinary prose.
    """
    token = text.strip()
    if not token:
        return None
    for spec in _SPECS:
        hit = spec.pattern.fullmatch(token)
        if hit:
            return f"{spec.prefix}:{hit.group('id')}"
    return None


def resolve_sequence_ontology(text: str) -> str | None:
    """Return the SO CURIE for ``text`` when it is one, else None.

    GFF3's type column admits only an SO term name or an SO accession, so the
    general :func:`resolve` is too wide there: it would accept a CURIE from any
    of the eleven other authorities and mint a feature type SO never issued.
    """
    curie = resolve(text)
    return curie if curie and curie.startswith(f"{_SO_PREFIX}:") else None


def authorities() -> dict[str, str]:
    """Return ``{CURIE prefix: issuing authority}`` for every grammar known here."""
    return {spec.prefix: spec.authority for spec in _SPECS}
