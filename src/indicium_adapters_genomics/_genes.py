"""Gene-symbol resolution against an HGNC-derived index.

A gene symbol has no distinguishing shape -- ``SET``, ``MAX`` and ``REST`` are
ordinary English words -- so it can only be resolved against a nomenclature
authority. This module reads one bundled table, ``data/gene_index.tsv``,
regenerated from the HGNC complete set by ``data/build_gene_index.py``.

**The table is not committed.** A hand-picked subset of symbols would resolve
the genes someone happened to look at and silently return None for every other
gene, which is the failure mode the ADR 0001 anti-patterns describe: a lookup
that manufactures confidence by looking healthy on the cases it was written
against. So the index is either the authority's full set or it is absent, and
absence is inspectable through :attr:`GeneIndex.is_available` and ``len()``
rather than being indistinguishable from "checked, no such gene".

Callers that already hold nomenclature data can bypass the file entirely by
constructing ``GeneIndex(rows=...)`` and passing it to ``GenomicsAdapter``.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping
from importlib.resources import files
from pathlib import Path
from typing import Final, NamedTuple

GENE_INDEX_FILENAME: Final[str] = "gene_index.tsv"

# Column in the index, and the CURIE prefix to mint from it. Order is the
# principal-CURIE preference: HGNC is the nomenclature authority, so it wins
# over the two database cross-references when a row carries all three.
_CURIE_COLUMNS: Final[tuple[tuple[str, str], ...]] = (
    ("hgnc", "HGNC"),
    ("ensembl", "ENSEMBL"),
    ("ncbigene", "NCBIGene"),
)


class GeneMatch(NamedTuple):
    """One resolved symbol: its principal CURIE and every cross-reference."""

    principal: str
    xrefs: dict[str, str]


def _to_curie(value: str, prefix: str) -> str:
    """Return ``value`` as a CURIE, prefixing it unless it already carries one."""
    return value if ":" in value else f"{prefix}:{value}"


def _match_from_row(row: Mapping[str, str]) -> GeneMatch | None:
    """Build a :class:`GeneMatch` from one index row, or None if it has no id."""
    xrefs = {
        prefix: _to_curie(row[column].strip(), prefix)
        for column, prefix in _CURIE_COLUMNS
        if row.get(column, "").strip()
    }
    if not xrefs:
        return None
    principal = next(iter(xrefs.values()))
    return GeneMatch(principal, xrefs)


def _rows_from_tsv(handle: Iterable[str]) -> dict[str, Mapping[str, str]]:
    """Read symbol-keyed rows from an open HGNC-derived TSV, first row winning."""
    rows: dict[str, Mapping[str, str]] = {}
    for row in csv.DictReader(handle, delimiter="\t"):
        symbol = (row.get("symbol") or "").strip()
        if symbol:
            rows.setdefault(symbol, row)
    return rows


class GeneIndex:
    """Symbol to HGNC / Ensembl / NCBI Gene lookup, loaded lazily and once.

    Construct with no arguments to read the bundled ``gene_index.tsv``, or pass
    ``rows`` (``{symbol: {column: value}}``) to supply nomenclature data from
    elsewhere and skip the file.
    """

    def __init__(self, rows: Mapping[str, Mapping[str, str]] | None = None) -> None:
        self._rows: dict[str, Mapping[str, str]] = dict(rows or {})
        self._is_loaded: bool = rows is not None

    @classmethod
    def from_tsv(cls, path: Path | str) -> GeneIndex:
        """Build an index from an HGNC-derived TSV anywhere on disk.

        Lets a deployment point at its own nomenclature dump instead of the
        bundled resource. ``path`` is the TSV; returns the loaded index.
        """
        with Path(path).open("r", encoding="utf-8", newline="") as handle:
            return cls(_rows_from_tsv(handle))

    def _load(self) -> None:
        """Read the bundled index on first use; a missing file leaves it empty."""
        if self._is_loaded:
            return
        self._is_loaded = True
        resource = files("indicium_adapters_genomics.data").joinpath(GENE_INDEX_FILENAME)
        if not resource.is_file():
            return
        with resource.open("r", encoding="utf-8", newline="") as handle:
            self._rows.update(_rows_from_tsv(handle))

    def __len__(self) -> int:
        """Return the number of symbols available for lookup."""
        self._load()
        return len(self._rows)

    @property
    def is_available(self) -> bool:
        """Whether any nomenclature data is loaded at all.

        False means "no index installed", which is a different state from a
        symbol that the authority does not list, and callers should say so.
        """
        return len(self) > 0

    def resolve(self, symbol: str) -> GeneMatch | None:
        """Return the :class:`GeneMatch` for ``symbol``, or None if unlisted.

        Lookup is exact, then retried upper-cased, because HGNC symbols are
        upper-case for human genes and papers vary the casing. No substring or
        fuzzy step: a wrong gene is worse than an unresolved one.
        """
        token = symbol.strip()
        if not token:
            return None
        self._load()
        row = self._rows.get(token) or self._rows.get(token.upper())
        return _match_from_row(row) if row else None


GENE_INDEX: Final[GeneIndex] = GeneIndex()
"""Module-level singleton reading the bundled index, shared across adapters."""
