"""VCF to asb:Claim projection (fixed columns + per-sample genotypes).

Implementation note: we parse the tab-separated body directly rather than take
a dependency on a VCF library. The projection only needs the eight fixed
columns and the GT subfield, all of which are positionally fixed by the spec,
and a parser dependency would pull htslib bindings into a pure-vocabulary
package.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Final, TextIO

from indicium_adapters_genomics._identifiers import resolve as resolve_identifier

HEADER_MARKER: Final[str] = "#CHROM"
FORMAT_COLUMN: Final[str] = "FORMAT"
GENOTYPE_KEY: Final[str] = "GT"
MISSING_VALUE: Final[str] = "."

_COMPOSITE: Final[str] = "semapv:CompositeMatching"


def _open_text(path: Path) -> TextIO:
    """Open a plain or gzip-compressed VCF for text reading."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _file_id(path: Path) -> str:
    """Return the evidence-id stem, with any ``.gz`` suffix removed first."""
    stem = path.name[: -len(".gz")] if path.name.endswith(".gz") else path.name
    return Path(stem).stem


def _read_records(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a VCF, returning its column names and one dict per data line.

    ``##`` meta-information lines are skipped; the ``#CHROM`` line supplies the
    column names. Rows shorter than the header keep only the columns present.
    """
    columns: list[str] = []
    records: list[dict[str, str]] = []
    with _open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith(HEADER_MARKER):
                columns = line.lstrip("#").split("\t")
            elif columns:
                records.append(dict(zip(columns, line.split("\t"), strict=False)))
    return columns, records


def _sample_columns(columns: list[str]) -> tuple[str, ...]:
    """Return the per-sample column names, i.e. every column after FORMAT."""
    if FORMAT_COLUMN not in columns:
        return ()
    return tuple(columns[columns.index(FORMAT_COLUMN) + 1 :])


def _variant_ids(record: dict[str, str]) -> tuple[str, ...]:
    """Return every identifier in the record's ID column, empty when missing.

    VCF 4.3 makes ID a semicolon-separated list of identifiers that all name the
    same variant, so every entry is kept rather than only the first.
    """
    identifier = record.get("ID", "").strip()
    if not identifier or identifier == MISSING_VALUE:
        return ()
    return tuple(part.strip() for part in identifier.split(";") if part.strip())


def _variant_label(record: dict[str, str]) -> str:
    """Return the record's first ID, or a positional ``CHROM:POSREF>ALT`` fallback.

    One of the identifiers has to be the subject the other claims hang off, and
    the first is the stable choice; the rest still reach ``identifies`` claims.
    """
    identifiers = _variant_ids(record)
    if identifiers:
        return identifiers[0]
    chrom, pos = record.get("CHROM", ""), record.get("POS", "")
    return f"{chrom}:{pos}{record.get('REF', '')}>{record.get('ALT', '')}"


def _evidence_id(record: dict[str, str], file_id: str) -> str:
    """Return the evidence id locating this record in its source file."""
    return f"vcf:{file_id}:{record.get('CHROM', '')}:{record.get('POS', '')}"


def _located_on_claim(record: dict[str, str], file_id: str) -> dict:
    """Build the ``located_on`` claim placing the variant on its chromosome."""
    return {
        "subject": _variant_label(record),
        "qualifier": "located_on",
        "object": record.get("CHROM", ""),
        "ontology_terms": {},
        "evidence_id": _evidence_id(record, file_id),
    }


def _identifies_claim(label: str, curie: str, evidence_id: str) -> dict:
    """Build one ``identifies`` claim tying ``label`` to the CURIE ``curie``."""
    return {
        "subject": label,
        "qualifier": "identifies",
        "object": curie,
        "ontology_terms": {"object": curie},
        "ontology_term_justifications": {"object": _COMPOSITE},
        "evidence_id": evidence_id,
    }


def _identifies_claims(record: dict[str, str], file_id: str) -> list[dict]:
    """Build one ``identifies`` claim per ID entry that names a known authority.

    The ID column is written by the calling pipeline, so a CURIE derived from it
    is a curated identifier taken verbatim rather than a lexical match. An entry
    no grammar recognises is skipped, not guessed at.
    """
    label = _variant_label(record)
    evidence_id = _evidence_id(record, file_id)
    curies = (resolve_identifier(entry) for entry in _variant_ids(record))
    return [_identifies_claim(label, curie, evidence_id) for curie in curies if curie]


def _subfield_at(value: str, position: int) -> str:
    """Return the colon-separated subfield of ``value`` at ``position``, or ''."""
    fields = value.split(":")
    return fields[position].strip() if position < len(fields) else ""


def _is_called(genotype: str) -> bool:
    """Whether a GT value names at least one called allele (``./.`` does not)."""
    return any(char.isdigit() for char in genotype)


def _genotype_claims(record: dict[str, str], file_id: str, samples: tuple[str, ...]) -> list[dict]:
    """Build one ``genotyped_as`` claim per sample with a called genotype."""
    keys = record.get(FORMAT_COLUMN, "").split(":")
    if GENOTYPE_KEY not in keys:
        return []
    position = keys.index(GENOTYPE_KEY)
    label = _variant_label(record)
    claims = []
    for sample in samples:
        genotype = _subfield_at(record.get(sample, ""), position)
        if not _is_called(genotype):
            continue
        claims.append({
            "subject": label,
            "qualifier": "genotyped_as",
            "object": genotype,
            "ontology_terms": {},
            "evidence_id": f"{_evidence_id(record, file_id)}:{sample}",
        })
    return claims


def parse_vcf_to_claims(vcf_path: Path | str) -> list[dict]:
    """Parse a VCF file and return its records as a list of asb:Claim dicts.

    Each claim has the 5-slot shape: subject, qualifier, object, ontology_terms,
    evidence_id. Accepts plain or gzip-compressed input. Raises FileNotFoundError
    if the file does not exist.
    """
    path = Path(vcf_path)
    if not path.exists():
        raise FileNotFoundError(f"VCF file not found: {path}")
    columns, records = _read_records(path)
    file_id = _file_id(path)
    samples = _sample_columns(columns)
    claims: list[dict] = []
    for record in records:
        claims.append(_located_on_claim(record, file_id))
        claims.extend(_identifies_claims(record, file_id))
        claims.extend(_genotype_claims(record, file_id, samples))
    return claims
