"""GFF3 to asb:Claim projection (feature location, provenance and nesting).

Implementation note: the nine columns are positional and the attribute column
is a simple percent-escaped ``tag=value`` list, so a direct reader is smaller
and clearer than a dependency on a genome-annotation toolkit.
"""

from __future__ import annotations

from pathlib import Path
from typing import Final
from urllib.parse import unquote

from indicium_adapters_genomics._identifiers import resolve_sequence_ontology

COLUMNS: Final[tuple[str, ...]] = (
    "seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes",
)

FASTA_DIRECTIVE: Final[str] = "##FASTA"
MISSING_VALUE: Final[str] = "."

_COMPOSITE: Final[str] = "semapv:CompositeMatching"


def _parse_attributes(column: str) -> dict[str, str]:
    """Parse the ninth column into ``{tag: value}``, values still percent-encoded.

    Values stay encoded because the separator for a multi-valued attribute is an
    *unescaped* comma: decoding first would turn a ``%2C`` inside one value into
    a separator and split a single identifier in two. Callers run :func:`_decoded`
    after any comma split. Tags that repeat keep their first occurrence.
    """
    attributes: dict[str, str] = {}
    for pair in column.split(";"):
        tag, separator, value = pair.strip().partition("=")
        if separator and tag:
            attributes.setdefault(unquote(tag), value)
    return attributes


def _decoded(value: str) -> str:
    """Percent-decode one attribute value, per the GFF3 escaping rules."""
    return unquote(value.strip())


def _read_features(path: Path) -> list[dict[str, str]]:
    """Read a GFF3 file and return one dict per feature line.

    Comment and directive lines are skipped, and ``##FASTA`` ends the body.
    Lines with fewer than nine columns are not features and are ignored.
    """
    features: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith(FASTA_DIRECTIVE):
                break
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= len(COLUMNS):
                features.append(dict(zip(COLUMNS, fields, strict=False)))
    return features


def _feature_label(feature: dict[str, str], attributes: dict[str, str]) -> str:
    """Return the feature's ID, else its Name, else a positional fallback."""
    for tag in ("ID", "Name"):
        value = _decoded(attributes.get(tag, ""))
        if value:
            return value
    return f"{feature['seqid']}:{feature['start']}-{feature['end']}"


def _evidence_id(feature: dict[str, str], file_id: str) -> str:
    """Return the evidence id locating this feature in its source file."""
    return f"gff3:{file_id}:{feature['seqid']}:{feature['start']}-{feature['end']}"


def _located_on_claim(feature: dict[str, str], file_id: str, label: str) -> dict:
    """Build the ``located_on`` claim placing the feature on its sequence.

    The ``type`` column may be written either as an SO term name or as an SO
    accession. Only the accession form is carried onto the subject slot: a name
    would need an SO label index to resolve, and guessing one manufactures a
    wrong CURIE. Only SO is accepted -- the spec admits no other authority in
    this column.
    """
    claim = {
        "subject": label,
        "qualifier": "located_on",
        "object": feature["seqid"],
        "ontology_terms": {},
        "evidence_id": _evidence_id(feature, file_id),
    }
    curie = resolve_sequence_ontology(feature["type"])
    if curie:
        claim["ontology_terms"]["subject"] = curie
        claim["ontology_term_justifications"] = {"subject": _COMPOSITE}
    return claim


def _annotated_by_claim(feature: dict[str, str], file_id: str, label: str) -> dict | None:
    """Build the ``annotated_by`` claim, or None when no source is recorded.

    Column two is the algorithm or operating procedure that generated the
    feature, which is exactly what ``annotated_by`` asserts.
    """
    source = feature["source"].strip()
    if not source or source == MISSING_VALUE:
        return None
    return {
        "subject": label,
        "qualifier": "annotated_by",
        "object": source,
        "ontology_terms": {},
        "evidence_id": _evidence_id(feature, file_id),
    }


def _parent_claims(parents: str, evidence_id: str, label: str) -> list[dict]:
    """Build one ``has_feature`` claim per comma-separated Parent identifier.

    ``parents`` is the raw, still-encoded attribute value; it is split on the
    unescaped comma first and each part decoded after.
    """
    claims = []
    for parent in parents.split(","):
        identifier = _decoded(parent)
        if not identifier:
            continue
        claims.append({
            "subject": identifier,
            "qualifier": "has_feature",
            "object": label,
            "ontology_terms": {},
            "evidence_id": evidence_id,
        })
    return claims


def _feature_claims(feature: dict[str, str], file_id: str) -> list[dict]:
    """Build every claim one GFF3 feature line yields."""
    attributes = _parse_attributes(feature["attributes"])
    label = _feature_label(feature, attributes)
    claims = [_located_on_claim(feature, file_id, label)]
    annotated_by = _annotated_by_claim(feature, file_id, label)
    if annotated_by:
        claims.append(annotated_by)
    parents = attributes.get("Parent", "")
    claims.extend(_parent_claims(parents, _evidence_id(feature, file_id), label))
    return claims


def parse_gff3_to_claims(gff3_path: Path | str) -> list[dict]:
    """Parse a GFF3 file and return its features as a list of asb:Claim dicts.

    Each claim has the 5-slot shape: subject, qualifier, object, ontology_terms,
    evidence_id. Raises FileNotFoundError if the file does not exist.
    """
    path = Path(gff3_path)
    if not path.exists():
        raise FileNotFoundError(f"GFF3 file not found: {path}")
    file_id = path.stem
    claims: list[dict] = []
    for feature in _read_features(path):
        claims.extend(_feature_claims(feature, file_id))
    return claims
