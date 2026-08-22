"""Hardcoded CURIE prefix to base-IRI table for the genomics adapter.

Nothing here is fetched at runtime: every expansion is pinned in this file, so
changing one changes how every genomics claim CURIE expands downstream.

Prefix spellings follow the sibling ``indicium-adapters-transcriptomics`` table
for every prefix the two share, so a CURIE minted by either adapter expands to
the same IRI. Two documented exceptions:

* ``edam`` is lower-case here and in ``indicium-adapters-proteomics`` because it
  is the prefix ``_edam.py`` mints; the transcriptomics table spells it
  ``EDAM``. A consumer merging both tables must carry both spellings.
* ``ClinVar`` has no counterpart there, and is split by accession family below.

Two families, and the distinction is load-bearing (ADR 0001):

* **Ontology prefixes** — the vocabularies claim text is resolved against.
  Expansions are OBO PURLs, or the ontology's own base IRI where it is not an
  OBO library member (EFO, EDAM).
* **Identifier-authority prefixes** — expansions for entity identifiers (gene,
  transcript, variant) and for the repository accessions that
  ``asb-domain-core``'s canonical ``accession_scanner`` mines. This table
  *expands* an accession that some other component already recognised; it never
  recognises one in prose. A second accession-detection table here would be the
  duplicate-then-drift failure ADR 0001 exists to prevent.

Database expansions use identifiers.org compact-URI resolution
(``https://identifiers.org/<registry-prefix>:``); the registry prefix inside the
URI is the authority's own string even where the dictionary key is spelled the
way the field writes it (``SRA`` -> ``insdc.sra``, ``ENA`` -> ``ena.embl``).
"""

from __future__ import annotations

from typing import Final

# --- Ontologies -------------------------------------------------------------

SO: Final[str] = "http://purl.obolibrary.org/obo/SO_"
GENO: Final[str] = "http://purl.obolibrary.org/obo/GENO_"
GO: Final[str] = "http://purl.obolibrary.org/obo/GO_"
OBI: Final[str] = "http://purl.obolibrary.org/obo/OBI_"
EFO: Final[str] = "http://www.ebi.ac.uk/efo/EFO_"
NCBITAXON: Final[str] = "http://purl.obolibrary.org/obo/NCBITaxon_"
UBERON: Final[str] = "http://purl.obolibrary.org/obo/UBERON_"
# Lower-case on purpose: it is the prefix `_edam.py` mints into every CURIE it
# returns, and `test_namespaces.py` binds the two together so they cannot drift.
EDAM: Final[str] = "http://edamontology.org/"

# --- Gene / transcript / variant identifier authorities ---------------------

HGNC: Final[str] = "https://identifiers.org/hgnc:"
ENSEMBL: Final[str] = "https://identifiers.org/ensembl:"
REFSEQ: Final[str] = "https://identifiers.org/refseq:"
NCBIGENE: Final[str] = "https://identifiers.org/ncbigene:"
DBSNP: Final[str] = "https://identifiers.org/dbsnp:"

# ClinVar issues three accession families and identifiers.org registers each
# separately, so one prefix cannot serve all three. `clinvar` accepts only the
# bare numeric variation id (`^\d+$`), never the `VCV` accession spelling, so
# the variation family expands through NCBI's own resolver instead.
CLINVAR_VARIATION: Final[str] = "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
CLINVAR_RECORD: Final[str] = "https://identifiers.org/clinvar.record:"
CLINVAR_SUBMISSION: Final[str] = "https://identifiers.org/clinvar.submission:"

# --- Sequence / expression repositories -------------------------------------

ENA: Final[str] = "https://identifiers.org/ena.embl:"
SRA: Final[str] = "https://identifiers.org/insdc.sra:"
GEO: Final[str] = "https://identifiers.org/geo:"
ARRAYEXPRESS: Final[str] = "https://identifiers.org/arrayexpress:"
BIOPROJECT: Final[str] = "https://identifiers.org/bioproject:"


GENOMICS_PREFIXES: Final[dict[str, str]] = {
    # Ontologies.
    "SO": SO,
    "GENO": GENO,
    "GO": GO,
    "OBI": OBI,
    "EFO": EFO,
    "NCBITaxon": NCBITAXON,
    "UBERON": UBERON,
    "edam": EDAM,
    # Gene / transcript / variant identifiers.
    "HGNC": HGNC,
    "ENSEMBL": ENSEMBL,
    "RefSeq": REFSEQ,
    "NCBIGene": NCBIGENE,
    "dbSNP": DBSNP,
    "ClinVar": CLINVAR_VARIATION,
    "ClinVar.record": CLINVAR_RECORD,
    "ClinVar.submission": CLINVAR_SUBMISSION,
    # Repositories.
    "ENA": ENA,
    "SRA": SRA,
    "GEO": GEO,
    "ArrayExpress": ARRAYEXPRESS,
    "BioProject": BIOPROJECT,
}
