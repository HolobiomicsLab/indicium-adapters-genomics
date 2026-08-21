"""GenomicsAdapter — the genomics DomainAdapter (api_version=1)."""

from __future__ import annotations

from typing import Final

from . import _edam

_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object", "qualifier")

_LEXICAL: Final[str] = "semapv:LexicalMatching"

# EDAM's Genomics topic DOES carry sub-topics, so unlike the sibling proteomics
# adapter this vocabulary is the field's own hierarchy rather than an editor's
# list. Read from EDAM at index-build time and frozen here so the runtime stays
# a pure function; regenerate with scripts/build_index.py on an EDAM release.
_EDAM_SUBTOPICS: Final[tuple[str, ...]] = (
    "comparative-genomics",
    "epigenomics",
    "epistasis",
    "functional-genomics",
    "metagenomics",
    "metatranscriptomics",
    "paleogenomics",
    "pharmacogenomics",
    "phylogenomics",
    "population-genomics",
    "proteogenomics",
    "structural-genomics",
    "transcriptomics",
)


class GenomicsAdapter:
    """Genome sequencing, assembly, variation and annotation.

    Owns the tier-3 concerns only (ADR 0001): classification vocabulary and
    claim/ontology enrichment. It deliberately does **not** recognise repository
    accessions (SRA, ENA, GEO, BioProject) — those belong to the canonical
    accession table in ``asb-domain-core``, and a second copy of those patterns
    here would be the duplicate-then-drift failure ADR 0001 exists to prevent.
    What this adapter contributes for identifiers is ``ontology_prefixes``, so a
    consumer can expand the CURIEs that tier produces.

    ``subdomains`` are read from EDAM's own Genomics branch rather than
    hand-listed, so the vocabulary is the field's and moves with an EDAM release.
    ``techniques`` are curated: EDAM models assays as operations and data rather
    than as a technique hierarchy, so there is nothing to derive.

    Every qualifier matches ``^[a-z][a-z0-9_]*$`` (snake_case), required because
    Perspicacité routes adapter qualifiers onto the open ``asb:domainQualifier``
    slot whose SHACL shape enforces that pattern.
    """

    domain_id = "genomics"
    api_version = 1

    primary_domain = "genomics"

    qualifiers: frozenset[str] = frozenset(
        {
            # --- Sequencing and assembly ---------------------------------------
            "sequenced_on",  # object is the sequencing platform
            "assembled_with",  # object is the assembler
            "aligned_to",  # object is the reference genome or index
            "basecalled_by",  # object is the basecaller (long-read)
            "polished_with",  # object is the consensus-polishing tool
            "scaffolded_by",  # object is the scaffolding evidence or tool
            # --- Variation ------------------------------------------------------
            "variant_called_by",  # object is the variant caller
            "genotyped_as",  # object is the genotype or allele
            "phased_by",  # object is the phasing method
            "imputed_from",  # object is the reference panel
            "associated_with",  # object is a trait or phenotype (GWAS)
            "segregates_with",  # object is a pedigree or population
            # --- Annotation and function ----------------------------------------
            "annotated_by",  # object is the annotation pipeline or database
            "predicted_to_encode",  # object is a gene product or feature
            "orthologous_to",
            "located_on",  # object is a chromosome, contig or scaffold
            "regulated_by",  # object is a regulatory element or factor
            # --- Validation ------------------------------------------------------
            "validated_by",  # object is an orthogonal assay
        }
    )

    ontology_prefixes: dict[str, str] = {
        "SO": "http://purl.obolibrary.org/obo/SO_",
        "GENO": "http://purl.obolibrary.org/obo/GENO_",
        "GO": "http://purl.obolibrary.org/obo/GO_",
        "OBI": "http://purl.obolibrary.org/obo/OBI_",
        "EFO": "http://www.ebi.ac.uk/efo/EFO_",
        "NCBITaxon": "http://purl.obolibrary.org/obo/NCBITaxon_",
        "UBERON": "http://purl.obolibrary.org/obo/UBERON_",
        "ensembl": "http://identifiers.org/ensembl/",
        "refseq": "http://identifiers.org/refseq/",
        "edam": "http://edamontology.org/",
    }

    @property
    def subdomains(self) -> tuple[str, ...]:
        """EDAM's own sub-topics of Genomics, lower-cased and hyphenated."""
        return _EDAM_SUBTOPICS

    # Curated — EDAM has no technique hierarchy to read.
    techniques: tuple[str, ...] = (
        "whole-genome-sequencing",
        "whole-exome-sequencing",
        "targeted-sequencing",
        "long-read-sequencing",
        "short-read-sequencing",
        "single-cell-sequencing",
        "hi-c",
        "optical-mapping",
        "genotyping-array",
        "amplicon-sequencing",
        "shotgun-metagenomics",
        "rna-seq",
        "chip-seq",
        "atac-seq",
        "bisulfite-sequencing",
        "crispr-screen",
    )

    # ---- protocol methods --------------------------------------------------

    def extraction_context(self) -> str:
        return (
            "Domain: genomics — sequencing, assembly, variation and annotation. "
            "Name subjects and objects with EDAM and Sequence Ontology "
            "terminology where one applies (for example 'sequence assembly', "
            "'variant calling'), using the singular preferred label rather than a "
            "paraphrase. Prefer a concrete tool, assembly or panel name over a "
            "generic word: 'GRCh38' rather than 'the reference', 'GATK "
            "HaplotypeCaller' rather than 'the caller'. Keep a gene, its "
            "transcript and its protein product distinct. Domain qualifiers "
            "available in addition to the base set: "
            + ", ".join(sorted(self.qualifiers))
            + "."
        )

    def enrich_claim(self, claim: dict) -> dict:
        """Attach EDAM CURIEs to the slots whose text names an EDAM concept.

        Pure: returns a new dict and never mutates the input. A slot that carries
        text but resolves to nothing is recorded under
        ``ontology_terms_unmatched`` rather than passed over, so "not an EDAM
        concept" stays distinguishable from "no attempt was made".
        """
        out = dict(claim)
        terms: dict[str, str] = dict(out.get("ontology_terms") or {})
        justifications: dict[str, str] = dict(out.get("ontology_term_justifications") or {})
        unmatched: dict[str, str] = dict(out.get("ontology_terms_unmatched") or {})

        for slot in _ENRICHABLE_SLOTS:
            if slot in terms:
                continue
            value = out.get(slot)
            if not isinstance(value, str) or not value.strip():
                continue
            curie = _edam.lookup(value)
            if curie:
                terms[slot] = curie
                justifications[slot] = _LEXICAL
                unmatched.pop(slot, None)
            else:
                unmatched[slot] = value

        if terms:
            out["ontology_terms"] = terms
            out["ontology_term_justifications"] = justifications
        if unmatched:
            out["ontology_terms_unmatched"] = unmatched
        return out
