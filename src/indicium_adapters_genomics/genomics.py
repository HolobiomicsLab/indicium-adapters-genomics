"""GenomicsAdapter -- the genomics DomainAdapter (api_version=1)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, NamedTuple

from . import _edam, _identifiers
from ._genes import GENE_INDEX, GeneIndex
from ._namespaces import GENOMICS_PREFIXES

_QUALIFIER_SLOT: Final[str] = "qualifier"

_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object", _QUALIFIER_SLOT)

_LEXICAL: Final[str] = "semapv:LexicalMatching"
_COMPOSITE: Final[str] = "semapv:CompositeMatching"


class _Resolution(NamedTuple):
    """One slot's resolution: the CURIE, its SEMAPV justification, its xrefs."""

    curie: str
    justification: str
    xrefs: dict[str, str]


@dataclass
class _Enrichment:
    """The four enrichment maps of a claim, accumulated slot by slot."""

    terms: dict[str, str]
    justifications: dict[str, str]
    unmatched: dict[str, str]
    xrefs: dict[str, dict[str, str]]

    @classmethod
    def from_claim(cls, claim: dict) -> _Enrichment:
        """Copy the maps out of ``claim`` so the caller's dicts are never shared."""
        return cls(
            terms=dict(claim.get("ontology_terms") or {}),
            justifications=dict(claim.get("ontology_term_justifications") or {}),
            unmatched=dict(claim.get("ontology_terms_unmatched") or {}),
            xrefs=dict(claim.get("xrefs") or {}),
        )

    def record(self, slot: str, resolution: _Resolution) -> None:
        """Attach ``resolution`` to ``slot``, clearing any earlier miss on it."""
        self.terms[slot] = resolution.curie
        self.justifications[slot] = resolution.justification
        self.unmatched.pop(slot, None)
        if resolution.xrefs:
            self.xrefs[slot] = resolution.xrefs

    def miss(self, slot: str, text: str) -> None:
        """Record that ``slot`` carried text which resolved to nothing."""
        self.unmatched[slot] = text

    def apply_to(self, claim: dict) -> None:
        """Write every non-empty map back onto ``claim``."""
        for key, value in (
            ("ontology_terms", self.terms),
            ("ontology_term_justifications", self.justifications),
            ("ontology_terms_unmatched", self.unmatched),
            ("xrefs", self.xrefs),
        ):
            if value:
                claim[key] = value


def _resolve_text(text: str, gene_index: GeneIndex) -> _Resolution | None:
    """Resolve one slot's text to a CURIE, or None when no resolver claims it.

    The three resolvers run in decreasing order of certainty: an identifier
    grammar first, because a string that *is* an rsID is a syntactic fact rather
    than a guess and is therefore justified as a composite match; then the gene
    nomenclature index; then the EDAM concept labels. Every step is exact --
    substring or fuzzy matching over an authority manufactures confident-looking
    wrong CURIEs, and a wrong identifier is worse than no identifier.
    """
    identifier = _identifiers.resolve(text)
    if identifier:
        return _Resolution(identifier, _COMPOSITE, {})
    gene = gene_index.resolve(text)
    if gene:
        return _Resolution(gene.principal, _LEXICAL, gene.xrefs)
    concept = _edam.lookup(text)
    if concept:
        return _Resolution(concept, _LEXICAL, {})
    return None


class GenomicsAdapter:
    """Genome sequencing, assembly, variation and annotation.

    Owns the tier-3 concerns only (ADR 0001): classification vocabulary and
    claim/ontology enrichment. It deliberately does **not** recognise repository
    accessions (SRA, ENA, GEO, BioProject) in article text -- those belong to the
    canonical accession table in ``asb-domain-core``, and a second copy of those
    patterns here would be the duplicate-then-drift failure ADR 0001 exists to
    prevent. What this adapter contributes for identifiers is ``ontology_prefixes``
    so a consumer can expand the CURIEs that tier produces, plus full-match
    recognition of *entity* identifiers (variant, gene, transcript, protein)
    already sitting in a claim slot.

    ``subdomains`` and ``techniques`` are both curated. Neither is derived: the
    EDAM walk in ``data/build_edam_index.py`` reaches a root and its direct
    subclasses only, and sorts what it finds, so it can produce neither the
    audited set nor its order. ``_edam.subtopics()`` still records what that
    walk yields, but it no longer feeds this declaration.

    Every qualifier matches ``^[a-z][a-z0-9_]*$`` (snake_case), required because
    Perspicacite routes adapter qualifiers onto the open ``asb:domainQualifier``
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
            "has_assembly_metric",  # object is a contiguity or completeness value
            # --- Variation ------------------------------------------------------
            "variant_called_by",  # object is the variant caller
            "genotyped_as",  # object is the genotype or allele
            "phased_by",  # object is the phasing method
            "imputed_from",  # object is the reference panel
            "has_consequence",  # object is a Sequence Ontology consequence term
            "associated_with",  # object is a trait or phenotype (GWAS)
            "segregates_with",  # object is a pedigree or population
            # --- Annotation and function ----------------------------------------
            "annotated_by",  # object is the annotation pipeline or database
            "identifies",  # object is the authority identifier for the subject
            "predicted_to_encode",  # object is a gene product or feature
            "orthologous_to",
            "located_on",  # object is a chromosome, contig or scaffold
            "has_feature",  # object is a nested child feature
            "regulated_by",  # object is a regulatory element or factor
            # --- Expression ------------------------------------------------------
            "differentially_expressed_in",  # object is the contrast or condition
            # --- Validation ------------------------------------------------------
            "validated_by",  # object is an orthogonal assay
        }
    )

    ontology_prefixes: dict[str, str] = dict(GENOMICS_PREFIXES)

    # Fixed by the 2026-08-25 cross-domain vocabulary audit: a term two
    # sciences both declare cannot tell them apart, so shared terms were dropped.
    subdomains: tuple[str, ...] = (
        "genetic-variation",
        "dna-polymorphism",
        "structural-variation",
        "gene-structure",
        "mobile-genetic-elements",
        "comparative-genomics",
        "population-genomics",
        # Issue #1: eight in-domain EDAM topics the harvest carries but this
        # adapter never uttered. Together with the four techniques added below,
        # they take the genomics silent share over genomics_verified.json from
        # 51.5% to 38.1% of judgeable rows (admitted 1,921 -> 2,459).
        "genetics",
        "phylogenetics",
        "population-genetics",
        "phylogeny",
        "gwas-study",
        "phylogenomics",
        "copy-number-variation",
        "molecular-genetics",
    )

    # Same audit, same rule -- each survivor is a form only genomics utters.
    techniques: tuple[str, ...] = (
        "genome-assembly",
        "genome-annotation",
        "genotyping-array",
        "nanopore-sequencing",
        "smrt-sequencing",
        # Restored after the first audit pass dropped them as collateral rather
        # than by measurement: neither is claimed by any other installed
        # science, and without them a whole-genome or exome tool matches
        # nothing its own field declares.
        "whole-genome-sequencing",
        "exome-sequencing",
        # Issue #1, same measurement: four more measured in-domain topics,
        # none claimed by any other installed science. The issue's fifth
        # candidate, the bare term "mapping", is deliberately left out: a
        # single common word is over-permissive under the containment matcher
        # domain_admission uses, and it buys only a further 38.1% -> 36.7%.
        "genotype-and-phenotype",
        "sequence-assembly",
        "oxford-nanopore",
        "pacbio",
    )

    def __init__(self, gene_index: GeneIndex | None = None) -> None:
        """Build the adapter, optionally over a caller-supplied gene index.

        ``gene_index`` defaults to the module-level index, which reads
        ``data/gene_index.tsv`` if a maintainer has generated one and is empty
        otherwise -- the file is deliberately not shipped (see :mod:`._genes`).
        Pass one to resolve symbols against nomenclature data the caller already
        holds instead.
        """
        self._gene_index = gene_index or GENE_INDEX

    @property
    def gene_index(self) -> GeneIndex:
        """The index backing gene-symbol resolution.

        Read ``.is_available`` before reporting a symbol as unresolved: an empty
        index means "no nomenclature installed", not "no such gene".
        """
        return self._gene_index

    # ---- protocol methods --------------------------------------------------

    def extraction_context(self) -> str:
        """Return the domain instructions handed to the extraction model."""
        return (
            "Domain: genomics -- sequencing, assembly, variation and annotation. "
            "Name subjects and objects with EDAM and Sequence Ontology "
            "terminology where one applies (for example 'sequence assembly', "
            "'variant calling'), using the singular preferred label rather than a "
            "paraphrase. Prefer a concrete tool, assembly or panel name over a "
            "generic word: 'GRCh38' rather than 'the reference', 'GATK "
            "HaplotypeCaller' rather than 'the caller'. Keep a gene, its "
            "transcript and its protein product distinct. Quote an identifier in "
            "the form its authority issues it -- rs334, ENSG00000139618, "
            "NM_000546.6, VCV000012345 -- rather than describing it. Domain "
            "qualifiers available in addition to the base set: "
            + ", ".join(sorted(self.qualifiers))
            + "."
        )

    def _is_own_vocabulary(self, slot: str, value: str) -> bool:
        """Whether ``value`` is one of this adapter's own declared qualifiers.

        A declared qualifier is vocabulary, not entity text. Resolving it can
        only ever miss, so recording the miss would put an unmatched entry on
        every single claim and leave the map useless as a signal.
        """
        return slot == _QUALIFIER_SLOT and value.strip() in self.qualifiers

    def enrich_claim(self, claim: dict) -> dict:
        """Attach CURIEs to the claim slots whose text names a resolvable entity.

        Pure: returns a new dict and never mutates the input or its nested maps.
        A slot that carries text but resolves to nothing is recorded under
        ``ontology_terms_unmatched`` rather than passed over, so "not a known
        entity" stays distinguishable from "no attempt was made" -- except for a
        qualifier drawn from this adapter's own set, which is vocabulary rather
        than entity text and is skipped outright. A gene symbol additionally
        populates ``xrefs[slot]`` with every cross-reference the nomenclature
        row carries.
        """
        out = dict(claim)
        enrichment = _Enrichment.from_claim(out)
        for slot in _ENRICHABLE_SLOTS:
            if slot in enrichment.terms:
                continue
            value = out.get(slot)
            if not isinstance(value, str) or not value.strip():
                continue
            if self._is_own_vocabulary(slot, value):
                continue
            resolution = _resolve_text(value, self._gene_index)
            if resolution:
                enrichment.record(slot, resolution)
            else:
                enrichment.miss(slot, value)
        enrichment.apply_to(out)
        return out
