# Changelog

All notable changes to this project will be documented in this file.

## [0.3.0] — 2026-08-22

An adversarial review of 0.2.0. Every item below is a defect that was found and
fixed, not a new feature; the two packaging ones made the package uninstallable
and the ontology ones minted CURIEs that resolve to nothing.

### Fixed (packaging)
- **`uv pip install -e .` — the command the README documents — could not resolve
  `indicium-adapters` at all**, because neither it nor `indicium-standard` is on
  PyPI and the repo declared no `[tool.uv.sources]`. Added the table the
  metabolomics sibling already carries.
- **A built wheel silently shipped no gene index.** `gene_index.tsv` is
  gitignored and hatchling excludes VCS-ignored files by default, so a
  maintainer could regenerate the index and still deploy a wheel that resolved
  no symbol. Declared it under `[tool.hatch.build.targets.wheel] artifacts`.
- **CI never installed `indicium`**, which supplies the `indicium-standard`
  distribution `indicium-adapters` requires; the resolver would have fallen
  through to PyPI. Added the install step, matching the sibling workflow.

### Fixed (ontology)
- **`ClinVar` expanded three different accession families to one IRI that
  resolves for none of them.** identifiers.org registers `clinvar` (numeric
  variation id, `^\d+$`), `clinvar.record` (`RCV…`) and `clinvar.submission`
  (`SCV…`) separately, and none of them accepts a `VCV` accession. Split into
  `ClinVar` (VCV, expanded through NCBI's own variation resolver),
  `ClinVar.record` and `ClinVar.submission`. **Breaking** for anything reading
  `ontology_prefixes["ClinVar"]` or matching on the old CURIE prefix.
- **The `NZ_` RefSeq alternative could not match a real accession.** It required
  an all-digit local part, but NZ accessions carry a WGS project code first
  (`NZ_CP012345`). Given its own grammar; `NZ_000546` no longer resolves.
- **Two `edam_index.json` entries were not EDAM.** `metatranscriptomics` mapped
  to `topic_3308` (Transcriptomics) although EDAM gives Metatranscriptomics its
  own `topic_3941`, and `environmental sequencing` was a label EDAM carries
  nowhere. Both were the wrong-CURIE failure the package documents as its own
  worst outcome. The index also omitted `topic_0199` (Genetic variation), a
  direct subtopic of Genomics.

### Added
- **`data/build_edam_index.py`** — the EDAM index generator. `_edam.py` and
  `genomics.py` both cited `scripts/build_index.py`, which does not exist in
  this repo, so the index was unmaintainable and its provenance unverifiable.
  The new script regenerates `edam_index.json` from EDAM via EBI OLS4 (stdlib
  only, never at install or import time) and reproduces the shipped file.

### Changed
- **`subdomains` is now read from the generated index instead of a parallel
  hand-written tuple**, which is what let the vocabulary claim "EDAM's own
  sub-topics of Genomics" while containing a term EDAM does not file there. It
  gains `genetic-variation`. A test binds every subdomain to a concept the
  index can resolve.
- **VCF: every entry of the `ID` column now yields an `identifies` claim.** The
  column is a semicolon-separated list of identifiers for the same variant, and
  all but the first were silently dropped.
- **GFF3: attribute values are decoded after the multi-valued split, not
  before.** An escaped `%2C` inside a `Parent` id was being turned into a
  separator, splitting one parent into two.
- **GFF3: the `type` column resolves against Sequence Ontology only.** It was
  using the general identifier resolver, which would have carried a CURIE from
  any of eleven other authorities onto the feature as its SO type.
- The `_namespaces.py` claim that prefix spellings follow the transcriptomics
  table "so a CURIE minted by either adapter expands to the same IRI" is
  narrowed to the truth: `edam` is lower-case here and in the proteomics
  sibling, `EDAM` there, and a consumer merging both must carry both.
- **A qualifier drawn from the adapter's own declared set is no longer reported
  as unresolved text.** `enrich_claim` walks the qualifier slot, so every claim
  came back carrying `ontology_terms_unmatched: {"qualifier": …}` for a term
  that was never entity text — a floor of one miss per claim, which made the
  map useless as the signal it exists to be. A qualifier outside the declared
  set is still resolved.
- `GenomicsAdapter.__init__` no longer describes `gene_index.tsv` as bundled; it
  is deliberately not shipped, as the rest of the package says.

## [0.2.0] — 2026-08-22

Turns the EDAM-only adapter into one that covers the claim shapes genomics
papers actually make, and adds the two file projections.

### Added
- **`_identifiers.py`** — full-match recognition of gene, transcript, variant
  and ontology identifiers carried in a claim slot: dbSNP `rs`, ClinVar
  `VCV`/`RCV`/`SCV`, Ensembl stable ids, RefSeq accessions, HGNC ids and OBO
  CURIEs. Every grammar records the authority that defines it, retrievable via
  `_identifiers.authorities()`. Nothing here scans article text: repository
  accessions mined from prose remain `asb-domain-core`'s concern (ADR 0001).
- **`_genes.py`** — gene-symbol resolution to HGNC, with Ensembl and NCBI Gene
  as cross-references on `claim["xrefs"][slot]`. Backed by an HGNC-derived index
  that is deliberately not committed; `GeneIndex.is_available` distinguishes "no
  nomenclature installed" from "no such gene".
- **`data/build_gene_index.py`** — maintainer script regenerating that index
  from the HGNC complete set (stdlib only, never run at install or import time).
- **`_namespaces.py`** — the prefix table, now covering the sequence and
  expression repositories (ENA, SRA, GEO, ArrayExpress, BioProject) and the
  identifier authorities (HGNC, Ensembl, RefSeq, NCBI Gene, dbSNP, ClinVar)
  alongside the ontologies. Spellings follow the sibling
  `indicium-adapters-transcriptomics` table so a CURIE minted by either adapter
  expands to the same IRI.
- **`vcf.parse_vcf_to_claims`** — projects VCF records to `located_on`,
  `identifies` and `genotyped_as` claims. Reads plain or gzip-compressed input,
  reads `GT` by its position in `FORMAT`, and skips uncalled genotypes.
- **`gff3.parse_gff3_to_claims`** — projects GFF3 features to `located_on`,
  `annotated_by` and `has_feature` claims, percent-decoding attributes and
  stopping at `##FASTA`.
- Five qualifiers: `identifies`, `has_feature`, `has_consequence`,
  `has_assembly_metric`, `differentially_expressed_in`.

### Changed (breaking)
- **`ontology_prefixes` keys `ensembl` and `refseq` are now `ENSEMBL` and
  `RefSeq`**, and their expansions moved to identifiers.org compact-URI form, to
  match the transcriptomics adapter. `edam` stays lower-case: it is the prefix
  `_edam.py` mints, and `test_namespaces.py` now binds the two together.
- `enrich_claim` copies the claim's nested maps instead of only its outer dict,
  so enrichment can no longer write back into a caller's `ontology_terms`.

## [0.1.0] — 2026-08-21

Initial release: `GenomicsAdapter` with EDAM-derived subdomains, curated
techniques and qualifiers, and exact-match EDAM concept enrichment.
