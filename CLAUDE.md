# CLAUDE.md — Dev Guide

## Read this before adding anything

This package sits in **tier 3** of the ASB domain architecture
(`AgenticScienceBuilder/docs/adding-a-domain.md`, implementing ADR 0001):
classification vocabulary and claim/ontology enrichment over *already extracted,
already validated* claim dicts. It never reads article text.

That boundary is the whole point. Ask of any addition: does it constrain the
extraction model, resolve an identifier, or validate output? If it would instead
mine free text, it belongs in the genomics JSON profile (tier 1) or in an
`asb-domain-core` enricher (tier 2), not here.

## What may go in a table here

Every term in a table must be a real, independently verifiable domain entity,
and the file must name the authority it came from:

| Table | Authority |
|---|---|
| `_identifiers.py` grammars | the body that issues the identifier (NCBI dbSNP, NCBI ClinVar, NCBI RefSeq, Ensembl, HGNC, OBO Foundry) |
| `_namespaces.py` expansions | OBO PURLs, or identifiers.org registry prefixes |
| `data/gene_index.tsv` | the HGNC complete set, in full |
| `data/edam_index.json` | EDAM, cut from a recorded set of roots |

**Do not** add a term because a particular paper needed it. A table grown from
the vocabulary of the papers that happened to fail is the documented overfit
failure (`docs/audits/2026-06-16_overfit-and-domain-audit.md` in the ASB repo):
it makes those papers pass while silently no-op'ing on every other one, and the
corpus then looks healthy while measuring nothing. If you cannot name the
vendor, standards body or database a term comes from, leave it out.

Two specific bans, both inherited from ADR 0001:

- **No repository-accession regex.** GEO / SRA / ENA / ArrayExpress /
  BioProject appear in `_namespaces.py` as *expansions* only. Detecting them in
  prose is `asb_domain_core.enrichers.accession_scanner`'s job, and it is the
  single source of truth for that.
- **No instrument, tool or organism keyword tables.** Those are extraction
  heuristics and live in `domain_profiles/genomics.json`.

## Project shape

One `DomainAdapter` plus two file projections in one wheel:

- `GenomicsAdapter` (`genomics.py`) — vocabulary and `enrich_claim`
- `_identifiers.py` — identifier grammars, full-match, one slot value at a time
- `_genes.py` — HGNC symbol index (not bundled; see below)
- `_edam.py` — vendored EDAM label index
- `_namespaces.py` — the CURIE prefix table
- `vcf/`, `gff3/` — path in, list of 5-slot claim dicts out

The adapter has `api_version = 1` matching `indicium_adapters.CURRENT_API_VERSION`
and is registered under the `indicium.adapters` entry-point group.

## Adding a qualifier

1. Add it to `GenomicsAdapter.qualifiers`. It must match `^[a-z][a-z0-9_]*$` —
   Perspicacité routes adapter qualifiers onto the open `asb:domainQualifier`
   slot whose SHACL shape enforces that pattern.
2. `extraction_context()` joins the set automatically, so there is nothing to
   update there.
3. If a projection should emit it, add the claim builder and a test.
   `test_protocol_conformance.py::test_the_projections_only_emit_declared_qualifiers`
   fails if a parser invents one.
4. If the claim shape it names is one of the four this package tracks (variant
   calls, differential expression, assembly metrics, annotation), add it to
   `test_the_claim_shapes_this_domain_needs_have_a_qualifier`.

## Adding an identifier grammar

1. Add an `_IdentifierSpec` naming the issuing authority.
2. Add its prefix to `_namespaces.GENOMICS_PREFIXES` — `test_namespaces.py`
   fails on a grammar that mints a CURIE the table cannot expand.
3. Add both halves of the pair to `test_identifiers.py`: a real identifier that
   resolves, **and** a near-miss that must not. The second is the one that stops
   the pattern being loosened until prose starts matching.

## Refreshing bundled data

- `data/build_gene_index.py` regenerates `gene_index.tsv` from the HGNC complete
  set. Stdlib only, maintainer-only, gitignored output:

      python -m indicium_adapters_genomics.data.build_gene_index

  The index is deliberately absent from the repo. Ship it whole or not at all —
  a partial symbol table is the overfit failure above. Absence is inspectable
  through `GeneIndex.is_available`, so it never reads as "checked, no such gene".

- `data/build_edam_index.py` regenerates `edam_index.json` from EDAM through
  EBI OLS4, and emits the `subtopics` slug list off the same walk. That list
  no longer feeds `GenomicsAdapter.subdomains`, which the 2026-08-25
  cross-domain vocabulary audit fixed by hand. Stdlib only, maintainer-only:

      python -m indicium_adapters_genomics.data.build_edam_index

  The roots in that file are the only vocabulary decision; everything else is
  read from EDAM. Do not hand-edit `edam_index.json` — the two entries that were
  hand-added before this script existed both mapped a label to a term EDAM never
  gave it.

Neither runs at install or import time.

## Test before commit

```bash
ruff check .
pytest
```

The `tests/test_entry_points.py` suite needs the package actually installed
(`pip install -e .`); it fails against a bare `PYTHONPATH=src` run because the
entry point is only registered by an install.

## Repo relationships

- The core Protocol lives in `~/git/indicium-adapters` (do not vendor it).
- `~/git/indicium-adapters-transcriptomics` owns the non-human nomenclature join
  (MGI, FlyBase, WormBase, ZFIN) and its own gene index. Do not duplicate it.
  Its prefix table is the reference spelling for every prefix the two share —
  except `edam`, lower-case here and in the proteomics sibling because that is
  what `_edam.py` mints, and `EDAM` there.
- Accession detection lives in `AgenticScienceBuilder/asb-domain-core`.
