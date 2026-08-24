# indicium-adapters-genomics

A `genomics` `DomainAdapter` for the [indicium](https://github.com/HolobiomicsLab/indicium)
ecosystem, plus VCF and GFF3 projections into ASB claim dicts.

It owns **tier-3 concerns only** (ADR 0001): the classification vocabulary and
claim/ontology enrichment. It does **not** scan article text for repository
accessions — that is the canonical `accession_scanner` in `asb-domain-core`, and
a second copy of those patterns here would be the duplicate-then-drift failure
ADR 0001 exists to prevent. What this package contributes for identifiers is the
prefix table that *expands* those accessions, plus full-match recognition of the
gene / transcript / variant identifiers that sit inside a claim slot.

## What it contributes

| Attribute | Source |
|---|---|
| `primary_domain` | `genomics` — without this, a genomics paper classifies as `other` in AgenticScienceBuilder |
| `subdomains` | curated — the set the 2026-08-25 cross-domain vocabulary audit fixed, no longer read off the EDAM walk |
| `techniques` | curated — same audit; EDAM models assays as operations and data, not as a technique hierarchy |
| `qualifiers` | curated domain claim qualifiers, all snake_case, covering variant calls, differential expression, assembly metrics and annotation |
| `ontology_prefixes` | prefix → base IRI (`_namespaces.py`): SO, GO, GENO, OBI, EFO, NCBITaxon, UBERON, EDAM, HGNC, Ensembl, RefSeq, NCBI Gene, dbSNP, ClinVar (VCV / `.record` / `.submission`), ENA, SRA, GEO, ArrayExpress, BioProject |
| `enrich_claim` | identifier grammars → gene nomenclature → EDAM labels, each exact-match |

## Enrichment is exact-match, in three passes

`enrich_claim` walks the `subject`, `object` and `qualifier` slots and tries
three resolvers in decreasing order of certainty:

1. **Identifier grammar** (`_identifiers.py`). `rs334`, `VCV000012345`,
   `ENSG00000139618`, `NM_000546.6`, `NZ_CP012345`, `HGNC:4827`, `SO:0000704`.
   A string that
   *is* an authority's identifier is a syntactic fact, not a guess, so it earns
   `semapv:CompositeMatching`. Every pattern full-matches one slot value; none
   of them scans prose.
2. **Gene nomenclature** (`_genes.py`). A symbol resolves to its HGNC id, with
   Ensembl and NCBI Gene carried in `claim["xrefs"][slot]`.
3. **EDAM labels** (`_edam.py`). The vendored concept index.

Substring and fuzzy matching over an ontology manufacture confident-looking
wrong CURIEs, and a wrong identifier is worse than no identifier. A slot that
carries text but resolves to nothing is recorded under
`ontology_terms_unmatched` rather than passed over, so "this text is not a known
entity" stays distinguishable from "no attempt was made".

## File projections

Both mirror the way the sibling metabolomics package projects mzTab-M rows: a
path in, a list of 5-slot claim dicts out, no network, no mutation.

```python
from indicium_adapters_genomics import parse_vcf_to_claims, parse_gff3_to_claims

parse_vcf_to_claims("cohort.vcf.gz")   # located_on / identifies / genotyped_as
parse_gff3_to_claims("genome.gff3")    # located_on / annotated_by / has_feature
```

`identifies` is emitted once per entry of the VCF `ID` column that a grammar
recognises — the column is a semicolon-separated list, so one record can yield
more than one. `located_on` carries an SO CURIE only when the GFF3 `type` column
is written as an SO accession rather than a term name, and only when that
accession is SO's: the spec admits no other authority in that column.

## The gene index is not bundled

`data/gene_index.tsv` is regenerated from the HGNC complete set by
`data/build_gene_index.py` and is **not committed**. A hand-picked subset of
symbols would resolve the genes someone happened to look at and silently return
nothing for every other gene — the failure mode the ADR 0001 anti-patterns
describe. So the index is either the authority's full set or absent, and absence
is inspectable:

```python
adapter = GenomicsAdapter()
adapter.gene_index.is_available   # False => no nomenclature installed
```

Deployments that already hold nomenclature data can bypass the file:

```python
GenomicsAdapter(GeneIndex.from_tsv("/data/hgnc_index.tsv"))
```

## Install

```bash
uv pip install -e .
```

The `indicium.adapters` entry point makes it discoverable via
`indicium_adapters.discover_adapters()` on install — no configuration.

## Regenerating the bundled indexes

- `data/edam_index.json` — `python -m
  indicium_adapters_genomics.data.build_edam_index`. Reads EDAM through EBI
  OLS4, records the roots it was cut from, and emits the `subtopics` list off
  the same walk. Regenerate on an EDAM release.
- `data/gene_index.tsv` — see above; `python -m
  indicium_adapters_genomics.data.build_gene_index`.

Neither runs at install or import time.

## Licence

MIT.
