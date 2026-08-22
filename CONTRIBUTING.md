# Contributing to indicium-adapters-genomics

Thanks for your interest!

## Setup

```bash
git clone https://github.com/HolobiomicsLab/indicium-adapters-genomics
cd indicium-adapters-genomics
pip install -e .[dev]
```

`indicium-adapters` is not on PyPI yet; install it from git first (see
`.github/workflows/ci.yml` for the exact commands CI uses).

## Development workflow

1. Create a feature branch from `main`.
2. Write tests first (TDD). The suite is `pytest`.
3. Run `ruff check .` and `pytest` before pushing.
4. Open a PR describing the change.

## Before you add a term

Read [CLAUDE.md § What may go in a table here](CLAUDE.md#what-may-go-in-a-table-here)
first. The short version: every term must be traceable to a named vendor,
standards body or database. If you cannot name the authority, leave it out.

## Adding a qualifier

See [CLAUDE.md § Adding a qualifier](CLAUDE.md#adding-a-qualifier).

## Refreshing bundled data

See [CLAUDE.md § Refreshing bundled data](CLAUDE.md#refreshing-bundled-data).
`data/build_gene_index.py` downloads from HGNC and is maintainer-only.

## Code of conduct

By participating you agree to abide by the
[Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
