"""Maintainer script: regenerate ``gene_index.tsv`` from the HGNC complete set.

The HUGO Gene Nomenclature Committee's complete set is the authority for human
gene symbols and carries the Ensembl and NCBI Gene cross-references in the same
row, so one download produces the whole table. Output schema::

    symbol  hgnc  ensembl  ncbigene  name

Only rows whose ``status`` is ``Approved`` are written -- withdrawn and entry-
level placeholder records would resolve a symbol HGNC no longer recognises.

Usage (never run at install or import time)::

    python -m indicium_adapters_genomics.data.build_gene_index
    python -m indicium_adapters_genomics.data.build_gene_index --hgnc-tsv local.txt

Re-run on an HGNC release. Non-human nomenclature (MGI, FlyBase, WormBase,
ZFIN) is deliberately out of scope here: the sibling
``indicium-adapters-transcriptomics`` already joins those authorities, and a
second copy of a nomenclature join is the duplicate-then-drift failure ADR 0001
exists to prevent.
"""

from __future__ import annotations

import argparse
import csv
import sys
import urllib.request
from pathlib import Path
from typing import Final

HGNC_URLS: Final[tuple[str, ...]] = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
    "hgnc_complete_set.txt",
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
)

OUT_COLUMNS: Final[tuple[str, ...]] = ("symbol", "hgnc", "ensembl", "ncbigene", "name")

APPROVED_STATUS: Final[str] = "Approved"

DOWNLOAD_TIMEOUT_SECONDS: Final[int] = 300


def _download_hgnc(target: Path) -> None:
    """Fetch the HGNC complete set into ``target``, trying each mirror in turn."""
    target.parent.mkdir(parents=True, exist_ok=True)
    for url in HGNC_URLS:
        print(f"[hgnc] downloading {url}", file=sys.stderr)
        try:
            with urllib.request.urlopen(url, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response:
                target.write_bytes(response.read())
            return
        except OSError as exc:
            print(f"[hgnc] {url} failed: {exc}", file=sys.stderr)
    raise RuntimeError("no HGNC mirror reachable; pass --hgnc-tsv with a local copy")


def _to_row(record: dict[str, str]) -> dict[str, str] | None:
    """Project one HGNC record onto the output schema, or None if unusable."""
    if record.get("status", "").strip() != APPROVED_STATUS:
        return None
    symbol = record.get("symbol", "").strip()
    if not symbol:
        return None
    return {
        "symbol": symbol,
        "hgnc": record.get("hgnc_id", "").strip(),
        "ensembl": record.get("ensembl_gene_id", "").strip(),
        "ncbigene": record.get("entrez_id", "").strip(),
        "name": record.get("name", "").strip(),
    }


def _write_index(source: Path, target: Path) -> int:
    """Read the HGNC TSV at ``source``, write the index to ``target``, return rows."""
    written = 0
    with source.open("r", encoding="utf-8", newline="") as reader_handle:
        with target.open("w", encoding="utf-8", newline="") as writer_handle:
            writer = csv.DictWriter(writer_handle, OUT_COLUMNS, delimiter="\t")
            writer.writeheader()
            for record in csv.DictReader(reader_handle, delimiter="\t"):
                row = _to_row(record)
                if row:
                    writer.writerow(row)
                    written += 1
    return written


def main() -> None:
    """Parse arguments, download the HGNC set if needed, and write the index."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hgnc-tsv", type=Path, help="local hgnc_complete_set.txt")
    args = parser.parse_args()
    here = Path(__file__).parent
    source = args.hgnc_tsv or here / "_hgnc_complete_set.txt"
    if not source.is_file():
        _download_hgnc(source)
    target = here / "gene_index.tsv"
    print(f"[hgnc] wrote {target} with {_write_index(source, target):,} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
