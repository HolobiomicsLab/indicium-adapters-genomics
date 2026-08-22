"""VCF parser submodule -- project variant records into asb:Claim dicts.

VCF is specified by the hts-specs group (samtools.github.io/hts-specs/VCFv4.3.pdf).
Three claim types are produced per data line, mirroring how the sibling
metabolomics package projects mzTab-M sections:

- ``located_on``   -- one per record: the variant and the CHROM it sits on.
- ``identifies``   -- one per entry of the ID column that a grammar recognises
                      (dbSNP ``rs``, ClinVar ``VCV``/``RCV``/``SCV``). The column
                      is a semicolon-separated list, so a record can yield more
                      than one. The CURIE is the file's own curated identifier,
                      so it is justified as ``semapv:CompositeMatching``.
- ``genotyped_as`` -- one per sample column carrying a called ``GT`` value.

Plain and gzip-compressed files are both accepted.
"""

from __future__ import annotations

from indicium_adapters_genomics.vcf.parser import parse_vcf_to_claims

__all__ = ["parse_vcf_to_claims"]
