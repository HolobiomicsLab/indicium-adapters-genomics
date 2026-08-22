"""GFF3 parser submodule -- project annotation features into asb:Claim dicts.

GFF3 is specified by the Sequence Ontology
(github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).
Three claim types are produced, mirroring how the sibling metabolomics package
projects mzTab-M sections:

- ``located_on``   -- one per feature: the feature and the sequence it sits on.
                      When the ``type`` column is written as an SO accession the
                      CURIE is carried on the subject slot.
- ``annotated_by`` -- one per feature whose ``source`` column names the
                      algorithm or procedure that produced it.
- ``has_feature``  -- one per ``Parent`` value, linking a parent feature to its
                      child (gene to mRNA, mRNA to exon).

The ``##FASTA`` section, if present, ends the annotation body and is not read.
"""

from __future__ import annotations

from indicium_adapters_genomics.gff3.parser import parse_gff3_to_claims

__all__ = ["parse_gff3_to_claims"]
