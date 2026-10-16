"""Codon / protein-level annotation tests using a tiny synthetic single-CDS reference."""
import io

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from annotate_variants import Cds, annotate


def make_reference(seq):
    rec = SeqRecord(Seq(seq), id="SYN_REF.1")
    rec.features.append(SeqFeature(FeatureLocation(0, len(seq)), type="CDS",
                                    qualifiers={"gene": ["TEST"], "codon_start": ["1"],
                                                "transl_table": ["1"],
                                                "translation": [str(Seq(seq).translate(to_stop=True))]}))
    return rec


REF_SEQ = "ATGCCAGAATTT"   # M P E F


def test_coding_substitution_proline_to_leucine():
    rec = make_reference(REF_SEQ)
    cds = [Cds(rec.features[0], rec)]
    variants = [{"sample_id": "S1", "segment": 4, "reference_accession": "SYN_REF.1",
                 "position": 5, "ref_allele": "C", "alt_allele": "T",
                 "variant_type": "substitution", "aln_column": 5, "note": ""}]
    out = annotate(variants, rec, cds, [])
    row = out.iloc[0]
    assert row.codon_change == "CCA>CTA"
    assert row.protein_change == "P2L"
    assert row.consequence == "missense"
    assert row.annotation_status == "annotated"


def test_synonymous_substitution():
    rec = make_reference(REF_SEQ)
    cds = [Cds(rec.features[0], rec)]
    # GAA (E) -> GAG (E) at the third codon position: synonymous
    variants = [{"sample_id": "S1", "segment": 4, "reference_accession": "SYN_REF.1",
                 "position": 9, "ref_allele": "A", "alt_allele": "G",
                 "variant_type": "substitution", "aln_column": 9, "note": ""}]
    out = annotate(variants, rec, cds, [])
    row = out.iloc[0]
    assert row.consequence == "synonymous"
    assert row.protein_change.endswith("=")


def test_stop_gained():
    rec = make_reference(REF_SEQ)
    cds = [Cds(rec.features[0], rec)]
    # GAA (E, codon 3) -> TAA (stop), single substitution at position 7
    variants = [{"sample_id": "S1", "segment": 4, "reference_accession": "SYN_REF.1",
                 "position": 7, "ref_allele": "G", "alt_allele": "T",
                 "variant_type": "substitution", "aln_column": 7, "note": ""}]
    out = annotate(variants, rec, cds, [])
    row = out.iloc[0]
    assert row.consequence == "stop_gained"
    assert row.protein_change.endswith("*")


def test_ambiguous_base_in_codon_withholds_protein_call():
    rec = make_reference(REF_SEQ)
    cds = [Cds(rec.features[0], rec)]
    variants = [{"sample_id": "S1", "segment": 4, "reference_accession": "SYN_REF.1",
                 "position": 5, "ref_allele": "C", "alt_allele": "T",
                 "variant_type": "substitution", "aln_column": 5, "note": ""}]
    coverage = {"S1": {"cov_start": 1, "cov_end": 12, "ambiguous_positions": [6]}}  # same codon, pos 6
    out = annotate(variants, rec, cds, [], coverage)
    row = out.iloc[0]
    assert row.annotation_status == "requires_review"
    assert row.protein_change == "" and row.consequence == ""
    assert "ambiguous_base_in_codon" in row.note
