"""Synthetic test cases with known expected results (project plan, section 7)."""
import pytest
from Bio.Align import PairwiseAligner

from align_sequences import align_pair, identity
from call_variants import call_variants

ALN_CFG = {"mode": "global", "match": 2, "mismatch": -3, "open_gap": -8, "extend_gap": -1}


def aligner():
    a = PairwiseAligner()
    a.mode = ALN_CFG["mode"]
    a.match_score = ALN_CFG["match"]
    a.mismatch_score = ALN_CFG["mismatch"]
    a.open_gap_score = ALN_CFG["open_gap"]
    a.extend_gap_score = ALN_CFG["extend_gap"]
    a.end_gap_score = 0.0
    return a


def align_and_call(ref, qry):
    r, q, _ = align_pair(ref, qry, aligner())
    return call_variants(r, q)


REF = "ATGCCAGAATTT"   # ATG CCA GAA TTT


def test_identical_reference_and_query():
    res = align_and_call(REF, REF)
    assert res["variants"] == []


def test_reference_compared_with_itself_has_zero_differences():
    res = align_and_call(REF, REF)
    assert len(res["variants"]) == 0
    assert res["ambiguous_positions"] == []


def test_one_known_substitution():
    qry = "ATGCTAGAATTT"   # position 5 (1-based): C -> T
    res = align_and_call(REF, qry)
    assert res["variants"] == [
        {"position": 5, "ref_allele": "C", "alt_allele": "T", "variant_type": "substitution",
         "aln_column": 5, "note": ""}
    ]


def test_one_deletion():
    # REF has "CC" at positions 4-5; deleting either gives the same query, so a
    # left-normalized (VCF-style) representation must anchor at position 3 (G).
    qry = "ATGCAGAATTT"    # one of the two C's deleted
    res = align_and_call(REF, qry)
    dels = [v for v in res["variants"] if v["variant_type"] == "deletion"]
    assert len(dels) == 1
    d = dels[0]
    assert d["position"] == 3
    assert d["ref_allele"] == "GC" and d["alt_allele"] == "G"   # anchor + deleted base


def test_one_insertion():
    qry = "ATGCCCAGAATTT"  # one extra C inserted after position 5
    res = align_and_call(REF, qry)
    ins = [v for v in res["variants"] if v["variant_type"] == "insertion"]
    assert len(ins) == 1
    assert ins[0]["alt_allele"].startswith(ins[0]["ref_allele"])
    assert len(ins[0]["alt_allele"]) - len(ins[0]["ref_allele"]) == 1


def test_ambiguous_nucleotide_is_flagged_not_substituted():
    qry = "ATGCNAGAATTT"   # N at position 5
    res = align_and_call(REF, qry)
    assert res["ambiguous_positions"] == [5]
    # an ambiguous base must never be reported as a substitution
    assert all(v["position"] != 5 for v in res["variants"] if v["variant_type"] == "substitution")


def test_incomplete_sequence_is_not_a_deletion():
    qry = "ATGCCAGA"       # truncated at the 3' end: missing coverage, not a deletion
    res = align_and_call(REF, qry)
    assert res["cov_end"] < len(REF)
    assert not any(v["variant_type"] == "deletion" for v in res["variants"])


def test_substitution_near_alignment_gap():
    # a deletion immediately followed by a substitution
    qry = "ATGCAGCATTT"    # delete C at pos5, then G->C at (shifted) position 9
    res = align_and_call(REF, qry)
    subs = [v for v in res["variants"] if v["variant_type"] == "substitution"]
    dels = [v for v in res["variants"] if v["variant_type"] == "deletion"]
    assert len(dels) == 1
    assert len(subs) == 1


def test_identity_is_one_for_identical_sequences():
    r, q, _ = align_pair(REF, REF, aligner())
    assert identity(r, q) == 1.0
