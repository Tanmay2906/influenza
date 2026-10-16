"""Stage 3a: pairwise alignment of each sample against the reference.

Semi-global alignment: end gaps are free so that truncated sequences are not penalised, and
terminal gaps are later treated as *missing coverage* rather than deletions. Internal gaps use affine
penalties from config.yaml. Alignments are saved as gapped strings for traceability.
"""
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner

from common import ROOT, load_config


def make_aligner(cfg):
    a = PairwiseAligner()
    a.mode = cfg["mode"]
    a.match_score = cfg["match"]
    a.mismatch_score = cfg["mismatch"]
    a.open_gap_score = cfg["open_gap"]
    a.extend_gap_score = cfg["extend_gap"]
    a.end_gap_score = 0.0           # free end gaps in both sequences
    return a


def align_pair(ref, qry, aligner):
    """Return (ref_aligned, qry_aligned, score) as gapped upper-case strings.
    Ties between equally good alignments are broken by Biopython's first alignment; indel
    placement is normalised afterwards in call_variants.left_normalize."""
    ref, qry = str(ref).upper(), str(qry).upper()
    aln = aligner.align(ref, qry)[0]
    return aln[0], aln[1], float(aln.score)


def identity(ref_aln, qry_aln):
    """Fraction of identical bases over columns where both sequences have a base."""
    both = [(a, b) for a, b in zip(ref_aln, qry_aln) if a != "-" and b != "-"]
    return sum(a == b for a, b in both) / len(both) if both else 0.0


def main():
    cfg = load_config()
    proc = ROOT / cfg["paths"]["processed"]
    ref = SeqIO.read(ROOT / cfg["paths"]["reference"] / "reference.fasta", "fasta")
    aligner = make_aligner(cfg["alignment"])
    rows = []
    for rec in SeqIO.parse(proc / "qc_passed.fasta", "fasta"):
        r, q, score = align_pair(ref.seq, rec.seq, aligner)
        ident = identity(r, q)
        rows.append({"sample_id": rec.id, "reference_accession": ref.id, "score": score,
                     "identity": ident,
                     "alignment_status": "ok" if ident >= cfg["alignment"]["min_identity"] else "lineage_mismatch",
                     "ref_aligned": r, "query_aligned": q})
    pd.DataFrame(rows).to_csv(proc / "alignments.tsv", sep="\t", index=False)
    bad = [r["sample_id"] for r in rows if r["alignment_status"] != "ok"]
    print(f"aligned {len(rows)} sequences to {ref.id}; {len(bad)} excluded as lineage_mismatch: {bad}")


if __name__ == "__main__":
    main()
