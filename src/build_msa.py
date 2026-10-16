"""Stage 6a: build a multiple-sequence alignment of reference + QC-passed, lineage-consistent samples.

We do not install an external MSA tool (MAFFT/MUSCLE); instead we build a progressive alignment with
Biopython by aligning every sequence to the reference (same engine/parameters as align_sequences.py)
and projecting all sequences onto the reference's coordinate system (a reference-anchored MSA). This
is a reasonable approximation for a small set of closely related sequences of one gene; it will not
place insertions relative to the reference as well as a true progressive aligner would, so insertions
are dropped from the MSA (every row keeps reference-length columns) and this limitation is documented.
"""
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from common import ROOT, load_config


def to_reference_columns(ref_aln, qry_aln):
    """Return the query's bases at every reference (non-gap) column; drop insertions."""
    return "".join(q for r, q in zip(ref_aln, qry_aln) if r != "-")


def main():
    cfg = load_config()
    proc = ROOT / cfg["paths"]["processed"]
    ref = SeqIO.read(ROOT / cfg["paths"]["reference"] / "reference.fasta", "fasta")
    aln = pd.read_csv(proc / "alignments.tsv", sep="\t")
    aln = aln[aln.alignment_status == "ok"] if "alignment_status" in aln.columns else aln

    records = [SeqRecord(ref.seq, id=ref.id, description="reference")]
    for row in aln.itertuples():
        proj = to_reference_columns(row.ref_aligned, row.query_aligned)
        records.append(SeqRecord(Seq(proj), id=row.sample_id, description=""))

    SeqIO.write(records, proc / "msa_reference_anchored.fasta", "fasta")
    print(f"wrote reference-anchored MSA: {len(records)} sequences x {len(ref.seq)} columns")


if __name__ == "__main__":
    main()
