"""Stage 2: sequence quality control.

Policy: ambiguous bases (IUPAC codes incl. N) are counted and kept as-is, never converted to A/C/G/T.
Records outside the length window, with invalid characters, or with too many ambiguous bases FAIL.
Records with any ambiguity, or whose sequence is identical to another record, are marked REVIEW/flagged.
"""
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

from common import ROOT, load_config

UNAMBIGUOUS = set("ACGT")
IUPAC_AMBIGUOUS = set("NRYSWKMBDHV")


def qc_record(seq, min_length, max_length, max_ambiguous_fraction, max_gap_fraction):
    seq = str(seq).upper()
    n = len(seq)
    invalid = sorted(set(seq) - UNAMBIGUOUS - IUPAC_AMBIGUOUS - {"-"})
    n_amb = sum(c in IUPAC_AMBIGUOUS for c in seq)
    n_gap = seq.count("-")
    flags, status = [], "pass"
    if invalid:
        flags.append("invalid_characters:" + "".join(invalid)); status = "fail"
    if n < min_length:
        flags.append("incomplete_too_short"); status = "fail"
    if n > max_length:
        flags.append("too_long"); status = "fail"
    if n and n_amb / n > max_ambiguous_fraction:
        flags.append("ambiguous_fraction_exceeds_threshold"); status = "fail"
    elif n_amb:
        flags.append("contains_ambiguous_bases")
        if status == "pass":
            status = "review"
    if n and n_gap / n > max_gap_fraction:
        flags.append("contains_gaps"); status = "fail"
    return {"length": n, "n_ambiguous": n_amb, "ambiguous_fraction": n_amb / n if n else 0.0,
            "n_gap": n_gap, "gap_fraction": n_gap / n if n else 0.0,
            "qc_status": status, "qc_flags": ";".join(flags)}


def run_qc(records, qc_cfg):
    rows, seen_ids, by_seq = [], defaultdict(int), defaultdict(list)
    for rec in records:
        seen_ids[rec.id] += 1
        by_seq[str(rec.seq).upper()].append(rec.id)
    for rec in records:
        r = qc_record(rec.seq, qc_cfg["min_length"], qc_cfg["max_length"],
                      qc_cfg["max_ambiguous_fraction"], qc_cfg["max_gap_fraction"])
        r["sample_id"] = rec.id
        flags = [f for f in r["qc_flags"].split(";") if f]
        if seen_ids[rec.id] > 1:
            flags.append("duplicate_accession"); r["qc_status"] = "fail"
        twins = [i for i in by_seq[str(rec.seq).upper()] if i != rec.id]
        r["identical_to"] = ";".join(twins)
        if twins:
            flags.append("identical_sequence_in_set")   # informational; does not change status
        r["qc_flags"] = ";".join(flags)
        rows.append(r)
    cols = ["sample_id", "length", "n_ambiguous", "ambiguous_fraction", "n_gap", "gap_fraction",
            "identical_to", "qc_status", "qc_flags"]
    return pd.DataFrame(rows)[cols]


def main():
    cfg = load_config()
    raw, proc, tables = (ROOT / cfg["paths"][k] for k in ("raw", "processed", "tables"))
    records = list(SeqIO.parse(raw / "samples.fasta", "fasta"))
    qc = run_qc(records, cfg["qc"])
    qc.to_csv(tables / "qc_summary.csv", index=False)
    keep = set(qc.loc[qc.qc_status.isin(["pass", "review"]), "sample_id"])
    SeqIO.write([r for r in records if r.id in keep], proc / "qc_passed.fasta", "fasta")
    print(qc.qc_status.value_counts().to_string())
    print(f"{len(keep)} of {len(records)} sequences carried forward")


if __name__ == "__main__":
    main()
