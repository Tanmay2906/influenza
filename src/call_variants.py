"""Stage 3b: call variants from a pairwise alignment.

Coordinates: 1-based positions on the reference. Indels use VCF-style anchored representation:
  deletion  -> position = anchor base (the base before the event), ref=anchor+deleted, alt=anchor
  insertion -> position = anchor base, ref=anchor, alt=anchor+inserted
Indels are left-normalised. Substitutions are reported one nucleotide at a time.

Missing data policy: terminal gaps / overhangs are *uncovered* (not deletions/insertions) and are
summarised as coverage; non-ACGT query bases are *ambiguous* positions (not substitutions).
"""
import pandas as pd

from common import ROOT, load_config

ACGT = set("ACGT")


def left_normalize(r, q):
    """Shift every gap run as far left as possible without changing the alignment score."""
    r, q = list(r), list(q)
    i, n = 0, len(r)
    while i < n:
        if (r[i] == "-") != (q[i] == "-"):
            gapped, other = (r, q) if r[i] == "-" else (q, r)
            j = i
            while j < n and gapped[j] == "-" and other[j] != "-":
                j += 1
            while (i > 0 and gapped[i - 1] != "-" and other[i - 1] != "-"
                   and gapped[i - 1] == other[i - 1] and gapped[i - 1] == other[j - 1]):
                gapped[j - 1], gapped[i - 1] = gapped[i - 1], "-"
                i -= 1
                j -= 1
            # run now starts at i; continue after it
            while j < n and gapped[j] == "-" and other[j] != "-":
                j += 1
            i = j
        else:
            i += 1
    return "".join(r), "".join(q)


def _nongap_span(s):
    idx = [k for k, c in enumerate(s) if c != "-"]
    return (idx[0], idx[-1]) if idx else (None, None)


def call_variants(ref_aln, qry_aln):
    """Return dict(variants, cov_start, cov_end, ambiguous_positions)."""
    ref_aln, qry_aln = ref_aln.upper(), qry_aln.upper()
    rs, re_ = _nongap_span(ref_aln)
    qs, qe = _nongap_span(qry_aln)
    if qs is None:
        return {"variants": [], "cov_start": None, "cov_end": None, "ambiguous_positions": []}
    lo, hi = max(rs, qs), min(re_, qe)               # columns where both sequences are present
    ref_ungapped = ref_aln.replace("-", "")
    pos = len(ref_aln[:lo].replace("-", ""))         # reference bases before column lo
    cov_start = pos + 1
    r, q = left_normalize(ref_aln[lo:hi + 1], qry_aln[lo:hi + 1])
    variants, ambiguous = [], []
    k, n = 0, len(r)
    while k < n:
        rc, qc = r[k], q[k]
        if rc != "-" and qc != "-":
            pos += 1
            if qc not in ACGT:
                ambiguous.append(pos)
            elif rc in ACGT and rc != qc:
                variants.append({"position": pos, "ref_allele": rc, "alt_allele": qc,
                                 "variant_type": "substitution", "aln_column": lo + k + 1, "note": ""})
            k += 1
        elif qc == "-":                              # deletion relative to reference
            j = k
            while j < n and q[j] == "-" and r[j] != "-":
                j += 1
            deleted = r[k:j]
            anchor = ref_ungapped[pos - 1]
            variants.append({"position": pos, "ref_allele": anchor + deleted, "alt_allele": anchor,
                             "variant_type": "deletion", "aln_column": lo + k + 1, "note": ""})
            pos += len(deleted)
            k = j
        else:                                        # insertion relative to reference
            j = k
            while j < n and r[j] == "-" and q[j] != "-":
                j += 1
            inserted = q[k:j]
            anchor = ref_ungapped[pos - 1]
            note = "contains_ambiguous_base" if set(inserted) - ACGT else ""
            variants.append({"position": pos, "ref_allele": anchor, "alt_allele": anchor + inserted,
                             "variant_type": "insertion", "aln_column": lo + k + 1, "note": note})
            k = j
    return {"variants": variants, "cov_start": cov_start, "cov_end": pos, "ambiguous_positions": ambiguous}


def main():
    cfg = load_config()
    proc, tables = ROOT / cfg["paths"]["processed"], ROOT / cfg["paths"]["tables"]
    aln = pd.read_csv(proc / "alignments.tsv", sep="\t")
    ref_len = len(aln.ref_aligned.iloc[0].replace("-", ""))
    var_rows, cov_rows = [], []
    excluded = aln[aln.alignment_status != "ok"]
    excluded[["sample_id", "identity", "alignment_status"]].to_csv(tables / "excluded_samples.csv", index=False)
    for row in aln[aln.alignment_status == "ok"].itertuples():
        res = call_variants(row.ref_aligned, row.query_aligned)
        for v in res["variants"]:
            var_rows.append({"sample_id": row.sample_id, "segment": cfg["study"]["segment"],
                             "reference_accession": row.reference_accession, **v})
        cov_rows.append({"sample_id": row.sample_id, "reference_length": ref_len,
                         "cov_start": res["cov_start"], "cov_end": res["cov_end"],
                         "covered_fraction": (res["cov_end"] - res["cov_start"] + 1) / ref_len,
                         "n_variants": len(res["variants"]),
                         "n_ambiguous_positions": len(res["ambiguous_positions"]),
                         "ambiguous_positions": ";".join(map(str, res["ambiguous_positions"]))})
    cols = ["sample_id", "segment", "reference_accession", "position", "ref_allele", "alt_allele",
            "variant_type", "aln_column", "note"]
    pd.DataFrame(var_rows, columns=cols).to_csv(tables / "mutation_calls.csv", index=False)
    pd.DataFrame(cov_rows).to_csv(tables / "coverage.csv", index=False)
    print(f"{len(var_rows)} variants across {len(cov_rows)} samples ({len(excluded)} excluded)")


if __name__ == "__main__":
    main()
