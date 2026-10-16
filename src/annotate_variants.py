"""Stage 4: gene / codon / protein annotation of called variants against the reference GenBank record.

Uses the record's real CDS coordinates and strand (compound locations supported). One output row per
(variant, overlapping CDS). Substitutions in the same codon of the same sample are combined when the
alternate codon is derived. Protein consequences are withheld (status requires_review) when the codon
holds an ambiguous/uncovered base or lies downstream of a frameshift in the same sample.
"""
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

from common import ROOT, load_config

COMP = str.maketrans("ACGT", "TGCA")


class Cds:
    def __init__(self, feat, rec):
        q = feat.qualifiers
        self.gene = (q.get("gene") or q.get("product") or ["CDS"])[0]
        self.strand = feat.location.strand or 1
        start = int(q.get("codon_start", ["1"])[0]) - 1
        self.positions = list(feat.location)[start:]            # 0-based ref indices, transcript order
        self.seq = "".join(str(rec.seq[p]).upper() for p in self.positions)
        if self.strand == -1:
            self.seq = self.seq.translate(COMP)
        self.table = CodonTable.unambiguous_dna_by_id[int(q.get("transl_table", ["1"])[0])]
        self.index = {p: i for i, p in enumerate(self.positions)}
        if "translation" in q:
            made = str(Seq(self.seq).translate(table=self.table.id, to_stop=True))
            if made != q["translation"][0]:
                raise ValueError(f"CDS {self.gene}: translation does not match annotation")

    def aa(self, codon):
        return "*" if codon in self.table.stop_codons else self.table.forward_table.get(codon, "X")


def load_reference(path):
    rec = SeqIO.read(path, "genbank")
    cds = [Cds(f, rec) for f in rec.features if f.type == "CDS"]
    domains = []
    for f in rec.features:
        if f.type in ("mat_peptide", "sig_peptide"):
            name = (f.qualifiers.get("product") or ["signal_peptide"])[0]
            domains.append((int(f.location.start), int(f.location.end), name))
    return rec, cds, domains


def _domain(p0, domains):
    names = [n for s, e, n in domains if s <= p0 < e]
    return ";".join(names)


def annotate(variants, ref_rec, cds_list, domains, coverage=None):
    """variants: list of dicts (needs sample_id, position, ref_allele, alt_allele, variant_type).
    coverage: dict sample_id -> {cov_start, cov_end, ambiguous_positions(list)}; optional."""
    coverage = coverage or {}
    out = []
    by_sample = {}
    for v in variants:
        by_sample.setdefault(v["sample_id"], []).append(v)
    for sid, vs in by_sample.items():
        cov = coverage.get(sid, {})
        amb0 = {p - 1 for p in cov.get("ambiguous_positions", [])}
        cs, ce = cov.get("cov_start"), cov.get("cov_end")
        # frameshift start per CDS (0-based transcript index)
        fs_from = {}
        for v in vs:
            if v["variant_type"] == "substitution":
                continue
            d = abs(len(v["alt_allele"]) - len(v["ref_allele"]))
            if d % 3:
                for c in cds_list:
                    i = c.index.get(v["position"] - 1)
                    if i is not None:
                        fs_from[c.gene] = min(fs_from.get(c.gene, 10**9), i)
        snvs = {}   # (gene, codon_idx) -> {pos0: alt}
        for v in vs:
            if v["variant_type"] == "substitution":
                for c in cds_list:
                    i = c.index.get(v["position"] - 1)
                    if i is not None:
                        alt = v["alt_allele"].translate(COMP) if c.strand == -1 else v["alt_allele"]
                        snvs.setdefault((c.gene, i // 3), {})[i] = alt
        for v in vs:
            p0 = v["position"] - 1
            base = {**{k: v.get(k, "") for k in ("sample_id", "segment", "reference_accession", "position",
                     "ref_allele", "alt_allele", "variant_type", "aln_column")}, "note": v.get("note", "")}
            if v["variant_type"] == "substitution":
                hit = [c for c in cds_list if p0 in c.index]
                touched = [p0]
            else:
                lo = p0 + (1 if v["variant_type"] == "deletion" else 0)
                n_del = len(v["ref_allele"]) - 1
                touched = list(range(p0 + 1, p0 + 1 + n_del)) if v["variant_type"] == "deletion" else [p0, p0 + 1]
                hit = [c for c in cds_list if all(t in c.index for t in touched)] or \
                      [c for c in cds_list if any(t in c.index for t in touched)]
            dom = _domain(p0, domains)
            if not hit:
                out.append({**base, "gene": "", "domain": dom, "codon_number": "", "codon_change": "",
                            "protein_change": "", "consequence": "non_coding", "annotation_status": "non_coding"})
                continue
            for c in hit:
                row = {**base, "gene": c.gene, "domain": dom, "codon_number": "", "codon_change": "",
                       "protein_change": "", "consequence": "", "annotation_status": "annotated"}
                notes = [base["note"]] if base["note"] else []
                if len(hit) > 1:
                    row["annotation_status"] = "ambiguous"
                    notes.append("overlapping_cds:" + ";".join(h.gene for h in hit))
                if v["variant_type"] == "substitution":
                    i = c.index[p0]
                    k = i // 3
                    row["codon_number"] = k + 1
                    codon_pos = c.positions[3 * k:3 * k + 3]
                    ref_codon = c.seq[3 * k:3 * k + 3]
                    alt_map = snvs[(c.gene, k)]
                    alt_codon = "".join(alt_map.get(3 * k + j, ref_codon[j]) for j in range(len(ref_codon)))
                    row["codon_change"] = f"{ref_codon}>{alt_codon}"
                    review = []
                    if len(codon_pos) < 3:
                        review.append("partial_codon_at_cds_end")
                    if any(x in amb0 for x in codon_pos):
                        review.append("ambiguous_base_in_codon")
                    if cs is not None and any(not (cs - 1 <= x <= ce - 1) for x in codon_pos):
                        review.append("codon_partially_uncovered")
                    if c.gene in fs_from and i > fs_from[c.gene]:
                        review.append("downstream_of_frameshift")
                    ra, aa = c.aa(ref_codon), c.aa(alt_codon)
                    row["protein_change"] = f"{ra}{k + 1}{'=' if ra == aa else aa}"
                    row["consequence"] = ("synonymous" if ra == aa else "stop_gained" if aa == "*"
                                          else "stop_lost" if ra == "*" else "missense")
                    if len(alt_map) > 1:
                        notes.append("multiple_substitutions_in_codon")
                    if review:
                        row["annotation_status"] = "requires_review"
                        row["protein_change"], row["consequence"] = "", ""
                        notes += review
                else:
                    delta = len(v["alt_allele"]) - len(v["ref_allele"])
                    i = c.index.get(p0 + 1, c.index.get(p0))
                    row["codon_number"] = i // 3 + 1 if i is not None else ""
                    kind = "deletion" if delta < 0 else "insertion"
                    if delta % 3:
                        row["consequence"] = "frameshift"
                    else:
                        row["consequence"] = f"inframe_{kind}"
                    if cs is not None and (p0 < cs - 1 or p0 > ce - 1):
                        pass
                    row["annotation_status"] = row["annotation_status"] if len(hit) > 1 else "annotated"
                    if row["consequence"] == "frameshift":
                        notes.append("frameshift_candidate_check_sequence_quality")
                        row["annotation_status"] = "requires_review"
                row["note"] = ";".join(n for n in notes if n)
                out.append(row)
    cols = ["sample_id", "segment", "reference_accession", "position", "ref_allele", "alt_allele",
            "variant_type", "gene", "domain", "codon_number", "codon_change", "protein_change",
            "consequence", "annotation_status", "note", "aln_column"]
    return pd.DataFrame(out, columns=cols)


def main():
    cfg = load_config()
    tables = ROOT / cfg["paths"]["tables"]
    rec, cds, domains = load_reference(ROOT / cfg["paths"]["reference"] / "reference.gb")
    calls = pd.read_csv(tables / "mutation_calls.csv").fillna({"note": ""})
    cov = pd.read_csv(tables / "coverage.csv").fillna({"ambiguous_positions": ""})
    coverage = {r.sample_id: {"cov_start": r.cov_start, "cov_end": r.cov_end,
                              "ambiguous_positions": [int(x) for x in str(r.ambiguous_positions).split(";") if x]}
                for r in cov.itertuples()}
    ann = annotate(calls.to_dict("records"), rec, cds, domains, coverage)
    ann.to_csv(tables / "annotated_variants.csv", index=False)
    print(ann.consequence.value_counts().to_string())
    print(ann.annotation_status.value_counts().to_string())


if __name__ == "__main__":
    main()
