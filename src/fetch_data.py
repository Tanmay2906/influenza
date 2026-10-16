"""Stage 1: download the reference and a year-stratified sample of H1N1 HA records from NCBI.

Outputs: data/reference/reference.gb + reference.fasta, data/raw/samples.fasta,
data/raw/metadata.csv, data/raw/retrieval_log.json.
"""
import datetime as dt
import json
import random
import re
import time

import pandas as pd
from Bio import Entrez, SeqIO
from io import StringIO

from common import ROOT, load_config, use_system_trust_store

QUERY = ('"Influenza A virus"[Organism] AND H1N1 AND hemagglutinin[Title] AND "complete cds"[Title] '
         'AND {lo}:{hi}[SLEN] AND {year}[PDAT]')


def efetch_gb(ids):
    out = []
    for i in range(0, len(ids), 100):
        chunk = ids[i:i + 100]
        for attempt in range(4):
            try:
                h = Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="gb", retmode="text")
                out.extend(SeqIO.parse(StringIO(h.read()), "genbank"))
                break
            except Exception:
                time.sleep(2 * (attempt + 1))
        time.sleep(0.4)
    return out


def record_metadata(rec):
    src = next((f for f in rec.features if f.type == "source"), None)
    q = src.qualifiers if src else {}
    g = lambda k: (q.get(k) or [""])[0]
    cds = [f for f in rec.features if f.type == "CDS"]
    date = g("collection_date")
    # Accept YYYY, YYYY-MM, YYYY-MM-DD, DD-Mon-YYYY and Mon-YYYY; anything else stays missing.
    m = re.match(r"^(\d{4})(?:-\d{2}){0,2}$", date) or re.match(r"^(?:\d{2}-)?[A-Za-z]{3}-(\d{4})$", date)
    return {
        "sample_id": rec.id,                      # accession.version
        "strain": g("strain") or g("isolate"),
        "collection_date": date,
        "collection_year": int(m.group(1)) if m else None,
        "country": g("country"),
        "host": g("host"),
        "segment": g("segment"),
        "serotype": g("serotype"),
        "length": len(rec.seq),
        "n_cds": len(cds),
        "description": rec.description,
    }


def main():
    cfg = load_config()
    use_system_trust_store()
    Entrez.email = cfg["ncbi"]["email"]
    s, qc = cfg["study"], cfg["qc"]
    rng = random.Random(s["seed"])
    raw, ref_dir = ROOT / cfg["paths"]["raw"], ROOT / cfg["paths"]["reference"]

    # Reference
    h = Entrez.efetch(db="nuccore", id=cfg["reference"]["accession"], rettype="gb", retmode="text")
    gb_text = h.read()
    (ref_dir / "reference.gb").write_text(gb_text)
    ref = SeqIO.read(StringIO(gb_text), "genbank")
    SeqIO.write(ref, ref_dir / "reference.fasta", "fasta")
    print(f"reference {ref.id}: {len(ref.seq)} bp")

    # Candidate pool by publication year, then assign collection year from the record itself
    rows, recs = [], {}
    for year in range(s["year_start"], s["year_end"] + 1 + 2):
        term = QUERY.format(lo=qc["min_length"], hi=qc["max_length"], year=year)
        r = Entrez.read(Entrez.esearch(db="nuccore", term=term, retmax=cfg["ncbi"]["candidates_per_year"]))
        ids = r["IdList"]
        print(f"pub year {year}: {r['Count']} hits, fetching {len(ids)}")
        for rec in efetch_gb(ids):
            if rec.id in recs or rec.id == ref.id:
                continue
            recs[rec.id] = rec
            rows.append(record_metadata(rec))

    df = pd.DataFrame(rows)
    # Selection: human host, strain-level H1N1, collection year within study window, complete collection year.
    ok = (df.collection_year.between(s["year_start"], s["year_end"])
          & df.host.str.lower().str.startswith("homo sapiens")
          & df.description.str.contains(r"\(H1N1\)", regex=True)
          & df.n_cds.ge(1))
    df["eligible"] = ok
    df["exclusion_reason"] = ""
    df.loc[~df.collection_year.between(s["year_start"], s["year_end"]), "exclusion_reason"] = "collection year missing/outside window"
    df.loc[df.collection_year.between(s["year_start"], s["year_end"]) & ~df.host.str.lower().str.startswith("homo sapiens"), "exclusion_reason"] = "host not human"
    df.loc[ok.eq(False) & df.exclusion_reason.eq(""), "exclusion_reason"] = "not H1N1 strain-level / no CDS"

    # Stratified random sample per year
    chosen = []
    for year, grp in df[df.eligible].groupby("collection_year"):
        ids = sorted(grp.sample_id)
        chosen += rng.sample(ids, min(s["per_year_target"], len(ids)))
    df["selected"] = df.sample_id.isin(chosen)
    df.loc[df.eligible & ~df.selected, "exclusion_reason"] = "not sampled (year quota reached)"
    df.sort_values(["collection_year", "sample_id"]).to_csv(raw / "metadata.csv", index=False)
    SeqIO.write([recs[i] for i in sorted(chosen)], raw / "samples.fasta", "fasta")

    log = {"retrieval_date": dt.date.today().isoformat(), "database": "NCBI nuccore (GenBank)",
           "reference": ref.id, "n_candidates": len(df), "n_eligible": int(df.eligible.sum()),
           "n_selected": len(chosen), "config_study": s}
    (raw / "retrieval_log.json").write_text(json.dumps(log, indent=2))
    print(json.dumps(log, indent=2))
    print(df[df.selected].groupby("collection_year").size().to_string())


if __name__ == "__main__":
    main()
