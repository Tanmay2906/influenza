# Viral Mutation Tracker

Comparative genomics pipeline for publicly available influenza A/H1N1 hemagglutinin (HA) sequences:
quality control, reference alignment, nucleotide variant calling, codon/protein annotation (Python),
and statistical summaries, figures and a phylogenetic tree (R).

**Research question:** how does the genetic variation of publicly available influenza A/H1N1 HA
sequences change across 2009-2019, and what nucleotide- and protein-level differences can be
identified relative to an annotated reference? See
[`docs/methods_and_limitations.md`](docs/methods_and_limitations.md) for the full methods, data
provenance and limitations, and [`docs/report.md`](docs/report.md) for the results write-up.

This analyzes already-deposited public sequences; it does not perform laboratory sequencing, and a
sequence comparison alone does not establish when a mutation arose or what phenotypic effect it has.

## Repository layout

```
viral-mutation-tracker/
├── config.yaml                 # single source of truth: reference, period, QC thresholds, etc.
├── data/
│   ├── raw/                    # downloaded sequences + metadata + retrieval log (stage 1)
│   ├── reference/              # the fixed reference record (GenBank + FASTA)
│   └── processed/              # QC-passed FASTA, pairwise alignments, reference-anchored MSA
├── src/                        # Python pipeline (stages 1-4, 6a)
├── R/                          # R analysis + tree building (stages 5, 6b)
├── tests/                      # pytest synthetic validation suite (stage 7)
├── results/
│   ├── tables/                 # qc_summary.csv, mutation_calls.csv, annotated_variants.csv, ...
│   └── figures/                # fig1-fig4 PNGs
└── docs/                       # methods, limitations, report, manual validation spot-check
```

## Setup

Requires Python 3.10+ and R 4.x.

```bash
git clone <this-repo> && cd viral-mutation-tracker
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
Rscript R/install_packages.R
```

## Running the pipeline

```bash
# Stage 1: download reference + samples from NCBI (needs network access)
python src/fetch_data.py

# Stages 2-4: QC, alignment, variant calling, annotation
python src/run_pipeline.py

# Stage 6a: project all samples onto reference coordinates (for the tree)
python src/build_msa.py

# Stage 5: R statistics and figures 1-3
Rscript R/analyze_variants.R

# Stage 6b: neighbor-joining tree and figure 4
Rscript R/build_tree.R

# Tests
pytest tests/ -v
```

Run the Python commands from the repository root (`src/` modules import `common.py` as a sibling
module when run as scripts) and the R commands from the repository root as well — both scripts
check for `config.yaml` in the current directory and fail with a clear message otherwise.

## Example outputs

| File | Contents |
| --- | --- |
| `results/tables/qc_summary.csv` | per-sample QC status (pass/review/fail) and flags |
| `results/tables/mutation_calls.csv` | raw nucleotide variant calls, reference-based coordinates |
| `results/tables/coverage.csv` | per-sample alignment coverage and ambiguous-base positions |
| `results/tables/annotated_variants.csv` | gene/codon/protein annotation per variant |
| `results/tables/excluded_samples.csv` | samples excluded as a different H1N1 lineage |
| `results/tables/variants_per_gene.csv` | Figure 1 source data |
| `results/tables/variant_frequency_by_period.csv` | Figure 2 source data, with denominators |
| `results/tables/sample_variant_matrix.csv` | Figure 3 source data |
| `results/tables/tree.nwk` | Newick tree (Figure 4) |
| `results/figures/fig1..fig4_*.png` | the four analysis figures |
| `docs/manual_validation_spotcheck.csv` | independent re-derivation of 4 sample variant calls |

## Current dataset (this run)

- Reference: `NC_026433.1`, A/California/07/2009(H1N1) segment 4 (HA)
- 88 sequences selected (8/year, 2009-2019); 85 carried forward after QC and a lineage-identity
  check excluded 3 non-pdm09 records (see `results/tables/excluded_samples.csv`)
- 2,994 substitution calls (1,806 synonymous, 1,188 missense), 318 distinct variant sites; no
  indels observed in this HA dataset

## Validation

- `tests/`: 13 pytest cases against synthetic sequences with known expected results (identical
  sequences, single substitution, single indel, ambiguous base, truncated sequence, substitution
  near a gap, coding substitution, synonymous substitution, stop-gain)
- `docs/manual_validation_spotcheck.csv`: 4 real variant calls independently re-derived from the
  saved alignment and reference CDS, outside the pipeline code, all matching

## Software versions

See `requirements.txt` for pinned Python package minimums and `R/install_packages.R` for the R
packages used (dplyr, tidyr, readr, ggplot2, stringr, ape).

## Citation

If you use this pipeline, please cite the reference sequence (`NC_026433.1`) and this repository.
