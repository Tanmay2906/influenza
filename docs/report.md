# Viral Mutation Tracker: Comparative Genomics of Influenza A/H1N1 HA (2009-2019)

*Independent bioinformatics project. See [`methods_and_limitations.md`](methods_and_limitations.md)
for full methodological detail, data provenance and known limitations; this report summarizes the
research question, approach and results.*

## Research question

How does the genetic variation of publicly available influenza A/H1N1 hemagglutinin (segment 4)
sequences change across a defined 2009-2019 sampling period, and what nucleotide- and protein-level
differences can be identified relative to an annotated reference?

## Approach

A stratified sample of up to 8 human-host H1N1 HA GenBank records per collection year (2009-2019,
88 total) was downloaded alongside the RefSeq reference `NC_026433.1`
(A/California/07/2009(H1N1), segment 4). Each sample was quality-controlled, pairwise-aligned to
the reference, and screened for lineage identity. Reference-based nucleotide variants were called
with explicit handling of missing/ambiguous data and left-normalized indel coordinates, then
annotated against the reference's CDS to predict codon and amino-acid consequences. R was used for
descriptive statistics, four figures, and a neighbor-joining tree. Full detail:
[`methods_and_limitations.md`](methods_and_limitations.md).

## Results

- **Dataset:** 88 sequences selected; 3 (all 2009 records) were excluded before variant calling
  because their identity to the pdm09 reference (~77%) was far below every other sample's (95-99%),
  indicating a different, seasonal H1N1 HA lineage rather than a comparable pdm09 sequence. 85
  samples were analyzed.
- **QC:** 79/88 passed cleanly; 9/88 were flagged `review` for containing ambiguous (IUPAC) bases,
  which were retained and tracked rather than resolved to a specific nucleotide; 0 failed QC.
- **Variants:** 2,994 substitution calls across 85 samples, corresponding to 318 distinct
  nucleotide sites relative to the reference. No insertions or deletions were observed in this HA
  dataset.
- **Consequences:** of the 2,994 substitution calls, 1,806 (60%) were synonymous and 1,188 (40%)
  were missense; none required withholding a protein call (no ambiguous, partially covered, or
  frameshift-downstream codons occurred in this sample set).
- **Figures** (`results/figures/`):
  - *Figure 1* — distinct variant sites per gene (all in `HA`, the only annotated CDS on this
    reference).
  - *Figure 2* — the proportion of eligible sequences per collection year carrying each of the
    6 most common missense variants, with the per-year sample count (denominator) shown explicitly.
  - *Figure 3* — a sample x variant-site presence/absence matrix for the top 20 variant sites,
    ordered by collection year, showing which variants are shared across years versus
    sample-specific.
  - *Figure 4* — a neighbor-joining tree of all 85 analyzed sequences plus the reference, colored
    by collection year, built from a reference-anchored multiple-sequence alignment.
- **Independent validation:** 4 real variant calls (3 missense, 1 synonymous) were re-derived from
  the saved alignment and reference CDS coordinates using code outside the pipeline itself;
  all 4 matched the pipeline's output exactly
  (`docs/manual_validation_spotcheck.csv`).

## Discussion

The dataset shows substantial nucleotide divergence from the 2009 pdm09 reference accumulating
over the sampling window, with missense substitutions comprising a sizeable minority (40%) of
observed changes — consistent with ongoing antigenic drift in HA, a gene known to be under
diversifying selection at antibody epitopes. The reference-anchored tree (Figure 4) shows
sequences broadly clustering by collection year, consistent with temporal accumulation of
divergence, though this is an NJ tree with pairwise-deletion handling of missing data, not a
maximum-likelihood phylogeny.

## What this analysis does *not* establish

- It does not determine *when* any individual mutation first arose — only that it differs from a
  fixed 2009 reference in this curated sample.
- It does not assess phenotypic effect (antigenicity, transmissibility, drug resistance) of any
  observed substitution.
- The tree shows inferred relationships from observed sequence similarity; it is not evidence of
  direct transmission between sampled isolates.
- With ~8 sequences per year from a non-systematic GenBank sample, year-to-year frequency
  differences in Figure 2 should be read as descriptive of this dataset, not as population-level
  prevalence estimates.

## Reproducing this analysis

See [`README.md`](../README.md) for exact commands. The full pipeline, from NCBI download through
figures and tree, runs with `python src/fetch_data.py && python src/run_pipeline.py &&
python src/build_msa.py && Rscript R/analyze_variants.R && Rscript R/build_tree.R`.
