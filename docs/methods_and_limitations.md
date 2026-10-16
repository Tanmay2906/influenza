# Methods, limitations and data provenance

## Research question

How does the genetic variation of publicly available influenza A/H1N1 hemagglutinin (HA, segment 4)
sequences change across 2009-2019, and what nucleotide- and protein-level differences can be
identified relative to an annotated reference?

This is a sequence-comparison study. It identifies differences between public GenBank records and a
fixed reference; it does not independently establish when a mutation arose, how it spread, or what
phenotypic effect it has. Those questions would require additional epidemiological, temporal and
experimental evidence beyond what is analyzed here.

## Data

- **Database:** NCBI GenBank (`nuccore`), accessed via the NCBI Entrez E-utilities.
- **Retrieval date:** recorded automatically in `data/raw/retrieval_log.json` at the time
  `src/fetch_data.py` was run.
- **Reference:** `NC_026433.1`, Influenza A virus (A/California/07/2009(H1N1)) segment 4 (HA),
  RefSeq, 1701 nt, annotated CDS (gene `HA`, with `sig_peptide`/`mat_peptide` HA1/HA2 subregions).
- **Segment:** only segment 4 (HA) is analyzed. Influenza A has 8 non-contiguous genome segments;
  this pipeline never concatenates segments or compares a sample's HA sequence against a
  non-HA reference.
- **Samples:** up to 8 records per collection year (2009-2019), selected by:
  - search term: GenBank title contains "H1N1" and "hemagglutinin" and "complete cds"
  - sequence length between 1650-1800 nt (near full-length HA CDS)
  - host recorded as human (`Homo sapiens`, including qualified forms such as
    `Homo sapiens; male`)
  - a parseable collection date within 2009-2019
  - stratified random sampling within each eligible year (fixed seed, see `config.yaml`)
  - **inclusion/exclusion criteria and counts** are recorded per-record in
    `data/raw/metadata.csv` (`eligible`, `selected`, `exclusion_reason` columns) and summarized in
    `data/raw/retrieval_log.json`.
- 88 records were selected (8/year x 11 years). Three (`GQ457548.1`, `GQ476049.1`, `GQ476108.1`,
  all 2009 "Saskatchewan"/"Massachusetts"/"New York" strains) were later excluded at the alignment
  stage: see Lineage-mismatch filter below.

## Quality control (`src/validate_sequences.py`)

- Nucleotide alphabet is checked against IUPAC codes; unrecognized characters fail QC.
- Ambiguous bases (`N` and other IUPAC ambiguity codes) are **counted, never silently resolved**
  to a specific nucleotide. A record with any ambiguous bases is flagged `review` (kept, with the
  positions recorded) rather than failed outright; a record above
  `qc.max_ambiguous_fraction` (1%, see `config.yaml`) fails.
- Length is checked against the expected HA CDS window (1650-1800 nt); shorter/longer records fail
  rather than being treated as containing a biological deletion.
- Duplicate accessions and identical sequences are both detected and recorded (duplicates fail;
  identical-but-distinct-accession sequences are flagged informationally, not failed, since
  biological identity between distinct isolates is possible).
- Output: `results/tables/qc_summary.csv` - 79 pass, 9 review (ambiguous bases present, otherwise
  clean), 0 fail, out of 88.

## Alignment (`src/align_sequences.py`)

- Method: Biopython `PairwiseAligner`, global mode, with **free end gaps** (`end_gap_score = 0`) so
  that a truncated/incomplete query is not penalized as if it contained a real deletion.
  Match/mismatch/gap-open/gap-extend scores are fixed in `config.yaml` and documented there.
- A single pairwise alignment per sample is a defensible starting point for ~90 closely related
  sequences of one gene; the project plan notes a multiple-sequence aligner is preferable for larger
  collections, and we use the reference as every other sequence's common alignment point.
- **Lineage-mismatch filter:** per-sample percent identity to the reference (over aligned,
  non-gap columns) is computed. Three 2009 records aligned at only ~77% identity versus ~95-99%
  for all other samples - consistent with a seasonal (pre-pandemic) H1N1 HA rather than the
  2009 pandemic (pdm09) lineage the reference represents. These are excluded from variant calling
  (`results/tables/excluded_samples.csv`) rather than forced through a comparison that would produce
  meaningless "variants" against an unrelated lineage.

## Variant calling (`src/call_variants.py`)

- Coordinates are **1-based, reference-relative**. The implementation deliberately keeps track of
  the 0-based Python string index separately from the reported 1-based biological position
  (see `call_variants.call_variants`).
- Indels use a left-normalized, VCF-style anchored representation: `position` is the reference base
  immediately before the event, `ref_allele`/`alt_allele` both include that anchor base. Indel
  placement is normalized to the leftmost equivalent alignment (e.g., a 1-nt deletion inside a
  homopolymer/repeat is always anchored at the same position regardless of which copy the aligner
  happened to delete).
- **Missing data is distinguished from true variants**: leading/trailing alignment columns outside
  both sequences' overlap (`cov_start`/`cov_end` in `results/tables/coverage.csv`) are *uncovered*,
  not deletions; query positions with a non-ACGT base are *ambiguous*, not substitutions, and are
  listed separately per sample.
- Result: 2994 substitution calls across the 85 lineage-consistent samples; no indels were observed
  in this HA dataset (consistent with HA being under strong purifying selection against frameshifts
  in a near-full-length CDS).

## Annotation (`src/annotate_variants.py`)

- Uses the reference GenBank record's actual CDS coordinates, strand and genetic code
  (`transl_table`), not an assumed reading frame. On load, the derived CDS translation is checked
  against the record's own `/translation` qualifier as a self-consistency check.
- For substitutions, the reference and alternate codon are reconstructed from all substitutions
  observed in that sample within the same codon (not just the single variant being annotated), so a
  double substitution in one codon is translated jointly rather than as two independent single-base
  changes.
- Protein consequences are **withheld** (`annotation_status = requires_review`, `protein_change`
  and `consequence` left blank) when: the codon is only partially covered by the alignment, the
  codon contains an ambiguous base, or the codon falls downstream of a frameshifting indel in that
  same sample. This prevents reporting a confident amino-acid call built partly from missing or
  ambiguous evidence.
- HA has one annotated CDS in this reference (`gene=HA`); no overlapping-CDS cases occur, so no
  variant needed multi-gene disambiguation in this run. The code supports it generically
  (`annotation_status = ambiguous`) for segments where this does occur (e.g. M1/M2, NS1/NEP), which
  is why segment expansion beyond HA will need this path re-tested.
- Result: of 2994 substitution calls, 1806 synonymous and 1188 missense; 0 required review (no
  ambiguous/partial/frameshift-downstream codons occurred among the 85 analyzed samples).

## Independent validation

- `tests/` (pytest): 13 unit tests covering the synthetic cases from the project plan - identical
  sequences, a single substitution, a single indel, an ambiguous base, a truncated/incomplete
  sequence, a substitution adjacent to an indel, a coding substitution with known codon/protein
  change, a synonymous substitution, and a stop-gain.
- `docs/manual_validation_spotcheck.csv`: four real variant calls (three missense, one synonymous)
  were independently re-derived directly from the saved pairwise alignment and the reference CDS
  coordinates, in code that does not reuse `call_variants`/`annotate_variants`, and compared against
  the pipeline's own output. All four matched exactly.

## R analysis (`R/analyze_variants.R`, `R/build_tree.R`)

- Figure 1: count of **distinct variant sites** (not calls, not samples) per gene.
- Figure 2: proportion of **eligible sequences per collection year** carrying each of the most
  common missense variants, with the per-year denominator (`n_eligible_samples`) reported alongside
  every frequency and shown as point size on the plot.
- Figure 3: presence/absence matrix of the top 20 variant sites across all analyzed samples,
  ordered by collection year.
- Figure 4: a neighbor-joining tree (`ape::nj`, raw/Hamming pairwise distance with pairwise
  deletion of gapped/ambiguous columns) built from a **reference-anchored** multiple-sequence
  alignment (`src/build_msa.py`): every sample is aligned to the reference with the same pairwise
  method used for variant calling, and insertions relative to the reference are dropped so every
  row has exactly the reference's 1701 columns. This is a simplification of a true multiple
  sequence alignment and should not be over-interpreted for insertion-containing lineages (there
  were none in this dataset). NJ is a distance method, not maximum likelihood; it does not evaluate
  alternative topologies or substitution-model fit, and the tree shows **inferred** relationships
  from observed similarity, not confirmed transmission chains, clinical severity or causality.

## Known limitations

- Sample size is modest (85 analyzed sequences, ~8/year) and was chosen by stratified random
  sampling from GenBank search results, not a systematic/representative epidemiological sample —
  frequencies in Figure 2 describe this specific curated dataset, not global H1N1 prevalence, and
  year-to-year comparisons with only 8 sequences per year are statistically underpowered to detect
  anything but large frequency differences.
- Only one segment (HA) and one reference lineage (pdm09) are analyzed in this version; the three
  excluded seasonal-H1N1 records show the dataset is not lineage-homogeneous at the source, and a
  real pre-pdm09/pdm09 comparison would need two reference lineages, not an identity filter.
- The reference-anchored MSA cannot represent insertions relative to the reference; a true
  progressive MSA tool (e.g. MAFFT) would be needed before extending this analysis to segments or
  lineages with indels.
- GenBank metadata (host, collection date, geography) is self-reported by submitters and of
  variable completeness/precision (e.g., year-only dates); records lacking a usable collection year
  were excluded rather than imputed.
