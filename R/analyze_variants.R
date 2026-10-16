# Stage 5: statistical summaries and figures for the annotated variant table.
#
# Inputs (from the Python pipeline):
#   results/tables/annotated_variants.csv
#   results/tables/coverage.csv
#   data/raw/metadata.csv
#
# Outputs:
#   results/tables/variants_per_gene.csv
#   results/tables/variant_frequency_by_period.csv
#   results/tables/sample_variant_matrix.csv
#   results/figures/fig1_variants_by_gene.png
#   results/figures/fig2_variant_frequency_over_time.png
#   results/figures/fig3_sample_variant_heatmap.png
#
# Design note: "number of variant sites" and "number of samples carrying a variant" are
# different quantities and are reported separately throughout (see project plan, stage 05
# quality check). Every frequency is reported with its denominator (n_eligible_samples).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# Run this script from the project root: `Rscript R/analyze_variants.R`
root <- getwd()
if (!file.exists(file.path(root, "config.yaml"))) {
  stop("Run this script from the project root (viral-mutation-tracker/), e.g.:\n",
       "  Rscript R/analyze_variants.R")
}

tables_dir  <- file.path(root, "results", "tables")
figures_dir <- file.path(root, "results", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

variants  <- read_csv(file.path(tables_dir, "annotated_variants.csv"), show_col_types = FALSE)
coverage  <- read_csv(file.path(tables_dir, "coverage.csv"), show_col_types = FALSE)
metadata  <- read_csv(file.path(root, "data", "raw", "metadata.csv"), show_col_types = FALSE)
excluded  <- tryCatch(read_csv(file.path(tables_dir, "excluded_samples.csv"), show_col_types = FALSE),
                       error = function(e) tibble(sample_id = character()))

# Eligible analysis set = samples that passed QC/alignment and were actually variant-called.
analyzed_samples <- coverage$sample_id
meta <- metadata %>%
  filter(sample_id %in% analyzed_samples) %>%
  mutate(collection_year = as.integer(collection_year))

n_eligible <- nrow(meta)
n_selected <- sum(metadata$selected, na.rm = TRUE)
cat(sprintf("Analyzed samples: %d (of %d selected; %d excluded as lineage mismatches: %s)\n",
            n_eligible, n_selected, nrow(excluded), paste(excluded$sample_id, collapse = ", ")))

substitutions <- variants %>% filter(variant_type == "substitution")

# ---------------------------------------------------------------------------
# Figure 1 — variant counts by gene (distinct variant SITES, not samples)
# ---------------------------------------------------------------------------
by_gene <- substitutions %>%
  mutate(gene = if_else(gene == "" | is.na(gene), "non_coding", gene)) %>%
  distinct(gene, position, ref_allele, alt_allele) %>%
  count(gene, name = "n_variant_sites") %>%
  arrange(desc(n_variant_sites))

write_csv(by_gene, file.path(tables_dir, "variants_per_gene.csv"))

fig1 <- ggplot(by_gene, aes(x = reorder(gene, -n_variant_sites), y = n_variant_sites)) +
  geom_col(fill = "#2c7fb8") +
  labs(title = "Figure 1. Distinct variant sites by gene",
       subtitle = sprintf("Influenza A/H1N1 HA segment, n = %d analyzed sequences", n_eligible),
       x = "Gene", y = "Number of distinct variant sites") +
  theme_minimal(base_size = 13)
ggsave(file.path(figures_dir, "fig1_variants_by_gene.png"), fig1, width = 7, height = 5, dpi = 150)

# ---------------------------------------------------------------------------
# Figure 2 — variant frequency over time, with explicit denominators
# ---------------------------------------------------------------------------
per_sample_gene <- substitutions %>%
  mutate(gene = if_else(gene == "" | is.na(gene), "non_coding", gene)) %>%
  distinct(sample_id, gene, position, ref_allele, alt_allele) %>%
  left_join(meta %>% select(sample_id, collection_year), by = "sample_id")

# The most frequent non-synonymous sites overall (site-level, sample-carrying counts)
top_sites <- substitutions %>%
  filter(consequence == "missense") %>%
  distinct(sample_id, position, protein_change) %>%
  count(position, protein_change, name = "n_samples_with_variant") %>%
  arrange(desc(n_samples_with_variant)) %>%
  slice_head(n = 6)

site_by_year <- substitutions %>%
  semi_join(top_sites, by = "position") %>%
  distinct(sample_id, position, protein_change) %>%
  left_join(meta %>% select(sample_id, collection_year), by = "sample_id") %>%
  count(collection_year, position, protein_change, name = "n_samples_with_variant")

denom_by_year <- meta %>% count(collection_year, name = "n_eligible_samples")

freq_by_year <- site_by_year %>%
  left_join(denom_by_year, by = "collection_year") %>%
  mutate(frequency = n_samples_with_variant / n_eligible_samples,
         label = sprintf("%s (pos %d)", protein_change, position))

write_csv(freq_by_year %>%
            select(collection_year, position, protein_change,
                   n_samples_with_variant, n_eligible_samples, frequency),
          file.path(tables_dir, "variant_frequency_by_period.csv"))

fig2 <- ggplot(freq_by_year, aes(x = collection_year, y = frequency, color = label)) +
  geom_line() + geom_point(aes(size = n_eligible_samples)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Figure 2. Observed frequency of the most common missense variants over time",
       subtitle = "Point size = number of eligible sequences analyzed that collection year (denominator)",
       x = "Collection year", y = "Proportion of eligible sequences carrying the variant",
       color = "Variant", size = "n eligible") +
  theme_minimal(base_size = 13)
ggsave(file.path(figures_dir, "fig2_variant_frequency_over_time.png"), fig2, width = 10, height = 5, dpi = 150)

# ---------------------------------------------------------------------------
# Figure 3 — sample-by-variant heatmap (top variant sites x all analyzed samples)
# ---------------------------------------------------------------------------
heat_sites <- substitutions %>%
  distinct(sample_id, position, protein_change) %>%
  count(position, protein_change, name = "n") %>%
  arrange(desc(n)) %>%
  slice_head(n = 20)

mat <- substitutions %>%
  semi_join(heat_sites, by = "position") %>%
  distinct(sample_id, position, protein_change) %>%
  mutate(present = 1L) %>%
  right_join(expand_grid(sample_id = meta$sample_id,
                         heat_sites %>% select(position, protein_change)),
             by = c("sample_id", "position", "protein_change")) %>%
  mutate(present = coalesce(present, 0L))

write_csv(mat %>% pivot_wider(id_cols = sample_id, names_from = c(position, protein_change),
                               values_from = present, values_fill = 0),
          file.path(tables_dir, "sample_variant_matrix.csv"))

mat_plot <- mat %>%
  left_join(meta %>% select(sample_id, collection_year), by = "sample_id") %>%
  mutate(label = sprintf("%d:%s", position, protein_change))

fig3 <- ggplot(mat_plot, aes(x = reorder(label, position),
                             y = reorder(sample_id, collection_year),
                             fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("0" = "grey90", "1" = "#d7301f"),
                     labels = c("absent", "present"), name = "") +
  labs(title = "Figure 3. Sample-by-variant matrix (top 20 variant sites)",
       x = "Variant (reference position : protein change)", y = "Sample (ordered by collection year)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5))
ggsave(file.path(figures_dir, "fig3_sample_variant_heatmap.png"), fig3, width = 9, height = 10, dpi = 150)

# ---------------------------------------------------------------------------
# Descriptive summary table
# ---------------------------------------------------------------------------
summary_tbl <- substitutions %>%
  group_by(consequence) %>%
  summarise(n_calls = n(), n_distinct_sites = n_distinct(position, ref_allele, alt_allele),
            n_samples = n_distinct(sample_id), .groups = "drop")
write_csv(summary_tbl, file.path(tables_dir, "consequence_summary.csv"))

cat("\n== Consequence summary ==\n"); print(summary_tbl)
cat("\nFigures written to", figures_dir, "\n")
cat("Tables written to", tables_dir, "\n")
