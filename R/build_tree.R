# Stage 6b: evolutionary context - a neighbor-joining tree from the reference-anchored MSA.
#
# Method and assumptions (documented per project plan, stage 06 quality check):
#   - Distances: raw (Hamming) pairwise distance between aligned sequences, as a fraction of
#     compared sites, ignoring columns where either sequence has a gap/ambiguous base.
#   - Tree: neighbor-joining (ape::nj), not maximum-likelihood. NJ is a reasonable, fast method for
#     a small set of closely related sequences but assumes a simple distance-additive model; it does
#     not test alternative topologies or evaluate substitution-model fit the way an ML method would.
#   - The input alignment is reference-anchored (see src/build_msa.py): insertions relative to the
#     reference are dropped, so this is an approximation of a true multiple-sequence alignment.
#   - This tree shows INFERRED relationships based on observed sequence similarity. It does not by
#     itself establish direct transmission, clinical severity, or the cause of any mutation.
#
# Inputs: data/processed/msa_reference_anchored.fasta, data/raw/metadata.csv
# Outputs: results/tables/tree.nwk, results/figures/fig4_phylogenetic_tree.png

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

root <- getwd()
if (!file.exists(file.path(root, "config.yaml"))) {
  stop("Run this script from the project root: Rscript R/build_tree.R")
}

msa <- read.dna(file.path(root, "data", "processed", "msa_reference_anchored.fasta"), format = "fasta")
meta <- read_csv(file.path(root, "data", "raw", "metadata.csv"), show_col_types = FALSE)

# pairwise.deletion = TRUE: for each pair, compare only columns where both sequences have a called
# base, so ambiguous/gap positions are excluded from that pair's distance rather than treated as
# a difference against an arbitrary resolved base.
d <- dist.dna(msa, model = "raw", pairwise.deletion = TRUE)
tree <- nj(d)
tree <- ladderize(tree)

write.tree(tree, file.path(root, "results", "tables", "tree.nwk"))

year_lookup <- setNames(meta$collection_year, meta$sample_id)
year_lookup[["NC_026433.1"]] <- 2009L  # reference: A/California/07/2009
tip_year <- year_lookup[tree$tip.label]

png(file.path(root, "results", "figures", "fig4_phylogenetic_tree.png"),
    width = 1600, height = 2000, res = 150)
pal <- colorRampPalette(c("#2c7fb8", "#d7301f"))(length(unique(na.omit(tip_year))))
year_levels <- sort(unique(na.omit(tip_year)))
tip_col <- ifelse(tree$tip.label == "NC_026433.1", "black",
                   pal[match(tip_year, year_levels)])
tip_label <- ifelse(tree$tip.label == "NC_026433.1",
                     paste0(tree$tip.label, " [reference]"),
                     paste0(tree$tip.label, " (", tip_year, ")"))
tree_for_plot <- tree
tree_for_plot$tip.label <- tip_label   # display labels only; results/tables/tree.nwk keeps accessions
plot(tree_for_plot, tip.color = tip_col, cex = 0.55, label.offset = 0.0005,
     main = "Figure 4. Neighbor-joining tree of H1N1 HA sequences\n(colored by collection year; black = reference)")
legend("bottomleft", legend = range(year_levels), fill = pal[c(1, length(pal))],
       title = "Collection year", bty = "n", cex = 0.7)
dev.off()

cat(sprintf("Tree with %d tips written to results/tables/tree.nwk and results/figures/fig4_phylogenetic_tree.png\n",
            length(tree$tip.label)))
