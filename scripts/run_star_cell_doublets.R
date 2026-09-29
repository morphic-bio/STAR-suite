#!/usr/bin/env Rscript

library(reticulate)
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)
library(anndata)

use_python("/usr/bin/python3", required = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: run_star_cell_doublets.R <counts.h5ad>")
}

h5ad_file <- args[1]
output_dir <- dirname(h5ad_file)
seed_text <- Sys.getenv("SCDBLFINDER_SEED", "1")
scdblfinder_seed <- suppressWarnings(as.integer(seed_text))
if (!grepl("^[+-]?[0-9]+$", seed_text) || is.na(scdblfinder_seed)) {
  stop("SCDBLFINDER_SEED must be a representable R integer")
}

ad <- read_h5ad(h5ad_file)
obs <- py_to_r(ad$obs)
if (class(ad$X)[1] != "dgCMatrix") {
  counts <- as(t(ad$X), "CsparseMatrix")
  counts <- as(counts, "dgCMatrix")
} else {
  counts <- as(t(ad$X), "dgCMatrix")
}

barcodes <- colnames(counts)
if (length(barcodes) == 0) {
  stop("counts.h5ad has no cell barcodes")
}

star_mask <- rep(TRUE, length(barcodes))
if ("is_cell" %in% colnames(obs)) {
  star_mask <- as.logical(obs[["is_cell"]])
} else if ("filter" %in% colnames(obs)) {
  star_mask <- as.logical(obs[["filter"]])
}

star_mask[is.na(star_mask)] <- FALSE
star_barcodes <- barcodes[star_mask]
if (length(star_barcodes) == 0) {
  stop("No STAR-called cells found in counts.h5ad")
}

writeLines(star_barcodes, file.path(output_dir, "non_empty_barcodes.txt"))

star_counts <- counts[, star_mask, drop = FALSE]
sce <- SingleCellExperiment(list(counts = star_counts))
barcode_order <- colnames(sce)

set.seed(scdblfinder_seed)
message("scDblFinder seed: ", scdblfinder_seed)
sce <- scDblFinder(sce, BPPARAM = BiocParallel::SerialParam(RNGseed = scdblfinder_seed))
result <- list(
  class = stats::setNames(as.character(sce$scDblFinder.class), colnames(sce)),
  score = stats::setNames(as.numeric(sce$scDblFinder.score), colnames(sce))
)
writeLines(as.character(scdblfinder_seed), file.path(output_dir, "scdblfinder_seed.txt"))

doublet_results <- data.frame(
  Barcode = barcode_order,
  Classification = unname(result$class[barcode_order]),
  Score = unname(result$score[barcode_order])
)
write.table(
  doublet_results,
  file.path(output_dir, "filtered_barcodes_with_scores.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

doublet_barcodes <- doublet_results$Barcode[doublet_results$Classification == "doublet"]
writeLines(doublet_barcodes, file.path(output_dir, "doublet_barcodes.txt"))

message("STAR-called cells: ", length(star_barcodes))
message("Doublets: ", length(doublet_barcodes))
