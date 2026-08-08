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

doublet_info <- rep("singlet", ncol(sce))
doublet_scores <- rep(NA_real_, ncol(sce))
names(doublet_info) <- barcode_order
names(doublet_scores) <- barcode_order

result <- tryCatch(
  {
    sce <- scDblFinder(sce)
    list(
      class = stats::setNames(as.character(sce$scDblFinder.class), colnames(sce)),
      score = stats::setNames(as.numeric(sce$scDblFinder.score), colnames(sce))
    )
  },
  error = function(err) {
    message("WARNING: scDblFinder failed, marking all STAR cells as singlets: ", err$message)
    list(class = doublet_info, score = doublet_scores)
  }
)

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
