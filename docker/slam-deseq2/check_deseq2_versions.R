#!/usr/bin/env Rscript

expected_deseq2 <- Sys.getenv("DESEQ2_VERSION", "1.52.0")
expected_tximport <- Sys.getenv("TXIMPORT_VERSION", "1.40.0")
expected_bioc <- Sys.getenv("BIOCONDUCTOR_VERSION", "3.23")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager is not installed")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is not installed")
}
if (!requireNamespace("tximport", quietly = TRUE)) {
  stop("tximport is not installed")
}

actual_bioc <- as.character(BiocManager::version())
actual_deseq2 <- as.character(utils::packageVersion("DESeq2"))
actual_tximport <- as.character(utils::packageVersion("tximport"))

cat("R_version\t", R.version.string, "\n", sep = "")
cat("Bioconductor_version\t", actual_bioc, "\n", sep = "")
cat("DESeq2_version\t", actual_deseq2, "\n", sep = "")
cat("tximport_version\t", actual_tximport, "\n", sep = "")

if (!identical(actual_bioc, expected_bioc)) {
  stop(sprintf("Bioconductor version mismatch: expected %s, got %s", expected_bioc, actual_bioc))
}
if (!identical(actual_deseq2, expected_deseq2)) {
  stop(sprintf("DESeq2 version mismatch: expected %s, got %s", expected_deseq2, actual_deseq2))
}
if (!identical(actual_tximport, expected_tximport)) {
  stop(sprintf("tximport version mismatch: expected %s, got %s", expected_tximport, actual_tximport))
}

optional <- c("apeglm", "ashr", "pheatmap", "RColorBrewer", "readr", "data.table")
for (pkg in optional) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "_version\t", as.character(utils::packageVersion(pkg)), "\n", sep = "")
  } else {
    cat(pkg, "_version\tMISSING\n", sep = "")
  }
}
