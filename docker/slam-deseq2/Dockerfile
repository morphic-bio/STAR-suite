# STAR-suite SLAM DESeq2 analysis runtime.
#
# Pin the Bioconductor/R image and assert the exact DESeq2 version during build.
# Sources checked 2026-05-11:
#   - Bioconductor 3.23 release uses R 4.6.
#   - DESeq2 release version is 1.52.0.

ARG BIOCONDUCTOR_IMAGE=bioconductor/bioconductor_docker:RELEASE_3_23-r-4.6.0
FROM ${BIOCONDUCTOR_IMAGE}

ARG BIOCONDUCTOR_VERSION=3.23
ARG DESEQ2_VERSION=1.52.0
ARG TXIMPORT_VERSION=1.40.0

ENV BIOCONDUCTOR_VERSION=${BIOCONDUCTOR_VERSION}
ENV DESEQ2_VERSION=${DESEQ2_VERSION}
ENV TXIMPORT_VERSION=${TXIMPORT_VERSION}

RUN Rscript -e ' \
    options(repos = BiocManager::repositories(version = Sys.getenv("BIOCONDUCTOR_VERSION"))); \
    pkgs <- c("DESeq2", "tximport", "apeglm", "ashr", "pheatmap", "RColorBrewer", "readr", "data.table"); \
    BiocManager::install(pkgs, ask = FALSE, update = FALSE, Ncpus = max(1L, parallel::detectCores() - 1L)); \
    actual_deseq2 <- as.character(packageVersion("DESeq2")); \
    expected_deseq2 <- Sys.getenv("DESEQ2_VERSION"); \
    if (!identical(actual_deseq2, expected_deseq2)) { \
      stop(sprintf("DESeq2 version mismatch: expected %s, got %s", expected_deseq2, actual_deseq2)); \
    }; \
    actual_tximport <- as.character(packageVersion("tximport")); \
    expected_tximport <- Sys.getenv("TXIMPORT_VERSION"); \
    if (!identical(actual_tximport, expected_tximport)) { \
      stop(sprintf("tximport version mismatch: expected %s, got %s", expected_tximport, actual_tximport)); \
    } \
  '

COPY docker/slam-deseq2/check_deseq2_versions.R /usr/local/bin/check_deseq2_versions.R
RUN chmod +x /usr/local/bin/check_deseq2_versions.R \
    && Rscript /usr/local/bin/check_deseq2_versions.R

WORKDIR /work
ENTRYPOINT ["Rscript"]
CMD ["/usr/local/bin/check_deseq2_versions.R"]
