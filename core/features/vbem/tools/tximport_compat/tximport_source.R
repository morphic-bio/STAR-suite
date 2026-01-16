function (files, type = c("none", "salmon", "sailfish", "alevin", 
    "piscem", "oarfish", "kallisto", "rsem", "stringtie"), txIn = TRUE, 
    txOut = FALSE, countsFromAbundance = c("no", "scaledTPM", 
        "lengthScaledTPM", "dtuScaledTPM"), tx2gene = NULL, varReduce = FALSE, 
    dropInfReps = FALSE, infRepStat = NULL, ignoreTxVersion = FALSE, 
    ignoreAfterBar = FALSE, geneIdCol, txIdCol, abundanceCol, 
    countsCol, lengthCol, importer = NULL, existenceOptional = FALSE, 
    sparse = FALSE, sparseThreshold = 1, readLength = 75, alevinArgs = NULL) 
{
    infRepImporter <- NULL
    type <- match.arg(type)
    countsFromAbundance <- match.arg(countsFromAbundance)
    if (countsFromAbundance == "dtuScaledTPM") {
        stopifnot(txOut)
        if (is.null(tx2gene)) 
            stop("'dtuScaledTPM' requires 'tx2gene' input")
    }
    if (!existenceOptional) 
        stopifnot(all(file.exists(files)))
    if (!txIn & txOut) 
        stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
    stopifnot(length(files) > 0)
    kallisto.h5 <- basename(files[1]) == "abundance.h5"
    if (type == "kallisto" & !kallisto.h5) {
        message("Note: importing `abundance.h5` is typically faster than `abundance.tsv`")
    }
    if (type == "rsem" & txIn & grepl("genes", files[1])) {
        message("It looks like you are importing RSEM genes.results files, setting txIn=FALSE")
        txIn <- FALSE
    }
    if (!is.null(infRepStat)) {
        if (dropInfReps) 
            stop("infRepStat requires infReps")
        if (type == "alevin") 
            stop("infRepStat does not currently work with alevin input")
        if (sparse) 
            stop("infRepStat does not currently work with sparse output")
    }
    if (type == "alevin") {
        if (is.null(alevinArgs)) {
            alevinArgs <- list(filterBarcodes = FALSE, tierImport = FALSE, 
                forceSlow = FALSE, dropMeanVar = FALSE)
        }
        stopifnot(is(alevinArgs, "list"))
        stopifnot(all(sapply(alevinArgs, is.logical)))
        alevinArgNms <- c("filterBarcodes", "tierImport", "forceSlow", 
            "dropMeanVar")
        stopifnot(all(names(alevinArgs) %in% alevinArgNms))
        filterBarcodes <- if (is.null(alevinArgs$filterBarcodes)) 
            FALSE
        else alevinArgs$filterBarcodes
        tierImport <- if (is.null(alevinArgs$tierImport)) 
            FALSE
        else alevinArgs$tierImport
        forceSlow <- if (is.null(alevinArgs$forceSlow)) 
            FALSE
        else alevinArgs$forceSlow
        dropMeanVar <- if (is.null(alevinArgs$dropMeanVar)) 
            FALSE
        else alevinArgs$dropMeanVar
        if (length(files) > 1) 
            stop("alevin import currently only supports a single experiment")
        if (compareVersion(getAlevinVersion(files), "0.14.0") == 
            -1) {
            stop("use of tximport version >= 1.18 requires alevin version >= 0.14")
        }
        mat <- readAlevin(files, dropInfReps, filterBarcodes, 
            tierImport, forceSlow, dropMeanVar)
        if (!is.list(mat)) {
            txi <- list(abundance = NULL, counts = mat)
        }
        else {
            txi <- list(abundance = NULL, counts = mat$counts, 
                mean = mat$mean, variance = mat$variance)
            if (tierImport) {
                txi$tier <- mat$tier
            }
            if ("infReps" %in% names(mat)) {
                txi$infReps <- mat$infReps
            }
        }
        txi$length <- NULL
        txi$countsFromAbundance = "no"
        return(txi)
    }
    readrStatus <- FALSE
    if (is.null(importer) & !kallisto.h5) {
        if (!requireNamespace("readr", quietly = TRUE)) {
            message("reading in files with read.delim (install 'readr' package for speed up)")
            importer <- read.delim
        }
        else {
            timeZoneMissing <- length(suppressWarnings(OlsonNames())) == 
                0 | is.na(Sys.timezone())
            if (timeZoneMissing) {
                message("reading in files with read.delim ('readr' installed but won't work w/o timezones)")
                importer <- read.delim
            }
            else {
                message("reading in files with read_tsv")
                readrStatus <- TRUE
            }
        }
    }
    if (type %in% c("salmon", "sailfish")) {
        txIdCol <- "Name"
        abundanceCol <- "TPM"
        countsCol <- "NumReads"
        lengthCol <- "EffectiveLength"
        if (readrStatus & is.null(importer)) {
            col.types <- readr::cols(readr::col_character(), 
                readr::col_integer(), readr::col_double(), readr::col_double(), 
                readr::col_double())
            importer <- function(x) readr::read_tsv(x, progress = FALSE, 
                col_types = col.types)
        }
        infRepImporter <- if (dropInfReps) {
            NULL
        }
        else {
            function(x) readInfRepFish(x, type)
        }
    }
    if (type == "piscem") {
        txIdCol <- "target_name"
        lengthCol <- "eelen"
        abundanceCol <- "tpm"
        countsCol <- "ecount"
        if (readrStatus & is.null(importer)) {
            col.types <- readr::cols(readr::col_character(), 
                readr::col_integer(), readr::col_double(), readr::col_double(), 
                readr::col_double())
            importer <- function(x) readr::read_tsv(x, progress = FALSE, 
                col_types = col.types)
        }
        infRepImporter <- if (dropInfReps) {
            NULL
        }
        else {
            readInfRepPiscem
        }
    }
    if (type == "oarfish") {
        txIdCol <- "tname"
        lengthCol <- "len"
        abundanceCol <- "num_reads"
        countsCol <- "num_reads"
        if (readrStatus & is.null(importer)) {
            col.types <- readr::cols(readr::col_character(), 
                readr::col_integer(), readr::col_double())
            importer <- function(x) readr::read_tsv(x, progress = FALSE, 
                col_types = col.types)
        }
        infRepImporter <- if (dropInfReps) {
            NULL
        }
        else {
            readInfRepPiscem
        }
    }
    if (type == "kallisto") {
        txIdCol <- "target_id"
        abundanceCol <- "tpm"
        countsCol <- "est_counts"
        lengthCol <- "eff_length"
        if (kallisto.h5) {
            importer <- read_kallisto_h5
        }
        else if (readrStatus & is.null(importer)) {
            col.types <- readr::cols(readr::col_character(), 
                readr::col_integer(), readr::col_double(), readr::col_double(), 
                readr::col_double())
            importer <- function(x) readr::read_tsv(x, progress = FALSE, 
                col_types = col.types)
        }
        infRepImporter <- if (dropInfReps) {
            NULL
        }
        else {
            readInfRepKallisto
        }
    }
    if (type == "rsem") {
        if (txIn) {
            txIdCol <- "transcript_id"
            abundanceCol <- "TPM"
            countsCol <- "expected_count"
            lengthCol <- "effective_length"
            if (readrStatus & is.null(importer)) {
                col.types <- readr::cols(readr::col_character(), 
                  readr::col_character(), readr::col_integer(), 
                  readr::col_double(), readr::col_double(), readr::col_double(), 
                  readr::col_double(), readr::col_double())
                importer <- function(x) readr::read_tsv(x, progress = FALSE, 
                  col_types = col.types)
            }
        }
        else {
            geneIdCol <- "gene_id"
            abundanceCol <- "TPM"
            countsCol <- "expected_count"
            lengthCol <- "effective_length"
            if (readrStatus & is.null(importer)) {
                col.types <- readr::cols(readr::col_character(), 
                  readr::col_character(), readr::col_double(), 
                  readr::col_double(), readr::col_double(), readr::col_double(), 
                  readr::col_double())
                importer <- function(x) readr::read_tsv(x, progress = FALSE, 
                  col_types = col.types)
            }
        }
    }
    if (type == c("stringtie")) {
        txIdCol <- "t_name"
        geneIdCol <- "gene_name"
        abundanceCol <- "FPKM"
        countsCol <- "cov"
        lengthCol <- "length"
        if (readrStatus & is.null(importer)) {
            col.types <- readr::cols(readr::col_character(), 
                readr::col_character(), readr::col_character(), 
                readr::col_integer(), readr::col_integer(), readr::col_character(), 
                readr::col_integer(), readr::col_integer(), readr::col_character(), 
                readr::col_character(), readr::col_double(), 
                readr::col_double())
            importer <- function(x) readr::read_tsv(x, progress = FALSE, 
                col_types = col.types)
        }
    }
    infRepType <- "none"
    if (type %in% c("salmon", "sailfish", "piscem", "oarfish", 
        "kallisto") & !dropInfReps) {
        infRepType <- if (varReduce & txOut) {
            "var"
        }
        else {
            "full"
        }
    }
    if (!txIn) {
        txi <- computeRsemGeneLevel(files, importer, geneIdCol, 
            abundanceCol, countsCol, lengthCol, countsFromAbundance)
        message("")
        return(txi)
    }
    if (is.null(tx2gene) & !txOut) {
        summarizeFail()
    }
    repInfo <- NULL
    if (infRepType != "none") {
        repInfo <- if (type %in% c("piscem", "oarfish")) {
            infRepImporter(files[1])
        }
        else {
            infRepImporter(dirname(files[1]))
        }
        if (is.null(repInfo)) {
            infRepType <- "none"
        }
    }
    if (sparse) {
        if (!requireNamespace("Matrix", quietly = TRUE)) {
            stop("sparse import requires core R package `Matrix`")
        }
        message("importing sparsely, only counts and abundances returned, support limited to\ntxOut=TRUE, CFA either 'no' or 'scaledTPM', and no inferential replicates")
        stopifnot(txOut)
        stopifnot(infRepType == "none")
        stopifnot(countsFromAbundance %in% c("no", "scaledTPM"))
    }
    for (i in seq_along(files)) {
        message(i, " ", appendLF = FALSE)
        raw <- as.data.frame(importer(files[i]))
        if (infRepType != "none") {
            repInfo <- if (type %in% c("piscem", "oarfish")) {
                infRepImporter(files[i])
            }
            else {
                infRepImporter(dirname(files[i]))
            }
        }
        else {
            repInfo <- NULL
        }
        stopifnot(all(c(abundanceCol, countsCol, lengthCol) %in% 
            names(raw)))
        if (i == 1) {
            txId <- raw[[txIdCol]]
        }
        else {
            stopifnot(all(txId == raw[[txIdCol]]))
        }
        if (!sparse) {
            if (i == 1) {
                mat <- matrix(nrow = nrow(raw), ncol = length(files))
                rownames(mat) <- raw[[txIdCol]]
                colnames(mat) <- names(files)
                abundanceMatTx <- mat
                countsMatTx <- mat
                lengthMatTx <- mat
                if (infRepType == "var") {
                  varMatTx <- mat
                }
                else if (infRepType == "full") {
                  infRepMatTx <- list()
                }
            }
            abundanceMatTx[, i] <- raw[[abundanceCol]]
            countsMatTx[, i] <- raw[[countsCol]]
            lengthMatTx[, i] <- raw[[lengthCol]]
            if (infRepType == "var") {
                varMatTx[, i] <- repInfo$vars
            }
            else if (infRepType == "full") {
                infRepMatTx[[i]] <- repInfo$reps
            }
            if (!is.null(infRepStat)) {
                countsMatTx[, i] <- infRepStat(repInfo$reps)
                tpm <- countsMatTx[, i]/lengthMatTx[, i]
                abundanceMatTx[, i] <- tpm * 1e+06/sum(tpm)
            }
        }
        else {
            if (i == 1) {
                txId <- raw[[txIdCol]]
                countsListI <- list()
                countsListX <- list()
                abundanceListX <- list()
                numNonzero <- c()
            }
            stopifnot(all(txId == raw[[txIdCol]]))
            sparse.idx <- which(raw[[countsCol]] >= sparseThreshold)
            countsListI <- c(countsListI, sparse.idx)
            countsListX <- c(countsListX, raw[[countsCol]][sparse.idx])
            numNonzero <- c(numNonzero, length(sparse.idx))
            if (countsFromAbundance == "scaledTPM") {
                abundanceListX <- c(abundanceListX, raw[[abundanceCol]][sparse.idx])
            }
        }
    }
    if (sparse) {
        countsMatTx <- Matrix::sparseMatrix(i = unlist(countsListI), 
            j = rep(seq_along(numNonzero), numNonzero), x = unlist(countsListX), 
            dims = c(length(txId), length(files)), dimnames = list(txId, 
                names(files)))
        if (countsFromAbundance == "scaledTPM") {
            abundanceMatTx <- Matrix::sparseMatrix(i = unlist(countsListI), 
                j = rep(seq_along(numNonzero), numNonzero), x = unlist(abundanceListX), 
                dims = c(length(txId), length(files)), dimnames = list(txId, 
                  names(files)))
        }
        else {
            abundanceMatTx <- NULL
        }
        lengthMatTx <- NULL
    }
    if (infRepType == "full") {
        if (length(infRepMatTx) != length(files)) {
            stop("Note: not all samples contain inferential replicates.\n  tximport can only import data when either all or no samples\n  contain inferential replicates. Instead first subset to the\n  set of samples that all contain inferential replicates.")
        }
        names(infRepMatTx) <- names(files)
    }
    message("")
    if (infRepType == "none") {
        txi <- list(abundance = abundanceMatTx, counts = countsMatTx, 
            length = lengthMatTx, countsFromAbundance = countsFromAbundance)
    }
    else if (infRepType == "var") {
        txi <- list(abundance = abundanceMatTx, counts = countsMatTx, 
            variance = varMatTx, length = lengthMatTx, countsFromAbundance = countsFromAbundance)
    }
    else if (infRepType == "full") {
        txi <- list(abundance = abundanceMatTx, counts = countsMatTx, 
            infReps = infRepMatTx, length = lengthMatTx, countsFromAbundance = countsFromAbundance)
    }
    if (type == "stringtie") {
        txi$counts <- txi$counts * txi$length/readLength
    }
    if (type == "rsem") {
        txi$length[txi$length < 1] <- 1
    }
    if (txOut) {
        if (countsFromAbundance != "no") {
            length4CFA <- txi$length
            if (countsFromAbundance == "dtuScaledTPM") {
                length4CFA <- medianLengthOverIsoform(length4CFA, 
                  tx2gene, ignoreTxVersion, ignoreAfterBar)
                countsFromAbundance <- "lengthScaledTPM"
            }
            txi$counts <- makeCountsFromAbundance(countsMat = txi$counts, 
                abundanceMat = txi$abundance, lengthMat = length4CFA, 
                countsFromAbundance = countsFromAbundance)
        }
        if (ignoreAfterBar) {
            for (matNm in c("counts", "abundance", "length")) {
                rowNms <- rownames(txi[[matNm]])
                rownames(txi[[matNm]]) <- sub("\\|.*", "", rowNms)
            }
        }
        return(txi)
    }
    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, varReduce, ignoreTxVersion, 
        ignoreAfterBar, countsFromAbundance)
    return(txiGene)
}
<bytecode: 0x623771bd8d98>
<environment: namespace:tximport>
