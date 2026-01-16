# Y-Chromosome BAM Split - Technical Documentation

## Overview

**Note**: This feature was developed for **Morphic requirements for KOLF cell lines**. It is **not connected to the Flex pipeline** and is a general-purpose feature usable with any STAR workflow (bulk RNA-seq, single-cell, etc.).

The Y-chromosome BAM split feature allows STAR to emit two additional BAM files during alignment:
- `<out>_Y.bam`: Contains all alignments for reads with any Y-chromosome alignment
- `<out>_noY.bam`: Contains all alignments for reads with no Y-chromosome alignments

This feature is useful for sex-specific analyses, separating male/female samples in mixed datasets, and meeting Morphic's requirements for KOLF cell line processing.

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              STAR Alignment                                  │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────┐     ┌──────────────────┐     ┌─────────────────────┐  │
│  │  Genome Load    │────▶│  Y-Contig        │────▶│  yTids set          │  │
│  │  (chrInfoLoad)  │     │  Detection       │     │  (chromosome IDs)   │  │
│  └─────────────────┘     └──────────────────┘     └─────────────────────┘  │
│                                                             │               │
│                                                             ▼               │
│  ┌─────────────────┐     ┌──────────────────┐     ┌─────────────────────┐  │
│  │  Read Alignment │────▶│  outputAlignments│────▶│  hasYAlignment_     │  │
│  │  (ReadAlign)    │     │  (per-read)      │     │  flag computed      │  │
│  └─────────────────┘     └──────────────────┘     └─────────────────────┘  │
│                                                             │               │
│                          ┌──────────────────────────────────┘               │
│                          ▼                                                  │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │                        BAM Output Path                                 │ │
│  │  ┌─────────────────────────┐    ┌─────────────────────────────────┐   │ │
│  │  │   Unsorted Path         │    │   Coordinate-Sorted Path        │   │ │
│  │  │   ─────────────         │    │   ─────────────────────         │   │ │
│  │  │   hasY → route to       │    │   hasY bit encoded in iRead     │   │ │
│  │  │   bgzfBAM_Y or          │    │   (bit 63 of uint64)            │   │ │
│  │  │   bgzfBAM_noY           │    │   │                             │   │ │
│  │  └─────────────────────────┘    │   ▼                             │   │ │
│  │                                 │   BAMbinSortByCoordinate         │   │ │
│  │                                 │   extracts bit, routes to        │   │ │
│  │                                 │   bgzfBin_Y or bgzfBin_noY       │   │ │
│  │                                 │   │                             │   │ │
│  │                                 │   ▼                             │   │ │
│  │                                 │   bamSortByCoordinate            │   │ │
│  │                                 │   concatenates temp files        │   │ │
│  │                                 └─────────────────────────────────┘   │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │  Output Files                                                         │ │
│  │  ├── Aligned.out_Y.bam          (unsorted Y)                          │ │
│  │  ├── Aligned.out_noY.bam        (unsorted noY)                        │ │
│  │  ├── Aligned.sortedByCoord.out_Y.bam    (sorted Y)                    │ │
│  │  └── Aligned.sortedByCoord.out_noY.bam  (sorted noY)                  │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Command-Line Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--emitNoYBAM` | `no` | Enable Y-chromosome BAM splitting (`yes`/`no`) |
| `--emitYReadNames` | `no` | Emit list of read names with any Y-chromosome alignment (one per line) |
| `--emitYNoYFastq` | `no` | Emit Y/noY FASTQ files directly during alignment (`yes`/`no`) |
| `--emitYNoYFastqCompression` | `gz` | Compression for FASTQ output (`gz`/`none`) |
| `--YFastqOutputPrefix` | - | Output prefix for Y FASTQ files (default: derived from input name, inserting `_Y` before the last `_R1`/`_R2`; if missing, falls back to `Y_reads.mateN`) |
| `--noYFastqOutputPrefix` | - | Output prefix for noY FASTQ files (default: derived from input name, inserting `_noY` before the last `_R1`/`_R2`; if missing, falls back to `noY_reads.mateN`) |
| `--keepBAM` | `no` | Keep primary BAM when splitting is enabled |
| `--noYOutput` | - | Override default path for noY BAM |
| `--YOutput` | - | Override default path for Y BAM |
| `--YReadNamesOutput` | - | Override output path for Y read names list (default: `<outFileNamePrefix>Aligned.out_Y.names.txt`) |

**Names-only mode**: You can emit Y-read names without producing Y/noY BAMs by using `--emitYReadNames yes` alone.

**FASTQ-only mode**: You can emit Y/noY FASTQ files without producing BAMs by using `--emitYNoYFastq yes` with `--outSAMtype None`.

**FASTQ edge cases**:
- If the reference has no Y contigs, Y FASTQs are empty and a warning is logged.
- FASTA inputs emit `.fa(.gz)` outputs with `>` headers and no `+`/quality lines.
- Multiple input files per mate derive output names from the first file for each mate.
- `--emitYNoYFastqCompression none` writes uncompressed outputs.
- Unmapped reads are routed to noY.

## Files Modified

### Header Files

| File | Changes |
|------|---------|
| `Parameters.h` | Added `emitNoYBAM`, `keepBAM` structs; `outBAMfileYName`, `outBAMfileNoYName` strings |
| `Genome.h` | Added `std::unordered_set<int> yTids` for Y-contig IDs |
| `ReadAlign.h` | Added `bool hasYAlignment_` flag (public) |
| `BAMoutput.h` | Added `BGZF *bgzfBAM_Y`, `bgzfBAM_noY`; `bool suppressPrimary_` |
| `InOutStreams.h` | Added `BGZF *outBAMfileY`, `*outBAMfileNoY` |

### Source Files

| File | Changes |
|------|---------|
| `parametersDefault` | Parameter definitions and help text |
| `Parameters.cpp` | Parameter parsing, path derivation, BGZF handle initialization |
| `Genome.cpp` | Y-contig detection in `chrInfoLoad()` |
| `ReadAlign_outputAlignments.cpp` | Mode-aware `hasYAlignment_` computation |
| `BAMoutput.cpp` | Routing logic for unsorted and coord-sorted paths |
| `BAMbinSortByCoordinate.cpp` | Y-bit extraction during merge, routing to temp files |
| `bamSortByCoordinate.cpp` | Final concatenation of Y/noY temp files |
| `samHeaders.cpp` | Header writing to Y/noY files |
| `STAR.cpp` | Cleanup and file closing |
| `InOutStreams.cpp` | Handle initialization and cleanup |
| `ChimericAlign_chimericBAMoutput.cpp` | Pass `hasYAlignment_` flag |
| `serviceFuns.cpp` | Added `funCompareArraysIgnoreYBit` comparison function |

## Implementation Details

### 1. Y-Contig Detection (Genome.cpp)

The `chrInfoLoad()` function populates `yTids` with chromosome IDs matching Y patterns:

```cpp
// In Genome::chrInfoLoad()
yTids.clear();
for (uint i = 0; i < nChrReal; i++) {
    string chrNameLower = chrName[i];
    transform(chrNameLower.begin(), chrNameLower.end(), chrNameLower.begin(), ::tolower);
    
    if (chrNameLower == "y" || chrNameLower == "chry" ||
        chrNameLower.rfind("chry_", 0) == 0 ||
        chrNameLower.rfind("y_", 0) == 0) {
        yTids.insert(i);
    }
}
```

**Patterns detected** (case-insensitive):
- `chrY`, `Y` - Standard Y chromosome
- `chrY_random`, `chrY_alt` - Patches and alternate assemblies
- Any contig starting with `chrY_` or `Y_`

### 2. Per-Read Y-Decision (ReadAlign_outputAlignments.cpp)

The `hasYAlignment_` flag is computed once per read, covering all alignments:

```cpp
// Mode detection
bool isSingleCellOrFlex = (P.pSolo.type != 0) || P.pSolo.flexMode;
bool hasMates = (P.readNmates == 2) && !isSingleCellOrFlex;

// For mapped reads
if (unmapType < 0) {
    for (uint iTr = 0; iTr < nTrCheck && !hasYAlignment_; iTr++) {
        if (mapGen.yTids.count(trCheck[iTr]->Chr) > 0) {
            hasYAlignment_ = true;
        }
    }
}

// For unmapped reads with mapped mates (bulk paired-end only)
if (unmapType >= 0 && hasMates) {
    for (uint iTr = 0; iTr < nTrCheck && !hasYAlignment_; iTr++) {
        if (mapGen.yTids.count(trCheck[iTr]->Chr) > 0) {
            hasYAlignment_ = true;
        }
    }
}
```

**Mode-specific behavior:**

| Mode | Detection Logic |
|------|-----------------|
| Single-cell / Flex | Only current read's alignments (R1/R2 are not mates) |
| Bulk paired-end | All transcripts (covers both mates' alignments) |
| Bulk single-end | Only current read's alignments |

### 3. Unsorted BAM Routing (BAMoutput.cpp)

```cpp
void BAMoutput::unsortedOneAlign(char *bamIn, uint bamSize, bool hasY) {
    if (P.emitNoYBAM.yes) {
        BGZF *targetHandle = hasY ? bgzfBAM_Y : bgzfBAM_noY;
        bgzf_write(targetHandle, bamIn, bamSize);
    }
    if (!suppressPrimary_) {
        bgzf_write(bgzfBAM, bamIn, bamSize);
    }
}
```

### 4. Coordinate-Sorted Path

#### 4.1 Y-Bit Encoding (BAMoutput.cpp)

The Y-flag is encoded as bit 63 of the `iRead` uint64:

```cpp
void BAMoutput::coordOneAlign(char *bamIn, uint bamSize, uint64 iRead, bool hasY) {
    // Encode hasY in bit 63
    if (hasY) {
        iRead |= (1ULL << 63);
    }
    // ... write to temp file
}
```

#### 4.2 Y-Bit Extraction (BAMbinSortByCoordinate.cpp)

During merge, the Y-bit is extracted and records are routed:

```cpp
// Custom comparison ignoring bit 63
inline int funCompareArraysIgnoreYBit(const void *a, const void *b) {
    uint64 *pa = (uint64*)a;
    uint64 *pb = (uint64*)b;
    uint64 a0 = pa[0] & ~(1ULL << 63);  // Mask out Y-bit
    uint64 b0 = pb[0] & ~(1ULL << 63);
    // ... comparison logic
}

// In BAMbinSortByCoordinate()
qsort((void*)startPos, binN, sizeof(uint64)*3, funCompareArraysIgnoreYBit);

for (uint64 iRead = 0; iRead < binN; iRead++) {
    bool hasY = (startPos[iRead * 3] >> 63) & 1;
    BGZF *targetHandle = hasY ? bgzfBin_Y : bgzfBin_noY;
    // ... write record
}
```

#### 4.3 Final Concatenation (bamSortByCoordinate.cpp)

Temporary files are collected and concatenated:

```cpp
// Collect Y and noY temp files
vector<string> yFiles, noYFiles;
for (int ib = 0; ib < nBins; ib++) {
    string yFile = tmpDir + "/b" + to_string(ib) + "_Y";
    string noYFile = tmpDir + "/b" + to_string(ib) + "_noY";
    if (stat(yFile.c_str(), &sb) == 0) yFiles.push_back(yFile);
    if (stat(noYFile.c_str(), &sb) == 0) noYFiles.push_back(noYFile);
}

// Concatenate using bam_cat
if (!yFiles.empty()) {
    bam_cat(yFiles, P.outBAMfileYName);
}
if (!noYFiles.empty()) {
    bam_cat(noYFiles, P.outBAMfileNoYName);
}
```

### 5. Path Derivation (Parameters.cpp)

Output paths are derived from the primary BAM path:

```cpp
if (emitNoYBAM.yes) {
    string baseName = outBAMfileCoordName;  // or outBAMfileUnsortedName
    size_t dotPos = baseName.rfind(".bam");
    if (dotPos != string::npos) {
        outBAMfileYName = baseName.substr(0, dotPos) + "_Y.bam";
        outBAMfileNoYName = baseName.substr(0, dotPos) + "_noY.bam";
    } else {
        outBAMfileYName = baseName + "_Y.bam";
        outBAMfileNoYName = baseName + "_noY.bam";
    }
}
```

## Testing

### Test Scripts

| Script | Purpose |
|--------|---------|
| `tests/run_ychrom_bam_split_test.sh` | Basic E2E test with synthetic data |
| `tests/run_ychrom_flex_validation_test.sh` | Flex mode validation |
| `tests/run_ychrom_bulk_pe_test.sh` | Bulk paired-end validation |
| `tests/run_ychrom_regression_test.sh` | Regression against STAR Solo |

### Validation Checks

1. **Count consistency**: Y + noY = baseline total
2. **noY exclusivity**: Zero chrY reads in `_noY.bam`
3. **Y exclusivity**: All reads in `_Y.bam` have chrY alignments
4. **Mate consistency**: No read names appear in both Y and noY files
5. **Regression**: Baseline matches upstream STAR Solo

### Test Reports

- `tests/TEST_REPORT_Y_SPLIT.md` - Basic test results
- `tests/TEST_REPORT_Y_SPLIT_FLEX.md` - Flex mode results
- `tests/TEST_REPORT_Y_SPLIT_BULK_PE.md` - Bulk paired-end results
- `tests/TEST_REPORT_Y_SPLIT_REGRESSION.md` - Regression results

## Usage Examples

### Single-Cell / Flex Mode

```bash
STAR \
  --genomeDir /path/to/flex_reference \
  --readFilesIn R2.fastq.gz R1.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --flex yes \
  --outSAMtype BAM SortedByCoordinate \
  --emitNoYBAM yes \
  --outFileNamePrefix output/
```

### Bulk Paired-End

```bash
STAR \
  --genomeDir /path/to/reference \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --emitNoYBAM yes \
  --outFileNamePrefix output/
```

### Keep Primary BAM

```bash
STAR \
  ... \
  --emitNoYBAM yes \
  --keepBAM yes \
  --outFileNamePrefix output/
```

### Custom Output Paths

```bash
STAR \
  ... \
  --emitNoYBAM yes \
  --YOutput /path/to/male_reads.bam \
  --noYOutput /path/to/female_reads.bam \
  --outFileNamePrefix output/
```

## Known Limitations

1. **Memory**: Y-bit encoding uses bit 63 of `iRead`, reducing available bits for read indexing (still supports up to 2^63 reads)
2. **Chimeric alignments**: Chimeric reads are routed based on `hasYAlignment_` from the main alignment
3. **Unmapped reads**: Unmapped reads with no mate information route to `_noY.bam` by default

---

# remove_y_reads - FASTQ Splitter Tool

## Overview

A companion C tool that splits FASTQ files based on a Y-only BAM file. Given the `_Y.bam` output from STAR's Y-chromosome split feature, this tool partitions original FASTQ files into Y/noY sets while preserving read order.

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         remove_y_reads Tool                                  │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  Phase 1: Build Y Hash Set                                           │   │
│  │  ┌─────────────────┐     ┌──────────────────┐     ┌──────────────┐  │   │
│  │  │  Load Y BAM     │────▶│  Normalize qname │────▶│  khash set   │  │   │
│  │  │  (htslib)       │     │  (strip @,/1,/2) │     │  (hash+len)  │  │   │
│  │  └─────────────────┘     └──────────────────┘     └──────────────┘  │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                      │                                      │
│                                      ▼                                      │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  Phase 2: Split FASTQs (file-level parallelism)                      │   │
│  │                                                                       │   │
│  │  ┌───────────┐  ┌───────────┐  ┌───────────┐                        │   │
│  │  │ Worker 1  │  │ Worker 2  │  │ Worker N  │  (bounded by semaphore)│   │
│  │  │  kseq.h   │  │  kseq.h   │  │  kseq.h   │                        │   │
│  │  └─────┬─────┘  └─────┬─────┘  └─────┬─────┘                        │   │
│  │        │              │              │                               │   │
│  │        ▼              ▼              ▼                               │   │
│  │  ┌───────────┐  ┌───────────┐  ┌───────────┐                        │   │
│  │  │ _Y.fastq  │  │ _Y.fastq  │  │ _Y.fastq  │                        │   │
│  │  │ _noY.fastq│  │ _noY.fastq│  │ _noY.fastq│                        │   │
│  │  └───────────┘  └───────────┘  └───────────┘                        │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Key Files

| File | Purpose |
|------|---------|
| [tools/remove_y_reads/remove_y_reads.c](../tools/remove_y_reads/remove_y_reads.c) | Main implementation (~400 lines) |
| [tools/remove_y_reads/Makefile](../tools/remove_y_reads/Makefile) | Build configuration |
| [tests/run_remove_y_reads_test.sh](../tests/run_remove_y_reads_test.sh) | Basic test script |
| [tests/run_y_removal_comprehensive_test.sh](../tests/run_y_removal_comprehensive_test.sh) | Comprehensive test suite |

## Command-Line Interface

```
Usage: remove_y_reads -y <Y.bam> [-o <out_dir>] [--threads N] [--gzip-level N] fastq1 [fastq2 ...]

Options:
  -y, --ybam FILE     Y-only BAM file (required)
  -o, --outdir DIR    Output directory (default: alongside input)
  -t, --threads N     Number of threads (default: 1)
  -z, --gzip-level N  Gzip compression level 1-9 (default: 6)
  -h, --help          Show this help
```

## Implementation Details

### Dependencies

- **htslib**: BAM reading via `sam_open()`, `sam_read1()`, `bam_get_qname()`
- **kseq.h**: Robust FASTQ parsing with dynamic allocation (no fixed buffers)
- **khash.h**: Efficient hash table for Y read name lookup
- **zlib**: Gzip compression/decompression with `gzbuffer()` for performance
- **pthreads**: File-level parallelism with semaphore-bounded concurrency

### Hash Collision Protection

Uses dual-hash approach to prevent misrouting:

```c
// Primary hash: FNV-1a 64-bit
static inline uint64_t hash_qname(const char *s, int len) {
    uint64_t h = 14695981039346656037ULL;
    for (int i = 0; i < len; i++) {
        h ^= (uint64_t)(unsigned char)s[i];
        h *= 1099511628211ULL;
    }
    return h;
}

// Secondary hash: djb2 64-bit
static inline uint64_t hash_qname2(const char *s, int len) {
    uint64_t h = 5381ULL;
    for (int i = 0; i < len; i++) {
        h = ((h << 5) + h) + (uint64_t)(unsigned char)s[i];
    }
    return h;
}

// Storage: hash1 (key) -> {len, hash2} (value)
typedef struct {
    uint32_t len;
    uint64_t hash2;
} qname_info_t;
```

Collision check requires: `hash1 == stored_hash1 && hash2 == stored_hash2 && length == stored_length`

### Name Normalization

Matches STAR's BAM qname format exactly:

1. Strip leading `@` (FASTQ format)
2. Truncate at first whitespace, tab, or comment
3. Strip trailing `/1` or `/2` (mate suffix)

```c
static int normalize_qname(const char *raw, int raw_len, char *out, int max_len) {
    const char *start = raw;
    int len = raw_len;
    
    // Skip leading '@'
    if (len > 0 && start[0] == '@') { start++; len--; }
    
    // Find end at whitespace
    int end = 0;
    while (end < len && start[end] != ' ' && start[end] != '\t' && 
           start[end] != '\n' && start[end] != '\r') { end++; }
    len = end;
    
    // Strip /1 or /2
    if (len >= 2 && start[len-2] == '/' && 
        (start[len-1] == '1' || start[len-1] == '2')) { len -= 2; }
    
    memcpy(out, start, len);
    out[len] = '\0';
    return len;
}
```

### Thread Management

Batched thread creation to avoid many parked threads:

```c
// Only creates up to num_threads at a time
while (next_file < num_files || active_count > 0) {
    // Start new threads up to the limit
    while (active_count < num_threads && next_file < num_files) {
        pthread_create(&threads[slot], NULL, worker_thread, &workers[next_file]);
        active_count++;
        next_file++;
    }
    
    // Join completed thread, start new one
    if (active_count > 0) {
        pthread_join(threads[0], NULL);
        active_count--;
        // Shift remaining threads
    }
}
```

### FASTQ Parsing with kseq.h

Handles arbitrarily long reads without buffer overflow:

```c
KSEQ_INIT(gzFile, gzread)

void process_fastq(const char *inpath, gzFile out_y, gzFile out_noy,
                   khash_t(qname_set) *y_set, ...) {
    gzFile fp = gzopen(inpath, "r");
    gzbuffer(fp, 128 * 1024);  // MANDATORY: large buffer for performance
    
    kseq_t *seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        // seq->name.s is dynamically allocated
        int norm_len = normalize_qname(seq->name.s, seq->name.l, norm, sizeof(norm));
        uint64_t h = hash_qname(norm, norm_len);
        uint64_t h2 = hash_qname2(norm, norm_len);
        
        khint_t k = kh_get(qname_set, y_set, h);
        bool is_y = (k != kh_end(y_set) && 
                     kh_val(y_set, k).len == norm_len &&
                     kh_val(y_set, k).hash2 == h2);
        
        gzFile out = is_y ? out_y : out_noy;
        gzprintf(out, "@%s", seq->name.s);
        if (seq->comment.l) gzprintf(out, " %s", seq->comment.s);
        gzprintf(out, "\n%s\n+\n%s\n", seq->seq.s, seq->qual.s);
    }
    kseq_destroy(seq);
    gzclose(fp);
}
```

## Usage Examples

### Basic Usage

```bash
# Split FASTQs based on Y BAM
./tools/remove_y_reads/remove_y_reads \
    -y Aligned.sortedByCoord.out_Y.bam \
    sample_R1.fastq.gz sample_R2.fastq.gz
```

Output: `sample_R1_Y.fastq.gz`, `sample_R1_noY.fastq.gz`, `sample_R2_Y.fastq.gz`, `sample_R2_noY.fastq.gz`

### With Output Directory and Threading

```bash
./tools/remove_y_reads/remove_y_reads \
    -y Aligned.sortedByCoord.out_Y.bam \
    --threads 4 \
    --gzip-level 6 \
    -o output_dir \
    sample_R1.fastq.gz sample_R2.fastq.gz
```

### Complete Workflow

```bash
# Step 1: Run STAR with Y/noY split
STAR \
  --genomeDir /path/to/reference \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --emitNoYBAM yes \
  --outFileNamePrefix output/

# Step 2: Split original FASTQs based on Y BAM
./tools/remove_y_reads/remove_y_reads \
    -y output/Aligned.sortedByCoord.out_Y.bam \
    --threads 4 \
    -o output/fastq_split \
    R1.fastq.gz R2.fastq.gz

# Result:
# output/fastq_split/R1_Y.fastq.gz     - Y-chromosome reads from R1
# output/fastq_split/R1_noY.fastq.gz   - Non-Y reads from R1
# output/fastq_split/R2_Y.fastq.gz     - Y-chromosome reads from R2
# output/fastq_split/R2_noY.fastq.gz   - Non-Y reads from R2
```

## Testing

### Test Scripts

| Script | Purpose |
|--------|---------|
| `tests/run_remove_y_reads_test.sh` | Basic self-contained test |
| `tests/run_y_removal_comprehensive_test.sh` | Comprehensive validation suite |

### Validation Checks

1. **Output existence**: `_Y.fastq.gz` and `_noY.fastq.gz` created for each input
2. **Count consistency**: Y + noY = original total reads
3. **Y hash membership**: All reads in `_Y.fastq` exist in Y BAM hash set
4. **noY hash membership**: No reads in `_noY.fastq` exist in Y BAM hash set
5. **Order preservation**: Relative read order maintained in outputs
6. **Compression round-trip**: Gzip handled correctly

### Test Report

Test results are written to `tests/TEST_REPORT_REMOVE_Y_FASTQ.md`.

### Running Tests

```bash
# Basic test
./tests/run_remove_y_reads_test.sh

# Comprehensive test (single-threaded, multithreaded, multiple files)
./tests/run_y_removal_comprehensive_test.sh
```

## Building

```bash
cd tools/remove_y_reads
make
```

Requires htslib to be built first:

```bash
cd source
make htslib/libhts.a
```

## Changelog

- **2025-12-12**: Initial implementation
  - Implemented Y BAM loader with dual-hash collision protection
  - Added kseq.h-based FASTQ parsing (no fixed buffers)
  - Implemented file-level threading with semaphore-bounded concurrency
  - Added comprehensive test suite with order preservation and hash membership validation

## References

- [README_flex.md](../README_flex.md) - User-facing documentation
- [plans/remove_y_tool_plan.md](../plans/remove_y_tool_plan.md) - Implementation plan
- [tests/TEST_REPORT_REMOVE_Y_FASTQ.md](../tests/TEST_REPORT_REMOVE_Y_FASTQ.md) - Test report
