# STAR-SLAM Binary Dump Format

This document specifies the **binary dump** and **weight sidecar** formats written by
STAR-SLAM when `--slamDumpBinary 1` and `--slamDumpWeights 1` are enabled.

The canonical implementation is in `slam/source/SlamDump.{h,cpp}`; this spec mirrors
those definitions.

## Conventions

- **Endianness:** little-endian (native on x86_64).
- **Alignment:** packed; fields are written back-to-back with no padding.
- **Booleans:** serialized as 1 byte (0 or 1), matching C++ `bool` writes.
- **Strings:** `u32 length` followed by raw bytes (no null terminator).
- **Base encoding:** STAR internal base encoding (`A=0, C=1, G=2, T=3`).

## Files

- Dump: `<sample>_slam_dump.bin`
- Weights: `<sample>_slam_weights.bin`

## Dump File (`SLAMDUMP`)

### Header Layout (in order)

```
char[8]   magic        = "SLAMDUMP"
u32       version      = 1
u32       flags        = weightMode in bits [1:0]
u32       nGenes
u32       nChrom
u64       nReads
f64       errorRate
f64       convRate

// Gene metadata (repeat nGenes)
string    geneId
string    geneName

// Chromosome metadata (repeat nChrom)
string    chrName
u64       chrStart
```

### Header Flags (bitwise)

- `flags & 0x3` → `weightMode`
  - `0` = **Alignments** (1 / nTr)
  - `1` = **Uniform** (1.0)
- All other bits are reserved and currently 0.

### Read Record Layout (repeat `nReads`)

```
string    readName
u32       readLength0          // mate 1 length
u32       readLength1          // mate 2 length (0 for SE)
u8        isMinus              // strand
u8        oppositeStrand
u8        isIntronic
u32       fileIndex            // source file index (batch/trimScope=per-file)
f64       weight               // assignment weight

u32       geneCount
u32[geneCount] geneIds

u32       posCount
for each position:
  u32     readPos              // concatenated read position
  u64     genomicPos           // STAR internal genomic coordinate
  u8      refBase              // 0..3
  u8      readBase             // 0..3
  u8      qual                 // Phred
  u8      secondMate           // 1 if mate2
  u8      overlap              // 1 if in PE overlap
```

Notes:
- `readPos` is in **concatenated** coordinates; `secondMate=1` indicates mate2.
- `genomicPos` is the internal STAR coordinate used during alignment (0-based in
  practice; do not assume BAM coordinate semantics).

## Weight Sidecar (`SLAMWGT1`)

This file stores read-level weights keyed by a stable hash of the buffered read.

### Header Layout

```
char[8]   magic        = "SLAMWGT1"
u32       version      = 1
u32       flags        = keyed + ordered + weightMode
u64       nReads
u32       weightMode   // same meaning as dump flags
```

### Flags (bitwise)

- bit 0 (`0x1`) = keyed records (always on)
- bit 1 (`0x2`) = ordered records (always on)
- bits [3:2] = `weightMode` (`0` Alignments, `1` Uniform)

### Record Layout (repeat `nReads`)

```
u64   key.h1
u64   key.h2
f64   weight
```

The key is an FNV64-based hash over the dump read contents (see
`computeSlamWeightKey` in `slam/source/SlamDump.cpp`).

## Versioning

- `version` currently fixed at 1 for both files.
- Backward-incompatible changes must bump the version and update readers.

