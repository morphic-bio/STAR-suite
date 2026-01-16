# STAR-Flex Changelog

This document tracks changes specific to the STAR-Flex fork. For upstream STAR changes, see [CHANGES.md](CHANGES.md).

## 2025-12-18: Y-Chromosome BAM Split (Morphic/KOLF)

**Feature**: Y-chromosome BAM splitting for **Morphic requirements for KOLF cell lines**. This is a general-purpose feature **not connected to the Flex pipeline**, usable with any STAR workflow.

**Parameters**:
- `--emitNoYBAM` (yes/no): Enable Y-chromosome BAM splitting
- `--keepBAM` (yes/no): Keep primary BAM when splitting is enabled
- `--noYOutput`, `--YOutput`: Optional custom output paths

**Use Case**: Developed specifically for Morphic's KOLF cell line processing requirements, but works with any bulk RNA-seq or single-cell workflow.

**Documentation**: See [docs/Y_CHROMOSOME_BAM_SPLIT.md](docs/Y_CHROMOSOME_BAM_SPLIT.md) for technical details.

## 2025-12-18: Cutadapt-Style Trimming (Bulk RNA-seq)

**Feature**: Added cutadapt-parity trimming for **bulk RNA-seq** workflows with perfect match to Trim Galore/cutadapt v5.1. This is a general-purpose trimming feature usable with any STAR workflow, not specific to the Flex pipeline.

**New Parameters**:
- `--trimCutadapt` (Yes/-): Enable cutadapt-style trimming
- `--trimCutadaptQuality` (default: 20): Quality threshold for 3' trimming
- `--trimCutadaptMinLength` (default: 20): Minimum read length after trimming
- `--trimCutadaptAdapter` (default: -): Custom adapter sequences or TruSeq defaults

**Algorithm**: Exact port of cutadapt v5.1 semiglobal alignment algorithm (`_align.pyx`, `Aligner.locate()` method) with edit distance (Levenshtein), scoring (MATCH=+1, MISMATCH=-1, INDEL=-2), and error threshold based on aligned adapter length.

**Parity**: Perfect parity achieved - 0 diff lines in integration tests, identical BAM outputs in end-to-end tests.

**Documentation**: See [docs/trimming.md](docs/trimming.md) for usage and algorithm details.

**Testing**: 
- FASTQ parity: `make test_trim_parity` (9/9 tests pass)
- End-to-end: `make test_trim_e2e` (alignment parity validated)
