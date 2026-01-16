# Trace-Reads Diff Harness (Alignment-Mode)

This directory contains tools to compare Salmon vs CLI alignment-mode decisions per read and verify decoy list parity and AS tag extraction.

**Note:** This CLI implements Salmon alignment-mode strictly (no pre-filtering). Compatibility gating and auxProb computation happen in `computeAuxProbs`.

## Tools

### 1. `trace_reads_diff.py`

Compares Salmon vs CLI trace outputs to pinpoint where filtering decisions diverge.

**Usage:**
```bash
python3 trace_reads_diff.py \
    --salmon-trace <salmon_trace.txt> \
    --cli-trace <cli_trace.txt> \
    --output <diff_report.txt>
```

**Input Format (Alignment-Mode):**

Both Salmon and CLI trace files use the same semicolon-delimited format:
```
<qname>\ttxpIDs=<ids>;as=<scores>;bestAS=<score>;logFragProb=<probs>;errLike=<likes>;logCompat=<compat>;orphan=<flags>;auxProb=<probs>;
```

Where:
- `txpIDs`: Comma-separated transcript IDs
- `as`: Comma-separated alignment scores
- `bestAS`: Best alignment score
- `logFragProb`: Comma-separated log fragment probabilities
- `errLike`: Comma-separated error likelihoods
- `logCompat`: Comma-separated log compatibility probabilities
- `orphan`: Comma-separated orphan flags (0/1)
- `auxProb`: Comma-separated auxiliary probabilities

**Output:**

The diff report includes:
- Summary of differences by type (bestAS, txpIDs, as, logCompat, auxProb mismatches)
- First divergence point
- Detailed differences (up to 100)

### 2. `verify_decoy_as_tags.py`

Verifies decoy list parity and AS tag extraction.

**Usage:**
```bash
python3 verify_decoy_as_tags.py \
    --bam <input.bam> \
    --salmon-decoys <salmon_decoys.txt> \
    --cli-decoys <cli_decoys.txt> \
    --sample-size <N> \
    --output <report.txt>
```

**Input Format:**

Decoy list files can contain either transcript IDs (numeric) or transcript names (one per line):
```
0
1
2
...
```

or:
```
ENST00000000001
ENST00000000002
...
```

The script automatically maps names to IDs using the BAM header.

**Output:**

The report includes:
- Decoy list comparison (Salmon vs CLI)
- AS tag presence rate
- List of reads missing AS tags (if any)

## Generating Trace Files

### CLI Trace

Enable trace output in the EC filter CLI:

```bash
tools/ec_filter_test/ec_filter_cli \
    --input <input.bam> \
    --transcripts <transcriptome.fa> \
    --output <eq_classes.txt> \
    --trace-reads <cli_trace.txt> \
    --trace-limit <N> \
    --lib-type <IU|ISF|...> \
    --threads <N>
```

The trace file will contain:
- `READ` lines with best_score, min_score_threshold, num_alignments
- `ALN` lines with tid, score, compat flag, prefilter_kept status
- `ALN_POST` lines with tid, kept status after filtering

### Salmon Trace

To generate Salmon trace output, you'll need to patch Salmon to output filtering decisions. See `salmon_trace_patch.patch` for reference.

## Example Workflow

1. **Generate CLI trace:**
```bash
tools/ec_filter_test/ec_filter_cli \
    --input /path/to/Aligned.toTranscriptome.name.bam \
    --transcripts /path/to/transcriptome.fa \
    --output /tmp/eq_classes.txt \
    --trace-reads /tmp/cli_trace.txt \
    --trace-limit 1000 \
    --lib-type A \
    --threads 16
```

2. **Generate Salmon trace** (requires patched Salmon):
```bash
# Set environment variable to enable tracing
export SALMON_TRACE_FILE=/tmp/salmon_trace.txt
salmon quant --alignments <input.bam> ...
```

3. **Compare traces:**
```bash
python3 trace_reads_diff.py \
    --salmon-trace /tmp/salmon_trace.txt \
    --cli-trace /tmp/cli_trace.txt \
    --output /tmp/diff_report.txt
```

4. **Verify decoy lists and AS tags:**
```bash
python3 verify_decoy_as_tags.py \
    --bam /path/to/Aligned.toTranscriptome.name.bam \
    --salmon-decoys /path/to/salmon_decoys.txt \
    --cli-decoys /path/to/cli_decoys.txt \
    --sample-size 1000 \
    --output /tmp/verification_report.txt
```

## What to Look For

### In Diff Report

- **bestAS_mismatch**: Best alignment score differs (indicates score extraction issue)
- **txpIDs_mismatch**: Transcript IDs differ (indicates alignment assignment issue)
- **as_mismatch**: Alignment scores differ (indicates score extraction issue)
- **logCompat_mismatch**: Log compatibility probabilities differ (indicates compatibility checking issue)
- **auxProb_mismatch**: Auxiliary probabilities differ (indicates auxProb computation issue)

### In Verification Report

- **Decoy list parity**: Should match exactly (same transcripts marked as decoys)
- **AS tag presence**: Should be 100% (all alignments should have AS tags)

## Troubleshooting

### No differences found but parity is low

- Check if trace files cover the same reads (same read order)
- Verify decoy lists match
- Verify AS tag extraction is correct

### First divergence is early

- Check decoy list parity first
- Verify AS tag extraction
- Compare best_score extraction logic

### logCompat mismatches

- Verify library format detection matches Salmon (use --lib-type A for auto-detect)
- Check paired-end compatibility implementation
- Verify observed format computation
- Check mate field population (mate_pos, mate_is_forward)