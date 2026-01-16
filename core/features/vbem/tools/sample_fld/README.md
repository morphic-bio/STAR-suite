# sample_fld - Fragment Length Distribution Sampler

Extracts fragment length distribution from paired-end BAM files, matching Salmon's FLD sampling approach.

## Usage

```bash
sample_fld --bam <input.bam> [--output <output.tsv>] [--min-len <int>] [--max-len <int>]
```

## Options

- `--bam FILE`: Input BAM file (required)
- `--output FILE`: Output TSV file (default: stdout)
- `--min-len INT`: Minimum fragment length to include (default: 0)
- `--max-len INT`: Maximum fragment length to include (default: 1000)

## Output Format

TSV file with three columns:
1. Fragment length (integer)
2. Count (number of fragments observed)
3. Probability (normalized frequency)

## Algorithm

- Reads paired-end BAM records
- Filters for properly paired reads (FLAG 0x2)
- Uses TLEN (template length) field for fragment size
- Only counts read1 to avoid double-counting pairs
- Builds histogram and normalizes to probabilities

## Building

```bash
make
```

Requires htslib (included in STAR-Flex source tree).
