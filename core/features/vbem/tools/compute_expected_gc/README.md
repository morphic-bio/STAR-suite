# compute_expected_gc

Compute expected GC distribution from transcriptome for GC bias correction.

## Usage

### From Transcriptome FASTA

```bash
compute_expected_gc --transcriptome transcripts.fa --output expected_gc.tsv
```

### From Genome + GTF

```bash
compute_expected_gc --genome genome.fa --gtf annotations.gtf --output expected_gc.tsv
```

## Options

- `--transcriptome FILE`: Transcriptome FASTA file
- `--genome FILE`: Genome FASTA file (requires `--gtf`)
- `--gtf FILE`: GTF annotation file (requires `--genome`)
- `--fld FILE`: Fragment length distribution file (default: built-in normal-like distribution)
- `--output FILE`: Output file (default: stdout)
- `--min-frag INT`: Minimum fragment length (default: 50)
- `--max-frag INT`: Maximum fragment length (default: 1000)
- `-h, --help`: Show help message

## Output Format

Tab-separated values with 101 rows (0-100% GC):
```
0	0.0001234567
1	0.0002345678
...
100	0.0003456789
```

## Fragment Length Distribution

If `--fld` is not provided, a default normal-like distribution is used:
- Mean: 200 bp
- Standard deviation: 80 bp
- Range: 50-1000 bp

FLD file format: one probability per line, indexed by fragment length (starting from length 0).

## Building

```bash
make
```

## Integration with STAR

This tool is automatically called during STAR genome generation if a GTF file is provided. The output `expected_gc.tsv` is stored in the genome directory and used during quantification for GC bias correction.
