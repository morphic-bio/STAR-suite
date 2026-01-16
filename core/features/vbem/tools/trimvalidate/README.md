# trimvalidate - Cutadapt Parity Trimming Validator

Standalone CLI tool for validating cutadapt-parity trimming implementation against Trim Galore/cutadapt outputs.

## Building

```bash
cd tools/trimvalidate
make
```

This will build `trimvalidate` and link it against `libtrim.a` from `source/libtrim/`.

## Usage

```bash
./trimvalidate -1 input_R1.fastq -2 input_R2.fastq \
               -o1 trimmed_R1.fastq -o2 trimmed_R2.fastq \
               [--quality 20] [--length 20] \
               [--adapter-r1 SEQ] [--adapter-r2 SEQ]
```

### Options

- `-1 <file>`: Input R1 FASTQ file (required)
- `-2 <file>`: Input R2 FASTQ file (required)
- `-o1 <file>`: Output R1 FASTQ file (required)
- `-o2 <file>`: Output R2 FASTQ file (required)
- `--quality <int>`: Quality threshold for trimming (default: 20)
- `--length <int>`: Minimum read length after trimming (default: 20)
- `--adapter-r1 <seq>`: Custom R1 adapter sequence (default: TruSeq R1)
- `--adapter-r2 <seq>`: Custom R2 adapter sequence (default: TruSeq R2)

## Output Format

- Standard 4-line FASTQ format
- Reads that pass min_length filter are written
- Dropped pairs (either mate < min_length) are omitted from output
- No read-name annotations (unlike cutadapt's `length=X` annotations)
- Same read order as input

## Parity Testing

### Running Tests

Run automated parity tests against Trim Galore expected outputs:

```bash
# From tools/trimvalidate directory
./run_parity_test.sh
# Or use Makefile:
make test

# From source directory (ensures libtrim is built)
cd ../../source && make test_trim_parity
```

### Test Coverage

The test suite includes:

1. **Synthetic Fixtures** (`test/fixtures/trim/`):
   - 7 test cases covering various trimming scenarios
   - Each fixture has `input_R*.fastq` and `expected_R*.fastq` files
   - Expected outputs generated with Trim Galore

2. **Integration Test** (`test/integration/trim/nfcore_smoke/`):
   - Real-world-like dataset derived from nf-core test data
   - Tests behavior on realistic read characteristics
   - Results stored in `test/integration/trim/nfcore_smoke/results/`:
     - `status.txt` - PASS/FAIL status
     - `diff_R1.txt` - Unified diff for R1 output
     - `diff_R2.txt` - Unified diff for R2 output

### Test Output

The test script reports:
- Individual fixture pass/fail status
- Summary counts (synthetic fixtures and integration test)
- Diff output for failed tests
- Overall pass/fail status (exit code)

See `test/integration/trim/nfcore_smoke/README.md` for details on the integration dataset and how to regenerate expected outputs.

## Regenerating Expected Outputs

To regenerate expected outputs for fixtures using Trim Galore:

```bash
cd test/fixtures/trim/<fixture_name>
trim_galore --paired --quality 20 --length 20 \
            --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            input_R1.fastq input_R2.fastq

mv input_R1_val_1.fq expected_R1.fastq
mv input_R2_val_2.fq expected_R2.fastq
```
