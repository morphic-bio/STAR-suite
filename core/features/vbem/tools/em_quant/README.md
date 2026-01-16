# em_quant - EM/VB Quantification Tool

A standalone C++ tool for quantifying transcript abundances from equivalence classes using Expectation-Maximization (EM) or Variational Bayes (VB) algorithms.

## Overview

`em_quant` implements EM and VB algorithms for transcript quantification from Salmon-style equivalence classes. It achieves parity with Salmon when bias correction is disabled.

## Building

```bash
cd tools/em_quant
make
```

This will:
1. Build `libem.a` library in `source/libem/`
2. Build `em_quant` CLI tool

## Dependencies

- C++11 compiler (g++ or clang++)
- OpenMP (for parallel E-step)
- OpenBLAS (optional, for future BLAS operations)

## Usage

```bash
em_quant [options] -e <eq_classes.txt> -o <output.tsv>
```

### Required Arguments

- `-e, --ec <file>`: Equivalence classes file (Salmon `eq_classes.txt` format)
- `-o, --output <file>`: Output TSV file

### Optional Arguments

- `-l, --lengths <file>`: Transcript lengths file (Salmon `quant.sf` format). If not provided, placeholder lengths (1.0) are used.
- `--vb`: Use Variational Bayes instead of EM
- `--vb-prior <float>`: Dirichlet prior concentration for VB (default: 0.01)
- `--max-iters <int>`: Maximum EM/VB iterations (default: 200 for VB, matches Salmon)
- `--tolerance <float>`: Convergence tolerance for relative change in log-likelihood/ELBO (default: 1e-8 for VB, matches Salmon)
- `--zero-threshold <float>`: Zero out NumReads below threshold (default: 0.001). Set to 0 to disable zeroing. This matches Salmon's behavior of zeroing very low-abundance transcripts.
- `--init-by-length`: Initialize abundances proportional to transcript length (default: uniform)
- `--uniform-init`: Use uniform initialization for VB (alpha = prior only). Default: prior + unique_counts (backward compatible). Salmon uses `txp.projectedCounts` + uniform mixing, not single-transcript EC counts. The unique_counts seeding can push to different local optima along flat/ambiguous directions.
- `--threads <int>`: Number of threads for parallel E-step (default: OMP_NUM_THREADS)
- `-v, --verbose`: Verbose output
- `-h, --help`: Show help message

## Input Format

### Equivalence Classes (`eq_classes.txt`)

Salmon format:
```
<num_transcripts>
<num_equivalence_classes>
<transcript_name_1>
<transcript_name_2>
...
<transcript_name_N>
<n_txp> <idx1> <idx2> ... <count>
...
```

Example:
```
3
3
T1
T2
T3
1 0 100
1 1 50
2 0 1 30
```

### Transcript Lengths (`quant.sf`)

Salmon TSV format with columns: `Name`, `Length`, `EffectiveLength`, `TPM`, `NumReads`

## Output Format

TSV file with columns:
- `Name`: Transcript identifier
- `Length`: Raw transcript length
- `EffectiveLength`: Effective length (equals Length for v1)
- `TPM`: Transcripts Per Million
- `NumReads`: Estimated fragment count

## Examples

### Basic EM quantification

```bash
em_quant -e eq_classes.txt -l quant.sf -o output.tsv
```

### VB quantification with custom prior

```bash
em_quant -e eq_classes.txt -l quant.sf -o output.tsv --vb --vb-prior 0.05
```

### Length-weighted initialization

```bash
em_quant -e eq_classes.txt -l quant.sf -o output.tsv --init-by-length
```

## Testing

### Unit Test (Synthetic Fixture)

```bash
cd tools/em_quant
./em_quant -e ../../test/fixtures/em_unit/eq_classes.txt -o /tmp/em_unit_output.tsv
```

### Parity Test (Salmon Fixture)

```bash
cd tools/em_quant
make test
# or manually:
./run_parity_test.sh
```

**Important:** The Salmon fixture must be generated with bias correction disabled:
```bash
salmon quant ... --noLengthCorrection --noFragLengthDist --noEffectiveLengthCorrection
```
This ensures Salmon uses raw transcript lengths (matching em_quant's v1 behavior) instead of FLD-adjusted effective lengths.

## Algorithm Details

### EM Algorithm

1. **Initialization**: Uniform or length-weighted abundances
2. **E-step**: Compute expected counts for each transcript from ECs (parallelized)
   - Responsibilities use effective-length weighting: `r_i = (abundance_i / eff_length_i) / sum_j(abundance_j / eff_length_j)`
   - Matches Salmon's `theta / effLen` formulation
3. **M-step**: Normalize expected counts to update abundances
4. **Convergence**: Check relative change in log-likelihood
5. **Zeroing**: Zero out transcripts with NumReads below threshold (default: 0.001) to match Salmon's behavior
6. **TPM**: Computed as `(count_i / eff_length_i) / sum_j(count_j / eff_length_j) * 1e6` (Salmon formula)

### VB Algorithm

Matches Salmon's VBEMOptimizer implementation:

1. **Initialization**: 
   - Default (backward compatible): `alpha_i = prior + unique_counts_i` where `unique_counts_i` = count of reads in single-transcript ECs for transcript i
   - With `--uniform-init`: `alpha_i = prior` only (matches Salmon's `txp.projectedCounts` + uniform mixing approach)
   - The unique_counts seeding can push to different local optima along flat/ambiguous directions
   - Initial abundances normalized from alpha

2. **E-step**: Compute expected counts using `exp(digamma(alpha_i)) / effLen_i` weights
   - Responsibilities: `r_i = (exp(digamma(alpha_i)) / effLen_i) / sum_j(exp(digamma(alpha_j)) / effLen_j)`
   - Matches Salmon's VBEMOptimizer E-step exactly

3. **M-step**: Update `alpha_i = prior + new_expected_counts_i`
   - Abundances normalized from alpha: `theta_i = alpha_i / sum_j(alpha_j)`

4. **Convergence**: ELBO (Evidence Lower BOund) with tolerance 1e-8, max iterations 200

5. **Zeroing**: Zero transcripts where `alpha_i < prior + 1e-8` AND no unique evidence
   - Matches Salmon's implicit zeroing behavior (not an external threshold)

6. **TPM**: Computed as `(count_i / eff_length_i) / sum_j(count_j / eff_length_j) * 1e6` (Salmon formula)

## Limitations (v1)

- No bias correction (sequence, GC, positional)
- Effective length = raw transcript length (no fragment length distribution)
- No BLAS operations (relying on compiler optimizations)

## Future Enhancements

- Fragment length distribution-based effective length
- Bias correction (sequence, GC, positional)
- xsimd SIMD optimizations
- Integration into STAR as `--quantMode TranscriptEM`

## License

GPL (compatible with STAR)
