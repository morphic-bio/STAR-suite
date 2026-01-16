# EM/VB Quantification - Questions for Planning Agent

## Algorithm Details

1. **Initialization Strategy**: What initialization should be used for transcript abundances?
   - Uniform across all transcripts?
   - Proportional to transcript length?
   - Small random perturbation to avoid symmetry?

2. **Convergence Criteria**: What are the default values for:
   - Maximum iterations? (e.g., 1000, 500?)
   - Relative tolerance threshold? (e.g., 1e-6, 1e-8?)
   - Should we use log-likelihood change or L1 norm on abundances as the primary convergence metric?

3. **Effective Length Handling**: The plan mentions "placeholder = length" for effective length. Should we:
   - Use raw transcript length throughout?
   - Implement fragment length distribution (FLD) based effective length correction (like Salmon)?
   - Make this configurable with a flag?

## BLAS/Vectorization

4. **BLAS Library Choice**: The plan specifies CBLAS/LAPACKE. Which implementation is preferred?
   - System CBLAS (may be slower)
   - OpenBLAS (faster, widely available)
   - Intel MKL (fastest on Intel, licensing considerations)
   - Should we support runtime detection/configuration?

5. **xsimd Requirement**: Is xsimd strictly required, or can we start with:
   - Compiler auto-vectorization with `-O3 -march=native`?
   - Add xsimd optimization in a later iteration?

## VB Engine Scope

6. **VB Priority**: The plan marks VB as "optional v1.1". Should the initial implementation:
   - Include VB infrastructure with a `--vb` flag from the start?
   - Defer VB entirely to a separate follow-up task?
   - What Dirichlet prior hyperparameters should be used? (Salmon uses alpha=1e-8 by default)

## Threading Model

7. **Parallelization Granularity**: What level of parallelism is expected?
   - Parallelize over ECs in the E-step (OMP parallel for)?
   - Parallelize over transcripts in the M-step?
   - Both?

8. **Thread Count**: Should the tool:
   - Use `OMP_NUM_THREADS` environment variable?
   - Accept a `--threads N` CLI argument?
   - Default to number of cores?

## Precision & Performance

9. **Floating Point Precision**: Should abundances and intermediate calculations use:
   - `double` (64-bit) throughout for numerical stability?
   - `float` (32-bit) for performance/memory with SIMD?
   - Configurable?

10. **Memory Layout**: For SIMD/cache efficiency, should we use:
    - Structure of Arrays (SoA) for transcript data?
    - Array of Structures (AoS) for ECs?
    - The plan mentions SoA - is this a hard requirement?

## Parity Testing

11. **Tolerance Threshold**: The plan specifies `relative diff < 1e-4` for parity with Salmon. Is this acceptable, or should we aim for tighter tolerance (e.g., 1e-6)?
    - Note: Salmon uses specific numerical methods (e.g., digamma approximations) that may cause minor differences.

12. **Which Salmon Output to Compare**: Should parity be measured against:
    - `NumReads` (estimated counts)?
    - `TPM` (normalized)?
    - Both?

## Directory Structure & Integration

13. **Source Location**: Should the EM library be placed at:
    - `source/libem/` (parallel to `libtrim`, `libflex`)?
    - `source/em/` (different naming)?
    - Something else?

14. **CLI Tool Location**: Should the CLI tool be placed at:
    - `tools/em_quant/` (parallel to existing tools)?
    - `tools/emquant/` (no underscore)?

15. **Future STAR Integration**: When EM is eventually integrated into STAR:
    - Will it replace existing `Transcriptome_quantAlign.cpp` logic?
    - Will it be an additional quantification mode (e.g., `--quantMode TranscriptEM`)?

## Deliverable Priority

16. **Implementation Order**: What's the priority order for deliverables?
    - Core EM engine first, then CLI, then tests?
    - Tests (TDD) first, then implementation?
    - CLI tool needed for manual testing before parity tests?

---

## Summary of Key Decisions Needed

| Topic | Options |
|-------|---------|
| Initialization | Uniform / Length-proportional / Random |
| BLAS library | System CBLAS / OpenBLAS / MKL / Configurable |
| xsimd | Required from start / Later optimization |
| VB scope | Include now / Defer to v1.1 |
| Precision | double / float / Configurable |
| Parity tolerance | 1e-4 / 1e-6 / Other |
| Source location | `source/libem/` / Other |
