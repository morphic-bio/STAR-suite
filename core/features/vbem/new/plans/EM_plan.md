# EM/VB Quantification Plan

## Objective
Implement a standalone C++ EM/VB quantifier over transcript equivalence classes (ECs), with strict parity tests against Salmon on small fixtures. For v1, parity will be measured against Salmon runs with bias correction disabled (no sequence/GC/positional bias), so we can match with raw-length effective lengths. Bias modeling will be added in a later stage to mirror Salmon’s bias corrections. Keep dependencies minimal (BLAS/CBLAS + xsimd/OpenMP; avoid Eigen).

## Scope
- Input: EC definitions and counts (from Salmon `--dumpEq` fixtures).
- Output: Per-transcript estimated counts/TPM-like values, matching Salmon within tolerance.
- Not in this phase: Bias correction, effective length adjustment, integration into STAR; focus on core EM/VB loop correctness and performance.

## Milestones & Tasks

### 1) Fixtures & I/O
- [ ] Define a compact EC input format (reuse Salmon `eq_classes.txt` + transcript map) and a small JSON/TSV for expected transcript abundances (from Salmon `quant.sf`).
- [ ] Add loader for `aux_info/eq_classes.txt` and transcript ID map from the generated fixture (chr22 slice, SRR6357070 subset).
- [ ] Add a tiny synthetic EC fixture (handcrafted 3–4 transcripts, few ECs) for deterministic unit tests.

### 2) Core Data Structures
- [ ] EC table: vector of {transcript_ids (sorted), count}.
- [ ] Transcript state: effective length (placeholder = length), current abundance, prior (for VB).
- [ ] Hash/interner for transcript-id tuples to EC IDs.

### 3) EM Engine (baseline, bias-free)
- [ ] Implement plain EM updates over ECs:
  - E-step: responsibility of each transcript within an EC (proportional to abundance).
  - M-step: update transcript abundances (normalize to counts/effective length).
- [ ] Convergence criteria: relative change in log-likelihood or L1 on abundances; max iters; tolerance.
- [ ] Multi-threading: parallelize over ECs with OpenMP; use SoA layout and xsimd/CBLAS for inner loops.

### 4) VB Engine (optional v1.1)
- [ ] Add Dirichlet priors; update variational parameters (gamma) per transcript.
- [ ] Use digamma-based updates (cephes/boost for psi) if needed; keep behind a flag.

### 5) CLIs & API
- [ ] CLI tool (e.g., `tools/em_quant`) that ingests EC fixture and emits per-transcript estimates (TSV).
- [ ] Library-style API (header) for future STAR integration: `run_em(const ECs&, Params) -> Estimates`.

### 6) Testing & Parity
- [ ] Unit tests on synthetic ECs (exact expected abundances).
- [ ] Parity test: run on Salmon chr22 fixture with bias correction disabled (no `--seqBias`, `--gcBias`, or `--posBias`), and compare estimates to Salmon `quant.sf` within tolerance (e.g., relative diff < 1e-4 for nonzero entries).
- [ ] Performance sanity: run on fixture and report runtime; ensure OMP/SIMD path works.

### 7) Bias Modeling (future phase)
- [ ] Add a bias layer to learn/apply sequence/GC/positional bias (mirroring Salmon) and effective-length corrections.
- [ ] Update parity tests to run against Salmon with bias enabled once the bias layer is implemented.

### 8) Documentation
- [ ] README for the EM tool describing input format, parameters, and how to regenerate fixtures.
- [ ] Note dependencies (CBLAS/LAPACKE or MKL/OpenBLAS, xsimd/OpenMP) and build flags.

## Dependencies & Constraints
- Prefer CBLAS/LAPACKE + xsimd; avoid Eigen/Armadillo.
- License compatibility with STAR (GPL) is required.
- Keep code self-contained; no Salmon linkage, only fixture-driven parity.

## Deliverables
- Source: `src/libem/` (or similar) with EM/VB implementation and CLI under `tools/em_quant/`.
- Tests: synthetic unit tests + Salmon fixture parity test wired into CI/Makefile.
- Docs: README describing usage, fixtures, and regeneration steps.
