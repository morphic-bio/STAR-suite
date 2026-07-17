# Runbook: molecule-first barcode assignment integration

Date: 2026-07-17  
Branch: `agent/molecule-first-barcode-assignment`  
Base: `v1.4.4` (`20efad80bd71ec9004b73c2f003fa8f3decdf20e`)  
Target release: `v1.5.0`

## Objective

Integrate the frozen, candidate-preserving barcode reference into STAR Suite
as a native, assay-neutral stage that runs after sequence/feature evidence is
available and before bin filtering, EmptyDrops, or any other cell call.

The operation order is **read-to-clique-to-UMI**:

```text
candidate-bearing reads
  -> globally candidate-consistent read cliques
  -> candidate-specific UMI correction and occupancy
  -> strict / soft expected / hard / gated-hard products
  -> spatial-bin or single-cell calling
```

The mandatory invariant is:

> Ambiguous molecules may refine a conservative broad field, but they may not
> manufacture unsupported fine-scale structure.

The current STAR Suite decoder and its existing ambiguity resolver remain the
default compatibility baseline. This tranche must not retune or silently
replace them.

## Release ordering

1. `v1.4.4` is the pre-integration patch baseline. It contains the feature
   barcode namespace-remapping bug fix made after `v1.4.3` and the associated
   workflow fixes.
2. This branch is created from the signed content of `v1.4.4`, not from a
   dirty working tree or an assay-specific development branch.
3. `v1.5.0` is cut only after native/reference conformance, default-off STAR
   regressions, and the applicable 100K assay smokes pass.

## Frozen evidence model

Allowed evidence in the first native implementation is limited to:

- the decoder's finite candidate set for each read;
- per-candidate sequence log likelihood from observed bases and Phred scores;
- PCR-inclusive exact-read candidate frequencies as a read prior;
- shared feature identity and raw UMI for read-clique construction;
- candidate-specific exact or directional one-mismatch UMI correction.

PCR amplification is intentionally retained in the read-frequency prior
because the input reads have not been deduplicated. The prior is added once
per read clique. It is never copied once per PCR read and is never reported as
UMI mass.

The frozen default parameters are:

```text
temperature=1
prior_alpha=1
prior_beta=1
spatial_lambda=0
umi_mode=1mm_cr
gate_min_posterior=0.95
gate_min_margin=0.90
```

Exact UMI is retained as a sensitivity output. No expression, cell-type,
neighborhood, graph, image, GPU, or spatial-density prior is permitted in this
release.

## Integration boundary

The native implementation has two layers:

1. A reusable C++ resolver library under `flex/source/solo/` owns deterministic
   clique construction, posterior normalization, candidate-specific UMI
   correction, weighted occupancy, and hard-call gates.
2. A STAR Suite executable under `flex/tools/molecule_first_resolver/` consumes
   a normalized candidate ledger and writes all four policy products in one
   invocation. The top-level build, Debian package, and release checks include
   this executable.

STAR 3-prime, STAR-Flex, and Visium HD adapters provide the same normalized
ledger. This separates feature assignment from barcode policy while placing
the policy before cell calling. The barcode universe is fixed by the supplied
candidate table; the resolver may neither add a candidate nor infer one from
expression.

The native tool is the release surface for `v1.5.0`. Existing STAR alignment,
Solo counting, Flex counting, and compatibility decoder code paths remain
byte-for-byte selectable by omitting the explicit molecule-first stage.

## Normalized input contract

The input is a UTF-8, tab-separated table with one row per read/candidate and
this header:

```text
read_id	feature_id	raw_umi	candidate	log_sequence_likelihood	exact_read_count
```

Rules:

- `(read_id, candidate)` is unique;
- all rows for a read have the same nonempty feature and raw UMI;
- log likelihoods are finite;
- exact-read counts are non-negative integers and agree for a candidate across
  the file;
- a read with no finite candidate is excluded before this stage and remains
  represented in upstream rejection accounting;
- input order has no semantic effect.

For composite Visium HD barcodes, an adapter may supply the frozen factorized
BC1/BC2 log prior as an explicit candidate prior column in a later schema.
The `v1.5.0` generic schema uses the whole-candidate exact-read count so the
same engine applies to Chromium and spatial barcode universes. The adapter
must record which prior representation it supplied.

## Clique and posterior contract

Reads are grouped by `(feature_id, raw_umi)`, sorted deterministically, and
greedily partitioned. A read joins only a clique whose global candidate
intersection remains nonempty. The selected destination is the
lexicographically smallest resulting intersection followed by the earliest
clique index.

This rule blocks an unrestricted `A-B, B-C, C-D` overlap chain from
manufacturing an unsupported `A-D` family.

For candidate `b` in clique `m`:

```text
log_evidence(m,b) =
    sum_read log_sequence_likelihood(read,b) / temperature
  + prior_beta * log(exact_read_count(b) + prior_alpha)
```

The candidate posterior is the deterministic log-sum-exp normalization over
the final clique intersection. Sequence evidence accumulates over PCR reads;
the exact-read prior is applied once to the shared clique coordinate.

## Candidate-specific UMI contract

For each `(feature, candidate)`, raw-UMI support is the sum of clique posterior
mass. Directional one-mismatch correction uses the existing frozen rule:

```text
support(parent) >= 2 * support(child) - 1
```

Roots are selected deterministically by descending support and then raw UMI.
For `(feature, corrected_umi, candidate)` containing clique probabilities
`p_1 ... p_n`, the expected molecule occupancy is:

```text
1 - product_i(1 - p_i)
```

Occupancy is not renormalized after deduplication.

## Required products

Every enabled run writes:

- `read_cliques.tsv`: member reads, candidate intersection, sequence score,
  prior, evidence, and posterior;
- `strict_molecules.tsv`: only single-candidate cliques, deduplicated by the
  candidate-specific corrected UMI;
- `soft_expected_molecules.tsv`: weighted candidate occupancies;
- `hard_molecules.tsv`: ungated post-clique posterior argmax, grouped after the
  same candidate-specific UMI correction;
- `gated_hard_molecules.tsv`: hard calls meeting both frozen gates;
- `hard_call_audit.tsv`: assigned/deferred state, top posterior, margin, and
  deterministic reason;
- `resolved_config.tsv` and `summary.tsv`: parameters, input/product counts,
  posterior mass, occupancy mass, and mass residuals.

Strict is the conservative fine-scale product. Ungated hard is the reference
integer recovery product. Soft expected count preserves calibrated
uncertainty. Gated hard is a coverage-changing safety view. None may overwrite
another.

## Determinism and failure behavior

- Use `LC_ALL=C` lexical ordering semantics implemented in C++, independent of
  hash iteration order.
- Stable tie breaking is smallest candidate identifier.
- Reject duplicate read/candidate rows, inconsistent read metadata,
  inconsistent candidate priors, non-finite values, invalid gates, and output
  creation failures.
- Write into a new or empty output directory. A failed run must not present a
  complete summary.
- Floating-point comparisons use an absolute tolerance of `1e-12` and a
  relative tolerance of `1e-10`; mass reconciliation uses
  `1e-9 * max(1, mass)`.

## Required tests

### Native unit tests

- global candidate intersection blocks overlap chains;
- clique construction is invariant to input permutation;
- likelihood accumulates per PCR read while the read prior is added once;
- posterior normalizes to one;
- exact duplicate PCR reads emit one clique mass;
- exact and directional one-mismatch UMI roots are deterministic;
- occupancy equals `1 - product(1 - p)` and is bounded by one;
- strict, hard, and gated products obey their definitions;
- unsupported candidates can never appear;
- malformed input fails closed.

### Reference conformance

Run the same tracked synthetic ledger through the native implementation and
the frozen Python reference in `visium-hd-processing` at its pinned commit.
Compare clique membership, candidates, posterior values, UMI roots, all four
products, and summaries within the declared tolerances.

Run the native implementation on the existing deterministic Chromium 3-prime,
Chromium Flex, Visium HD 3-prime, and Visium HD Flex 100K ledgers. Compare
native product hashes or canonical numeric tables with the already sealed
reference products.

### STAR Suite regressions

From a clean build, with no molecule-first invocation:

- `tests/run_cbub_regression_test.sh`;
- standard STARsolo smoke;
- Flex smoke and Flex configuration smoke;
- bridge determinism harness;
- Tier A container CI;
- release packaging/install checks.

The default-off outputs must match the `v1.4.4` baseline for the selected
fixtures. Any intentional packaging-only difference must be listed by path.

### Enabled integration checks

- top-level and Debian builds install the native executable;
- the executable reports the STAR Suite version;
- all four policies materialize before a fixed cell/bin call;
- clean-read controls do not regress;
- challenged ambiguous subsets reproduce the frozen direction of effect;
- no Space Ranger or Cell Ranger assignment field is required by the primary
  run.

## Release gates

`v1.5.0` requires all of the following:

1. this runbook is committed before implementation;
2. native unit and Python-reference conformance tests pass;
3. four 100K conformance smokes pass or any unavailable external fixture is
   recorded explicitly rather than silently skipped;
4. clean STAR core build and applicable default-off regressions pass;
5. GitHub PR CI passes on the integration commit;
6. a release-preparation commit updates `VERSION`, Debian metadata, binary
   distribution documentation, release defaults, and `v1.5.0` notes;
7. the annotated `v1.5.0` tag points to the merged, tested master commit;
8. the GitHub release workflow and published artifacts complete successfully.

Space Ranger and Cell Ranger are final compatibility diagnostics only. Their
source code is outside this work's clean-room boundary.

## Execution ledger

Populate this section during implementation. Do not replace planned commands
with claims until their exit status and artifact paths are recorded.

| Gate | Command or artifact | Result |
| --- | --- | --- |
| v1.4.4 local build | clean legacy STAR build | passed before branch creation |
| v1.4.4 focused regression | namespace plus workflow-contract tests | passed before tag |
| v1.4.4 Tier A | GitHub Actions run `29560701174` | passed |
| Native unit tests | `make -C flex/tools/molecule_first_resolver clean test`; ASan/UBSan build | passed |
| Python conformance | `tests/run_molecule_first_reference_conformance.sh` | passed; tracked ledger exact within tolerance |
| Four 100K conformance smokes | `tests/molecule_first/run_100k_conformance.py`; reports under `20260717_star_suite_molecule_first_integration_v1` | passed; all integer residuals zero |
| Default-off STAR regressions | STARsolo, Flex on/off, Flex config, and 100K bridge determinism; v1.4.4 tree/hash comparisons | passed; CB/UB fixture test explicitly unavailable |
| Committed package snapshot | Ubuntu 22.04 release-build and runtime-check containers; installed native smoke | passed |
| Integration PR CI | GitHub Actions run `29563199754`; PR #5 merge `b03930adff707f97c009d4d6116a677ff38168dc` | passed |
| v1.5.0 release workflow | pending | pending |

## Recovery after interruption

1. Verify `v1.4.4` resolves to `20efad80bd71ec9004b73c2f003fa8f3decdf20e`.
2. Verify the current branch descends from that tag and the working tree is
   clean except for the intended tranche.
3. Read this runbook and the implementation report before editing code.
4. Re-run native unit tests before restarting an external 100K conformance
   job.
5. Never reuse a partially written output directory as a passed artifact.
6. Do not modify or re-estimate frozen oracle parameters during integration.
