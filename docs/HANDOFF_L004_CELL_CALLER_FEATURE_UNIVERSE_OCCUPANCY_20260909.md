# L004 cell caller: feature universe and GEM occupancy

Written 2026-09-09. Continues
[HANDOFF_L004_CELL_CALLER_DIVERGENCE_REVIEW_20260908.md](HANDOFF_L004_CELL_CALLER_DIVERGENCE_REVIEW_20260908.md),
which left the divergence unexplained and named OrdMag calibration as the next priority.

**No production source was changed. No production default was adopted.** Every result below
is a diagnostic run of saved binaries on saved matrices, or a post-hoc filter on saved
callsets. Clean-room boundary respected throughout: public 10x documentation, permitted
Cell Ranger logs, CR output matrices and callsets, and our own tools only. No Cell Ranger
source was inspected, and the excluded `cr_native_mex_calling_20260908` material was not
opened.

## Summary

Two changes account for most of the L004 gap, and two previously suspected causes are now
ruled out by measurement.

| Change | Lymph-node Jaccard | Scope of evidence |
|---|---|---|
| Baseline, adaptive ambient window | 0.833996 | measured |
| Fixed ambient window `[90,000, 180,000)` | 0.892115 | measured, prior work |
| Model universe 18,129 + 375 included=FALSE genes, on CR stage sets | 0.934369 | measured, LN only |
| Plus GEM-occupancy filter | **0.945550** | measured, LN only |

Across all eight L004 samples the occupancy filter alone moves the pooled Jaccard from
**0.943858 to 0.955637**.

Ruled out: the Monte Carlo simulation count, and STAR's counting itself. On a matched
feature universe our caller on STAR counts and on CR counts agree at Jaccard **0.967682**,
with STAR marginally closer to the official callset. Raising simulations from 10,000 to
100,000 flips 134 of 31,492 tail decisions and slightly *lowers* agreement.

## 1. GEM-occupancy post-filter

Cell Ranger's public Flex algorithm page states that after per-sample cell calling it
estimates the mean of the Poisson distribution describing probe barcodes observed per GEM
"under optimal chip loading", estimates the 99.9th percentile of that distribution, and
removes all cells belonging to GEMs beyond it. The current tag-aware path has no such step;
`--soloFlexCellCaller legacy` had one, and the 2026-09-08 integration deliberately dropped it.

On L004, CR's final calls never exceed **7** called tags per GEM (749 GEMs at 7, zero at 8,
a hard cliff). Ours reach **13**.

| Lambda variant | lambda | threshold k | GEMs flagged | calls removed | official lost | extras removed |
|---|---:|---:|---:|---:|---:|---:|
| Plain mean over occupied GEMs | 2.147 | 8 | 144 | 1,358 | 7 | 1,351 |
| Zero-truncated Poisson MLE | 1.788 | **7** | 404 | 3,438 | 119 | 3,319 |
| Fixed `--max-occupancy 7` | n/a | 7 | 404 | 3,438 | 119 | 3,319 |

**The zero-truncated MLE reproduces CR's observed ceiling; the plain mean does not.** The
boundary between k=7 and k=8 sits at lambda 1.970. Fitting to CR's own callsets gives the
same split (mean 2.163 to k=8, ZTP 1.809 to k=7).

Eight-sample effect, applied to the saved adaptive-window callsets, k=7:

| Sample | Jaccard before | Jaccard after | Removed extras / official |
|---|---:|---:|---|
| Glioblastoma | 0.980929 | 0.992911 | 414 / 14 |
| Kidney | 0.976075 | 0.991747 | 457 / 15 |
| Breast cancer | 0.968894 | 0.987574 | 256 / 9 |
| Skin melanoma | 0.966886 | 0.978650 | 468 / 15 |
| Endo | 0.958807 | 0.970734 | 435 / 12 |
| Lung cancer | 0.950730 | 0.958645 | 491 / 22 |
| Colorectal | 0.886268 | 0.899881 | 594 / 25 |
| Lymph node | 0.833996 | 0.842743 | 204 / 7 |
| **Pooled** | **0.943858** | **0.955637** | 3,319 / 119 |

The cross-tab is a sharper test than the Jaccards. Of the **404** GEMs where we call more
than 7 tags, CR has **zero cells in 387 and exactly 7 in the other 17, never 1 to 6**. The
119 official cells we lose are exactly those 17 GEMs where CR sits at its ceiling and we
call one or two tags more. Conversely, 54 GEMs where we call only 6 or 7 tags are CR-empty
because CR's own pre-filter occupancy exceeded 7; they hold 374 extras this rule cannot
reach from our side.

CR's LN log shows 18,423 initial-plus-rescued calls against 18,173 final, a 250-cell
reduction. Our LN calls inside CR-emptied high-occupancy GEMs number 219, and the k=7 filter
removes 211 to 233 LN calls depending on the arm. Same magnitude; this is supporting
evidence, not an exact reconstruction.

The effect is close to orthogonal to the ambient and OrdMag work, adding **+0.009 to +0.011**
to every lymph-node arm tested, including the fixed-window baseline, the 5,773-primary
restriction, the CR-rank ambient arm and the 56,000 budget arm.

**Why this looked harmless before.** The 2026-03-19 benchmark measured the legacy filter
removing 21 cells on the 4-sample JAX library, because mean GEM occupancy there is about
0.2. On L004 it is 2.2. The step is load-dependent, not inert.

## 2. Feature universe

CR called cells on its full **39,186**-feature raw matrix. Its logged 38,749 eligible
barcodes and 1,765-UMI primary boundary match raw all-feature totals exactly, not totals
over the 18,129 included probes. STAR exports only the 18,129.

Excluded UMIs in the CR LN raw matrix, from the output H5 and the public probe CSV:

| Class | Features | UMIs |
|---|---:|---:|
| `DEPRECATED_*` rows | 516 with counts, 542 present | 1,078,671 |
| Genes with active `included=FALSE` probes | 375 with counts, 397 present | 637,322 |
| Non-probe genes | 20,118 | 0 |
| **Total outside the 18,129** | | **1,715,993 of 55,770,908 (3.08%)** |

The loss is differential, which is the missing scaling factor: **4.45%** of UMIs in the top
5,773 barcodes versus **2.91%** in the ambient window. Largest contributors in the top
barcodes are `DEPRECATED_IGKC` (384,930), `IGLC2` (76,472), `IGLC3` (56,712), `IGHG2`
(45,737) and `DEPRECATED_IGLC1` (29,053). The probe CSV carries 1,121 `included=FALSE` rows,
498 active plus 623 deprecated.

CLI arms, all lymph node, fixed window, FDR 0.01, against the 18,173 official CR cells:

| Arm | Matrix / universe | sims | Primary | Tail passes | Calls | Missed | Extra | Jaccard |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| A0 | STAR, 18,129 | 10K | 6,912 | 10,886 | 17,798 | 1,213 | 838 | 0.892115 |
| A1 | STAR, 18,129 | 100K | 6,912 | 10,814 | 17,726 | 1,273 | 826 | 0.889521 |
| B1 | CR, all 39,186 | 10K | 5,222 | 15,412 | 20,634 | 134 | 2,595 | 0.868596 |
| B2 | CR, 18,129 | 10K | 6,513 | 10,759 | 17,272 | 1,521 | 620 | 0.886075 |
| B3 | CR, all 39,186 | 100K | 5,222 | 15,618 | 20,840 | 133 | 2,800 | 0.860154 |
| B4 | CR, nonzero-total 18,926 | 10K | 5,222 | 16,119 | 21,341 | 119 | 3,287 | 0.841286 |
| B5 | CR, 18,129 + 375 incl=FALSE | 10K | 6,334 | 11,853 | 18,187 | 704 | 718 | **0.924726** |
| B6 | CR, 18,129 + 516 DEPRECATED | 10K | 5,450 | 14,877 | 20,327 | 154 | 2,308 | 0.879791 |
| CR native log | all 39,186 | | 5,773 | 12,650 | 18,423 | | | |

A0 reproduces the saved library replay exactly: identical 17,798-barcode callset, symmetric
difference zero, ledger identical at six significant digits.

Key pairwise Jaccards: A0 versus B2 **0.967682** (same universe, different counts); B1
versus B2 **0.803416** (same counts, different universe). The universe dominates.

### Which feature class matters

The two excluded classes act in opposite directions. Adding the 375 genes with active
`included=FALSE` probes improves every stage and brings the rescue count closest to CR's
logged 12,650 from below. Adding the deprecated rows drives large over-rescue in every
combination tested. Dropping the 20,118 zero-count genes changes nothing, so a
Simple Good-Turing unseen-mass dilution is **not** the mechanism: P0 is 0.00230 to 0.00247
in every arm and the per-unseen log probability is -14.24 to -14.29 except in the
full-39,186 arms.

### Harness arms with CR's own stage sets

To separate the model universe from stage selection, the stage-restriction replay harness
was rerun with CR's reconstructed primary set (5,773), eligible pool (38,749) and CR-count
ranked ambient window, varying only the universe used for ambient modelling and likelihoods.

| Arm | Model universe | Tail passes | vs CR 12,650 | Missed | Extra | Jaccard | With occupancy |
|---|---|---:|---:|---:|---:|---:|---:|
| C2 | 18,129 included | 11,082 | -1,568 | 1,662 | 344 | 0.891667 | 0.901672 |
| C1 | 18,129 + 375 incl=FALSE | **12,159** | **-491** | 733 | 492 | **0.934369** | **0.945550** |
| C3 | C1 + 516 DEPRECATED rows | 15,690 | +3,040 | 0 | 3,290 | 0.846713 | 0.857089 |
| C4 | C1 with DEPRECATED merged into parent genes | 16,519 | +3,869 | 0 | 4,119 | 0.815225 | 0.825245 |

C1 is the closest reconstruction of CR's tail behaviour obtained so far, within 3.9% of the
logged rescue count. Its common-tail cross-tab against the 12,495 official tail cells is
11,762 passing in both, 397 ours only, 733 CR only.

C3 and C4 recover every official cell, including all 733 C1 misses, but admit 3,000 to 4,000
barcodes CR rejects. The deprecated counts therefore carry real signal that CR uses in some
way we have not identified. The C4 merge mapped 538 of 542 deprecated rows to parent genes
(1,072,747 of 1,078,671 UMIs; largest `IGKC` +521,481), so the negative result is not a
mapping artefact.

**Inference, not measurement:** CR appears to gate barcodes on raw all-feature totals, which
include deprecated counts and match its logged stage numbers, while excluding deprecated rows
from the ambient model and likelihood. What distinguishes the 733 cells cannot be determined
from the evidence we are permitted to use.

## 3. Ruled out

- **Simulation count.** A0 to A1 changes 134 of 31,492 tail decisions, 103 lost and 31
  gained, with observed log-likelihoods bit-identical and every flip within 0.0011 of the
  BH boundary. Median absolute p change is 0.000527, the binomial standard error at 10,000
  simulations. Jaccard falls from 0.892115 to 0.889521. On the CR full-feature input, 100,000
  simulations add 229 passes and one official cell. The earlier inference from CR's logged
  minimum p-values does not hold; 10,000 simulations are adequate.
- **STAR's counting.** A0 versus B2 reach Jaccard 0.967682 on the same universe, with STAR
  marginally closer to the official callset than CR's own counts. Combined with the earlier
  read-level audit, count differences are not the source of the divergence.
- **Dead features and SGT dilution.** B4 drops all zero-count features and moves further from
  CR, not closer.
- **BH context.** Saved q-values already include primaries at p=0, which is the lenient
  choice; the tail-only alternative passes fewer barcodes.

## 4. Recommendations

Ranked by measured effect. None has been implemented; all require the usual validation.

1. **Reinstate the GEM-occupancy filter in the tag-aware path.** Fit lambda by zero-truncated
   Poisson MLE over GEMs with at least one call, take k as the smallest integer with Poisson
   CDF at least 0.999, and drop every call in GEMs with occupancy above k. Compute occupancy
   across all tags of the run, preserving the 2026-03-19 fix that restricts the denominator to
   observed tags. Worth +0.012 pooled across the eight L004 samples. Re-validate on the 4-sample
   JAX library, where the historical effect was 21 cells, to confirm no regression at low
   occupancy. A reference implementation with all three lambda variants is in
   `apply_occupancy_filter.py` in the analysis directory below.
2. **Widen the calling universe to the 375 genes with active `included=FALSE` probes** while
   continuing to export the 18,129-feature filtered matrix as CR does. This is the largest
   single lever on lymph node, +0.042 with CR's stage sets. It requires quantification to stop
   discarding `included=FALSE` probes.
3. **Exclude deprecated probe features from the ambient model and likelihood**, whether as
   separate rows or merged into parents. Both alternatives over-rescue badly.
4. **Keep the fixed ambient window** `[90,000, 180,000)`, scaled by tag count and never
   extended to a UMI target. Independently corroborated: CR logs a Simple Good-Turing slope of
   -2.1749 for lymph node; the fixed window on CR-count ranks gives -2.1765 and on STAR ranks
   -2.1670, while the adaptive 1,000,000-UMI window gives -2.0208.
5. **Do not raise the simulation count and do not tune OrdMag toward N=5,718.** The lymph-node
   loss surface is multimodal, and with the universe corrected the estimator should reach CR's
   basin without a fitted constant.

### Open questions and the next experiment

The decisive untested arm is a **two-universe design driven from STAR's own counts**: gate
primaries and candidates on totals that include the deprecated-probe equivalents, model
ambient and likelihoods on the 18,504-feature universe, then apply the occupancy filter.
Every arm above either used CR's reconstructed stage sets or a single universe for both roles.
This arm would also need a decision on whether STAR can or should count deprecated probes at
all, which is a quantification question, not a calling question.

Secondary: identify what distinguishes the 733 official cells that only deprecated-probe
signal recovers. This may not be answerable within the clean-room boundary.

Scope limits to respect. The universe and harness arms cover **lymph node only**, one library,
one depth. The occupancy result covers all eight samples but is a post-hoc filter on saved
callsets, and CR applies the same rule to its own pre-filter calls rather than to ours. No
result here justifies a production default on its own.

## 5. Artifacts and provenance

Local analysis root `/mnt/pikachu/star_suite_paper/analysis/`. That tree is untracked by
convention and is 2.4 GB; the reports below carry every hash, argv and completion record.

| Directory | Contents |
|---|---|
| `l004_gem_occupancy_filter_20260909/` | `README.md`, `apply_occupancy_filter.py`, `crosstab_occupancy.py`, `run_all.sh`, `summarize.py`, `results/` with per-arm filtered callsets and tables |
| `l004_feature_universe_sims_20260909/` | `README.md` in four sections (CLI arms, follow-up classes, harness universes, C4 merge), `results_table.tsv`, `pairwise.json`, `flip_summary.json`, `harness_table.tsv`, all SSM launch/execution JSON, runners, `SHA256SUMS` over 367 files |

Remote execution used AWS instance `i-06de289faa5d78117`, region `us-west-2`, profile `uw`.
Every arm ran once, serially, in a fresh directory, launched detached so no SSM timeout could
kill a caller, with exit status 0 recorded.

- Roots: `/scratch/l004_feature_universe_sims_20260909_v1` (SSM
  `7cf5bbd6-6af8-4db0-9f12-5b5a9454643d`), `..._followup_v1` (`6919ccd7-71c2-4794-9537-4386cb04d7d4`),
  `/scratch/l004_harness_model_universe_20260909_v1` (`4b978c49-d967-4f98-80fc-3a70d6b03be0`),
  `..._c4_v1` (`5fe9644a-ecc0-4a8e-a4a7-936611d66168`).
- Caller: `/scratch/tools/scrna_simpleed-internal-quality-dd14a4e1906e`, SHA-256
  `dd14a4e1906e6f6dac0bac0709158d4463c9684f0be7e2d8cc8ac6a06554725f`, verified before each use.
- Harness: `/scratch/l004_stage_restriction_replay_20260909_v2/replay_stages`, SHA-256
  `e760e2af8799bcd0363ebcbed1b85db8c3cc03166332d47202827d855e483c94`, unchanged and not rebuilt.
  Its CR stage sets were reused byte for byte.
- Restricted matrices were verified per barcode against the audit vectors in
  `l004_ordmag_parameter_audit_20260909/corrected_feature_totals/`, with zero mismatches over
  726,968 barcodes.
- Result archives are on the approved bucket's `analysis-tools/` prefix, downloaded and
  hash-verified locally; the URIs and SHA-256 values are in each report.

## 6. Boundaries

Do not inspect Cell Ranger source, and do not use the excluded `cr_native_mex_calling_20260908`
material. Preserve full `CB16+TAG8` identity and normalize only a terminal `-1`; CB16
truncation conflates droplets across tags. Read `AGENTS.md` and `AGENTS.local.md` before any
new execution, serialize jobs, use fresh output directories, and take completion from the
wrapper's explicit exit status. An absolutely identical rerun requires explicit user
authorization. Do not commit, reset or overwrite unrelated dirty work in either repository.
