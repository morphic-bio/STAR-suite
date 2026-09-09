# Independent review handoff: L004 cell-caller divergence

## Current OrdMag findings and next investigation — 2026-09-09

**The new aggregation diagnostic does not resolve the STAR/CR discrepancy.**
For the worst-case L004 lymph-node sample (`LNReactive_BC9-10`), it estimates
**N=7,930 on STAR counts versus N=5,127 on CR counts**, despite small optimized
losses. Native CR logs **N=5,718 and loss=0.004181068568834495**. Reproducing
that native N/loss pair remains the calibration target; it has not been met.

Keep four quantities separate: **N** is the expected-cell parameter;
**F(N)** is a direct threshold count; the internal caller can subsequently
estimate a primary count by another bootstrap stage; the final count also
includes tail rescue and subsequent filtering.

| Lymph-node quantity | STAR counts, new diagnostic | CR counts, same diagnostic | Native CR |
|---|---:|---:|---:|
| Expected-cell estimate N | **7,930** | **5,127** | **5,718** |
| Mean minimized diagnostic loss / native logged loss | 0.003390498 | 0.003740071 | 0.004181069 |
| Standard deviation of replicate N | 6,491.26 | 2,822.49 | unavailable |
| Direct F(N) on original counts at the reported N | **7,873** | **5,219** | not established as the logged quantity |
| Primary count from native CR execution | — | — | **5,773** |

The earlier internal STAR caller estimated **N=6,659** and retained **6,912
primaries**. At the new N=7,930, direct STAR thresholding uses at least **1,455
UMIs** and retains 7,873 barcodes. That number is not a new internal-caller or
final-cell result. The new STAR diagnostic covers **LN only**; the CR diagnostic
covers all eight samples, with the other seven N estimates within 0.7% of native.

### What was tested

The public [CR9 OrdMag documentation](https://www.10xgenomics.com/support/software/cell-ranger/9.0/algorithms-overview/cr-gex-algorithm)
describes minimizing `(F(N)-N)^2/N`, where F counts barcodes above one tenth
of the 99th-percentile UMI count among the top N barcodes. The
[Flex documentation](https://www.10xgenomics.com/support/software/cell-ranger/9.0/algorithms-overview/cr-flex-frp-algorithm)
says cell calling is performed separately per sample. Those pages do not
fully specify the numerical grid, percentile convention, resampling, or
meaning of the logged estimate/loss pair.

A direct, unresampled test used every integer N=2..45,000, linear percentile
interpolation and strict UMI > m/10. Exact zero-loss solutions were:

- **CR counts:** 4,682; 5,057; 5,574; 5,608; 6,576; 6,616; 7,033.
- **STAR counts:** 4,615; 4,805; 4,806; 7,832; 7,842.

Neither list includes 5,718. At that N, the direct CR calculation gives
F=6,091 and loss=24.331759; STAR gives F=6,087 and loss=23.812697.
The user correctly required treating native CR's near-zero logged loss as
the reference and auditing our reconstruction before judging the estimator.
**These direct curves do not establish that CR selected an invalid solution
or that its full estimator has the same instability.**

The subsequent aggregation diagnostic uses all nonzero barcode totals,
100 multinomial bootstrap replicates, seed 20260909, and a rounded 2,000-point
logarithmic grid over 2..45,000 (1,279 unique integers). Each replicate uses
the same linear 99th percentile and strict threshold; its minimizing N and
minimum loss are then averaged separately. This is an explicit diagnostic
model, not a verified CR implementation. Bootstrap concerns estimation of N;
the agreed deterministic barcode quality-tie policy is unchanged.

### What the parameter audit ruled out, and what remains uncertain

- The raw CR cache exactly matches a fresh read of the output H5. Raw totals
  reproduce the native **38,749 eligible barcodes** and **1,765-UMI primary
  boundary**, favoring them for the observed calling stages.
- The **18,129 included-feature mask** agrees independently between probe
  CSV, filtered H5 feature IDs and raw H5 target-set indices. Applying it
  gives direct F=6,142 and loss=31.440364 at N=5,718. It does not fix the metric.
- Thirteen percentile conventions and three threshold conventions on seven
  feature masks do not recover the native pair. Dividing by N-squared looks
  close for LN alone but fails to explain the other seven logged losses.
- Logged `N * loss = 23.907350...` is not an integer square, so the displayed
  N/loss cannot both describe one exact integer evaluation as previously
  assumed. Aggregation or unreported numerical/logging details remain possible.
- Generic resampling demonstrates the distinction: one CR control has mean
  optimized loss **0.004374817**, yet direct loss on original counts at its
  averaged N is **2.622414**. This supports aggregation as an explanation for
  the metric mismatch, but does not identify native CR's exact procedure.

**Next priority:** establish the search/aggregation/logging conventions on
the CR counts using the logged N/loss pairs as controls, then validate the
same procedure on STAR counts. Check explanations across all eight CR samples;
do not tune a constant to LN alone or treat a small optimized loss as proof
of robust cell calling. Primary target estimation is upstream of gene/MT
barcode tiebreaks. Exact native LN N/loss and improved final-cell behavior
remain unresolved; no new production default is justified by these tests.

### Artifacts and execution status

- [Direct formula sweeps and objective plot](/mnt/pikachu/star_suite_paper/analysis/l004_direct_ordmag_formula_20260909/README.md):
  sixteen completed CR/STAR evaluations, full curves and all exact minima.
- [Parameter and aggregation audit](/mnt/pikachu/star_suite_paper/analysis/l004_ordmag_parameter_audit_20260909/README.md):
  native logs, input checks, explicit conventions and complete replicate
  results. Use `corrected_feature_totals/`; earlier target-set metadata was
  corrected and superseded, with the recovery documented in that report.
- [Latest STAR LN estimate](/mnt/pikachu/star_suite_paper/analysis/l004_star_aggregation_estimate_20260909/README.md):
  one new completed STAR estimate on all **707,278 nonzero barcodes** and
  **54,293,016 UMIs**; saved CR comparator reused. Includes compact `.npy`
  totals, all replicate results, direct primary identities and input hashes.

All completed diagnostic controls have successful wrapper evidence. These
latest experiments made **no production caller changes and no tail runs**;
final-cell Jaccards below remain the latest measured results. Retain the
accepted fixed ambient window `[90,000,180,000)` for subsequent tail controls.
Preserve full **24bp CB16+TAG8** identities. **Do not inspect CR source** or
the excluded `cr_native_mex_calling_20260908` investigation; use public docs,
reference data, output matrices and permitted logs. Do not repeat an absolutely
identical caller/benchmark execution without explicit user authorization.

## Earlier controls and final-cell results

**2026-09-09 matched native expected-input control:** [Result and interpretation](/mnt/pikachu/star_suite_paper/analysis/l004_ordmag_native_expected_control_20260909/README.md)
runs our unchanged OrdMag on the same CR LN totals, supplying CR's logged
expected input **5,718**. Our target is **5,889**, versus native CR **5,773**:
**+2.01%**, compared with **-19.05%** under the 48,000 median budget (4,673).
This identifies expected-input calibration as the leading explanation for
the large LN primary deficit; divisor 10 and percentile 0.99 were unchanged.
Mapping logged input 5,718 through the CR candidate-median weights corresponds
to budget 61,409.48, or 3,838.09/tag. This is a local calibration point, not an
inferred universal CR constant. Bootstrap/reference construction remains a
candidate for the much smaller residual. One new serial OrdMag-only control
completed with exit 0; no tail run, production edit or baseline rerun.

**2026-09-09 3,500-per-tag / 56,000-budget follow-up:** [Results and provenance](/mnt/pikachu/star_suite_paper/analysis/l004_median_budget_56000_20260909/README.md)
tests the user's higher OrdMag expected-input budget with the same candidate
median weights. LN final official-CR Jaccard improves **0.889545 -> 0.903060**
relative to the 48,000 budget (original automatic OrdMag with fixed ambient:
0.892115; earlier STAR-quality 5,773 restriction: 0.905032). New LN primary
counts are STAR **5,468**, CR-count **5,308**; final STAR calls 17,281, shared
16,824, missed 1,349, extras 457. Compared with 48,000 this adds 258 official
cells and 7 extras. Primary STAR/CR-count Jaccard is 0.941622, slightly below
48,000's 0.950714 but much above original 0.747156. All-eight new primary counts
are STAR 88,232 / CR-count 88,857; pooled Jaccard bounds 0.978714–0.979400.
Sixteen new OrdMag arms and one LN tail arm completed once, serially, exit 0.
Counts, ambient, candidate pool, likelihoods and common-tail raw p-values
remain unchanged. Final-cell scope is LN only. This is user-defined expected
input scaling, not proof of a CR normalization mechanism or a 3,500 constant.
No production default was adopted or source changed.

**2026-09-09 final-cell follow-up to median-budget OrdMag:** [Final results and transition audit](/mnt/pikachu/star_suite_paper/analysis/l004_median_budget_final_20260909/README.md)
tests the new 4,747 STAR primaries on the cached STAR LN matrix with the accepted
fixed ambient and unchanged tail policy. **Final official-CR Jaccard slips
0.892115 -> 0.889545** (17,016 calls; 16,566 shared; 1,607 missed; 450 extras).
The new set removes 394 official cells and 388 extras. Of 2,165 demoted
primaries, 1,446 are rescued and 719 rejected (338 official). Another 63
existing tail calls fail after BH (56 official) despite identical raw p-values.
All counts, ambient profile, likelihoods and 31,492 common-tail raw p-values
are unchanged. The earlier STAR-quality 5,773-primary restriction retains
higher final Jaccard 0.905032; moving from that restriction to 4,747 loses
321 official cells and 36 extras. This follow-up covers **LN only**, not all
eight final callsets. One new tail execution completed with exit 0; no OrdMag
or baseline rerun and no production source change. Improved cross-matrix
primary consistency does not establish improved final-cell accuracy.

**2026-09-09 median-budget OrdMag experiment:** [Results and provenance](/mnt/pikachu/star_suite_paper/analysis/l004_median_budget_ordmag_20260909/README.md)
applies the user's 48,000 expected-cell budget proportionally to each sample's
median UMI among all >=500-UMI barcodes, **for OrdMag only**. Sixteen distinct
serial runs completed with exit 0 using the unchanged STAR OrdMag library:
eight STAR quality inputs and eight CR-total controls. Lymph-node primary
Jaccard between our caller on STAR and CR counts improves **0.747156 -> 0.950714**;
primary counts change from 6,912/5,222 to **4,747/4,673**. Across all eight samples,
new primary Jaccards are at least 0.924600; pooled boundary-tie bounds are
0.980525–0.980927. Total STAR primaries change **127,324 -> 82,563**, versus
83,306 on CR counts. The 48,000 sum is an input budget, not a cap on primary
calls or a final-cell budget. Native CR remains a separate comparator: its
logged LN primary count is 5,773, and the reconstructed native-primary overlap
is essentially unchanged. No tail rescue or production caller change was
performed. Do not infer final-cell gains or losses from primary demotions.

**2026-09-09 parameter/count clarification:** [CR distribution audit and stage reconciliation](/mnt/pikachu/star_suite_paper/analysis/l004_umi_expected_cells_20260909/README.md)
verifies **253,873 final CR cells for L004** directly from all eight filtered
H5 matrices. The sum **128,053** is the sum of logged OrdMag `recovered_cells`
parameters, not the final cell total. The full-data comparator remains
325,410 cells. Keep the OrdMag parameter, initial calls, tail rescues and final
counts separate. CR raw total UMI shares and candidate-median shares do not
match its OrdMag parameter distribution; median-based weights are closer to
final-cell proportions than raw total-UMI weights by a descriptive eight-sample
distribution comparison. This audit itself was not a calling validation. The
subsequent user-selected median-budget OrdMag experiment is linked above;
no production source change has been made. Complete CR barcode
total caches for all eight samples, permitted log evidence and public-documentation
links are retained with the report.

**2026-09-09 OrdMag follow-up:** [OrdMag with ambient held fixed](/mnt/pikachu/star_suite_paper/analysis/l004_ordmag_fixed_ambient_20260909/README.md)
separates primary count from primary identity. With the same STAR counts,
38,404 eligible candidates, fixed [90,000,180,000) ambient window and FDR 0.01,
reducing primaries from 6,912 to 5,773 under STAR quality ranking improves
official-CR Jaccard 0.892115 -> 0.905032, removing 352 extras but also 73 official
cells. CR-ranked identities at the same count reach 0.906398. The 1,139
demoted STAR primaries yield 741 rescues and 398 rejections; another 27 existing
tail calls fail BH with unchanged raw p-values. The new projected-CR-count
OrdMag trace gives recovered estimate 6,210 and retained target 6,513, between
the saved full-feature CR-input target 5,222 and STAR-input target 6,912.
This is a count-estimation trace, not projected-matrix final calling. The
estimated primary count remains sensitive before quality tiebreaking. The
5,773 restriction is diagnostic, not an adopted target. No production source
changed. All new executions completed once, serially, with exit 0 and archived
ledgers/provenance. The user has retained the base ambient window and set
OrdMag as the next priority; candidate/FDR policy is unchanged.

**2026-09-09 follow-up:** [Stage restrictions on unchanged STAR counts](/mnt/pikachu/star_suite_paper/analysis/l004_stage_restriction_replay_20260909/README.md)
isolates a large ambient-selection effect using CR logs/results and the existing
STAR library. The baseline replay is exact. Changing only the ambient window
raises official-CR Jaccard from 0.833996 to 0.892115; cumulative reconstructed CR
primary/candidate restrictions plus the logged ambient window reach 0.911886.
A CR-count-ranked ambient reconstruction reaches 0.919928, with explicitly
unknown native CR ambient tie choices. No production caller change was made.
The native CR log reports 5,773 primaries; the 5,222 below refers to our caller
on CR counts. The linked report includes reconstruction limits, per-stage
ledgers, completion records, and a verified raw-count H5AD cache.

Written 2026-09-08 at the user's request for another agent to examine the
problem. This document was assembled by a subagent from saved results and
STAR-suite source. No caller, STAR execution, benchmark, or count correction was launched
to write it. The user wants to understand why very similar STAR/Cell Ranger
counts can produce substantially different cell calls, especially in lymph
node, and requested Jaccards at each stage.

The internal caller integration is complete and matches the external caller
exactly on the current STAR L004 counts. **The unresolved problem is caller
sensitivity to the input counts/features, including the OrdMag target and
EmptyDrops rescue.** The latest full-depth stage audit establishes that most
missing reference cells reach EmptyDrops and are rejected there. It does not
yet separate feature policy, ambient composition, likelihood calculation, and
multiple-testing effects as causes.

## Strongest findings to start from

- The full-depth L004 lymph-node candidate sets have Jaccard **0.991097**.
  Every STAR candidate is also a candidate when our caller uses the CR MEX.
  All 18,173 official CR cells occur in both nonzero input sets and both
  180,000-barcode OrdMag input sets.
- The primary sets have Jaccard **0.747156**, but this largely reflects
  different target sizes: 5,222 on CR counts versus 6,912 on STAR counts.
  They share 5,189 cells, or 99.37% of the CR-arm primary set. The first
  logged divergence is already in the bootstrapped `recovered_cells`
  estimate, before quality tiebreaking.
- Among **31,459 barcodes actually tail-tested in both arms**, 9,418 pass
  both at FDR 0.01, 2,969 pass only the CR-count arm, and 44 pass only the
  STAR-count arm. This remains a large rescue difference after removing the
  confounding effect of different primary memberships.
- All **2,969** CR-only common-tail rescues already have larger STAR **raw
  p-values**, with median per-barcode STAR/CR ratio **6.222233**. Their median
  raw p-values are 0.00239976 on CR counts and 0.0123988 on STAR counts. The
  difference therefore precedes BH adjustment; BH alone cannot explain it.
- Of the **2,448 official CR cells missed by STAR**, only **120** are outside
  the STAR candidate set; **2,328** are tail-tested and rejected. The 2,969
  CR-arm-only common-tail rescues contain **2,208 official CR cells**.
- For those 2,969 discordant rescues, median totals are 780 UMIs on all CR
  features, 766 on CR counts projected to STAR feature IDs, and 769 on STAR
  counts. Median STAR-versus-projected-CR relative total difference is
  **+0.3653%**. This supports the user's observation that totals are close.
  It does **not** establish identical per-gene vectors or ambient models.
- Projecting CR totals to STAR feature IDs makes initial barcode agreement
  even closer: nonzero-input Jaccard **0.999878**, top-180K Jaccard
  **0.993367**, and >=500-UMI Jaccard **0.998229**. This was a saved-count/rank
  calculation, **not a projected-matrix caller rerun**.

Primary evidence is the [stage audit summary](/mnt/pikachu/star_suite_paper/analysis/l004_barcode_stage_audit_20260908/results/summary.json),
[stage table](/mnt/pikachu/star_suite_paper/analysis/l004_barcode_stage_audit_20260908/results/stage_jaccards.tsv),
[discordant-tail count context](/mnt/pikachu/star_suite_paper/analysis/l004_barcode_stage_audit_20260908/tail_count_context.json),
[raw-p-value context](/mnt/pikachu/star_suite_paper/analysis/l004_barcode_stage_audit_20260908/discordant_probability_context.json),
and the root agent's [stage-audit report](/mnt/pikachu/star_suite_paper/analysis/l004_barcode_stage_audit_20260908/README.md).

## Exact comparators and scope

The stage audit concerns **full-depth L004, LNReactive_BC9-10**, a biological
sample with TAG8 sequences `ACAGTCTG` and `AGTGAGTG`. It compares:

| Arm | Counts | Caller |
|---|---|---|
| CR-count arm | Saved L004 CR MEX with its full feature content | Our external quality caller, SHA256 prefix `73d0dd68ab7f` |
| STAR-count arm | Current H1X2/no-align L004 STAR MEX, 18,129 exported feature IDs | Our integrated native/shared quality caller, STAR SHA256 prefix `3e896a3769e3` |
| Official reference | 18,173 saved CR L004 lymph-node cell barcodes | Actual CR output, used only as the reference callset |

The CR-count arm is **our caller on CR counts**, not CR's own internal OrdMag
or EmptyDrops diagnostics. In every stage table below, `CR` means this
CR-count arm unless explicitly described as the official reference.

Both caller arms use the new quality-ranking policy, 48 bootstrap workers,
8 Monte Carlo workers, 10,000 simulations, and FDR 0.01. For two tags, OrdMag
receives the top 180,000 nonzero barcodes; maximum expected cells is 45,000;
the ambient window starts at zero-based rank 90,000 and adapts to a target of
1,000,000 UMIs; the tested tail floor is 500 UMIs. The complete recorded argv
and diagnostics take precedence over this abbreviated description.

**Preserve full CB16+TAG8 identity.** Normalize only terminal `-1`. CB16
truncation conflates droplets across tags. Older parity tools and a historical
note in `AGENTS.local.md` suggest truncation; that is inappropriate for this
experiment and contrary to the established full-barcode comparison policy.

This full-depth **barcode-stage** audit is separate from the **read-level**
audit, which covers a shallow 400K-pair control (100K from each lane). There
is no current full-depth L004 per-read decision ledger. Do not extrapolate the
control's read concordance or molecule residual to all of L004.

## Stage-set results

Jaccard is intersection divided by union. These are our-caller-versus-our-caller
sets on different count matrices. The ambient set is a separate modeling
branch, not another subset of the candidate set.

| Stage | CR | STAR | Shared | CR only | STAR only | Jaccard |
|---|---:|---:|---:|---:|---:|---:|
| Nonzero input | 717,686 | 707,278 | 707,227 | 10,459 | 51 | 0.985357 |
| OrdMag input, top 180K | 180,000 | 180,000 | 176,548 | 3,452 | 3,452 | 0.962366 |
| Candidates | 38,749 | 38,404 | 38,404 | 345 | 0 | 0.991097 |
| OrdMag primary | 5,222 | 6,912 | 5,189 | 33 | 1,723 | 0.747156 |
| Tail tested | 33,527 | 31,492 | 31,459 | 2,068 | 33 | 0.937396 |
| Ambient droplets | 272,744 | 285,755 | 265,651 | 7,093 | 20,104 | 0.907129 |
| Tail passes, FDR 0.01 | 13,689 | 9,495 | 9,418 | 4,271 | 77 | 0.684149 |
| Final calls, FDR 0.01 | 18,911 | 16,407 | 15,820 | 3,091 | 587 | 0.811365 |

The candidate/primary cross-tab is particularly useful: 5,189 are primary in
both; 33 are CR primary/STAR tail; 1,723 are CR tail/STAR primary; 31,459 are
tail in both; and 345 are CR tail/STAR input but below candidate eligibility.
The marginal tail-pass Jaccard of 0.684149 includes this stage switching.
Restricting to the common tested tail gives passing-set Jaccard **0.757622**:

| Decision on the 31,459 common tail barcodes | Barcodes | Official CR cells |
|---|---:|---:|
| Both pass | 9,418 | 9,270 |
| CR-count arm only passes | 2,969 | 2,208 |
| STAR-count arm only passes | 44 | 0 |
| Both fail | 19,028 | 120 |

All official reference cells are present initially. Candidate eligibility
retains all 18,173 on CR counts and 18,053 on STAR counts. Consequently, the
large STAR deficit cannot be attributed primarily to absent starting barcodes
or the 500-UMI floor.

The previous headline **0.833996** compares STAR's 16,407 calls against the
18,173 official CR cells: shared 15,725, missed 2,448, extra 682. Our caller
on CR counts has **0.931660** against that same reference: 18,911 calls,
shared 17,886, missed 287, extra 1,025. These are different comparators from
the **0.811365** final stage-set Jaccard above.

### Where OrdMag first diverges

| Logged quantity | CR-count arm | STAR-count arm |
|---|---:|---:|
| Input barcodes | 180,000 | 180,000 |
| Estimated `recovered_cells` | 5,082 | 6,659 |
| Bootstrap mean primary count | 5,221.67 | 6,912.39 |
| Bootstrap SD | 660.937 | 644.628 |
| Exact rounded primary target | 5,222 | 6,912 |
| Total UMIs at primary cutoff | 1,865 | 1,569 |

The first bootstrap estimates `recovered_cells` from count-distribution fits;
a second bootstrap estimates retained count at the resulting baseline index.
Only then does the deterministic ranking select exactly the rounded target.
The tiebreak fix does not force these estimators to agree on different
matrices. Bootstrap sensitivity to count perturbations and worker count is
still unresolved. These mean/SD values are preserved in the
[reconciliation report](/mnt/pikachu/star_suite_paper/analysis/flex_internal_l004_20260908/CELL_JACCARD_RECONCILIATION.md)
with the underlying manifest/diagnostic paths.

### Feature projection and discordant totals

The stage audit streamed the saved 604 MB CR MEX once and cached per-barcode
totals both over all CR features and over exact STAR feature IDs. It reused
STAR's saved rank/candidate diagnostics. It did not reread the huge STAR MEX
or run a caller.

| CR totals projected to STAR feature IDs versus STAR | CR | STAR | Shared | Jaccard |
|---|---:|---:|---:|---:|
| Nonzero input | 707,250 | 707,278 | 707,221 | 0.999878 |
| Top 180K by total UMI/full barcode | 180,000 | 180,000 | 179,401 | 0.993367 |
| At least 500 UMIs | 38,336 | 38,404 | 38,336 | 0.998229 |

For the 2,969 CR-only common-tail rescues, 195 barcodes have **exactly equal
projected-CR and STAR total UMIs**. The median fraction of their CR counts
outside STAR feature IDs is 1.8519%. Equal totals can conceal different gene
composition, and an identical cell vector can receive different likelihoods
under different ambient profiles. Neither equality is yet checked here.

## What has already been implemented and validated

The accepted primary ranking is total UMIs descending, non-mitochondrial UMIs
descending when an annotation mask is supplied, detected non-MT genes
descending, then full barcode ascending. It retains the exact rounded target.
**No bootstrap is used for tiebreaks.** Original counts remain in bootstrap,
ambient modeling, likelihoods, and output MEX; MT counts are excluded only
from the optional quality scores. Equal-UMI ambient ranks use barcode identity
without quality scores. Bootstrap input counts are canonicalized before
sampling to remove dependence on matrix barcode order at fixed counts and
worker configuration.

Native Flex and external `scrna_simpleed --flex-tag-aware` now share
`runSimpleEDWithAmbient`. Native Flex groups tags sharing a biological sample
label, preserves CB24, exports one sample MEX, and saves bootstrap/rank/ambient
and candidate evidence. STAR defaults to `--soloFlexCellCaller tag-aware`
when Flex filtering runs; `legacy` remains an explicit option. Shared flags
include `--soloCellFilterMitochondrialGenes` and
`--soloCellFilterBootstrapThreads`.

Completed checks:

- Clean isolated default STAR and CLI builds; all five EmptyDrops unit tests
  passed. The working checkout executable was **not** replaced.
- Integration fixture v4 covers single- and two-tag biological samples,
  CB16 collisions across tags, MT/gene ties, and the 500-UMI boundary. Internal
  and external ledgers agree exactly (175 candidates/137 calls and
  350 candidates/274 calls); paired results also match the frozen quality
  caller. Duplicate tag assignment is rejected.
- Full L004 processed **1,823,648,323 read pairs**, STAR exit 0, eight sample
  MEX outputs. Every output callset reconstructs exactly from diagnostics;
  every quality tie at the primary cutoff has size one.
- External replay of current STAR lymph-node counts agrees with native
  results across all **38,404 candidates and 16,407 calls**, with no changed
  identities, numerical fields, primary flags, or indices; diagnostic JSON
  also agrees. Integration is not a demonstrated source of divergence.
- Earlier quality tests made barcode-axis reversal change **zero** decisions;
  the old caller changed 9,490 under that permutation. However, the shallow
  colorectal count perturbation still changes 4,561 decisions through the
  estimator. The tiebreak fix did not resolve that separate non-robustness.

The initial full-L004 wrapper wrongly required every sample to have a nonempty
tail and failed after STAR succeeded. Three samples validly have primary-only
calls. A corrected **read-only** validator and postprocessing succeeded; STAR
ran once. Earlier synthetic fixtures v1-v3 hit an existing Simple Good-Turing
requirement for five distinct nonzero count frequencies; v4 exercises the
sampler successfully. Neither issue was hidden or interpreted as a successful
failed workload.

See the [integration report](/mnt/pikachu/star_suite_paper/analysis/flex_internal_l004_20260908/README.md),
[quality-ranking report](/mnt/pikachu/star_suite_paper/analysis/ordmag_quality_ties_20260908/README.md),
and [protocol source documentation](/mnt/pikachu/STAR-suite/docs/ORDMAG_QUALITY_RANKING.md).

## Earlier favorable results and the read audit: keep these separate

The remembered >0.9 current STAR lymph-node result is real: the full STAR MEX
has Jaccard **0.903207 at FDR 0.04** (31,174 calls), versus **0.782993 at 0.01**
(26,183 calls). The quality change preserved these calls. Current L004 has
**0.893341 at 0.04** (19,016 calls) versus **0.833996 at 0.01** (16,407).
These FDR 0.04 results are thresholding saved q-values, not new caller runs or
an adopted production default. Comparing >0.9 full-depth/0.04 directly with
L004/0.01 changes both depth and threshold. On L004 CR counts, increasing to
0.04 overcalls (23,458 calls; reference Jaccard 0.765971), so the effect is
matrix-dependent.
At matched FDR 0.04, final arm-to-arm Jaccard is only **0.778196**; raising
FDR does not restore equivalence between the two count inputs.

Historical controls with the older favorable prototype (`0ca77536eba2...`,
two bootstrap workers, FDR 0.01) changed **only CR feature content** and gave
32,104 calls/Jaccard 0.959803 with full CR features; 24,934/0.738070 after
restricting to STAR feature IDs; and 30,253/0.901575 after adding deprecated
features back. This establishes feature-content sensitivity without any STAR
read processing. It does not quantify its exact current L004 effect or prove
that exporting deprecated features is the correct fix. That older prototype's
ambient-feature and tail-only BH changes were not restored by integration.
The [reconciliation report](/mnt/pikachu/star_suite_paper/analysis/flex_internal_l004_20260908/CELL_JACCARD_RECONCILIATION.md)
records the historical evidence locations and the precise comparison caveats.

The independent [read-to-molecule audit](/mnt/pikachu/star_suite_paper/analysis/l004_read_count_audit_20260908/README.md)
compares the first 100K pairs per lane, 400K total, on included nondeprecated
2024-A features. Assigned-read Jaccard is 0.965851 overall and 0.963927 for the
100K L004 subset. All comparable full barcodes agree; all retained STAR genes
are compatible with CR assignment (one CR assignment has two gene IDs).
Ambiguous/missing CR molecule keys are explicitly excluded from exact-key
comparisons. Read assignment is distinguished from final counted molecules.

Production-source UMI replay reproduces every control STAR MEX coordinate and
count. It explains 327,320 STAR versus 326,998 CR molecules as +343 from tied
UMI read-count components, +33 other STAR-only molecules, and -54 unsupported
CR molecules with STAR `UNMATCHED_TAG`. The active Flex bucket path invokes
`UMICorrector::correctClique` despite the `--soloUMIdedup 1MM_CR` flag. All
404 multi-raw CR families have the same raw read support in both assignments;
63 merge in STAR and 341 tied components remain split. This correction issue
is **not** the primary-cell ranking tiebreak policy. It remains unfixed, and
its effect on full-depth cell calling has not been measured.

High read-set Jaccard and low cell-set Jaccard use different units and stages.
The former weights individual assigned reads, while the latter weights each
barcode's final binary decision. The latest audit gives a more specific
finding than that general observation: nearly identical candidate membership
and close shared-feature totals coexist with systematically different tail
scores. The exact source of that score shift remains the next question.

## Source entry points for review

Use the frozen source snapshot for the tested behavior; the main checkout is
dirty and includes unrelated pre-existing work.

| Question | Entry point |
|---|---|
| Recovered-cell search, two bootstrap stages, exact retained target | [OrdMagStage.cpp](/mnt/pikachu/STAR-suite/core/features/libscrna/src/OrdMagStage.cpp:55) |
| Deterministic quality comparator | [OrdMagRank.h](/mnt/pikachu/STAR-suite/core/features/libscrna/include/OrdMagRank.h:44) |
| Retain ranks, adaptive ambient window, matrix assembly, caller invocation | [SimpleEDCaller.cpp](/mnt/pikachu/STAR-suite/core/features/libscrna/src/SimpleEDCaller.cpp:120) |
| Ambient profile, likelihood, CR Monte Carlo sampler, BH adjustment | [EmptyDropsMultinomial.cpp](/mnt/pikachu/STAR-suite/core/features/libscrna/src/EmptyDropsMultinomial.cpp:26), adjacent `EmptyDropsCRSampler.cpp`, and [CRLogProb.h](/mnt/pikachu/STAR-suite/core/features/libscrna/include/CRLogProb.h) |
| Biological sample grouping and tag-scaled settings | [FlexFilterTagAware.cpp](/mnt/pikachu/STAR-suite/flex/source/libflex/FlexFilterTagAware.cpp:86) |
| Active UMI correction dispatch and tied-component rule | [bucket collapse](/mnt/pikachu/STAR-suite/flex/source/SoloFeature_collapseUMI_fromBuckets.cpp:278), [UMICorrector.cpp](/mnt/pikachu/STAR-suite/flex/source/UMICorrector.cpp:68) |

The current `applyFDR` operates on the supplied result vector, including primary
entries assigned p=0; the older favorable prototype's tail-only BH policy is
a distinct historical experiment. Inspect raw p-values as well as q-values
before attributing differences to BH. Different primary sets alter that
multiple-testing context as well as which cells are simulated.

## Artifacts, provenance, and reconstruction details

Local analysis root is `/mnt/pikachu/star_suite_paper/analysis/`.

| Directory | Contents to inspect first |
|---|---|
| `l004_barcode_stage_audit_20260908/` | `results/summary.json`, `stage_jaccards.tsv`, `barcode_membership.tsv.gz`, `candidate_decisions.tsv.gz`, `cr_barcode_totals.tsv.gz`, `manifest.json`; `tail_count_context.json`, `summarize_discordant_tail.py`, `run_stage_audit.py`, SSM payload/execution records, `archive.json` |
| `flex_internal_l004_20260908/` | `README.md`, `CELL_JACCARD_RECONCILIATION.md`, `results/diagnostics/LNReactive_BC9-10/`, `results/external_ln/parity.json`, `source/`, `source_sha256.json`, `integration.patch`, `binaries.json`, `validated_build_status.json`, `BENCHMARK_TODO.md` |
| `ordmag_quality_ties_20260908/` | Previous external quality-caller source/build; `results/l004/mt_genes/manifest.json`, caller summaries, `results/l004_comparison.json`, `results/full_star_comparison.json`, annotation provenance |
| `l004_read_count_audit_20260908/` | `README.md`, exact per-read ledger, molecule-family decomposition, production UMI replay, coordinate validation and manifests |
| `flex_caller_recovery_20260908/` | Preserved historical prototype and feature-content controls, with original command/time/comparison files |

The stage audit's older CR ledger has an all-zero `is_simple_cell` diagnostic
bug. Its primary set was reconstructed as the p=q=0 prefix and independently
checked: exactly 5,222 records, exactly all CR barcodes with >=1,865 UMIs, no
simulated tail p=0, and exact reconstruction of the saved final callset.
Native diagnostic flags are fixed. CR retain/ambient ranks use total UMI
descending and full barcode ascending; the saved adaptive CR ambient interval
is zero-based `[90000, 362744)`. STAR uses its saved explicit `is_ambient`
flags. Membership masks use the bit map in `summary.json`; do not assume bit
positions without reading it.

STAR nonzero input membership was reconstructed from its observed-barcode MEX
axis for the two tags and checked against the native `input_cells=707278`
diagnostic. Candidate identities/counts/primary flags were cross-checked
against the rank ledger. The stage audit records hashes of its exact matrix,
feature/barcode, rank, caller-ledger and configuration inputs in `manifest.json`.

Remote provenance:

- AWS instance `i-06de289faa5d78117`, profile `uw`, region `us-west-2`, 48 CPUs.
- Full L004 root `/scratch/flex_internal_l004_20260908_v1`; STAR raw MEX under
  `star/Solo.out/Gene/raw`; diagnostics under `diagnostics/<sample>`.
- CR LN raw MEX `/scratch/l004_diff_20260906_v1/cr_mex_recall/LNReactive_BC9-10/raw`.
- Saved external CR caller `/scratch/ordmag_quality_ties_20260908_v1/l004/mt_genes`.
- Stage audit `/scratch/l004_barcode_stage_audit_20260908_v1`; SSM
  `f7f18a58-06a1-4fa6-9fe0-cd057b7e790a`, **Success, response code 0**.
- Stage archive: 9,228,471 bytes, SHA256
  `999132fa85f19e049ceb70f3127c97885b49a278ffefc7f82d645491ed19ea3d`,
  downloaded and verified locally. S3 URI is in `archive.json`, under the
  existing benchmark bucket's `analysis-tools/` prefix. Later local
  tail/probability context JSON files and their analysis scripts are
  supplemental to that archive.
- Base source HEAD `f13cb19d9e4f40c2ee3a657b0706ac4be08dcd9d` plus saved changes.
  Clean isolated build `/tmp/flex_internal_l004_20260908_v1`.
- Validated STAR SHA256
  `3e896a3769e3874c7e556ebb8ac1241a4fd263ed9d398d3506c033a9cf101cbc`;
  remote `/scratch/tools/STAR-internal-quality-3e896a3769e3`.
- Shared external CLI SHA256
  `dd14a4e1906e6f6dac0bac0709158d4463c9684f0be7e2d8cc8ac6a06554725f`;
  prior CR quality caller SHA256
  `73d0dd68ab7f97e0e94fe2f7cdc6e8b6031b90ff4c6a7829dff0339396172e24`.
- Full L004 SSM `2c1f4d6a-20a8-4cb8-89ef-3f1bb1e8ab0f`; corrected
  postprocessing `2fc8141d-43a2-4c06-bbdf-b0e133fbe6a4` succeeded.
  `completion.json` records actual STAR exit 0. Diagnostic archive SHA256
  `121d6524b37460bc911cccb5f77ce87167d4652cd68c1e1ba02cbff032c70c4a`;
  `diagnostics_archive.json` preserves the full transfer manifest.

The current full-L004 run already uses `noGenome`/`noAlign`: it dispatches
before `Genome::genomeLoad()` and logs 48 fused no-genome threads. An earlier
statement attributing startup to genome loading was corrected. Before a future
performance benchmark, verify this path stays active and separate cache/setup,
read processing, counting and calling times. This diagnostic run was not a
performance benchmark.

## Read-only inspection commands

These read downloaded artifacts, recompute set metrics, and print a few
discordant examples. They do not rerun a caller or access the cloud.

```bash
python3 - <<'PY'
from pathlib import Path
import csv, gzip, json
p = Path('/mnt/pikachu/star_suite_paper/analysis/l004_barcode_stage_audit_20260908')
s = json.loads((p/'results/summary.json').read_text())
sets = {name: (set(), set()) for name in s['stage_mask_bits']}
with gzip.open(p/'results/barcode_membership.tsv.gz', 'rt') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        for name, bit in s['stage_mask_bits'].items():
            for i, key in enumerate(('cr_stage_mask', 'star_stage_mask')):
                if int(r[key]) & (1 << bit):
                    sets[name][i].add(r['barcode'])
for name, (a, b) in sets.items():
    m = s['stages'][name]
    assert (len(a), len(b), len(a & b)) == (m['cr'], m['star'], m['shared'])
    assert abs(len(a & b)/len(a | b) - m['jaccard']) < 1e-12
    print(name, len(a), len(b), len(a & b), f"J={m['jaccard']:.9f}")
n = 0
with gzip.open(p/'results/candidate_decisions.tsv.gz', 'rt') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['cr_stage'] == r['star_stage'] == 'tail':
            if float(r['cr_q']) <= .01 < float(r['star_q']):
                if n < 5:
                    print('CR-only common-tail rescue:', r)
                n += 1
assert n == 2969
print('CR-only common-tail rescues:', n)
PY
```

`run_stage_audit.py` preserves the original reconstruction but creates a fixed
remote output directory and uploads an archive. It is provenance, **not a
command to rerun blindly**. The compact gzip TSVs already avoid repeated large
MEX reads for barcode/total/p-value questions. Persistent sparse binary matrix
storage (H5AD or equivalent) was suggested by the user but is still deferred.

## Discriminating next analyses — proposed, not executed here

1. Join the 2,969 CR-only common-tail rescues to full gene vectors on identical
   feature IDs. Quantify exact vector equality and feature-wise residuals,
   including the 195 barcodes with identical projected totals. The larger
   STAR raw p-values are established; locate their source without treating
   total-UMI equality as vector equality.
2. Compare ambient profiles on identical feature IDs: gene support, ambient
   mass, smoothing/unseen-feature mass, and probability changes for genes
   carried by discordant cells. Ambient membership Jaccard alone does not
   measure profile similarity. Use the saved sets to separate changes in
   which droplets enter from changes in counts within the same droplets.
3. Evaluate a matched-feature CR caller arm with the current shared source and
   identical depth/settings/worker counts. Save both bootstrap estimates,
   primary membership, ambient profile, raw p and adjusted p. The historical
   feature intervention and the current projected totals are not substitutes
   for this still-unperformed current L004 caller experiment.
4. Use controlled component substitutions to distinguish count-vector effects
   from ambient-model and BH effects: identical common-tail candidates and
   feature axis, then fixed primary membership and/or fixed ambient profile.
   Such interventions would be diagnostic experiments, not production changes.
5. Trace the recovered-cell objective/loss curve and bootstrap sample outcomes
   on matched counts. The quality rule removes rank ambiguity; it leaves the
   estimator's response to small count changes and worker-dependent sampling
   as separate robustness questions.

Before any new execution, read repository `AGENTS.md` and `AGENTS.local.md`.
Serialize benchmark jobs, use fresh output directories, and obtain completion
from the wrapper's explicit result/exit. An absolutely identical execution
requires explicit user repeat authorization; changing only logging or output
paths does not make a new run. A changed source/binary validation within the
authorized fix scope is distinct. Clean rebuild before debugging a crash or
regression. Do not commit, reset, revert, or overwrite unrelated dirty work.
The historical AGENTS technical claims that Flex always uses fixed 3,000
expected cells/raw p-values are superseded by this validated tag-aware path;
they do not describe the current experiment.
