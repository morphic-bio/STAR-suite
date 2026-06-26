# Handoff: nf-core parity harness + tightening TranscriptVB toward Salmon

**Date:** 2026-06-25. **Audience:** an agent tasked with closing the small remaining
gap between STAR Suite's integrated quant (TranscriptVB) and Salmon, by reading
Salmon's source. **Status:** the comparison harness is built and validated; the gap is
precisely isolated to the quant algorithm. This doc hands off both the harness and the
concrete improvement task.

## TL;DR

- Harness: `STAR-suite/benchmarks/nfcore_rnaseq_compare/` — `./run_compare.sh --mode chr22`
  (fast, cached chr22 fixture) compares STAR Suite vs the nf-core `star_salmon` chain
  on the same reads/reference and emits/compares everything. README has usage.
- **The alignment is already identical.** "Exact mode" (Salmon run on STAR Suite's
  transcriptome BAM) reproduces the reference Salmon at **Pearson 1.0000** transcript
  and gene. So nothing in alignment/trimming needs work.
- **The whole gap is TranscriptVB vs Salmon on the *same* BAM:** NumReads
  Pearson **0.9916** (transcript), **0.9999** (gene). The task is to make TranscriptVB's
  quant match Salmon's, given identical input alignments.
- This is an **enhancement, not a blocker**: the harness already offers a pluggable
  "BAM → your Salmon" mode that is exact, so nf-core adoption does not depend on closing
  this. Closing it just lets the integrated one-binary path be numerically Salmon-equal.

## The harness (so you can iterate)

- `run_compare.sh --mode chr22|full` — builds reference, builds/【caches】a chr22-filtered
  read fixture (`fixtures/chr22/`), runs (A) Trim Galore→STAR 2.7.11b→Salmon, (B) STAR
  Suite integrated (`--trimCutadapt Yes --quantMode TranscriptomeSAM TranscriptVB
  --quantVBLibType A`), and (B') Salmon on STAR Suite's transcriptome BAM. Version guard
  asserts upstream STAR == STAR Suite `--upstream-version` (2.7.11b == nf-core's pin).
- `compare_outputs.py` — parity-first report (`compare/report.md`).
- **Iterate fast:** the chr22 fixture is cached, so re-runs are quick. After changing
  TranscriptVB, re-run and watch **integrated transcript Pearson climb toward 1.0**
  (gene is already 0.9999; do not regress it, and do not touch exact mode = the
  alignment proof).

## The problem, precisely isolated

`B_starsuite/Aligned.toTranscriptome.out.bam` is fed to BOTH:
- STAR Suite TranscriptVB → `B_starsuite/quant.sf`  (to improve)
- Salmon → `B2_starsuite_bam_salmon/quant.sf`       (the target; == reference Salmon)

So comparing those two `quant.sf` files is a **pure quant-algorithm diff on identical
alignments** — no alignment confound. Validated chr22 result:

```
TranscriptVB vs Salmon on the IDENTICAL transcriptome BAM:
  transcripts shared 4803, expressed 796
  NumReads Pearson           0.9916
  EffectiveLength Pearson    1.0000  but mean |EL_vb - EL_sa| = 83.4 bp  <-- systematic offset
  total assigned reads       TranscriptVB 12580  vs  Salmon 12924  (-344)  <-- TVB drops reads
```

Diagnostic of the worst-divergent transcripts:
`runs/chr22_20260625_154221/compare/transcriptvb_vs_salmon_sameBAM_topdiff.tsv`
(regenerate for any run with the snippet at the bottom of this doc).

## The integrated gap is the KNOWN auto-detect detection-window tradeoff — already decided

**Do NOT send an agent to read Salmon's code to "fix the algorithm." The algorithm
already matches Salmon; this was characterized and decided previously.** See
`STAR-suite/core/features/vbem/docs/STAR_Salmon_Quantification_Parity_Report.md`.

What the integrated 0.99 actually is: with `--quantVBLibType A`, STAR Suite runs a
separate library-format detection pass over the first `quantVBAutoDetectWindow` reads
(default **1000**, `Parameters.h:483`; `STAR.cpp:1817-1868`) that votes on the format but
does not count those reads (no rewind afterward). Salmon infers the type without
discarding fragments, so it counts them. That is the entire gap — a **small-data sampling
artifact of the detection window**, not a quant-algorithm difference. The prior report's
findings:

- **Fixed `--quantVBLibType ISF` → read gap = 0%**, and multi-mapper resolution then
  matches Salmon exactly (reads ∝ 1/effLen). I.e. the algorithm is already correct.
- The auto-detect gap is **negligible at production scale**: 48k→6.2%, 1M→0.3%,
  **10M+→<0.1% (correlation >0.999)**.
- Documented recommendation: small datasets (<100k reads) → use fixed `--quantVBLibType
  ISF`; large datasets (>1M) → auto-detect is fine, gap negligible.

So my earlier "−344 dropped reads / ~83 bp EL offset / read Salmon's code" framing was
**re-discovering a settled, accepted tradeoff.** The chr22 fixture just amplifies it (1000
detection reads is a big fraction of ~13k chr22 fragments). At full scale, or with fixed
libtype, integrated parity is ~1.0. Exact mode (BAM→Salmon) is already 1.0 and gene-level
is 0.9999 — those reflect the true parity.

**Implication for this harness (the only real action item):** the chr22 integrated number
is honest-but-amplified. To present true parity on the small fixture, either run STAR
Suite with the detected fixed libtype (`--quantVBLibType ISF`) — expect ~zero gap — or
add a note citing the decision above. Consider adding a `--libtype {A,ISF}` switch to
`run_compare.sh` and reporting both so the detection-window effect is visible and
explained rather than mistaken for an algorithm gap.

**If (and only if) someone wants auto-detect itself to be exact even on tiny data:** the
change is to rewind the read files after detection (`closeReadsFiles()/openReadsFiles()`
after `STAR.cpp:1867`, suppressing detection-pass output to avoid the documented duplicate
output) so the window reads are re-counted. But the team explicitly judged this **not
worth it** vs. fixed-libtype / production-scale negligibility — so this is a deliberate
non-goal, recorded here only for completeness.

## Where TranscriptVB lives (STAR Suite side)

In the STAR-suite repo (`core/legacy/source/` and the suite's quant code). Flags:
`--quantMode TranscriptVB`, `--quantVBLibType A`, `--quantVBgcBias 1`. Find the TranscriptVB
implementation (grep `TranscriptVB`, `quantVB`) — that's the file to align to Salmon's
effective-length + acceptance + VBEM behaviour.

## How to validate progress

Re-run `./run_compare.sh --mode chr22` (reuses the cached fixture). Target:
- integrated transcript NumReads Pearson 0.9916 → approaching exact mode's 1.0;
- keep gene Pearson ≈ 0.9999+ and total-reads gap → 0;
- **do not change** exact mode (it must stay 1.0 — it's the alignment proof).
Then confirm on `--mode full` before claiming parity.

## Honesty / scope

- Perfect numerical equivalence with Salmon may be unattainable (Salmon has specific
  heuristics); "match to within Salmon's own run-to-run/version variation" is the realistic
  bar. The earlier note in the nf-core plan stays true: the integrated quant is ~0.99 and
  the **BAM → Salmon path is exact** — so this work improves the integrated path, it does
  not gate nf-core adoption.
- Keep changes scoped to TranscriptVB quant; don't perturb alignment (exact mode guards it).

## Appendix — regenerate the same-BAM diagnostic for any run dir

```python
import sys, numpy as np
R = sys.argv[1]
def load(p):
    d={}
    for i,l in enumerate(open(p)):
        x=l.rstrip().split('\t')
        if i==0: continue
        d[x[0]]=(float(x[3]),float(x[4]),float(x[2]))  # TPM, NumReads, EffLen
    return d
vb=load(f"{R}/B_starsuite/quant.sf"); sa=load(f"{R}/B2_starsuite_bam_salmon/quant.sf")
k=sorted(set(vb)&set(sa)); nv=np.array([vb[x][1] for x in k]); ns=np.array([sa[x][1] for x in k])
m=(nv>0)|(ns>0); print("NumReads Pearson", np.corrcoef(nv[m],ns[m])[0,1])
print("total VB", nv.sum(), "Salmon", ns.sum())
```
Run: `python3 thissnippet.py runs/<dir>`.
