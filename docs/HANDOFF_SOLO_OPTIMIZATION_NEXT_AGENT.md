# Solo optimization — handoff for next agent

**Date:** 2026-03-24  
**Branch:** `feature/solo-optimization-20260324-mastermerge`  
**Workspace tip (verify):** `git rev-parse HEAD` — was `735ed6e` when this file was written.

**Deep dive / history / benchmark tables:** [HANDOFF_SOLO_OPTIMIZATION_20260324.md](./HANDOFF_SOLO_OPTIMIZATION_20260324.md) (long document; start with sections on direct bridge v5–v8, **v9**, **v10**).  
**Bridge tuple + ambiguous aggregate redesign (2026-03-25):** [HANDOFF_SOLO_BRIDGE_TUPLE_AMBIG_REDESIGN_NEXT_AGENT_20260325.md](./HANDOFF_SOLO_BRIDGE_TUPLE_AMBIG_REDESIGN_NEXT_AGENT_20260325.md) — packed `(wlCb, umi24, gene16)` exact hash, per-key ambiguous accounting (no pin-replay rows), Bayes-resolved read attribution policy, worktree `/tmp/star-suite-v10-redesign-20260325`, full EBs2_2 benchmark root under `/storage/100K/ucsf_solo_bridge_redesign_20260325/`.  
**Repo agent rules:** [AGENTS.md](../AGENTS.md) at repo root (clean rebuild before debugging; do not touch `README.md` if it is dirty unless asked).

---

## What this branch is for

Solo / STARsolo **post-map** optimization (counting, UMI collapse, MEX), especially the **non-Flex** path. Flex mapping-phase work belongs in other handoffs.

---

## Binary spool downstream HASH_MEX (non-Flex, mapping unchanged)

**Environment:** `STAR_SOLO_BINARY_SPOOL=1` and `STAR_SOLO_BINARY_SPOOL_HASH_MEX=1`. Legacy control for A/B: `env -u STAR_SOLO_BINARY_SPOOL_HASH_MEX`. Do **not** set `STAR_SOLO_NONFLEX_HASH_BRIDGE` or `--soloInlineHashMode` for this track.

**Code:** `SoloFeature_collapseUMI_fromBinarySpoolHash.cpp` (spool scan + global hash) → `solo_bridge_collapse::runFromGlobalAggregates` in `SoloFeature_collapse_bridge_global.cpp` (shared with inline-hash bridge drain). `SoloFeature_countCBgeneUMI.cpp` selects the branch.

**Validated 2M pair (CSR v4, 2026-03-25):** exact `Summary.csv` + raw/filtered MEX triple parity vs legacy binary spool. Roots: `…/iPSC2_1_GEX_2M_binaryspool_{legacy,hashmex}_v4_csrval/` under `/storage/100K/ucsf_solo_spool_hash_20260325/`. Details and timing: long handoff **“Binary spool HASH_MEX: CSR global sort in `runFromGlobalAggregates` (v4, 2026-03-25)”** and `tests/ARTIFACTS.md` (`ucsf_solo_spool_hash_20260325`).

---

## What is already landed (direct-hash bridge, non-Flex)

**Environment:** `STAR_SOLO_NONFLEX_HASH_BRIDGE=1` enables the direct bridge from packed inline hash into collapse without the legacy `rGeneUMI` / `rCBp` replay.

**Constraints the path respects (do not break casually):**

- `Unique` + `1MM_CR` + `MultiGeneUMI_CR`, optional CR multimap rescue at alignment
- CB shard ownership: `wlCb % nShards`
- Bayesian ambiguous-CB resolution unchanged in intent
- No `countMatMult` on this path
- `outSAMtype None` is common in benchmarks

**v10 (current MultiGene hot path):** same outer architecture as v9; inside **`collapseOneBarcodeRows`**
(`SoloFeature_collapse_bridge_global.cpp`) nested `umiGeneMapCount` / `umiGeneMapCount0` maps and
**`umiCorrected`** are replaced by flat `corrUmiGeneCount` + `corrPackedKeys` + sparse
**`origByOrigUmi`** (see long handoff “Direct bridge v10”). The **inline-hash bridge** path uses the
same MultiGene **semantics** via per-CB **`MgRow`** sorting in
`SoloFeature_collapseUMI_fromBridgeHash.cpp` (no nested maps there either after 2026-03-25).

**v9:** same as v8 **except** no unused `bridgeSlotsByCompactCb_`; collapse uses
`unordered_map<wl, vector<slotId>>` for **observed** barcodes only (not `resize(cbWLsize)`). Tuple
hash → stable slot id; global fold + per-CB sort + **`collapseOneBarcodeRows`** (superseded in v10 for MultiGene aggregates only).

**v8 (prior checkpoint):** packed slots + (unused) per-compact-CB slot list; global fold into
`vector<vector>` sized to full whitelist.

**v7a (prior):** per-shard **`khash_t(cg_agg)`** + sorted **`vector<khiter_t>`**; no full-shard **`vector<ShardRow>`** before collapse.

**Key source files:**

| File | Role |
|------|------|
| `core/legacy/source/SoloFeature_collapseUMI_fromBridgeHash.cpp` | Direct bridge: CSR drain + per-CB sort + **flat `MgRow` MultiGeneUMI_CR** (no nested umi→gene maps; `umiArrayCorrect_CR` with `readInfoRec=false` when BAM tags unused) |
| `core/legacy/source/SoloFeature_collapse_bridge_global.cpp` | Shared `runFromGlobalAggregates` (bridge + binary-spool HASH_MEX) |
| `core/legacy/source/SoloFeature_collapseUMI_fromBinarySpoolHash.cpp` | `STAR_SOLO_BINARY_SPOOL_HASH_MEX` spool path |
| `core/legacy/source/SoloFeature_countCBgeneUMI.cpp` | Chooses bridge vs legacy vs HASH_MEX path |
| `core/legacy/source/SoloFeature_sumThreads.cpp` | Drains thread hashes + ambiguous merge |
| `core/legacy/source/SoloFeature.cpp` / `SoloFeature.h` | Wiring, declarations |
| `tests/ARTIFACTS.md` | Canonical benchmark output roots (v5, v6, v7a, legacy) |

---

## Validation you should trust

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

**2M UCSF fixture (documented):** see `tests/ARTIFACTS.md` — **v10** root (current):

`/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_postmerge_v10/`

**v9** (nested MultiGene maps, same outer bridge):  
`/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_postmerge_v9/`

**Clean paired validation roots (use these first for v9 vs v10 claims):**

- v9 paired control:
  `/storage/100K/ucsf_solo_optimization_20260324/validation_pair_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_paired_v9/`
- v10 paired candidate:
  `/storage/100K/ucsf_solo_optimization_20260324/validation_pair_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_paired_v10/`

These paired runs were rebuilt cleanly from the same commit with only the v9/v10 MultiGene hot-path
difference. They show that v10 is **not parity-identical** to v9:

- `Summary.csv`: different
- `Features.stats`: different (`yesCellBarcodes 41319 -> 41320`, `yesUMIs 1352870 -> 1352810`)
- `raw/matrix.mtx`: different, raw nnz `1255640 -> 1255585`
- `filtered/barcodes.tsv`: `7225 -> 7214`, with `11` v9-only barcodes and `0` v10-only
- `filtered/matrix.mtx`: different, filtered nnz `1206518 -> 1206089`

Compare timing, RSS, and outputs vs **v8**, **v5**, **v7a**, and **legacy** paths in the same doc. The
stored v9/v10 roots remain useful historically, but the paired roots above are the ones to trust for the
current v9/v10 validation question.

**Parity expectation:** direct-hash vs **stored** v5/v6 artifacts shows **small** raw-nnz / header drift (`cmp` often non-zero), but the clean paired v9/v10 rerun demonstrates that the v10 MultiGene rewrite itself currently changes outputs. **Legacy** matrix parity vs direct bridge is still an **open** product question, and **v9 vs v10 parity is also still open**.

**Important bugfix already applied after that investigation:** the ambiguous-CB resolver no longer depends on
an arbitrary representative `cbQual` string. Shared ambiguous entries now accumulate per-position
log-likelihood evidence (`cbLogLikMatch`, `cbLogLikMismatch`, `cbEvidenceReads`) and `CbBayesianResolver`
consumes that aggregated evidence when present. The same-binary 2M determinism check is now stable; see the
long handoff section “Ambiguous CB quality aggregation fix (2026-03-25)”.

**Flex note:** the same aggregated-evidence fix was also carried into `flex/source/InlineCBCorrection.*`
for merged ambiguous-CB handling, so do not reintroduce single-representative `cbQual` behavior there.

---

## Suggested next work (pick one track)

1. **Fix v10 vs v9 parity before a full benchmark:** the clean paired A/B still shows the v10 MultiGene rewrite is genuinely changing `Summary.csv`, `Features.stats`, raw MEX, and filtered cell selection vs v9. The ambiguous-CB quality-order bug is no longer the right explanation, so the next agent should debug the remaining v10-vs-v9 collapse delta directly.
2. **Memory:** v7a holds **khash + sorted iterator vectors** concurrently; if RSS spikes on larger runs, consider streaming destroy policies or a more compact CB index (see original design notes in the long handoff).
3. **Legacy parity:** if the goal is a **drop-in** replacement for non-bridge Solo, need a defined golden run and systematic `matrix.mtx` / filtered matrix comparison vs **legacy** (not only vs v5/v9).
4. **MEX / output:** `MexWriter` / `SoloFeature_outputResults` remain broad-cost surfaces per the long handoff “What matters most” section.

---

## Local git note

At handoff time, **`README.md` was modified** in the working tree; avoid committing or reverting it unless the task explicitly includes it.
