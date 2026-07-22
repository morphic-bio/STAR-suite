# Handoff: non-Flex bridge tuple + ambiguous-CB aggregate redesign

**Date:** 2026-03-25  
**Audience:** next coding agent on Solo **non-Flex** `STAR_SOLO_NONFLEX_HASH_BRIDGE` + `--soloInlineHashMode yes`.

---

## Where the code lives (important)

Primary implementation was done in a **detached worktree**, not necessarily merged to your current branch:

- **Worktree:** `/tmp/star-suite-v10-redesign-20260325`  
- **Stated base:** `735ed6e` (clean v9/v10 validation base)  
- **Before debugging crashes:** `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR` ([AGENTS.md](../AGENTS.md)).

If `/mnt/pikachu/STAR-suite` does not yet contain these edits, **diff or merge from the worktree** before continuing.

---

## What problem this slice solved

1. **Exact-read bridge path:** use a single packed key `(wlCb, umi24, gene16)` in `inlineHash_` with count as value; **do not** use per-thread CB/gene compaction maps or `bridgePackedSlots_`-style slot payloads for that path.

2. **Ambiguous-CB path:** stop paying for **per-read pin-replay** vectors (`bridgeAmbigPinFlat_`, `bridgePinOrphans_`, row-at-a-time replay in `applyAmbiguousBridgePinReplay`). Replace with **per ambiguous key**:
   - `(umi,gene) -> count` in `bridgeAmbigUmiGene_` (kept),
   - aggregated gene **unique / multi-feature** counts + **sample** `readFlag` templates for read-stats bulk update,
   - aggregated **readInfo** counts + sample flag (either folded into `pendingAmbiguous_` or held in **`bridgeAmbigReadInfoOrphan_`** until first gene hit — preserves legacy “orphan until splice” lifecycle **without** per-read rows),
   - **`bridgeAmbigPinCandQuals_`** aligned with `candidateIdx` for **one** key-level pin evaluation (refreshed when representative `cbQual` wins).

3. **Merge correctness:** `mergePendingAmbiguous` must **sum** aggregates, merge `bridgeAmbigUmiGene_`, apply deterministic best-`cbQual` + pin-qual sync, and merge **`bridgeAmbigReadInfoOrphan_`** (including folding orphan into pending when the same key exists post-merge). No order-sensitive orphan vector concatenation bug.

4. **Resolution accounting (product decision, locked for recent benchmarks):** when Bayesian resolution **succeeds** and pin passes, attribute ambiguous **read / nRead / readStats** aggregates to the **resolved whitelist CB**, not only the pin winner. This is **intentionally not** legacy pin-only parity.

---

## Key source files (worktree)

| File | Role |
|------|------|
| `flex/source/hash_shims_cpp_compat.h` | `packBridgeWlUmiGeneKey` / unpack / guards |
| `core/legacy/source/SoloReadFeature.h` | `ExtendedAmbiguousEntry` aggregate fields; `bridgeAmbigReadInfoOrphan_` |
| `core/legacy/source/SoloReadFeature.cpp` | `mergePendingAmbiguous`; `applyBridgeAmbiguousAggregatedReadAccounting`; pin eval helper |
| `core/legacy/source/SoloReadFeature_record_base.cpp` | Exact hash insert; `accumulateBridgeAmbiguousCB`; gene/readInfo aggregates; orphan merge on first gene |
| `core/legacy/source/SoloFeature.cpp` | `resolvePendingAmbiguousToHash`: accounting call + orphan epilogue |
| `core/legacy/source/SoloFeature_sumThreads.cpp` | Release/clear orphan map on thread state |
| `flex/source/solo/CbBayesianResolver.cpp` | Skip UMI term when histogram empty (bridge uses empty UMI map for resolve) |

---

## Design docs and artifacts (update targets)

- **Runbook / design intent:** [RUNBOOK_SOLO_BRIDGE_AMBIGUOUS_CB_DESIGN_20260325.md](./RUNBOOK_SOLO_BRIDGE_AMBIGUOUS_CB_DESIGN_20260325.md) — extend with new benchmark subsection if missing on your branch.
- **Long history / tables:** [HANDOFF_SOLO_OPTIMIZATION_20260324.md](./HANDOFF_SOLO_OPTIMIZATION_20260324.md) — search for **“tuple redesign”**, **“Full EBs2_2 direct bridge — tuple redesign”**.
- **Output roots:** [tests/ARTIFACTS.md](../tests/ARTIFACTS.md) — sections **Solo non-Flex bridge v10 tuple redesign** and **Full UCSF EBs2_2 GEX bridge tuple redesign**.

---

## Validation already run (trust but re-check after merge)

### 2M UCSF `iPSC2_1` GEX (paired v9 control family)

- **Root (example):** `/tmp/v10_redesign_2m_parity3/` (or paths in `tests/ARTIFACTS.md` on worktree).
- **Env:** `STAR_SOLO_NONFLEX_HASH_BRIDGE=1`, `--soloInlineHashMode yes`, fresh `--outTmpDir`.
- **Parity:** `pending` / `added_to_hash` matched archived v9 log line; **Summary / matrices not byte-identical** (expected under Bayes-resolved accounting + key-level pin aggregation).

### Full corrected UCSF `EBs2_2` (32 threads, same CLI family as direct bridge v4)

- **Root:** `/storage/100K/ucsf_solo_bridge_redesign_20260325/ucsf_ebs2_2_gexonly_bridge_redesign_v1/`
- **Wall:** ~**20:08** (`/usr/bin/time -v`); **max RSS** ~**70 985 716 kB**.
- **Mapping log interval:** **531 s** (vs **678 s** on archived bridge **v4** log stamps on a prior day — treat as strong signal, not a controlled micro-benchmark).
- **Solo (original nested-map collapse on this artifact):** `collapseUMIall_fromBridgeHash` ~**278 s** (unchanged vs v4 — dominated wall until CSR + flat MultiGene rewrite; see below).
- **Summary:** GeneFull unique **0.85468**, cells **13715**, UMIs **254181669**, Total GeneFull **33750**.

**Serialize** full-sample benchmarks; do not overlap with other heavy STAR jobs on the same machine ([AGENTS.md](../AGENTS.md) benchmark hygiene).

---

## Known limitations / residual risk

1. **Legacy parity:** not a goal for the chosen accounting policy; 2M and full-sample **Summary** will drift from pin-replay / pin-only semantics.

2. **Key-level pin:** one `bridgeAmbigPinCandQuals_` per key vs per-read candidate quals → small statistical difference vs old replay for `noTooManyWLmatches` / pin-attributed counts.

3. **Read-stats bulk:** one **sample** flag per (unique vs multi) class × count — same class of approximation as discussed in the runbook.

4. **Collapse CPU (updated 2026-03-25):** post-map `collapseUMIall_fromBridgeHash` was optimized in worktree `/tmp/star-suite-v10-redesign-20260325`: **CSR + per-CB `(gene,umi)` sort** (no global tuple sort), then **flat `MgRow` MultiGeneUMI_CR** (sort by corrected UMI, vector merges — no nested `umi→gene→count` maps). **`umiArrayCorrect_CR` uses `readInfoRec=false`** on this path (`outSAMtype None` — no per-read `recordReadInfo` replay). Full EBs2_2 (32 threads): **`collapseUMIall_fromBridgeHash` ~48 s**, `processRecords` ~**278 s**, wall ~**16:22**, RSS ~**70.4 GiB** — artifact `/storage/100K/ucsf_solo_bridge_redesign_20260325/ucsf_ebs2_2_gexonly_flat_multigene_20260325/`. Prior CSR-only step-down: `…/ucsf_ebs2_2_gexonly_csr_groupedcollapse_20260325_v2/` (~**208 s** collapse). Compare tuple-redesign control `ucsf_ebs2_2_gexonly_bridge_redesign_v1` (~**278 s** collapse).

---

## Suggested next steps (prioritized)

1. **Merge worktree → integration branch** with `git merge --no-ff` if shared files overlap (see safe-merge policy in [AGENTS.md](../AGENTS.md) for large Solo files).

2. **Reproduce** full EBs2_2 wall + mapping interval on **two** serialized runs (same host, quiet machine) to confirm mapping win vs v4 is stable.

3. **HASH / collapse phase:** further gains likely in **parallel per-CB collapse** or shared patterns with `SoloFeature_collapse_bridge_global.cpp` / binary-spool HASH_MEX — see [HANDOFF_SOLO_OPTIMIZATION_NEXT_AGENT.md](./HANDOFF_SOLO_OPTIMIZATION_NEXT_AGENT.md).

4. **Optional parity mode:** only if product demands legacy stats — gate “pin-only ambiguous nRead” behind a flag (explicitly **not** the current default).

5. **2M regression:** after any collapse change, rerun the 2M bridge fixture and compare `Summary.csv` + logs to the last known good redesign artifact.

---

## Related handoffs

- [HANDOFF_SOLO_OPTIMIZATION_NEXT_AGENT.md](./HANDOFF_SOLO_OPTIMIZATION_NEXT_AGENT.md) — branch-wide Solo optimization index (direct bridge v9/v10, spool HASH_MEX).  
- [HANDOFF_SOLO_OPTIMIZATION_20260324.md](./HANDOFF_SOLO_OPTIMIZATION_20260324.md) — full benchmark narrative.
