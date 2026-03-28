# Handoff: Full Velocyto Baseline Triage (2026-03-28)

## Summary

The original handoff scope was wrong.

There are **two different historical failures** in the preserved March 27 runs:

1. a **full PF + Velocyto** run on corrected `EBs2_2` that stopped during
   mapping at 8 threads
2. a **GEX-only Velocyto** run that finished mapping and `GeneFull`, then
   stopped at the start of `Gene` post-map counting

After rerunning the full PF + Velocyto surface on March 28:

- the old PF + Velocyto case is **no longer strong evidence of a hang**
- it is more consistent with a **prematurely terminated 8-thread run**
- the stronger remaining bug candidate is the **GEX-only Velocyto `Gene`
  post-map stall**

## What Was Reproduced On 2026-03-28

### Full PF + Velocyto, 8 threads

Command surface:

- [run_star_velocyto_canonical.sh](/mnt/pikachu/STAR-suite/scripts/run_star_velocyto_canonical.sh)
- `--profile 2m`
- `--threads 8`
- `--soloFeatures Gene GeneFull Velocyto`
- `--pfMultiConfig /storage/ucsf-velocyto-validation/pf_multi_config_ebs2_2_full.csv`
- `--outSAMtype BAM Unsorted`

New repro run:

- [/tmp/full_ebs2_2_pf_velocyto_repro_20260328_031211](/tmp/full_ebs2_2_pf_velocyto_repro_20260328_031211)

Observed:

- the run was still healthy in mapping when manually stopped
- it passed the historical old stop point of `153,300,815` reads
- therefore the old PF + Velocyto artifact did **not** reproduce as a hang

### Full PF + Velocyto, 32 threads

New repro run:

- [/tmp/full_ebs2_2_pf_velocyto_repro32_20260328_033441](/tmp/full_ebs2_2_pf_velocyto_repro32_20260328_033441)
- [Log.out](/tmp/full_ebs2_2_pf_velocyto_repro32_20260328_033441/stream_t32/Log.out)
- [Log.final.out](/tmp/full_ebs2_2_pf_velocyto_repro32_20260328_033441/stream_t32/Log.final.out)

Observed:

- mapping completed
- `GeneFull`, `Gene`, and `Velocyto` completed
- pf-multi merge/finalize completed
- `ALL DONE!`

Conclusion:

- the full PF + Velocyto baseline currently **completes cleanly**

## Reinterpreted Historical Runs

### A. Full PF + Velocyto run: likely stopped early, not proven hung

Historical artifact:

- [Log.out](/storage/ucsf-velocyto-validation/full_ebs2_2_velocyto_stream_vs_hash_20260327_055805/stream_tN/Log.out)
- [Log.progress.out](/storage/ucsf-velocyto-validation/full_ebs2_2_velocyto_stream_vs_hash_20260327_055805/stream_tN/Log.progress.out)

What the logs show:

- command used `--runThreadN 8`
- command included `--pfMultiConfig` and `--soloFeatures Gene GeneFull Velocyto`
- `Log.progress.out` was still advancing through mapping
- last progress line was `153,300,815` reads at `2026-03-27 06:11:01 UTC`
- it never reached Solo post-map

Current interpretation:

- this looks more like a run that was **stopped before completion**
- it is **not strong evidence of a shutdown hang**

### B. Full GEX-only Velocyto run: still the main unresolved suspect

Historical artifact:

- [Log.out](/storage/ucsf-velocyto-validation/full_gex_velocyto_stream_vs_hash_20260327_061129/stream_tN/Log.out)
- [Log.progress.out](/storage/ucsf-velocyto-validation/full_gex_velocyto_stream_vs_hash_20260327_061129/stream_tN/Log.progress.out)

What the logs show:

- no `--pfMultiConfig`
- no guide FASTQs
- `--soloFeatures Gene GeneFull Velocyto`
- mapping finished
- `GeneFull` completed
- then:
  - `Mar 27 06:38:37 ... Starting Solo post-map for Gene`
  - `Mar 27 06:38:38 ... Allocated and initialized readInfo array, nReadsInput = 444896732`
  - no further progress
- no `Log.final.out`
- directory contains only:
  - `Log.out`
  - `Log.progress.out`
  - `SJ.out.tab`
  - `Chimeric.out.junction`

Current interpretation:

- this still looks like a **real post-map stall candidate**
- however, it is **not proven** from artifacts alone
- there is no preserved wrapper console log, PID capture, or exit status
- so manual termination because it looked too slow cannot be ruled out

## Best Current Answer To "Did It Stall Or Was It Stopped?"

### Full PF + Velocyto case

Best answer:

- most likely **stopped early because it was actually an 8-thread run**
- not a confirmed hang

### GEX-only Velocyto case

Best answer:

- **looks more like a stall than a simple "too slow" stop**
- because it reached a very specific late phase boundary and then stopped
- but this is still **not proven**

Why the GEX-only case looks more suspicious:

- it finished mapping
- it finished `GeneFull`
- it stopped at the exact start of `Gene` post-map
- a later rerun on the same general surface did complete:
  - [stream_t32](/storage/ucsf-velocyto-validation/full_gex_velocyto_stream_vs_hash_t32_20260327_064213/stream_t32/Log.out)
  - [hash_t32](/storage/ucsf-velocyto-validation/full_gex_velocyto_stream_vs_hash_t32_20260327_064213/hash_t32/Log.out)

That pattern fits an intermittent late-phase issue better than the PF run does.

## Debug Instrumentation Added

These changes are already in the working tree and are opt-in only:

- [assignBarcodes.c](/mnt/pikachu/STAR-suite/core/features/process_features/src/assignBarcodes.c)
  - `PF_TRACE_BOOTSTRAP=1`
- [SoloFeature_processRecords.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_processRecords.cpp)
  - `STAR_SOLO_PHASE_DEBUG=1`
- [SoloFeature_countCBgeneUMI.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countCBgeneUMI.cpp)
  - `STAR_SOLO_PHASE_DEBUG=1`

These were validated on successful runs and should be reused for any new repro.

## Most Relevant Next Step

If this bug is pursued further, the best target is:

- the **full GEX-only Velocyto stream** surface
- same binary
- fresh outdir
- debug envs enabled

Command family:

- [run_star_velocyto_gexonly_canonical.sh](/mnt/pikachu/STAR-suite/scripts/run_star_velocyto_gexonly_canonical.sh)

Goal:

- determine whether the old `Gene` post-map stop is reproducible
- if it reproduces, capture PID and stacks before killing

## Bottom Line

The earlier claim of a **combined GEX + feature-assignment hang** was too broad.

Current state:

- **full PF + Velocyto baseline is not currently reproducing a hang**
- the old PF + Velocyto artifact is likely a prematurely terminated 8-thread run
- the remaining unresolved bug candidate is the **GEX-only Velocyto `Gene`
  post-map stall**
