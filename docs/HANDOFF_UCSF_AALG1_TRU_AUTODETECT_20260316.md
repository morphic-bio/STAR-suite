# Handoff: UCSF `AALG1` TRU Autodetect Debug (2026-03-16)

## Resolution / Status Update (2026-07-03)

**Verified status: PARTIAL / SUPERSEDED, not confirmed as a poly-G-driven
AALG1 autodetect bug.**

Current evidence from code, commits, and tracked benchmark artifacts supports
the following:

- `AALG1` is a UCSF dataset label, not an A375 benchmark label. The corrected
  UCSF organization treats `AALG1` as the TRU GEX arm and `AALG2` as the NXT
  guide arm. See `docs/HANDOFF_UCSF_FASTQ_MISLABEL_INVESTIGATION_20260316.md`.
- The old failing `iPSC2_1_AALG1` auto run is superseded by the corrected UCSF
  organization and by current benchmark/production commands that explicitly use
  a TRU Solo whitelist for GEX and a two-column NXT whitelist for guide
  assignment.
- `--clip3pPolyG` is implemented and wired into the current STAR clipping path
  by commit `fce35365f73eedeae931897bd3e466d10fee4179`:
  `Parameters.cpp` registers the parameter, `ParametersClip_initialize.cpp`
  resolves `yes|no|auto`, and `ClipCR4::polyTail3p()` trims poly-A and,
  when enabled, poly-G 3' tails. This plausibly fixes the UCSF LINC00486
  gene-count artifact documented in
  `docs/HANDOFF_UCSF_FULL_GENE_PEARSON_ANOMALY_20260220.md`.
- I did **not** find code evidence that poly-G trimming directly changes the
  16 bp barcode namespace decision. The NXT/TRU autodetect path samples raw
  barcode strings and compares raw versus translated whitelist hits in
  `core/features/process_features/src/assignBarcodes.c`, then composes the
  namespace decision in `core/legacy/source/PfMultiProcess.cpp`.
- The STAR-side mixed-chemistry fixes that are relevant to current UCSF
  processing are separate namespace/autodetect fixes:
  - `8c402fc` (`2026-02-25`): NXT/TRU chemistry autodetection and per-library
    override infrastructure.
  - `f692f45` (`2026-03-03`): source-namespace-aware normalization and shallow
    bootstrap OrdMag behavior.
  - `6452600` (`2026-03-17`): NXT-whitelist autodetect inversion and filtered
    barcode normalization fixes.
- Tracked benchmark evidence that the corrected UCSF surface is healthy is in
  commit `207581e84355adce437fe536d4a2e7a18f2ecf67`
  (`comparisons/paper_benchmarks_20260318/ucsf_ebs2_2/`). The recorded command
  uses AALG1 GEX FASTQs, `--soloCBwhitelist ...3M-february-2018_TRU.txt`,
  `--crChemistry auto`, `--crWhitelist ...translation/3M-february-2018_NXT.txt`,
  and `--clip3pPolyG yes`. The tracked parity report records STAR/CR cells
  `13,721 / 13,760`, barcode Jaccard `0.976`, gene Pearson `0.995`, and
  CRISPR set-equivalent calls `98.9%`.

What remains unconfirmed:

- I did not find a tracked rerun of the original `iPSC2_1_AALG1` auto command
  proving that the exact old failure now resolves to TRU without the corrected
  UCSF organization and current command surface.
- I did not find a dedicated automated regression test for "AALG1 auto chooses
  TRU." Current protection is the corrected UCSF benchmark/script surface plus
  mixed NXT/TRU autodetect and namespace smokes.

Release-note check:

- `docs/RELEASE_NOTES_*` do not currently call out the AALG1/AALG2 UCSF
  organization correction, `--clip3pPolyG`, or the NXT/TRU autodetect fixes by
  name. `v1.0.0` is the natural release-note bucket because these fixes landed
  before the first production release.

## Problem

UCSF `AALG1` samples showed a severe STAR vs Cell Ranger cell-calling mismatch
while `AALG2` parity was previously strong.

The concrete failure case was `iPSC2_1_AALG1`:

- STAR auto-detect run:
  `/mnt/pikachu/ucsf-perturb-yremove_all_velocyto_globus_20260313_093723/samples/iPSC2_1_AALG1/run/`
- Cell Ranger run:
  `/mnt/pikachu/cellranger_iPSC2_1_AALG1_full_20260316_001/`

Observed mismatch:

- STAR auto: `2613` filtered cells
- Cell Ranger: `8011` filtered cells
- raw barcode overlap with suffix-normalization: `4`
- Jaccard: `0.000377`

This is far outside normal UCSF CR-compat parity.

## Key Finding

Forcing `TRU` on the same STAR command fixes the cell-calling parity.

Forced-TRU rerun:

- `/mnt/pikachu/star_iPSC2_1_AALG1_force_tru_20260316_0327/`
- command:
  [RUN_COMMAND.sh](/mnt/pikachu/star_iPSC2_1_AALG1_force_tru_20260316_0327/RUN_COMMAND.sh)

The rerun explicitly sets:

- `--soloCBwhitelist /home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt`
- `--crChemistry TRU`
- `--crOutputChemistry TRU`
- `--crWhitelist /home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt`

`GeneFull` result from
[Log.out](/mnt/pikachu/star_iPSC2_1_AALG1_force_tru_20260316_0327/run/Log.out):

- `finished emptyDrops_CR filtering (libscrna): total pass=8010`

That is essentially identical to Cell Ranger:

- STAR forced `TRU`: `8010`
- Cell Ranger: `8011`

After normalizing Cell Ranger barcodes by stripping `-1`, overlap is:

- shared cells: `7992`
- Jaccard: `0.995392`
- STAR-only: `18`
- CR-only: `19`

So the main `AALG1` cell-calling mismatch is not an EmptyDrops issue. It is
caused by the chemistry/autodetect path selecting the wrong barcode regime.

## Important Detail: Suffix vs Namespace

One misleading intermediate result was raw zero overlap between forced-TRU STAR
and Cell Ranger. That was only because STAR writes bare barcodes while Cell
Ranger writes the same strings with `-1`.

Example:

- STAR forced TRU:
  `AAACCCAAGCAGAAAG`
- Cell Ranger:
  `AAACCCAAGCAGAAAG-1`

This was not a namespace mismatch.

## Evidence Collected

### 1. STAR auto-detect claims `NXT` for both `AALG1` and `AALG2`

`iPSC2_1_AALG1`:

- [Log.out](/mnt/pikachu/ucsf-perturb-yremove_all_velocyto_globus_20260313_093723/samples/iPSC2_1_AALG1/run/Log.out)
- reports:
  `pf-multi chemistry: requested=auto inferred=NXT ... effective=NXT ... soloOutputNamespace=TRU output_effective=TRU`

`iPSC2_1_AALG2`:

- [Log.out](/mnt/pikachu/ucsf-perturb-yremove_all_velocyto_globus_20260313_093723/samples/iPSC2_1_AALG2/run/Log.out)
- also reports `inferred=NXT ... output_effective=TRU`

### 2. Raw barcode sequence evidence disagrees for `AALG1`

Using the first 200k GEX `R1` barcodes from one lane, exact matches against the
two columns of
`/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt`
showed:

`iPSC2_1_AALG1`

- matches whitelist col1: `29,687`
- matches whitelist col2: `190,355`

`iPSC2_1_AALG2`

- matches whitelist col1: `189,878`
- matches whitelist col2: `31,133`

Interpretation:

- `AALG1` barcode reads look strongly `TRU`-like
- `AALG2` barcode reads look strongly `NXT`-like

This does **not** necessarily mean different wet-lab chemistries. It means the
barcode namespace in the reads looks different.

### 3. Cell Ranger chemistry detection also differs

`AALG1` Cell Ranger detection:

- `/mnt/pikachu/cellranger_iPSC2_1_AALG1_full_20260316_001/iPSC2_1_AALG1_full_20260316_001/SC_MULTI_CS/SC_MULTI_CORE/MULTI_CHEMISTRY_DETECTOR/DETECT_COUNT_CHEMISTRY/fork0/join/_outs`
- detected:
  - `SC3Pv3-polyA`
  - whitelist `3M-february-2018_TRU`
  - `translation: false`

Earlier matched `AALG2` Cell Ranger run:

- `/mnt/pikachu/ucsf-perturb-yremove-smoke_1m_integrated_20260311_221028/cellranger_iPSC2_1_AALG2_1m_match/iPSC2_1_AALG2_1m_match/SC_MULTI_CS/SC_MULTI_CORE/MULTI_CHEMISTRY_DETECTOR/DETECT_COUNT_CHEMISTRY/fork0/join/_outs`
- detected:
  - `SC3Pv3-CS1`
  - whitelist `3M-february-2018_NXT`
  - `translation: true`

So this is not just a STAR-only guess problem. Cell Ranger itself treats
`AALG1` and `AALG2` differently.

## What This Means

For `AALG1`, STAR auto-detect is landing on the wrong chemistry path for CR
parity. The forced-`TRU` run restores near-perfect cell-set parity with Cell
Ranger.

For `AALG2`, the prior good-parity benchmark appears to be correct under `NXT`
matching chemistry with `TRU` output namespace.

So the likely bug is:

- not a generic `NXT`/`TRU` merge bug
- not EmptyDrops
- not GEX quantification itself
- but specifically the auto-chemistry detection or chemistry-resolution logic
  for `AALG1`-like barcode inputs

## Minimal Repro

Compare these three runs:

1. Cell Ranger `AALG1`
   - `/mnt/pikachu/cellranger_iPSC2_1_AALG1_full_20260316_001/`
2. STAR auto `AALG1`
   - `/mnt/pikachu/ucsf-perturb-yremove_all_velocyto_globus_20260313_093723/samples/iPSC2_1_AALG1/run/`
3. STAR forced-TRU `AALG1`
   - `/mnt/pikachu/star_iPSC2_1_AALG1_force_tru_20260316_0327/`

Filtered barcode files:

- Cell Ranger:
  `/mnt/pikachu/cellranger_iPSC2_1_AALG1_full_20260316_001/iPSC2_1_AALG1_full_20260316_001/outs/per_sample_outs/iPSC2_1_AALG1_full_20260316_001/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz`
- STAR auto:
  `/mnt/pikachu/ucsf-perturb-yremove_all_velocyto_globus_20260313_093723/samples/iPSC2_1_AALG1/run/outs/filtered_feature_bc_matrix/barcodes.tsv.gz`
- STAR forced TRU:
  `/mnt/pikachu/star_iPSC2_1_AALG1_force_tru_20260316_0327/run/Solo.out/GeneFull/filtered/barcodes.tsv`

## Next Steps For The Next Agent

1. Finish the forced-TRU run and confirm final artifacts.
   - At handoff time it was still running later stages after `GeneFull`.
   - The parity result is already established from `GeneFull`.

2. Debug the auto-detect path, not the downstream comparison code.
   Focus on:
   - `detectChemistryFromWhitelistPath`
   - any barcode-content-driven chemistry inference
   - the interaction between `pf-multi chemistry: requested=auto inferred=...`
     and CR-compat chemistry resolution

3. Reproduce the decision point that makes `AALG1` resolve to `NXT`.
   - Compare the raw `R1` barcode signal for `AALG1` vs `AALG2`
   - Compare this to Cell Ranger chemistry detection outputs

4. Determine whether the fix should be:
   - better autodetect logic
   - a stronger confidence gate for auto mode
   - or explicit override behavior when Cell Ranger-style evidence is clearly
     `TRU`-like

5. After fixing autodetect, rerun `iPSC2_1_AALG1` in STAR auto mode and verify:
   - filtered cells near `8011`
   - Jaccard near `0.995`
   - no need for manual `TRU` override

## Non-Issues / False Leads Already Eliminated

- EmptyDrops parameters: not the main cause
- barcode suffix mismatch (`-1`): only cosmetic
- generic `NXT`/`TRU` output namespace confusion: not sufficient to explain the
  `AALG1` failure
- old `AALG2` parity benchmark being a pure `TRU` chemistry benchmark: false;
  it was an `NXT`-matching run with `TRU` output namespace

## Build Note

Per repo policy, a clean rebuild was done before this debugging session:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
```
