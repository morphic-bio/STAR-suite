# STAR Perturb Handoff (2026-02-17)

## Branch scope

This `noBCfix` branch focused on CRISPR feature-call parity/stability and
integration cleanup for STAR perturb-mode workflows, not full end-to-end GEX
parity.

## Main changes completed

1. CRISPR call-only parity harness and results bundle
- Added reusable comparison bundle under:
  `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/`
- Captures run scripts, thresholds, and reproducible STAR vs Cell Ranger
  call-only parity outputs for UCSF iPSC2/AALG2.

2. Internal naming cleanup (`CrMulti` -> `PfMulti`)
- Renamed CR-compat internal wiring classes/files from `CrMulti*` to
  `PfMulti*` to avoid implying Cell Ranger code reuse and to align with
  `process_features` naming.

3. `process_features` integration harmonization
- Reduced interface drift between external `process_features` behavior and
  in-STAR integration paths.
- Consolidated parameter flow so in-proc feature assignment and standalone
  behavior are closer under matched settings.

4. `noBC` edge-case/performance work
- Added handling for no-BC edge behavior in perturb workflows and improved
  matching path behavior under high guide-cardinality settings.
- Resulted in large speedup vs earlier linear bootstrap behavior.

5. 2M downsample parity experiments
- Built and used UCSF downsampled set in `/storage/ucsf-2M`.
- Re-ran CR baseline with `--min-crispr-umi=3` and repeated STAR/CR parity
  comparisons for feature UMIs and call-level agreement.

## Current parity status (feature-focused)

Using the current CR rerun and STAR NXT run:
- CRISPR UMI-per-barcode Pearson/Spearman are effectively at parity
  (~0.999 / ~0.999 on mapped overlap).
- Call-level agreement (`feature_call`) is ~99% on overlapping mapped cells.

Important caveat:
- The STAR run used for one GEX comparison ingested mixed GEX+guide FASTQs
  (`3,000,000` reads total), which inflates STAR GEX UMIs.
- This is not feature-to-GEX matrix double counting; it is read-set mixing in
  the GEX pass.

## What still needs to be fixed (core compatibility)

1. Barcode whitelist compatibility
- Add first-class support for 2-column translation whitelists (NXT/TRU map).
- Persist and use both barcode namespaces consistently (input correction vs
  emitted/output namespace).
- Add chemistry-aware whitelist defaults and explicit override controls.

2. Auto-detection and pattern robustness
- Harden CB/UMI/pattern auto-detection for perturb libraries.
- Keep deterministic fallback behavior for small datasets and noisy starts.

3. Cell-calling mode parity (small sets)
- Add explicit compatibility path for ORDMAG/simple EmptyDrops behavior on
  small datasets.
- Ensure feature EM/GMM can be restricted to the intended filtered cell set.

4. Prevent feature-library leakage into GEX quant
- In perturb multi-library mode, ensure GEX quant only consumes GEX libraries
  (unless explicitly requested otherwise).
- Add regression test asserting stable GEX UMIs when guide FASTQs are present.

5. Regression matrix to gate merges
- A375 parity (existing baseline).
- UCSF 2M parity (feature UMIs, calls, runtime).
- One tiny small-set test for ORDMAG behavior and whitelist mapping edge cases.

## Recommended next branch

- Create from updated `master`:
  `core-compatibility-fixes`
- First tasks:
  1. Whitelist namespace mapping (2-column support).
  2. GEX-library-only quant path in perturb multi mode.
  3. ORDMAG/small-set compatibility switch and tests.
