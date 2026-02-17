# UCSF Call-Only GMM Parity Bundle

This directory captures a reproducible comparison between Cell Ranger CRISPR
calls and STAR helper call-only GMM calls on the UCSF iPSC2_1_AALG2 sample.

Contents:
- `run_callonly_gmm_parity.sh`: end-to-end reproducible script
- `RESULTS.md`: summary from the 2026-02-17 runs (`--min-umi 10` and `--min-umi 3`)

The workflow mirrors the earlier A375 call-only parity process, but uses the
UCSF sample and 548-guide feature set.
