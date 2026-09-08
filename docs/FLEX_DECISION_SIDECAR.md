# Flex decision sidecar

`--soloFlexDecisionSidecar FILE` enables a diagnostic, fixed-width binary
sidecar for the fused Flex inline-hash classifier and the ordinary
BAM-producing Flex path. It is disabled by default (`-`). It does not change
classification or counting policy.

The file uses schema `FLXDEC2`, version 2: a 512-byte little-endian header and
one 48-byte record at `512 + 48 * global_ordinal`. The final file is published
only after every input pair has a record and the header is marked complete.
Interrupted runs retain only `FILE.tmp`.

Each record contains:

- the run's zero-based global read ordinal (implicit in its slot);
- input lane, lane-local record ordinal when the fused reader provides it, and
  a 64-bit FNV-1a hash of the normalized read name;
- cache action, cache class, gene, negative code, probe offset and probe
  region;
- whether the conservative exactly-one-N retry ran and resolved, including
  the underlying matching tier before its runtime H1 normalization;
- sample-tag checked/matched/rejected state and token;
- residual-alignment handoff and its resolved/rejected probe/genomic result.
- H1X2 probe seed outcome (absent, ambiguous, split-probe, or score failure);
- legacy residual half-anchor/alignment fields retained for sidecars written by
  the earlier experimental alignment-gated implementation.

Fused lane work stealing can assign global ordinals in a different interleaving
when the thread count changes. `lane + lane_ordinal` and the normalized-name
hash are therefore the cross-run join keys. Global ordinals remain the exact
join key to STAR's BAM/read-ID state from the same run.

The diagnostic uses random-access output and a read-modify-write operation for
residual-alignment updates. It is intended for focused audits, not timed
benchmarks. With the option omitted the writer is never opened and the hot path
only checks a null pointer. Storage is exactly 48 bytes per read plus the
512-byte header: about 19.2 MB for 400,000 reads, but about 350 GB for 7.3
billion reads. Full-dataset diagnostics therefore require explicit capacity
planning; downsampled read-level audits are the intended use.

Build and dump a completed sidecar with:

```sh
make -C flex/tools/molecule_first_resolver flex_decision_sidecar_dump
flex/tools/molecule_first_resolver/flex_decision_sidecar_dump FILE
```

The dump is a global-ordinal-ordered TSV suitable for joining to BAM ledgers or
cache-replay output.

For a BAM-independent read/molecule audit, set `STAR_FLEX_DECISION_LEDGER` to an
output TSV path. Its TRIAGE and FINAL events contain the bounded read name,
lane/local ordinal, sample tag, corrected cell barcode, UMI, chosen gene,
cache tier, final molecule key, rejection reason, and (for seed-and-extend)
the full-probe Hamming distance. It is diagnostic and disabled by default.
