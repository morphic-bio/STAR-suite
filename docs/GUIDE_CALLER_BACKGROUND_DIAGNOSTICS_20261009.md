# Guide-caller background diagnostics

The standalone `core/features/process_features/call_features` tool now writes
`guide_ambient_model.json` alongside its existing ambient-FDR outputs. The
receipt identifies the all-ambient composition null, states that no ambient
fraction was fitted, and counts calls that depend on a zero estimated ambient
rate. Numerical tail underflow for a positive rate is not classified as a
zero-background estimate.

The optional `--guide-fdr-pseudocount X` adds a finite nonnegative count per guide
to the ambient profile. Its default is zero, preserving existing calls. A value
of 0.5 was evaluated on public guide-only count fixtures; it removes zero-rate
estimates but does **not** repair the severe Zhao CRISPRa coverage failure.
It is an experimental diagnostic, not a recommended replacement model.

The normalization remains `(ambient_count + X) / (ambient_total + X * guides)`.
Cell-specific null means multiply that proportion by the entire observed cell
guide depth. This option does not estimate a contamination fraction or claim
calibration across different assays.

Caller-only cloud validation used the Zhao CRISPRa 2M-read fixture, CAT-ATAC
2M-read fixture, and a fixed 1,024-cell Zhao count fixture retaining original
guide counts. Default outputs matched the archived implementation; smoothed
outputs matched an independent SciPy calculation. Nonfinite, negative and
malformed pseudocount arguments are rejected. No GEX or ATAC processing is
required to exercise these checks.

Workflow acceptance and explicit fallback selection belong to the orchestration
layer. GMM output must remain identified as GMM, with no copied ambient q-values.
These changes do not modify the default caller, add a STAR processing parameter,
or authorize a release.
