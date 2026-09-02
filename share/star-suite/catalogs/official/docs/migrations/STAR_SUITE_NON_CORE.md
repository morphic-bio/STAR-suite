# STAR Suite non-core migration

Source snapshot: `morphic-bio/STAR-suite` commit
`c12ad552c846eeb16412058a06573cb6c4dc6164`.

The SLAM operational workflows, helper, and container recipe were copied from
STAR Suite into this catalog. Host-specific input, reference, binary, output,
and transfer defaults were removed. The algorithms and output formats remain in
STAR Suite; these files only orchestrate them.

STAR Suite 1.7 retains the original built-in workflow IDs as compatibility
surfaces. The official IDs are namespace-qualified and carry the former ID in
`logical_id`, so clients can keep variants separate, prompt users, or select a
version through the configured reconciliation policy.

MorPhiC-specific UCSF, JAX, and MSK production wrappers remain canonical in
`morphic-recipes`; they are not copied here. Their reusable engines and public
reporting/building recipes are imported separately.
