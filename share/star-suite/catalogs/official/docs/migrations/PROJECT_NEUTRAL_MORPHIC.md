# Project-neutral MorPhiC recipe migration

Source snapshot: `morphic-bio/morphic-recipes` commit
`d68fc54be5c0f8f27a25373541407706067847ae`.

The reusable multiome engine, scRNA downstream, MuData builders, QC reports,
and public multiomics smoke wrappers moved to the official STAR Suite catalog.
Their scripts remain open and self-contained within the catalog root.

The import removed MorPhiC provenance assumptions and workstation defaults.
Dataset, reference, transfer, SSH, and output locations are caller-supplied.
Remote execution remains opt-in and requires explicit host and staging-root
parameters. The official multiome profile defaults to the minimal local
matrices-and-peaks surface.

Project-specific JAX, MSK, and UCSF wrappers, transfer destinations, sample
manifests, and handoff automation remain canonical in `morphic-recipes`.
