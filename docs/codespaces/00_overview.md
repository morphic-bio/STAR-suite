# STAR-suite Codespaces Demos

This Codespaces surface is a guided entry point for the stable STAR-suite demo flows.

Current status:
- Stable now: setup + mini-reference, bulk, SLAM, shared single-cell fixture derivation
- Scaffold only for now: perturb, Flex

Design rules:
- In Codespaces, `.devcontainer/postCreate.sh` packages and installs the repo-local walkthrough extension on a best-effort basis. The fallback entry point is Command Palette -> `STAR-suite: Open Codespaces Walkthrough`.
- The walkthroughs use helper scripts under [`scripts/codespaces/`](../../scripts/codespaces/).
- The mini-reference is a public human `chr22 + chrY` `FASTA + GTF` source surface built directly with STAR.
- Advanced users can swap in their own indices, references, or FASTQs at any step.
- Perturb and Flex are intentionally present as placeholders so the module structure is stable even while their public demo fixtures are still being tightened.

Fast delta guide:
- [Guidance for STAR and Cell Ranger Users](./07_star_cellranger_users.md)

Chapter model:
- Optional shared foundation:
  - [Setup + Mini Reference](./01_setup_reference.md)
- Independent module chapters:
  - [Bulk](./02_bulk.md)
  - [SLAM](./03_slam.md)
  - [Perturb Placeholder](./05_perturb.md)
  - [Flex Placeholder](./06_flex.md)
- Shared utility chapter for perturb and Flex only:
  - [Single-cell Fixture Derivation](./04_single_cell_fixture.md)

Use this non-linearly:
- skip the foundation chapter if you already have a suitable reference or index
- go directly to Bulk or SLAM if that is all you need
- use the single-cell fixture chapter only if you want the bounded public input used by perturb and Flex

Starting points by experience:
- Experienced users:
  - skip directly to [Bulk](./02_bulk.md) or [SLAM](./03_slam.md) if you already have a suitable index
  - skip to [Single-cell Fixture Derivation](./04_single_cell_fixture.md) only if you specifically want the bounded public input used by perturb or Flex
  - use [Perturb Placeholder](./05_perturb.md) or [Flex Placeholder](./06_flex.md) as command-shape scaffolds if you already understand the assay-specific inputs
- STAR / Cell Ranger users:
  - start with [Guidance for STAR and Cell Ranger Users](./07_star_cellranger_users.md) if you want the walkthrough to focus on deltas instead of basics
  - then jump directly to the relevant module chapter
- Novice users:
  - start with [Setup + Mini Reference](./01_setup_reference.md) if you want to see how the demo reference is built
  - then try [Bulk](./02_bulk.md) first, because it is the smallest stable end-to-end chapter
  - then move to [SLAM](./03_slam.md)
  - leave perturb and Flex for later unless you specifically want the prototype surfaces
