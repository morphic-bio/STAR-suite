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

Suggested order:
1. [Setup + Mini Reference](./01_setup_reference.md)
2. [Bulk](./02_bulk.md)
3. [SLAM](./03_slam.md)
4. [Single-cell Fixture Derivation](./04_single_cell_fixture.md)
5. [Perturb Placeholder](./05_perturb.md)
6. [Flex Placeholder](./06_flex.md)
