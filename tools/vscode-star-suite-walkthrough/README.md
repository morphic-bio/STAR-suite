# STAR-suite Demo Walkthrough

This repo-local VS Code extension adds guided Codespaces walkthroughs for the stable STAR-suite demo surfaces.

Current split:
- stable now: setup + mini-reference, bulk, SLAM, shared single-cell fixture derivation
- placeholder for now: perturb, Flex

The extension intentionally keeps all execution logic in the existing shell scripts under `scripts/codespaces/`. The extension only opens module guides and launches those scripts in an integrated terminal.
