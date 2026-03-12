# Perturb Demo

Current status: scaffold only.

The repository already has a runnable prototype wrapper for a public `chr22 + chrY` perturb demo surface, but it is not yet strong enough to present as a faithful public assay walkthrough. The structure is still exposed in Codespaces so the module slot is stable and easy to swap later.

Current prototype preview:
```bash
bash scripts/codespaces/run_perturb_public_demo.sh --dry-run
```

Current prototype run surface:
```bash
bash scripts/codespaces/run_perturb_public_demo.sh
```

What is stable here:
- shared mini-reference creation
- shared GEX fixture derivation
- CR-config based perturb wrapper integration

What is still placeholder quality:
- the generated guide companion does not yet produce a strong public demo signal

Planned replacement:
- swap in a bounded real public perturb assay fixture while keeping the same walkthrough slot and command shape
