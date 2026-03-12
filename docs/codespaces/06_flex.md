# Flex Demo

Current status: scaffold only.

The repository already has a runnable prototype wrapper for a public `chr22 + chrY` Flex demo surface, but it is not yet strong enough to present as a faithful public assay walkthrough. The structure is still exposed in Codespaces so the module slot is stable and easy to swap later.

Current prototype preview:
```bash
bash scripts/codespaces/run_flex_public_demo.sh --dry-run
```

Current prototype run surface:
```bash
bash scripts/codespaces/run_flex_public_demo.sh
```

What is stable here:
- shared mini-reference creation
- shared GEX fixture derivation
- FLEX wrapper/config integration

What is still placeholder quality:
- the synthetic probe/sample-tag overlay is too weak for a proper public demo

Planned replacement:
- swap in a bounded real public Flex assay fixture while keeping the same walkthrough slot and command shape

## If You Already Know STAR or Cell Ranger

- STAR users: this chapter is about the STAR-suite Flex wrapper surface rather than raw STAR alone.
- Cell Ranger users: treat this as a prototype config-driven scaffold for now, not a finished public parity walkthrough.
