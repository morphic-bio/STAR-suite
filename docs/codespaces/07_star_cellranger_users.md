# Guidance for STAR and Cell Ranger Users

This note is for users who already know normal STAR and/or Cell Ranger and want the walkthroughs to focus on what is different in STAR-suite rather than re-explaining the basic run model.

## How to use the walkthroughs

Use the chapter model non-linearly:
- skip directly to [Bulk](./02_bulk.md) or [SLAM](./03_slam.md) if that is your target surface
- use [Setup + Mini Reference](./01_setup_reference.md) only if you want the demo reference build or do not already have a suitable index
- use [Single-cell Fixture Derivation](./04_single_cell_fixture.md) only for the bounded public perturb/Flex input surface

## What is different from normal STAR

Across these walkthroughs, the main differences from a normal STAR session are:
- the demo wrappers fetch or derive bounded public fixtures for you
- the single-cell demos use a shared public `chr22 + chrY` mini-reference to keep Codespaces runs small
- the walkthrough scripts write `RUN_COMMAND.sh` so users can inspect and then replace the defaults
- perturb and Flex use STAR-suite wrappers around STAR-compatible inputs rather than raw STAR alone

## What is different from Cell Ranger expectations

The main Cell Ranger deltas in this walkthrough surface are:
- no `mkref` path in Codespaces; the contract is the processed `FASTA + GTF` source surface and STAR builds the runnable index directly
- barcode geometry is treated as assay-specific input metadata, not a global constant
- CR-style configs are used where appropriate, but STAR-suite may pin or supersede some behavior with STAR-side parameters
- the current perturb and Flex chapters are scaffold/prototype surfaces, not final public parity demos yet

## Chapter-specific deltas

### Setup + Mini Reference
- For STAR users: this is just a small public source reference plus `genomeGenerate`
- For Cell Ranger users: the important parity target is the processed `genome.fa` and `genes.gtf`, not a `mkref` artifact

### Bulk
- For STAR users: this is closest to a normal STAR run; the main extra value is fixture download and `RUN_COMMAND.sh` generation
- For Cell Ranger users: there is no CR-specific abstraction here; this is simply the bulk STAR command shape

### SLAM
- For STAR users: the key deltas are `--slamQuantMode`, `--slamGrandSlamOut`, `--autoTrim variance`, and `--slamQcReport`
- For Cell Ranger users: this is STAR-SLAM specific and not a CR-like flow

### Single-cell Fixture Derivation
- For STAR users: this is a preprocessing utility to keep the public single-cell demo bounded in Codespaces
- For Cell Ranger users: this is where the walkthrough deliberately constrains the public input surface rather than trying to mirror a full CR-scale ingest

### Perturb
- For STAR users: the chapter is mainly about CR-config compatible wrapper inputs, not a plain `soloType` tutorial
- For Cell Ranger users: the important part is the CR-config ingestion surface and the STAR-side parameters that remain pinned by the wrapper

### Flex
- For STAR users: the chapter is about the STAR-suite Flex wrapper surface, not just raw STAR indexing/alignment
- For Cell Ranger users: the wrapper is config-driven, but the current public demo is still a prototype scaffold and not the final parity surface
