# CR Config Compatibility Matrix

This note describes the current wrapper-level compatibility surface for using
Cell Ranger `config.csv` inputs with STAR-suite. It is not a claim of full
Cell Ranger config parity. The goal is narrower:

- use the same library inputs that Cell Ranger used
- reuse the same feature references where practical
- document exactly which fields STAR-suite consumes, extends, ignores, or
  supersedes with pinned STAR wrapper parameters

The active wrappers are:

- perturb:
  - [run_ucsf_full_compat_forward_rescue_guides.sh](/mnt/pikachu/STAR-suite/scripts/run_ucsf_full_compat_forward_rescue_guides.sh)
  - [run_ucsf_perturb_yremove_batch.sh](/mnt/pikachu/STAR-suite/scripts/run_ucsf_perturb_yremove_batch.sh)
- Flex:
  - [run_flex_cr_config.sh](/mnt/pikachu/STAR-suite/scripts/run_flex_cr_config.sh)

## General Rule

The Cell Ranger config is the source of truth for input libraries.

It is not the source of truth for all STAR execution parameters.

STAR-suite still pins core runtime choices such as:

- `--genomeDir`
- CB whitelist paths
- CR-compat Solo settings
- Flex settings
- Y-removal behavior

Those are controlled by the STAR wrapper, not imported blindly from the Cell
Ranger config.

## Perturb Support

### Consumed from standard CR sections

`[libraries]`
- `fastqs`
  - consumed
  - exact FASTQ source for GEX and feature libraries
- `feature_types` / `feature_type` / `library_type`
  - consumed
  - used to classify GEX vs non-GEX and identify guide-capture libraries
- `sample` / `fastq_id` / `library_id`
  - partially consumed
  - used only for sample ID inference if a single logical sample can be
    inferred

`[feature]`
- `reference` / `ref`
  - consumed
  - used as the global feature reference when per-library `star_feature_ref`
    is absent

`[gene-expression]`
- `reference`
  - recorded in the run manifest only
  - not used to set STAR `--genomeDir`
- `chemistry`
  - recorded in the run manifest only
  - not used to choose STAR chemistry settings automatically

### STAR-suite perturb extensions

`[libraries]`
- `star_chemistry`
  - consumed
  - passed through into generated `pf_multi_config.from_cr.csv`
- `star_feature_ref`
  - consumed
  - passed through into generated `pf_multi_config.from_cr.csv`
- `star_library_id`
  - consumed
  - passed through into generated `pf_multi_config.from_cr.csv`

`[star]`
- `sample-id` / `sample_id`
  - STAR extension
  - sets the logical STAR sample name when the CR config contains multiple
    library IDs that should collapse to one STAR sample

CLI override
- `--cr-sample-id`
  - takes precedence over `[star] sample-id`

### Ignored or superseded for perturb

Ignored if present in the config:
- `[gene-expression] create-bam`
- `[gene-expression] check-library-compatibility`
- other Cell Ranger-specific run-control keys not explicitly listed above

Superseded by pinned STAR wrapper parameters:
- `--genomeDir`
  - local pinned STAR reference policy, not Cell Ranger bundle path
- `--clipAdapterType CellRanger4`
- `--soloCellFilter EmptyDrops_CR`
- `--soloUMIfiltering MultiGeneUMI_CR`
- `--soloUMIdedup 1MM_CR`
- `--soloCrGexFeature genefull`
- `--soloCrMultimapRescue yes`
- `--crMinUmi`
- `--crWhitelist`
- all Y-removal flags in the batch Y-removal runner

## Flex Support

### Consumed from standard CR sections

`[libraries]`
- `fastqs`
  - consumed
  - exact GEX FASTQ source
- `feature_types` / `feature_type` / `library_type`
  - consumed
  - used only to identify GEX libraries

`[samples]`
- `sample_id`
  - consumed
  - becomes STAR Flex sample label
- `probe_barcode_ids`
  - consumed
  - mapped through the sample probe catalog into whitelist and sample-probe TSVs

`[gene-expression]`
- `probe-set` / `probeset`
  - consumed
  - used to derive `probe_list.from_cr.txt`
- `reference`
  - recorded in the run manifest only
  - not used to set STAR `--genomeDir`

### STAR-suite Flex extensions

`[star]`
- `sample-probe-catalog` / `sample_probe_catalog`
  - STAR extension
  - path to the sample probe variant-to-canonical mapping catalog
- `sample-probe-offset` / `sample_probe_offset`
  - STAR extension
  - used to set `--soloSampleProbeOffset`

CLI overrides
- `--sample-probe-catalog`
  - fallback if the `[star]` section does not provide one
- `--sample-probe-offset`
  - overrides the config offset when passed explicitly

### Ignored or superseded for Flex

Ignored if present in the config:
- non-GEX library rows
  - currently not used by the Flex wrapper
- extra sample columns beyond `sample_id` and `probe_barcode_ids`
- other Cell Ranger run-control keys not explicitly listed above

Superseded by pinned STAR wrapper parameters:
- `--genomeDir`
- Flex whitelist path
- `--clipAdapterType CellRanger4`
- `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`
- `--soloUMIfiltering MultiGeneUMI_CR`
- `--soloUMIdedup 1MM_CR`
- `--soloMultiMappers Rescue`
- `--soloCellFilter None`
- `--flex yes`
- `--soloFlexExpectedCellsPerTag 3000`

## Intentional Non-Goals

These wrappers do not attempt to import all Cell Ranger semantics from the
config file.

Current non-goals:

- auto-switching STAR `--genomeDir` from the Cell Ranger reference path
- blindly inheriting Cell Ranger chemistry or run-mode flags
- full generic support for every Cell Ranger config field
- generic multisample orchestration beyond what a single config already
  expresses

## Validation Coverage

Perturb:
- [test_ucsf_cr_config_input.sh](/mnt/pikachu/STAR-suite/tests/test_ucsf_cr_config_input.sh)
- [test_ucsf_batch_cr_config_input.sh](/mnt/pikachu/STAR-suite/tests/test_ucsf_batch_cr_config_input.sh)
- [test_ucsf_cr_config_multifeature_input.sh](/mnt/pikachu/STAR-suite/tests/test_ucsf_cr_config_multifeature_input.sh)
- [test_ucsf_batch_cr_config_tiny_smoke.sh](/mnt/pikachu/STAR-suite/tests/test_ucsf_batch_cr_config_tiny_smoke.sh)
- [test_ucsf_batch_cr_config_multifeature_tiny_smoke.sh](/mnt/pikachu/STAR-suite/tests/test_ucsf_batch_cr_config_multifeature_tiny_smoke.sh)
- [run_ucsf_cr_config_1m_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_ucsf_cr_config_1m_smoke.sh)

Flex:
- [test_flex_cr_config_input.sh](/mnt/pikachu/STAR-suite/tests/test_flex_cr_config_input.sh)
- [run_flex_cr_config_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_flex_cr_config_smoke.sh)
- [run_flex_cr_config_star_section_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_flex_cr_config_star_section_smoke.sh)
