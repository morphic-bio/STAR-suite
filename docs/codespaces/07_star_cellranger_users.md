# If you already use STAR or Cell Ranger

This page is for users who already know the tools, but still want simple setup instructions.

## You can usually skip straight to a guide

- If you already have a STAR index:
  - skip setup
  - go straight to [Bulk](./02_bulk.md) or [SLAM](./03_slam.md)
- If you want to try perturb or Flex:
  - build the small shared input with [Single-cell input for perturb and Flex](./04_single_cell_fixture.md)
  - then open the perturb or Flex guide

## What is different here

- the guides download small public datasets for you
- the guides write `RUN_COMMAND.sh` so you can see the exact command
- the small single-cell demos use a tiny human reference built from `chr22` and `chrY`
- perturb and Flex are still being finished, so those guides are there mainly to show the current command shape

## Cell Ranger note

These guides do not use `cellranger mkref`.
They build the small demo reference directly with STAR.
