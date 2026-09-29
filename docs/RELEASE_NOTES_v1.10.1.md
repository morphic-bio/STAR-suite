# STAR Suite v1.10.1 Release Notes

Unreleased. These changes were merged into `dev-release-v1.10.0` after the
`v1.10.0` source was fixed and are not part of the `v1.10.0` tag.

## Launchpad fixes merged from 1.9.5.b

The portable config/forms, relocatable installed launcher, bundled browser
assets, runtime capability report and release test gate are included. The
pinned catalog's legacy Multiome recipe gains input/runtime validation and
managed jobs/logs/cancellation for an explicitly selected external 1.9.5.b
runtime. Standalone STAR 1.10 reports `chromap_atac:false` and cannot execute it.
The deleted built-in `morphic_multiome` schema and all removed native Chromap
integration remain deleted. Current Multiome ownership stays in Multiomics
Suite. The 1.10 version and existing immutable rc tags are unchanged.

## Multiome compatibility guide

The tested STAR 1.9 Multiome setup and build commands are documented in
[Multiome ownership and compatibility](LAUNCHPAD_MULTIOME.md).
