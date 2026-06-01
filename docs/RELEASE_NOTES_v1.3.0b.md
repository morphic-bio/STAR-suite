# STAR Suite v1.3.0b Release Notes

Date: 2026-06-01

`v1.3.0b` is a minor feature-label refresh of `v1.3.0`, not a bug-fix
release. The artifact version is `1.3.0b`, and `STAR --version` reports
`1.3.0b`. The upstream STAR base remains `2.7.11b`, and the genome index
compatibility string remains `2.7.4a`.

## Scope

- Keep the `v1.3.0` parallel CBQ integration feature set.
- Add an explicit STAR-Flex Fixed RNA CBQ MCP recipe:
  `star_flex_fixed_rna_cbq`.
- Update recipe metadata so FLEX FASTQ recipes default to STAR internal gzip
  rather than the external `zcat` path, with legacy mode still available as an
  explicit opt-in.
- Document and expose native multiome CBQ recipe inputs for both GEX and ATAC:
  `--gex-cbq`, `--atac-read-pair-cbq`, and `--atac-barcode-cbq`.
- Update OCM CBQ recipe metadata for CBQ Y/noY sidecar support through
  `--emitYNoY yes --emitYNoYFormat cbq`.

## Validation

- Parsed updated MCP workflow YAMLs.
- Ran shell syntax checks for updated recipe wrappers.
- Verified Debian changelog version parsing for `1.3.0b-1`.
- Rebuilt STAR and verified `STAR --version` reports `1.3.0b`.
