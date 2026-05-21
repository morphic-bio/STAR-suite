# MSK 30KO GemX Production Draft Runbook

Date: 2026-05-13
Status: planning only; do not launch runs from this document until the sample
set and output locations are confirmed.

## Goal

Process the MSK 30KO single-cell perturb datasets in a UCSF-production-style
workflow, with MSK-specific changes:

- no Y-chromosome removal
- keep BAM output, but upload BAMs by Globus because they are large
- do not upload FASTQs
- run CellBender on the remote GPU server, as in UCSF production
- include LARRY barcode outputs and compare against provider annotations
- use provider cell types from `/mnt/pikachu/df.meta.rds` to seed the
  Scimilarity/cell-typing path

## Primary Inputs

| Purpose | Path | Notes |
| --- | --- | --- |
| Provider analysis metadata | `/mnt/pikachu/df.meta.rds` | `data.frame`, 132087 x 40. Authoritative reference for dataset identities, cell types, provider calls, and comparison groupings. |
| FASTQ manifest | `docs/MSK_30KO_FASTQ_MANIFEST.tsv` | FASTQ-only manifest for the nine provider `orig.ident` groups and GEX, LARRY, CRISPR PolyIII libraries. |
| New DE GemX FASTQs | `/mnt/pikachu/scRNAseq_30KO_DE_GEM_X` | Corrected update for the old `MSK_30_KO` DE GemX surface, supplied because the previous GemX set was missing LARRY barcodes. Size: ~80G. |
| gRNA feature reference, STAR generic | `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv` | 29/30 guide-style features, currently `Custom`. |
| gRNA feature reference, CRISPR type | `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC_crispr.csv` | Use for `CRISPR Guide Capture` routing. |
| LARRY feature reference | `/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv` | 245,979 LARRY features. |
| STAR genome | `/storage/autoindex_110_44/bulk_index` | Same MSK benchmark STAR index. |
| Solo TRU whitelist, non-GemX | `/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt` | Non-GemX GEX and LARRY use the February-2018 TRU whitelist family. |
| NXT whitelist, non-GemX | `/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt` | Non-GemX PolyIII/gRNA uses the February-2018 NXT whitelist family. |
| Solo TRU whitelist, GemX | `/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt` | GemX GEX and LARRY use the May-2023/GEM-X TRU whitelist family. Staged from Cell Ranger 9.0.1. |
| NXT whitelist, GemX | `/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_NXT.txt` | GemX PolyIII/gRNA uses the May-2023/GEM-X NXT whitelist family. Staged from Cell Ranger 9.0.1 translation whitelist column 1. |
| STAR binary | `core/legacy/source/STAR` or release-pinned binary | Pin before production. Clean rebuild first if using source binary. |

Historical-village outputs under `/mnt/pikachu/MSK_KO_village_cell_typing` are
not primary inputs for this run. They can be used for context only.

Historical `.h5ad` files are not production inputs for this run. Build the
production surface from FASTQs and use new downstream `.h5ad` files only as
outputs from the planned workflow.

When local historical outputs and the new FASTQ/update surface disagree, treat
the new FASTQ/update surface as the newer correct version.

## Existing Code To Reuse

| Purpose | Path |
| --- | --- |
| MSK 3-library STAR precedent | `scripts/paper/run_msk_30polyko_benchmark.sh` |
| MSK FASTQ whitelist-family preflight | `scripts/run_msk_30ko_fastq_preflight.sh` |
| Generic whitelist-family preflight | `scripts/preflight_whitelist_family.py` |
| GEX/feature pairing preflight | `scripts/preflight_library_pairing.py` |
| UCSF production orchestration precedent | `scripts/paper/run_ucsf_corrected_production_workflow.sh` |
| UCSF per-sample STAR/downstream/Globus precedent | `scripts/run_ucsf_perturb_yremove_batch.sh` |
| Downstream h5ad + feature integration | `scripts/run_scrna_downstream_gene_full_velocyto.sh` |
| Remote GPU CellBender watcher | `scripts/run_remote_cellbender_scan.sh` |
| Remote CellBender one-shot helper | `scripts/run_remote_cellbender_rsync.sh` |
| Scimilarity/cell typing repo | `/mnt/pikachu/ucsf-perturb-analysis` |
| Scimilarity GPU container | `/mnt/pikachu/ucsf-perturb-analysis/containers/scimilarity-gpu` |
| MSK scRNA-seq widget precedent | `/mnt/pikachu/scRNA-seq/workflows/scRNA_seq_features_msk_7_16` |

The current MSK benchmark wrapper does not yet implement UCSF-style per-sample
downstream, remote CellBender, or BAM-only Globus cleanup. Plan to create/adapt
an MSK production wrapper rather than running the paper benchmark wrapper as-is.

## New DE GemX Grouping

New source directory:

`/mnt/pikachu/scRNAseq_30KO_DE_GEM_X`

Observed FASTQ groups:

| Logical library | FASTQ sample id | Lanes | File roles |
| --- | --- | --- | --- |
| GEX | `mRNA_DE_30KO_ZN_IGO_16692_C_5_S23` | L001-L008 | I1/I2/R1/R2 |
| GEX | `mRNA_DE_30KO_ZN_IGO_16692_C_5_S33` | L001-L008 | I1/I2/R1/R2 |
| GEX | `mRNA_DE_30KO_ZN_IGO_16692_D_5_S4` | L001-L008 | I1/I2/R1/R2 |
| GEX | `mRNA_DE_30KO_ZN_IGO_16692_D_5_S53` | L005-L008 | I1/I2/R1/R2 |
| gRNA / PolyIII | `PolyIII_DE_30KO_ZN_IGO_16692_C_8_S54` | L001-L008 | I1/I2/R1/R2 |
| LARRY | `LARRY_DE_30KO_ZN_IGO_16692_C_6_S37` | L001-L008 | I1/I2/R1/R2 |

Proposed provider label mapping:

- production sample label: `30_KO_DE_GemX`
- provider metadata `orig.ident`: `DE_GemX`
- provider stage: `S1`

This is the replacement/update for the older `/mnt/pikachu/MSK_30_KO` GemX
surface. It was provided because the earlier GemX material was missing LARRY
barcodes. It should be treated as the correct GemX source when there is a
conflict with older local outputs.

There are two DE datasets in the provider reference:

- `DE`: older DE/XM-style surface
- `DE_GemX`: separate GemX run, produced by a different technician and kit

The provider reports no noticeable batch effect for the GemX surface, so the
two DE groups should remain separate provider-reference datasets while still
being eligible for joint downstream interpretation where appropriate.

## Authoritative Dataset Set

The production datasets should match the nine `orig.ident` groups in
`/mnt/pikachu/df.meta.rds`. Treat the `.rds` as the source of truth for
dataset identities and comparison groups.

| Provider group | Stage | Cells in `df.meta.rds` | Notes |
| --- | --- | ---: | --- |
| `ES` | S0 | 22445 | Older ES raw source exists. |
| `DE` | S1 | 17960 | Older DE/XM raw source; distinct from GemX. |
| `DE_GemX` | S1 | 28496 | New `/mnt/pikachu/scRNAseq_30KO_DE_GEM_X` source is primary/corrected. |
| `PP1` | S3 | 6310 | Older PP source exists. |
| `PP2` | S4 | 12072 | Older PP source exists. |
| `S5_1` | S5 | 13397 | Older S5/S6 source exists. |
| `S5_2` | S5 | 13601 | Older S5/S6 source exists. |
| `S6_1` | S6 | 8823 | Older S5/S6 source exists. |
| `S6_2` | S6 | 8983 | Older S5/S6 source exists. |

Do not collapse these to six biological/timepoint groups for production. Any
merged or timepoint-level summaries can be derived after per-`orig.ident`
processing and provider comparison.

## FASTQ Manifest And GemX Inspection

Use `docs/MSK_30KO_FASTQ_MANIFEST.tsv` as the working FASTQ manifest. It has 27
rows: nine provider groups x three libraries (`GEX`, `CRISPR_PolyIII`, `LARRY`).
Each row records the FASTQ root, sample IDs, lanes, file roles, chemistry,
feature reference, and whitelist.

Run the whitelist-family preflight before generating any STAR command:

```bash
scripts/run_msk_30ko_fastq_preflight.sh \
  --manifest docs/MSK_30KO_FASTQ_MANIFEST.tsv \
  --outdir plans/artifacts/msk_30ko_fastq_preflight_YYYYMMDD
```

Expected result: all 27 rows pass, with `DE_GemX` rows resolving to
`may2023_gemx:*` and all other rows resolving to `feb2018:*`.

FASTQ inspection notes:

- GemX has no sidecar metadata files in `/mnt/pikachu/scRNAseq_30KO_DE_GEM_X`
  from the inspected directory contents; the available primary evidence is the
  FASTQ names and headers.
- GemX, older DE/XM, PP, and S5/S6 representative FASTQs all show four file
  roles per lane: `I1`, `I2`, `R1`, `R2`.
- Representative read lengths are consistent across inspected sources:
  `I1=10`, `I2=10`, `R1=29`, `R2=89`.
- GemX differs in whitelist family. Sampling the first 200,000 R1 reads per
  representative FASTQ showed GemX GEX and LARRY match the May-2023/GEM-X TRU
  whitelist, not the February-2018 TRU whitelist. GemX PolyIII matches the
  May-2023/GEM-X NXT whitelist.
- Representative GemX exact whitelist hit rates from R1 first 16 bp:
  GEX `90.22-91.01%` for May-2023 TRU and about `0.43-0.47%` for February-2018
  TRU; PolyIII `97.59%` for May-2023 NXT; LARRY `96.55%` for May-2023 TRU.
- Older DE/XM controls showed the expected February-2018 family behavior:
  GEX `92.56%` February-2018 TRU, PolyIII `97.69%` February-2018 NXT, and
  LARRY `97.43%` February-2018 TRU.
- GemX is still distinct operationally: it uses the `ZN` sample naming and the
  `/mnt/pikachu/scRNAseq_30KO_DE_GEM_X` corrected FASTQ root. Header inspection
  shows GemX GEX split across NovaSeq X flowcells including `22FL2LLT3`,
  `22FLFMLT3`, and `22FL2FLT3`; GemX LARRY and PolyIII are on `22FLFMLT3` in
  the inspected lanes.
- Older ES/DE/XM FASTQs also use the `LH00288` / `IGO_16692` run family, while
  PP/S5/S6 FASTQs use `LH00834` / `IGO_17014` run families. This looks like a
  run/source difference, not a different FASTQ layout.
- The PP1, PP2, and S6_1 LARRY rows have two sample IDs each on distinct
  flowcells with matching sample-index pairs. Treat these as multi-flowcell
  LARRY inputs unless the demultiplexing sheet contradicts this.

## Older Raw Sources For Mapping

Use these only after the final dataset list is confirmed:

| Raw source | Candidate groups |
| --- | --- |
| `/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM` | `DE`, `ES`; each has GEX + PolyIII + LARRY. |
| `/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_PP1_PP2` | `PP1`, `PP2`; GEX + PolyIII + LARRY. |
| `/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_S5_S6` | `S5_1`, `S5_2`, `S6_1`, `S6_2`; directory also contains some PP/DE carryover FASTQs, so select by sample ID, not by whole-directory glob. |
| `/mnt/pikachu/scRNAseq_30KO_DE_GEM_X` | `DE_GemX`; corrected update with LARRY barcodes. |
| `/mnt/pikachu/MSK_30_KO` | Historical `.h5ad` outputs only; do not use for this FASTQ-first production run. |

## STAR Production Surface

Start from the MSK paper benchmark parameter surface, with production changes:

- `--soloFeatures GeneFull Velocyto` if we want the UCSF downstream wrapper
  unchanged; otherwise add a GeneFull-only downstream wrapper variant.
- `--outSAMtype BAM Unsorted` so BAMs can be transferred.
- no `--emitNoYBAM`, no `--emitYNoYFastq`, no Y/noY FASTQs.
- `--pfMultiConfig` per dataset with three logical libraries:
  - GEX: `Gene Expression`, chemistry `TRU`, no feature ref
  - PolyIII/gRNA: `CRISPR Guide Capture`, chemistry `NXT`,
    `ref_feature_geneBC_crispr.csv`, NXT whitelist
  - LARRY: `Custom`, chemistry `TRU`, `ref_feature_larryBC.csv`, TRU whitelist
- Use the correct whitelist family per dataset:
  - non-GemX groups: February-2018 TRU/NXT whitelists
  - `DE_GemX`: May-2023/GEM-X TRU/NXT whitelists
- keep paper-critical perturb flags:
  - `--soloCrMultimapRescue yes`
  - `--dynamicThreadInterface 1`
  - `--crAssignConsumerThreads -1`
  - `--crAssignSearchThreads 1`
  - `--soloCrGexFeature genefull`
  - `--clip3pPolyG yes`

Run each logical dataset into a fresh output directory. Do not run benchmark or
production jobs in parallel on this host.

## Downstream, CellBender, And Uploads

Per dataset, target local output layout should mirror UCSF production:

```text
<run_root>/samples/<dataset>/
  run/
    Aligned.out.bam
    outs/raw_feature_bc_matrix/
    outs/filtered_feature_bc_matrix/
    outs/crispr_analysis/
    cr_assign/
  downstream_genefull_velocyto_cellbender/
    counts.h5ad
    unfiltered_counts.h5ad
    filtered_counts.h5ad
    default_singlet_filtered_counts.h5ad
    final_counts.h5ad
    cellbender/cellbender_counts.h5
    feature_libraries/
  globus_batch.tsv
  RUN_GLOBUS_TRANSFER.sh
```

Globus policy:

- upload `run/Aligned.out.bam` only
- do not upload source FASTQs
- do not upload Y/noY FASTQs because they are not produced
- after a successful Globus task, local BAM cleanup is allowed if the manifest
  and task id are preserved

Remote CellBender:

- use `scripts/run_remote_cellbender_scan.sh` or
  `scripts/run_remote_cellbender_rsync.sh`
- UCSF precedent remote host: `10.159.4.53`
- remote staging root and GPU slot policy need confirmation before launch
- keep adaptive QC enabled in downstream prep

## Provider Metadata Comparison

Use `/mnt/pikachu/df.meta.rds` as the provider reference. Relevant columns:

- `orig.ident`
- `stage`
- `celltype_20250623`
- `subcelltype`
- `larryBC`, `larryUMI`, `larryBCCategory`, `larryBC_type`,
  `larryBCpresent`
- `geneBC`, `geneBCUMI`, `geneBCCategory`, `geneBC_type`, `geneBCpresent`
- `Qualified`, `Qualified_BC`, `Qualified_BC_group`

Provider row names look like:

```text
AAACCAAAGATTCCAG-DE_GemX
```

For STAR outputs, derive join keys as:

- `canonical_cb`: strip trailing `-1` from STAR barcodes
- `provider_group`: map production dataset label to provider `orig.ident`
- provider key: split row name into `canonical_cb` and provider group

Initial comparisons to plan:

- STAR filtered cell set vs provider `Qualified`/`singlet`/group-specific cells
- STAR gRNA calls vs provider `geneBC`, `geneBCUMI`, `geneBCpresent`
- STAR LARRY calls vs provider `larryBC`, `larryUMI`, `larryBCpresent`
- downstream cell types seeded/transferred from `celltype_20250623` and
  `subcelltype`

## Scimilarity Seeding Plan

Use provider labels as the seed/reference surface, not the village metadata.

Suggested reference seed:

- input: `/mnt/pikachu/df.meta.rds`
- primary label: `celltype_20250623`
- secondary label: `subcelltype`
- stage/time label: `stage`
- batch/sample label: `orig.ident`

Candidate implementation paths:

- `/mnt/pikachu/ucsf-perturb-analysis/pipeline/cell_typing/`
- `/mnt/pikachu/ucsf-perturb-analysis/pipeline/reference_prep/`
- `/mnt/pikachu/ucsf-perturb-analysis/containers/scimilarity-gpu/`
- `/mnt/pikachu/scRNA-seq/workflows/scRNA_seq_features_msk_7_16/widgets/scRNA_seq_features_msk_7_16/Assign_cell_type/Dockerfiles/create_reference_from_subset.py`

Before running Scimilarity, confirm whether the provider metadata alone is
sufficient or whether the provider expression object is also required. Only
`df.meta.rds` was identified as the compact provider reference in this pass.

## Preflight Checklist Before Any Run

1. Use the nine `orig.ident` groups in `/mnt/pikachu/df.meta.rds` as the final
   production dataset list.
2. Keep `DE` and `DE_GemX` as separate datasets; `DE_GemX` is a different
   technician/kit GemX run with no noticeable batch effect and supersedes older
   local GemX outputs when conflicts appear.
3. Confirm GemX uses the May-2023/GEM-X whitelist family in generated
   `pfMultiConfig` and `--soloCBwhitelist`; do not run GemX against the
   February-2018 TRU/NXT whitelists.
4. Run `scripts/run_msk_30ko_fastq_preflight.sh` and require all 27 manifest
   rows to pass before launching STAR.
5. Confirm STAR binary pin and do a clean rebuild if using `core/legacy/source/STAR`.
6. Decide whether MSK production emits Velocyto so the UCSF downstream wrapper
   can be reused unchanged.
7. Add or adapt an MSK production wrapper for:
   - no Y-removal
   - BAM-only Globus transfer
   - remote GPU CellBender
   - provider metadata comparison hooks
8. Confirm Globus source endpoint, destination endpoint, and destination root.
9. Confirm remote CellBender host, remote root, GPU slots, and image policy.
10. Run one small dry-run/smoke on DE GemX before full production.
