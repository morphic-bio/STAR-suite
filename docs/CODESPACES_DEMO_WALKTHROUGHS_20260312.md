# Codespaces Demo Walkthroughs

This document defines the Codespaces demo surface for the four major STAR-suite modules:

1. STAR-core bulk
2. STAR-perturb
3. STAR-Flex
4. STAR-SLAM

The goal is not to hide STAR. The demos assume the user already understands the normal STAR run model and can replace the default public inputs with their own indices and FASTQs at any point.

## Scope

- Demos live under `.codespaces-demo/` and are intentionally untracked.
- Bulk and SLAM use small public SRA subsets.
- Perturb and Flex now share a real public human single-cell source rather than the earlier `universc` tiny placeholder.
- The canonical demo reference is the source `FASTA + GTF` pair under `.codespaces-demo/data/public_human_chr22y_ref/`. STAR builds the runnable demo index from that source surface directly. We do not introduce `cellranger mkref` into the Codespaces walkthrough path.
- Barcode geometry is treated as assay-specific source metadata, not a global constant. The pinned public 10x source currently uses the defaults in `scripts/codespaces/public_demo_sources.sh`, and the walkthrough wrappers expose overrides if a different public assay is substituted.
- The shared single-cell demo path is:
  1. fetch public male 10x GEX FASTQs
  2. fetch public `chr22` and `chrY` reference sequences plus GTF
  3. build a tiny `chr22 + chrY` STAR index in Codespaces
  4. align a bounded number of public reads to that mini-reference
  5. retain only raw read pairs whose cDNA mate maps to `chr22` or `chrY`
  6. generate the perturb guide companion or Flex probe/sample-tag overlay on top of that real public barcode/read-name surface
- This keeps the walkthroughs small enough for Codespaces while avoiding a fully synthetic barcode surface.

## Current Validation Status

- The shared public `chr22 + chrY` reference fetch/build path is working.
- The public GEX fixture derivation from real 10x FASTQs is working.
- The bulk and SLAM walkthrough surfaces are structurally sound.
- The current synthetic perturb and Flex overlays are useful as wiring/prototype scaffolds, but they are not yet strong enough to present as faithful assay demos:
  - perturb emits the expected CR-style outputs, but the generated guide companion does not produce meaningful guide counts
  - Flex runs, but the synthetic probe overlay is too weak and does not yield a useful `per_sample` result on the current mini-reference surface

That means the clean next step for a polished Codespaces experience is:

1. keep the shared public `chr22 + chrY` reference/index-creation walkthrough
2. keep the public GEX fixture derivation walkthrough
3. keep bulk and SLAM on real public datasets
4. replace the synthetic perturb and Flex overlays with small real public assay-specific datasets before calling those two demos complete

The first concrete step toward that replacement is already pinned:

- perturb candidate: official 10x K562 CRISPR Guide Capture dataset
- Flex candidate: official 10x K562 Fixed RNA Profiling dataset

Small metadata fetchers for those sources are available:

- `scripts/codespaces/fetch_public_perturb_assets.sh`
- `scripts/codespaces/fetch_public_flex_assets.sh`

These download the official `config.csv` and feature/probe reference files and pin the raw FASTQ tar URLs. The remaining open problem is writing a bounded stream-extract path for those large public raw tar archives so Codespaces can derive small runnable fixtures without downloading the full tarballs.

## VS Code Walkthrough Surface

- In Codespaces, `.devcontainer/postCreate.sh` packages and installs the repo-local walkthrough extension on a best-effort basis. The fallback entry point is Command Palette -> `STAR-suite: Open Codespaces Walkthrough`.

Codespaces now includes a repo-local VS Code walkthrough extension under `tools/vscode-star-suite-walkthrough/`. The extension presents the demo surface as independent chapters:

- one optional shared foundation chapter for mini-reference creation
- standalone module chapters for bulk and SLAM
- a shared utility chapter for the single-cell fixture used by perturb and Flex
- placeholder perturb and Flex chapters that keep the walkthrough shape stable while their real public assay fixtures are being tightened

Structured module guides also live under `docs/codespaces/`:

- `00_overview.md`
- `01_setup_reference.md`
- `02_bulk.md`
- `03_slam.md`
- `04_single_cell_fixture.md`
- `05_perturb.md`
- `06_flex.md`

## One-Time Setup

Open the repo in Codespaces. The devcontainer installs:

- `git-lfs`
- `sra-tools`
- `samtools`
- `pigz`
- build toolchain
- Python + pandas

If the STAR binary is missing, the walkthrough scripts build it automatically.

The intended user flow is non-linear:

- start with the foundation chapter only if you need the shared mini-reference
- go directly to Bulk or SLAM otherwise
- use the single-cell fixture chapter only when you want the bounded public input for perturb or Flex

Recommended chapter picks:

- users who already know STAR or Cell Ranger can start with `docs/codespaces/07_star_cellranger_users.md` and then jump straight to the relevant chapter
- Experienced users:
  - go directly to Bulk or SLAM if you already have an index and only want the command shape
  - use the single-cell fixture chapter only if you need the bounded public input used by perturb or Flex
  - treat perturb and Flex as prototype chapters for now
- Novice users:
  - start with the optional foundation chapter if you want to see the reference build once
  - then use Bulk as the first stable end-to-end demo
  - then try SLAM

## Demo Matrix

| Module | Data source | Fully runnable in Codespaces | Main helper |
| --- | --- | --- | --- |
| STAR-core bulk | GEO/SRA `GSE88509 / GSM2344101 / SRR4422207` | Yes, after the user points to a human STAR index | `scripts/codespaces/run_bulk_public_demo.sh` |
| STAR-perturb | Public male 10x GEX + derived `chr22 + chrY` fixture + generated guide companion | Yes | `scripts/codespaces/run_perturb_public_demo.sh` |
| STAR-Flex | Public male 10x GEX + derived `chr22 + chrY` fixture + generated chr22/chrY probe/sample-tag overlay | Yes | `scripts/codespaces/run_flex_public_demo.sh` |
| STAR-SLAM | SRA `SRR32576116` subset | Yes, after the user points to a human STAR index | `scripts/codespaces/run_slam_public_demo.sh` |

## Shared Public Inputs

### Public male 10x GEX surface

Pinned by:

- Dataset page: `https://www.10xgenomics.com/datasets/human-b-cells-from-a-healthy-donor-1-k-cells-2-standard-6-0-0`
- Raw FASTQs: `https://cf.10xgenomics.com/samples/cell-vdj/6.0.0/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar`

Fetched and reduced by:

```bash
bash scripts/codespaces/derive_public_chr22y_gex_fixture.sh --dry-run
```

Outputs:

- `.codespaces-demo/data/public_10x_male_bcell/raw_fastqs`
- `.codespaces-demo/data/public_chr22y_gex_fixture/gex/public_chr22y_demo`

Current pinned barcode geometry for this public source:

- `soloCBstart=1`
- `soloCBlen=16`
- `soloUMIstart=17`
- `soloUMIlen=10`

These values are demo-source metadata for the chosen public assay. They are not treated as global STAR-suite defaults.

### Public human `chr22 + chrY` mini-reference

This is the single source reference surface for the single-cell demos. The same `genome.fa` and `genes.gtf` files can be treated as the CR-compatible source inputs, while STAR builds the executable demo `genomeDir` directly from them.

Pinned by:

- Ensembl release `110`
- `chr22` FASTA
- `chrY` FASTA
- GRCh38 GTF

Fetched by:

```bash
bash scripts/codespaces/fetch_public_chr22y_reference.sh
```

Index built by:

```bash
bash scripts/codespaces/build_public_chr22y_index.sh
```

Outputs:

- `.codespaces-demo/data/public_human_chr22y_ref/fasta/genome.fa`
- `.codespaces-demo/data/public_human_chr22y_ref/genes/genes.gtf`
- `.codespaces-demo/indices/public_human_chr22y_star`

### Public bulk fixture

Pinned by the existing downloader:

- GEO series: `GSE88509`
- GEO sample: `GSM2344101`
- SRA run: `SRR4422207`

### Public SLAM fixture

Pinned by the SLAM helper:

- SRA run: `SRR32576116`

### Real public perturb replacement candidate

Pinned by:

- official dataset page: `10k K562 Transduced with Small Guide Library, 5 HT v2, Chromium X`
- official raw FASTQ tar
- official `config.csv`
- official feature-reference CSV

Fetch the small metadata assets with:

```bash
bash scripts/codespaces/fetch_public_perturb_assets.sh
```

### Real public Flex replacement candidate

Pinned by:

- official dataset page: `10k Human K562-r Cells, Singleplex Sample 1, Standard`
- official raw FASTQ tar
- official `config.csv`
- official probe-set CSV

Fetch the small metadata assets with:

```bash
bash scripts/codespaces/fetch_public_flex_assets.sh
```

## Walkthrough 1: STAR-core Bulk

This is the smallest real public bulk demo in the repo.

Prepare the fixture and emit the command:

```bash
bash scripts/codespaces/run_bulk_public_demo.sh \
  --dry-run \
  --genome-dir /path/to/human_star_index
```

Run it for real:

```bash
bash scripts/codespaces/run_bulk_public_demo.sh \
  --run \
  --genome-dir /path/to/human_star_index
```

What it does:

- downloads `SRR4422207`
- trims with native `--trimCutadapt Yes`
- runs a standard bulk alignment with `GeneCounts`

Outputs:

- `.codespaces-demo/runs/bulk_public_demo/RUN_COMMAND.sh`
- `.codespaces-demo/runs/bulk_public_demo/Aligned.sortedByCoord.out.bam`
- `.codespaces-demo/runs/bulk_public_demo/ReadsPerGene.out.tab`

## Walkthrough 2: STAR-perturb

This demo uses a real public male 10x GEX barcode/read-name surface, but reduces it to the reads whose cDNA mates map to `chr22` or `chrY`. It then generates a small CRISPR guide-capture companion library on top of that surface.

Run it:

```bash
bash scripts/codespaces/run_perturb_public_demo.sh
```

Dry-run only:

```bash
bash scripts/codespaces/run_perturb_public_demo.sh --dry-run
```

What it does:

- fetches the public `chr22 + chrY` mini-reference
- builds the corresponding STAR index in Codespaces
- derives a small public GEX fixture by aligning a bounded number of reads from the public male 10x source and keeping only read pairs with `chr22/chrY` cDNA mappings
- generates:
  - a small CRISPR feature ref
  - a small guide-capture companion library
  - a whitelist from the public barcode reads
  - a CR-style `config.csv`
- runs the hardened `--cr-config` perturb wrapper

Outputs:

- `.codespaces-demo/runs/perturb_public_demo/assets/config.csv`
- `.codespaces-demo/runs/perturb_public_demo/run/outs/filtered_feature_bc_matrix`
- `.codespaces-demo/runs/perturb_public_demo/run/crispr_analysis`

Scale-up path:

- replace the generated `config.csv` with a real public or local Cell Ranger `config.csv`
- point the same wrapper at a full STAR index and real FASTQ directories

## Walkthrough 3: STAR-Flex

This demo shares the same public `chr22 + chrY` GEX fixture as perturb, but overlays a tiny probe/sample-tag surface using a small probe set chosen from `chr22` and `chrY` genes. That keeps the demo small and reproducible while making index creation part of the walkthrough.

Run it:

```bash
bash scripts/codespaces/run_flex_public_demo.sh
```

Dry-run only:

```bash
bash scripts/codespaces/run_flex_public_demo.sh --dry-run
```

What it does:

- fetches the public `chr22 + chrY` mini-reference
- derives the shared public GEX fixture from the male 10x source
- chooses a tiny probe set from `chr22` and `chrY`
- generates:
  - a small sample-probe catalog
  - a small synthetic Flex GEX/read2 overlay using the real public barcode/read-name surface
  - a CR-compatible Flex `config.csv`
- builds a tiny hybrid Flex reference with the existing Flex helper scripts
- runs the hardened `--cr-config` Flex wrapper

Outputs:

- `.codespaces-demo/runs/flex_public_demo/assets/config.csv`
- `.codespaces-demo/runs/flex_public_demo/flex_reference/filtered_reference`
- `.codespaces-demo/runs/flex_public_demo/run/per_sample/flexfilter_summary.tsv`

Scale-up path:

- keep the same wrapper shape
- swap in a real public or local Flex `config.csv`, probe set, and STAR index

## Walkthrough 4: STAR-SLAM

This is the smallest public SLAM entry point.

Prepare the fixture and emit the command:

```bash
bash scripts/codespaces/run_slam_public_demo.sh \
  --dry-run \
  --genome-dir /path/to/human_star_index
```

Run it for real:

```bash
bash scripts/codespaces/run_slam_public_demo.sh \
  --run \
  --genome-dir /path/to/human_star_index
```

What it does:

- downloads `SRR32576116`
- runs `--slamQuantMode 1`
- writes GRAND-SLAM compatible output and a QC report

Outputs:

- `.codespaces-demo/runs/slam_public_demo/RUN_COMMAND.sh`
- `.codespaces-demo/runs/slam_public_demo/SlamQuant.grandslam.tsv`
- `.codespaces-demo/runs/slam_public_demo/slam_qc.html`

## Updatable Parts

The walkthroughs are designed so the public source pins are easy to update:

- `scripts/codespaces/public_demo_sources.sh`
- `scripts/codespaces/fetch_public_chr22y_reference.sh`
- `scripts/codespaces/fetch_public_10x_male_fastqs.sh`
- `scripts/codespaces/derive_public_chr22y_gex_fixture.sh`
- `scripts/codespaces/fetch_public_slam_fixture.sh`
- `scripts/download_public_bulk_fixture.sh`

## Boundaries

- Bulk and SLAM still expect a human STAR index to be supplied or imported.
- Perturb and Flex are intentionally Codespaces-sized, not paper-sized.
- The public single-cell fixture builder is bounded by `--source-read-limit`; it is a demo surface, not a parity benchmark.
- The perturb and Flex demos are module walkthroughs, not parity benchmarks.
- Users who want biological parity should replace the demo inputs with their real public or local datasets.
