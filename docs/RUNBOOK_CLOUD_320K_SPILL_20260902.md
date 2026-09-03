# Runbook: Cloud 320k matrix — spill-free benchmarks on public data (r6id.8xlarge)

Date drafted: 2026-09-02 (rev 2 — expanded from the two-run spill experiment to the
full matrix per user direction). Status: GATED on the bucket-tail release — STAR
runs use the **released precompiled binary**, never a dev build. Executes on the
operator's AWS account.

## Purpose

One instance session produces the definitive spill-free 320k table — every tool,
both input formats, identical hardware, fully public data — plus two bonus claims:

1. **The complete matrix** (7 runs + 1 optional): cyto and STAR-Flex on both
   fastq.gz and CBQ; STAR-Flex with BGZF detection off (serial gzip control);
   STAR-Flex hash-then-align-residual (the accuracy-maximal mode, timed); and
   Cell Ranger 9.0.1. All STAR runs in `--soloBucketMode ram` (256 GiB holds the
   ~190 GB working set) — a spill-free environment, as none of the lab rows can be.
2. **Worked reproducibility example** for the paper/letter: fresh machine, the
   released Debian package installed with a single apt command, data streamed
   from 10x's own storage
   (feeds `[[SYNC-CLOUD-TIME]]` = the STAR-Flex BGZF RAM run, and `[[SYNC-CLOUD-COST]]`).
3. **Optional spill-cost run** (#8): one spill-mode repeat gives the same-hardware
   RAM-vs-spill ratio for the bounded-memory sentence.

## Shape and cost

| Item | Value |
| --- | --- |
| Instance | **r6id.8xlarge**, us-west-2 (32 vCPU = lab thread count, 256 GiB, 1.9 TB NVMe) |
| AMI | Ubuntu 24.04 LTS (matches the glibc239 assets and the release's own container checks; enables the one-command apt install) |
| Price | ~$2.42/h on-demand |
| Session | setup ~45 min + runs (~5 h) + Cell Ranger (~8–12 h) + CBQ encode (~70 min) ≈ **16–20 h ≈ $40–50** |
| Transfer | 10x's tar pulls in-region, free and fast; lab staging is ~55 GB once |

## One-time staging from the lab (~55 GB)

```bash
aws s3 mb s3://<staging-bucket> --region us-west-2
# STAR-Flex references (~36 GB)
aws s3 sync /storage/flex_filtered_reference_2024/star_index  s3://<staging-bucket>/refs/star_index
aws s3 cp /storage/downsampled_100K/SC2300771/results/flex_h01_2024_20260320_081246/h01_cache.bin s3://<staging-bucket>/refs/
aws s3 cp /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt      s3://<staging-bucket>/refs/
aws s3 cp /mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv                 s3://<staging-bucket>/refs/
aws s3 cp /mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt s3://<staging-bucket>/refs/
# Cell Ranger: your own licensed install + reference (private staging of your own
# copy for your own use — not redistribution) (~18.5 GB)
aws s3 sync /home/lhhung/cellranger-9.0.1                    s3://<staging-bucket>/cellranger-9.0.1
aws s3 sync /mnt/pikachu/CR-references/refdata-gex-GRCh38-2024-A s3://<staging-bucket>/refs/refdata-gex-GRCh38-2024-A
aws s3 cp /mnt/pikachu/benchmark_cr9_320k_config.csv         s3://<staging-bucket>/refs/   # from the lab chain; paths rewritten on-instance
aws s3 cp /mnt/pikachu/process_features/tables/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv s3://<staging-bucket>/refs/   # REQUIRED by the CR config rewrite
# cyto (~20 MB): binary, its whitelists, the gene-probe table
aws s3 cp /home/lhhung/.cargo/bin/cyto                       s3://<staging-bucket>/tools/
aws s3 sync /home/lhhung/.cyto                               s3://<staging-bucket>/tools/dot-cyto
aws s3 cp /mnt/pikachu/cyto_headtohead_20260901/gene_probes_v1.1.0.tsv s3://<staging-bucket>/refs/
# CBQ encoder: NOT staged and NOT shipped in the package — the script builds it on
# the instance from the released source tarball before the (untimed) encode step.
```

## Instance setup

```bash
sudo mkfs.xfs /dev/nvme1n1 && sudo mkdir /scratch && sudo mount /dev/nvme1n1 /scratch
sudo chown $USER /scratch && mkdir -p /scratch/{fastqs,cbq,refs,tools,runs}
sudo apt-get update -qq && sudo apt-get install -y -qq awscli   # Ubuntu AMI ships without it

# Released package — the drop-in demonstration starts here: ONE apt command
curl -L -O https://github.com/morphic-bio/STAR-suite/releases/download/<RELEASE_TAG>/star-suite_<VERSION>-1_amd64.deb
sudo apt-get install -y ./star-suite_<VERSION>-1_amd64.deb
STAR --version   # record — installed system-wide by the released Debian package

aws s3 sync s3://<staging-bucket>/refs  /scratch/refs
aws s3 sync s3://<staging-bucket>/tools /scratch/tools && chmod +x /scratch/tools/{cyto,cbq_ordered_encoder}
aws s3 sync s3://<staging-bucket>/cellranger-9.0.1 /scratch/cellranger-9.0.1 && chmod -R +x /scratch/cellranger-9.0.1/bin /scratch/cellranger-9.0.1/cellranger
curl -LsSf https://astral.sh/uv/install.sh | sh   # cyto's cell-filter shells to uvx (accommodation, note in protocol)

# FASTQs streamed from 10x's public bucket; index reads dropped in flight (~15–30 min)
aws s3 cp --no-sign-request \
  s3://10x.files/samples/cell-exp/9.0.0/320k_scFFPE_16-plex_GEM-X_FLEX_Multiplex/320k_scFFPE_16-plex_GEM-X_FLEX_Multiplex_fastqs.tar - \
  | tar -x -C /scratch/fastqs --strip-components=1 --wildcards '*_R1_001.fastq.gz' '*_R2_001.fastq.gz'

{ date -u; curl -s http://169.254.169.254/latest/meta-data/instance-type; uname -a; } > /scratch/runs/preflight.txt
```

Shared STAR argument block: as in the lab run `scffpe320k_bgzf_20260902T173415Z.command.txt`
with `/scratch/refs/...` paths — define `$common` exactly as rev-1 of this runbook did
(threads 32, whitelists, probe list/offset 68, hash cache, `--soloBucketCount 256`,
flexPipeline, `--flexNoAlign 1`, dynamicThreadInterface). `$R1`/`$R2` from
`/scratch/fastqs`. Before EVERY timed run: `sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'`.
Every run under `/usr/bin/time -v`, its own dir under `/scratch/runs/<name>/`.

## The run matrix (execute in this order — disk plan depends on it)

| # | Name | Command delta | Expected |
| --- | --- | --- | --- |
| 1 | `flex_bgzf_ram` | `$common --soloBucketMode ram --readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1` | ~30–35 min, ~190 GiB — **the quoted [[SYNC-CLOUD-TIME]] run** |
| 2 | `flex_gzip_ram` | as #1 but `--readFilesBgzfMode off` (serial-decode control) | ~45–60 min |
| 3 | `flex_align_ram` | as #1 but `--flexNoAlign 0` (hash first, align the residual — accuracy-maximal; genome loads, RAM fine) | ~1.5–2.5 h; record the accuracy-mode time cost |
| 4 | `cyto_fastq` | `/scratch/tools/cyto workflow gex --gex /scratch/refs/gene_probes_v1.1.0.tsv -w /scratch/tools/dot-cyto/737K-fixed-rna-profiling.txt.gz -p /scratch/tools/dot-cyto/probe-barcodes-fixed-rna-profiling-rna.txt --probe-regex 'BC0(0[1-9]\|1[0-6])' --preset gex-v1 -T 32 -o /scratch/runs/cyto_fastq <8 fastq files>` with `TMPDIR=/scratch/tmp` | ~50–80 min, ~100 GiB |
| 5 | `cellranger` | rewrite the staged config's paths to `/scratch` (reference, probe-set CSV, fastqs), then `/scratch/cellranger-9.0.1/cellranger multi --id=cr9_320k --csv=... --localcores=32 --localmem=200` from `/scratch/runs/` | ~8–12 h — start before stopping for the day |
| 6 | encode CBQ (untimed) | after CR completes, free its scratch (`rm -rf /scratch/runs/cr9_320k/SC_MULTI_CS` keeping `outs/`), then 4 parallel `cbq_ordered_encoder --readFilesIn <R2> <R1> --outFile /scratch/cbq/lane_00N.cbq` (mate0 = R2) | ~70 min, ~200 GB |
| 7 | `flex_cbq_ram` | `$common --soloBucketMode ram --readFilesType Binseq PE --readFilesCbqRangeMode range --readFilesIn <4 cbq comma-joined>` (no BGZF flags) | ~30 min |
| 8 | `cyto_cbq` | as #4 but `-g '[gex][:18][probe] \| [barcode][umi:12]'` replacing the preset, inputs = 4 CBQ files | ~40–60 min |
| 9 (opt) | `flex_bgzf_spill` | as #1 but `--soloBucketMode spill --soloBucketMemGB 32 --soloBucketSpillDir /scratch/runs/flex_bgzf_spill/bucket_spill` | +35 min, $1.40 — same-hardware spill-cost ratio |

Disk peaks: FASTQs 465 + refs 52 + CR scratch (large) before step 6; CBQ 200 GB
arrives only after CR scratch is freed. Monitor `df -h /scratch` before steps 5–7.

## Verify, retrieve, tear down

Per run: exit status, `Elapsed`/`Max RSS`/`CPU%` from time output, 16 sample outputs
present. Output-identity checks: `flex_bgzf_ram` vs `flex_cbq_ram` (md5 manifests of
per_sample) and, if run, vs `flex_bgzf_spill` — all should be byte-identical.
`flex_align_ram` differs by design (the residual-alignment recoveries) — record its
counter deltas instead. Ship logs/times/md5s/Log.out/preflight (not matrices) to
`s3://<staging-bucket>/results-<UTC>/`, then **terminate the instance**.

## Print obligations this session feeds

- `[[SYNC-CLOUD-TIME]]` / `[[SYNC-CLOUD-COST]]` in paper + letter quote run #1
  (flex_bgzf_ram) and the whole-session compute cost, rounded.
- The session's command script must be COMMITTED to the repository at release so
  "the commands are in the repository" is true before any referee reads it.
- Run #9 feeds the spill-cost sentence added at resync with the bounded-memory
  Methods text.
- The full matrix populates the 320k table in Supplementary Note 1 (hardware column:
  r6id.8xlarge, us-west-2; input locus: instance NVMe; cyto accommodations noted:
  uv installed, TMPDIR on scratch).

## Executable companion (the Zenodo artifact)

The whole instance session is scripted at
`docs/benchmarks/cloud_320k/run_cloud_320k_matrix.sh` — parameterized by
`STAGING_BUCKET`, `RELEASE_TAG`, `BINARY_ASSET` (plus `RUN_SPILL=1`,
`SKIP_CELLRANGER=1`, `RETRY_FAILED=1`), resumable via step markers, results
accumulated in a TSV, output-identity checks and result shipping included. It is
fail-closed: every arm verifies its own outputs (16 per-sample results), a failed
arm blocks whatever depends on it and is not retried automatically, a failed Cell
Ranger run keeps its pipestance so the logs survive, partial CBQ encodes are
deleted rather than passed downstream, and the script exits nonzero unless every
arm is verified. Run it inside tmux or nohup so an SSH drop cannot kill a run. That script is what gets
committed at release and archived on Zenodo: the paper's "commands are in the
repository" points at it, and the reproducibility demonstration is literally
`bash run_cloud_320k_matrix.sh` on a fresh instance. This document is its
narrative; the script is authoritative where they differ.

## Recording convention

Back at the lab: pull the results prefix, place logs under `docs/benchmarks/`
following the existing naming (`scffpe320k_cloud_<name>_*`), record binary version
string and instance identity with every row.
