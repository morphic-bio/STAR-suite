# Run 10x Multiome v1 with STAR Launchpad

**Version scope:** STAR Suite 1.9.x is the last line that hosts the integrated
RNA + ATAC engine. From **STAR Suite 1.10 onward**, Multiome support is maintained
in the separate **Multiomics Suite** repository, which builds the combined
executable using STAR Suite and Chromap Suite as components. The commands below
are specifically for STAR Suite **1.9.5.b**. They do not apply to the 1.10 STAR
binary. Follow Multiomics Suite's installation instructions for that newer
architecture; this 1.9.5.b guide does not require access to that repository.

STAR Suite 1.9.5.b adds local Multiome execution to Launchpad. It uses the
official recipe bundled with STAR Suite; no separate recipe checkout is needed.
Multiome requires STAR compiled with Chromap support and the native ATAC
peak-matrix helper. The portable tarball and Debian STAR binary omit Chromap.
Installing the standalone `chromap` command does not add ATAC to that binary.

## Build the Multiome engine from source

Use the STAR Suite 1.9.5.b source tree for this Launchpad flow. For the command
line recipe alone, 1.9.5.a also supports the external HTSlib discovery used here.
Both STAR Suite and Chromap Suite must be source checkouts. Run the following in
the same Bash terminal on an Ubuntu 22.04/24.04 x86-64 processing host. The clone
commands assume the destination directories do not already exist and pin the
versions independently of future changes to `master`.

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential xxd cmake pkg-config git ca-certificates \
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev \
  libglib2.0-dev libhts-dev libdeflate-dev python3 python3-venv

mkdir -p "$HOME/software/star-multiome"
cd "$HOME/software/star-multiome"
git clone https://github.com/morphic-bio/STAR-suite.git
git -C STAR-suite checkout --detach 4f444062fe4e168553e4a8e242604197e140a933
git clone --branch v1.1.0 --single-branch \
  https://github.com/morphic-bio/Chromap-suite.git

export STAR_SRC="$PWD/STAR-suite"
export CHROMAP_SRC="$PWD/Chromap-suite"

# Use HTTPS for the RapidMACS submodule, including on machines without SSH keys.
git -C "$CHROMAP_SRC" \
  -c submodule.third_party/rapidmacs.url=https://github.com/morphic-bio/rapidmacs.git \
  submodule update --init --recursive
make -C "$CHROMAP_SRC" -j8

make -C "$STAR_SRC" core-clean
make -C "$STAR_SRC" -j8 core \
  WITH_CHROMAP=1 CHROMAP_SUITE_DIR="$CHROMAP_SRC"
make -C "$STAR_SRC/core/features/libchromap_contract" -j8 \
  star_multiome_atac_peak_mex CHROMAP_DIR="$CHROMAP_SRC"

export STAR_BIN="$STAR_SRC/core/legacy/source/STAR"
export BUILD_ATAC_MEX_NATIVE="$STAR_SRC/core/features/libchromap_contract/star_multiome_atac_peak_mex"
"$STAR_BIN" --version
"$STAR_BIN" --build-features
```

The commands should report `1.9.5.b` and
`{"suite_version":"1.9.5.b","chromap_atac":true}`. `libchromap.a` is linked
into STAR; the ATAC alignment arm runs inside STAR. The separate peak-matrix
helper produces the downstream ATAC matrix. For custom HTSlib installations,
see [compile instructions](compile_instructions.md).

## Start the UI

From the source tree, create the user-owned Python environment once, then start:

```bash
python3 "$STAR_SRC/scripts/launchpad_cli.py" --setup
python3 "$STAR_SRC/scripts/launchpad_cli.py"
```

For a 1.9.5.b installed package, use `star-suite-launchpad --setup` followed by
`star-suite-launchpad`. Keep the `STAR_BIN` and `BUILD_ATAC_MEX_NATIVE` exports
above, or select those executables in the form's **Installed runtime** fields.
The package's portable STAR cannot run Multiome itself.

Open <http://127.0.0.1:8765/launchpad/> and choose **10x Multiome v1 — RNA + ATAC**.
For a remote server, forward the port with
`ssh -L 8765:127.0.0.1:8765 user@server` and open the same local URL.
Run/cancel operations require a loopback connection. Paths always refer to files
on the server running Launchpad.

The launcher works from any working directory. Python dependencies live under
`~/.local/share/star-suite/launchpad-venv` (or `XDG_DATA_HOME`); `--env-dir` overrides
this. The browser assets are bundled and need no CDN. Initial `--setup` needs
access to Python packages. `--port` and `--config` accept site overrides.

## Select inputs

Use the raw demultiplexed sequencing FASTQs, without I1/I2 sample-index files:

| Form field | Sequencing file | Meaning |
| --- | --- | --- |
| RNA R1 | RNA `_R1_` | Cell barcode and UMI |
| RNA R2 | RNA `_R2_` | cDNA |
| ATAC R1 | ATAC `_R1_` | Genomic mate 1 |
| ATAC R2 | ATAC `_R2_` | Cell barcode |
| ATAC R3 | ATAC `_R3_` | Genomic mate 2 |

Enter absolute paths, comma-separated or one per line, with matching lane order
within each library. RNA and ATAC may have different lane counts. The recipe
maps RNA R2/R1 to STAR's read order and extracts/reverse-complements bases 9–24
of the raw ARC v1 ATAC barcode read. CBQ is an alternative input mode.

Supply a STAR genome index, the ARC v1 RNA and ATAC whitelists, reference FASTA,
Chromap index, and ATAC-to-RNA barcode translation file. The two alignment arms
must use the same reference assembly. Choose a new output directory. No dataset
or reference paths are preselected. For mounts outside the installation, home
directory and `/tmp`, set `STAR_SUITE_DATA_ROOTS` to colon-separated trusted
directories before starting the server, or list them in a site config.

This recipe takes explicit paired RNA/ATAC inputs; it does not require an OCM
`--ocmMultiConfig` file. It retains the recipe's RNA `CellRanger4` clipping mode.

## Validate and run

Click **Validate**, **Generate command**, then **Run sample**. Run always repeats
server-side path and runtime checks. A portable or older STAR binary is rejected
before creating the output directory. Launchpad runs one Multiome job at a time,
shows status and logs, and provides **Cancel run**. Browser reload reconnects to
the current server session; restarting the server does not resume jobs.

Success requires the recipe to exit zero and write `LOCAL_MEX_READY.txt`. Outputs
include RNA matrices, ATAC peaks and the peak-by-cell matrix. This UI flow stops
at local matrices; it does not start remote postprocessing or CellBender.

Each run writes `job.json` and `run.log` under the configured `artifact_log_root`.
The UI shows their paths. The record includes arguments, status, exit code,
compiled feature report, and SHA-256 hashes of both executables. Changing the
server's configuration or executables while a run is active is not supported.

## Command-line example: matrices and peaks

This runs the same bundled recipe as Launchpad. Replace the example paths with
your server paths. Each FASTQ variable can also hold a comma-separated lane list;
keep the lanes in matching order within each library. The reference files and
indexes must already exist and use the same assembly.

```bash
# Use the STAR_SRC, STAR_BIN and BUILD_ATAC_MEX_NATIVE exports from the build.
GEX_R1="/path/to/rna/sample_R1.fastq.gz"
GEX_R2="/path/to/rna/sample_R2.fastq.gz"
ATAC_R1="/path/to/atac/sample_R1.fastq.gz"
ATAC_R2="/path/to/atac/sample_R2.fastq.gz"  # Raw barcode read
ATAC_R3="/path/to/atac/sample_R3.fastq.gz"  # Second genomic mate

GENOME_DIR="/path/to/references/star_index"
GEX_WHITELIST="/path/to/references/arc_v1_rna_whitelist.txt"
CHROMAP_REF="/path/to/references/genome.fa"
CHROMAP_INDEX="/path/to/references/genome.chromap.index"
ATAC_WHITELIST="/path/to/references/arc_v1_atac_whitelist.txt"
ATAC_TO_GEX="/path/to/references/atac_to_rna.tsv"
OUT_DIR="$HOME/multiome-results/sample"

bash "$STAR_SRC/share/star-suite/catalogs/official/scripts/run_star_multiome_lane_smoke.sh" \
  --profile matrices-peaks \
  --gex-input-format fastq --atac-input-format fastq \
  --gex-r1 "$GEX_R1" --gex-r2 "$GEX_R2" \
  --atac-r1 "$ATAC_R1" --atac-barcode "$ATAC_R2" --atac-r2 "$ATAC_R3" \
  --genome-dir "$GENOME_DIR" --gex-whitelist "$GEX_WHITELIST" \
  --chromap-ref "$CHROMAP_REF" --chromap-index "$CHROMAP_INDEX" \
  --atac-whitelist "$ATAC_WHITELIST" --atac-to-gex "$ATAC_TO_GEX" \
  --threads 16 --chromap-threads 8 \
  --chromap-low-mem --chromap-macs3-frag-low-mem \
  --skip-build --stop-after-local-mex \
  --out-dir "$OUT_DIR"
```

To inspect the command without processing reads, add `--dry-run` and choose a
separate preview output directory. The recipe writes `RUN_STAR_MULTIOME.sh` and
`DRY_RUN_PREVIEW.txt`. The former contains the exact STAR invocation for those
inputs. No `--ocmMultiConfig` or I1/I2 FASTQs are needed.

## Direct STAR example: RNA + ATAC alignment

The command below shows the STAR stage of the matrices-peaks recipe, using the
input/reference variables above and a separate output directory. It processes
RNA and ATAC together, writing GeneFull RNA outputs, an ATAC BAM and a binary
fragment sidecar. **The direct STAR command does not perform the subsequent
ATAC peak-MEX step.** Use the recipe above for the complete local matrix/peak
workflow, which also runs `star_multiome_atac_peak_mex` and packages the RNA MEX.

```bash
DIRECT_OUT="$HOME/multiome-results/sample-star-only"
mkdir -p "$DIRECT_OUT/run" "$DIRECT_OUT/chromap_tmp"

"$STAR_BIN" \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir "$GENOME_DIR" \
  --readFilesIn "$GEX_R2" "$GEX_R1" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$DIRECT_OUT/run/" \
  --outTmpDir "$DIRECT_OUT/star_tmp" \
  --outSAMtype None \
  --clipAdapterType CellRanger4 --clip3pPolyG yes \
  --alignEndsType Local --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist "$GEX_WHITELIST" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no --soloStrand Forward \
  --soloFeatures GeneFull --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes --soloInlineHashMode no \
  --chromapAtacEnable 1 --chromapAtacStartMode concurrent \
  --chromapAtacReferenceFasta "$CHROMAP_REF" \
  --chromapAtacIndex "$CHROMAP_INDEX" \
  --chromapAtacRead1 "$ATAC_R1" --chromapAtacRead2 "$ATAC_R3" \
  --chromapAtacBarcode "$ATAC_R2" --chromapAtacReadFormat "bc:8:23:-" \
  --chromapAtacBarcodeWhitelist "$ATAC_WHITELIST" \
  --chromapAtacBarcodeTranslate "$ATAC_TO_GEX" \
  --chromapAtacBarcodeTranslateFromFirst 1 \
  --chromapAtacOutputFormat BAM \
  --chromapAtacOutputFragments "$DIRECT_OUT/run/atac_possorted.bam" \
  --chromapAtacSecondaryFragments "$DIRECT_OUT/run/atac_fragments.bin" \
  --chromapAtacSortBam 1 \
  --chromapAtacSummary "$DIRECT_OUT/run/chromap_summary.csv" \
  --chromapAtacThreads 8 \
  --chromapAtacLowMem 1 --chromapAtacLowMemRam 0 \
  --chromapAtacMacs3FragLowMem 1 \
  --chromapAtacTempDir "$DIRECT_OUT/chromap_tmp" \
  --chromapAtacTn5ShiftMode classical
```

`--readFilesIn` is **RNA R2 followed by RNA R1**. ATAC genomic reads are R1/R3;
the raw ATAC R2 barcode read is supplied separately. `bc:8:23:-` performs the
ARC v1 barcode extraction/reverse-complement; do not pre-trim that read for this
example. `--outSAMtype None` suppresses the RNA alignment BAM while retaining
RNA counting; the separately configured ATAC BAM is still written.

## Validation of these instructions

Checked on **Ubuntu 22.04.5 x86-64, 29 September 2026**:

- Fresh HTTPS clones of STAR commit `4f44406` (1.9.5.b) and Chromap Suite
  `v1.1.0` (`a47f077`), including the pinned RapidMACS submodule (`34df448`).
- Chromap, Chromap-enabled STAR and the ATAC peak-matrix helper compiled
  successfully. STAR reported `1.9.5.b` and `chromap_atac:true`.
- Launchpad dependency installation succeeded in a new Python 3.10 environment.
  The source launcher started from a different working directory.
- Chromium drove the actual Launchpad form and launched the freshly compiled
  executables on synthetic data: a 120-kb reference, two barcodes, 4,000 RNA
  read pairs and 4,000 ATAC read pairs, with two threads per arm.
- The run exited zero and wrote `LOCAL_MEX_READY.txt`, RNA raw/filtered matrices
  (one gene, two barcode columns), and an ATAC peak matrix (one peak, two barcode
  columns). Browser reload restored the completed job; no browser errors occurred.
- All Bash examples pass syntax checking. The direct STAR options above match
  the recipe command actually executed, allowing for thread counts, temporary
  directory names and the explicitly written default `--runMode alignReads`.

This verifies compilation and a small end-to-end synthetic run. It does not
establish biological accuracy or performance on real datasets, and Ubuntu
24.04 was not separately tested in this check. Artifact locations are recorded
in [tests/ARTIFACTS.md](../tests/ARTIFACTS.md).
