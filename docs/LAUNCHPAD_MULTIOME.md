# Run 10x Multiome v1 with STAR Launchpad

STAR Suite 1.9.5.b adds local Multiome execution to Launchpad. It uses the
official recipe bundled with STAR Suite; no separate recipe checkout is needed.
Multiome requires STAR compiled with Chromap support and the native ATAC
peak-matrix helper. The portable tarball and Debian STAR binary omit Chromap.
Installing the standalone `chromap` command does not add ATAC to that binary.

## Build the Multiome engine from source

Use the STAR Suite 1.9.5.b source tree for this Launchpad flow. For the command
line recipe alone, 1.9.5.a also supports the external HTSlib discovery used here.
Both STAR Suite and Chromap Suite must be source checkouts. On Ubuntu 22.04/24.04:

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential xxd cmake pkg-config git ca-certificates \
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev \
  libglib2.0-dev libhts-dev libdeflate-dev python3 python3-venv

export STAR_SRC="$(realpath /path/to/STAR-suite)"
export CHROMAP_SRC="$(realpath /path/to/Chromap-suite)"

# Use HTTPS for the RapidMACS submodule, including on machines without SSH keys.
git -C "$CHROMAP_SRC" \
  -c submodule.third_party/rapidmacs.url=https://github.com/morphic-bio/rapidmacs.git \
  submodule update --init --recursive
make -C "$CHROMAP_SRC" -j8 libchromap.a librapidmacs

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

The last command must report `"chromap_atac":true`. `libchromap.a` is linked
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
