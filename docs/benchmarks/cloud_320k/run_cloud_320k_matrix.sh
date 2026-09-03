#!/bin/bash
# run_cloud_320k_matrix.sh — reproduce the STAR Suite 320k scFFPE benchmark matrix
# on a fresh AWS r6id.8xlarge (us-west-2), entirely from public data and the
# released STAR Suite package. Companion narrative:
# docs/RUNBOOK_CLOUD_320K_SPILL_20260902.md
#
# Usage (on the instance, inside tmux/nohup so an SSH drop cannot kill a run):
#   export STAGING_BUCKET=s3://<your-staging-bucket>   # refs/tools staged per runbook
#   export RELEASE_TAG=v1.8.1 BINARY_ASSET=star-suite_1.8.1-1_amd64.deb
#   bash run_cloud_320k_matrix.sh
# Options:
#   INSTALL_MODE=ppa PPA=ppa:<owner>/star-suite   install from Launchpad instead
#   RUN_SPILL=1        add the optional spill-mode row (RAM-vs-spill cost ratio)
#   SKIP_CELLRANGER=1  omit the Cell Ranger row
#   RETRY_FAILED=1     re-attempt steps that failed in a previous invocation
#   FORCE_MKFS=1       allow formatting the scratch device when it holds a filesystem
#
# Behavior: fail-closed. Every arm verifies its own outputs; a failed arm is
# recorded, blocks anything that depends on it, and makes the whole matrix exit
# nonzero. Completed steps leave markers in /scratch/runs/.done and are skipped on
# re-invocation; failed steps leave markers in /scratch/runs/.failed and are NOT
# retried automatically (set RETRY_FAILED=1 to override).
set -uo pipefail

: "${STAGING_BUCKET:?set STAGING_BUCKET}"
INSTALL_MODE="${INSTALL_MODE:-asset}"   # asset (GitHub release) | ppa (Launchpad)
if [ "$INSTALL_MODE" = ppa ]; then : "${PPA:?set PPA=ppa:<owner>/star-suite}";
else : "${RELEASE_TAG:?set RELEASE_TAG}"; : "${BINARY_ASSET:?set BINARY_ASSET}"; fi
RUN_SPILL="${RUN_SPILL:-0}"; SKIP_CELLRANGER="${SKIP_CELLRANGER:-0}"
RETRY_FAILED="${RETRY_FAILED:-0}"; FORCE_MKFS="${FORCE_MKFS:-0}"
SRC_TAG="${SRC_TAG:-v1.8.1}"; SRC_TARBALL="${SRC_TARBALL:-star-suite_1.8.1.orig.tar.gz}"
EXPECT_SAMPLES="${EXPECT_SAMPLES:-16}"
SCRATCH_DEV="${SCRATCH_DEV:-/dev/nvme1n1}"

S=/scratch; R=$S/runs; D=$R/.done; F=$R/.failed; RESULTS=$R/results.tsv
TENX_TAR="s3://10x.files/samples/cell-exp/9.0.0/320k_scFFPE_16-plex_GEM-X_FLEX_Multiplex/320k_scFFPE_16-plex_GEM-X_FLEX_Multiplex_fastqs.tar"
FAILED=()

log()    { echo "[$(date -u +%FT%TZ)] $*"; }
mark()   { touch "$D/$1"; rm -f "$F/$1"; }
done_p() { [ -e "$D/$1" ]; }
fail()   { FAILED+=("$1"); touch "$F/$1" 2>/dev/null; log "FAIL $1: $2"; }
die()    { log "FATAL: $*"; exit 1; }
cold()   { sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' || log "WARN: cache drop failed"; }

# A step is skipped when already complete; refused when it failed previously
# (unless RETRY_FAILED=1) so a rerun cannot silently build on bad state.
runnable() {
  local name=$1
  done_p "$name" && { log "skip $name (already complete)"; return 1; }
  if [ -e "$F/$name" ] && [ "$RETRY_FAILED" != 1 ]; then
    FAILED+=("$name")
    log "REFUSE $name: failed in a previous invocation; inspect $R/$name.stderr, then rerun with RETRY_FAILED=1"
    return 1
  fi
  return 0
}

require() { [ -e "$1" ] || die "required input is absent: $1"; }

record() { # name status timefile
  local wall rss cpu
  wall=$(grep -oP 'Elapsed \(wall clock\).*: \K.*' "$3" 2>/dev/null || echo NA)
  rss=$(grep -oP 'Maximum resident set size \(kbytes\): \K\d+' "$3" 2>/dev/null || echo NA)
  cpu=$(grep -oP 'Percent of CPU this job got: \K.*' "$3" 2>/dev/null || echo NA)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$wall" "$rss" "$cpu" "$(date -u +%FT%TZ)" >> "$RESULTS"
  log "RESULT $1: status=$2 wall=$wall rss_kb=$rss cpu=$cpu"
  # Durable progress: scratch is ephemeral, so mirror the table after every row.
  aws s3 cp "$RESULTS" "$STAGING_BUCKET/progress/results.tsv" --only-show-errors 2>/dev/null || true
}

# Scratch high-water sampling: cheap df deltas, so a run's disk appetite is
# measured rather than inferred. Sampling costs no I/O against the data.
disk_sample_start() { # name -> pid on stdout
  local name=$1
  ( local base cur
    base=$(df --output=used -B1 "$S" 2>/dev/null | tail -1 | tr -d ' ')
    printf 'utc\tused_bytes\tdelta_bytes\n' > "$R/$name.disk.tsv"
    while :; do
      cur=$(df --output=used -B1 "$S" 2>/dev/null | tail -1 | tr -d ' ')
      printf '%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$cur" "$((cur-base))" >> "$R/$name.disk.tsv"
      sleep 300
    done ) >/dev/null 2>&1 &
  echo $!
}
disk_sample_stop() { # pid name
  kill "$1" 2>/dev/null
  local gib
  gib=$(awk -F'\t' 'NR>1 && $3>m {m=$3} END {printf "%.1f", m/1073741824}' "$R/$2.disk.tsv" 2>/dev/null)
  echo "${gib:-NA}" > "$R/$2.peak_disk_gib"
  log "$2 scratch high-water: ${gib:-NA} GiB above baseline"
}

count_dirs() { find "$1" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | wc -l; }
count_h5ad() { find "$1" -maxdepth 1 -name '*.filt.h5ad' 2>/dev/null | wc -l; }

# ---------- Step 0: scratch filesystem, install, staged inputs ----------
if ! done_p setup; then
  if ! mountpoint -q $S; then
    [ -b "$SCRATCH_DEV" ] || die "$SCRATCH_DEV is not a block device"
    if lsblk -no MOUNTPOINT "$SCRATCH_DEV" 2>/dev/null | grep -qE '[^[:space:]]' ; then
      die "$SCRATCH_DEV is mounted elsewhere; refusing to format"
    fi
    # Test the reported type, not blkid's exit status: an unprivileged blkid can
    # exit 0 with no output, which would falsely read as "holds a filesystem".
    existing_fs=$(sudo blkid -o value -s TYPE "$SCRATCH_DEV" 2>/dev/null)
    if [ -n "$existing_fs" ] && [ "$FORCE_MKFS" != 1 ]; then
      die "$SCRATCH_DEV already holds a $existing_fs filesystem; set FORCE_MKFS=1 only if it is the ephemeral instance store"
    fi
    sudo mkfs.xfs -f "$SCRATCH_DEV" || die "mkfs failed on $SCRATCH_DEV"
    sudo mkdir -p $S && sudo mount "$SCRATCH_DEV" $S || die "mount failed"
  fi
  sudo chown "$(id -u)" $S || die "cannot take ownership of $S"
  mkdir -p $S/{fastqs,cbq,refs,tools,tmp} $R $D $F || die "cannot create scratch tree"
  [ -f "$RESULTS" ] || printf 'run\tstatus\twall\tmax_rss_kb\tcpu\tutc\n' > "$RESULTS"

  command -v aws >/dev/null || { sudo apt-get update -qq && sudo apt-get install -y -qq awscli; }
  command -v aws >/dev/null || die "awscli unavailable"

  # ---- install the released suite ----
  if [ "$INSTALL_MODE" = ppa ]; then
    # The demonstration's opening line: the suite installs from a public apt repository.
    sudo apt-get update -qq && sudo apt-get install -y -qq software-properties-common \
      || die "apt prerequisites failed"
    sudo add-apt-repository -y "$PPA" || die "add-apt-repository failed"
    sudo apt-get install -y star-suite || die "apt install star-suite failed"
    command -v STAR >/dev/null || die "STAR not on PATH after apt install"
    ln -sf "$(command -v STAR)" $S/tools/STAR
  else
    curl -fL -o "$S/tools/release_asset.${BINARY_ASSET##*.}" \
      "https://github.com/morphic-bio/STAR-suite/releases/download/${RELEASE_TAG}/${BINARY_ASSET}" \
      || die "release asset download failed"
    case "$BINARY_ASSET" in
      *.deb)    sudo apt-get install -y "$S/tools/release_asset.deb" || die "apt install of the released deb failed"
                command -v STAR >/dev/null || die "STAR not on PATH after deb install"
                ln -sf "$(command -v STAR)" $S/tools/STAR ;;
      *.tar.gz) mkdir -p $S/tools/release && tar -xzf "$S/tools/release_asset.gz" -C $S/tools/release \
                  || die "release tarball extraction failed"
                bin=$(find $S/tools/release -type f -name STAR | head -1)
                [ -n "$bin" ] || die "no STAR binary inside the release tarball"
                cp "$bin" $S/tools/STAR && chmod +x $S/tools/STAR ;;
      *)        cp "$S/tools/release_asset.${BINARY_ASSET##*.}" $S/tools/STAR && chmod +x $S/tools/STAR ;;
    esac
  fi
  $S/tools/STAR --version > $R/star_version.txt 2>&1 || die "the installed STAR does not run"
  log "installed STAR version: $(cat $R/star_version.txt)"

  # ---- staged references and tools ----
  aws s3 sync "$STAGING_BUCKET/refs"  $S/refs  --only-show-errors || die "refs sync failed"
  aws s3 sync "$STAGING_BUCKET/tools" $S/tools --only-show-errors || die "tools sync failed"
  chmod +x $S/tools/cyto 2>/dev/null || true
  for f in star_index 737K-fixed-rna-profiling.txt sample_whitelist_full_16.tsv \
           probe-barcodes-fixed-rna-profiling-rna.txt h01_cache.bin gene_probes_v1.1.0.tsv; do
    require "$S/refs/$f"
  done
  require "$S/refs/star_index/flex_probe_artifacts/probe_list.txt"
  require "$S/tools/cyto"
  require "$S/tools/dot-cyto/737K-fixed-rna-profiling.txt.gz"
  if [ "$SKIP_CELLRANGER" != 1 ]; then
    aws s3 sync "$STAGING_BUCKET/cellranger-9.0.1" $S/cellranger-9.0.1 --only-show-errors \
      || die "Cell Ranger sync failed"
    chmod +x $S/cellranger-9.0.1/cellranger 2>/dev/null
    chmod -R +x $S/cellranger-9.0.1/bin 2>/dev/null || true
    require "$S/cellranger-9.0.1/cellranger"
    require "$S/refs/refdata-gex-GRCh38-2024-A"
    require "$S/refs/benchmark_cr9_320k_config.csv"
    require "$S/refs/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv"
  fi
  command -v uvx >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh   # cyto's cell-filter step
  { date -u; curl -s http://169.254.169.254/latest/meta-data/instance-type; echo; uname -a; } > $R/preflight.txt
  mark setup   # only after every required artifact is verified present
fi
export PATH="$HOME/.local/bin:$PATH"

# ---------- Step 0b: public inputs ----------
if ! done_p fastqs; then
  log "streaming 320k FASTQs from 10x's public bucket (index reads dropped in flight)"
  aws s3 cp --no-sign-request "$TENX_TAR" - \
    | tar -x -C $S/fastqs --wildcards '*_R1_001.fastq.gz' '*_R2_001.fastq.gz' \
    || die "FASTQ stream/extract failed"
  # The archive stores members as ./<dir>/<file>, so the reads land nested. Flatten
  # whatever depth was used instead of assuming a --strip-components level.
  find $S/fastqs -mindepth 2 -name '*.fastq.gz' -exec mv -t $S/fastqs {} + 2>/dev/null
  find $S/fastqs -mindepth 1 -type d -empty -delete 2>/dev/null
  n=$(find $S/fastqs -maxdepth 1 -name '*.fastq.gz' | wc -l)
  [ "$n" -eq 8 ] || die "expected 8 FASTQs after extraction, found $n"
  mark fastqs
fi

R2=$(ls $S/fastqs/*_R2_001.fastq.gz | sort | paste -sd,)
R1=$(ls $S/fastqs/*_R1_001.fastq.gz | sort | paste -sd,)
common=(--runThreadN 32 --genomeDir $S/refs/star_index
  --soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17 --soloCBlen 16 --soloUMIlen 12
  --soloBarcodeReadLength 0 --soloCBwhitelist $S/refs/737K-fixed-rna-profiling.txt
  --flex yes --soloFlexExpectedCellsPerTag 3000
  --soloSampleWhitelist $S/refs/sample_whitelist_full_16.tsv
  --soloProbeList $S/refs/star_index/flex_probe_artifacts/probe_list.txt
  --soloSampleProbes $S/refs/probe-barcodes-fixed-rna-profiling-rna.txt
  --soloSampleProbeOffset 68 --soloFlexAllowedTags $S/refs/sample_whitelist_full_16.tsv
  --limitIObufferSize 50000000 50000000 --outSJtype None --outSAMtype None --outSAMattributes None
  --soloFeatures Gene --soloCellFilter None --soloMultiMappers Rescue
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR --soloStrand Unstranded --clipAdapterType CellRanger4
  --alignEndsType Local --chimSegmentMin 0 --soloKeysCompat cr --soloSampleSearchNearby no
  --soloHashScreenFile $S/refs/h01_cache.bin --soloBucketCount 256
  --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0
  --dynamicThreadInterface 1 --crAssignConsumerThreads -1 --crAssignSearchThreads 1)

star_run() { # name extra-args...
  local name=$1 st n; shift
  runnable "$name" || return 0
  mkdir -p $R/$name; cold
  local sampler; sampler=$(disk_sample_start "$name")
  /usr/bin/time -v -o $R/$name.time $S/tools/STAR "${common[@]}" "$@" \
    --soloFlexOutputPrefix $R/$name/per_sample --outFileNamePrefix $R/$name/ \
    > $R/$name.stdout 2> $R/$name.stderr
  st=$?
  disk_sample_stop "$sampler" "$name"
  if [ $st -ne 0 ]; then
    record "$name" "exit=$st" "$R/$name.time"
    fail "$name" "STAR exited $st (see $R/$name.stderr)"; return 0
  fi
  n=$(count_dirs "$R/$name/per_sample")
  if [ "$n" -ne "$EXPECT_SAMPLES" ]; then
    record "$name" "bad-output($n/$EXPECT_SAMPLES)" "$R/$name.time"
    fail "$name" "produced $n per-sample outputs, expected $EXPECT_SAMPLES"; return 0
  fi
  record "$name" ok "$R/$name.time"; mark "$name"
}

cyto_run() { # name args...
  local name=$1 st n; shift
  runnable "$name" || return 0
  mkdir -p $S/tmp; cold
  TMPDIR=$S/tmp /usr/bin/time -v -o $R/$name.time $S/tools/cyto workflow gex \
    --gex $S/refs/gene_probes_v1.1.0.tsv \
    -w $S/tools/dot-cyto/737K-fixed-rna-profiling.txt.gz \
    -p $S/tools/dot-cyto/probe-barcodes-fixed-rna-profiling-rna.txt \
    --probe-regex 'BC0(0[1-9]|1[0-6])' -T 32 -o $R/$name "$@" \
    > $R/$name.stdout 2> $R/$name.stderr
  st=$?
  if [ $st -ne 0 ]; then
    record "$name" "exit=$st" "$R/$name.time"
    fail "$name" "cyto exited $st (see $R/$name.stderr)"; return 0
  fi
  n=$(count_h5ad "$R/$name/counts")
  if [ "$n" -ne "$EXPECT_SAMPLES" ]; then
    record "$name" "bad-output($n/$EXPECT_SAMPLES)" "$R/$name.time"
    fail "$name" "produced $n sample matrices, expected $EXPECT_SAMPLES"; return 0
  fi
  record "$name" ok "$R/$name.time"; mark "$name"
}

# ---------- The matrix ----------
star_run flex_bgzf_ram  --soloBucketMode ram --flexNoAlign 1 \
  --readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1 --readFilesIn "$R2" "$R1"
star_run flex_gzip_ram  --soloBucketMode ram --flexNoAlign 1 \
  --readFilesBgzfMode off --readFilesIn "$R2" "$R1"
star_run flex_align_ram --soloBucketMode ram --flexNoAlign 0 \
  --readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1 --readFilesIn "$R2" "$R1"
cyto_run cyto_fastq --preset gex-v1 $(ls $S/fastqs/*.fastq.gz | sort)

if [ "$SKIP_CELLRANGER" != 1 ] && runnable cellranger; then
  sed -e "s|reference,.*|reference,$S/refs/refdata-gex-GRCh38-2024-A|" \
      -e "s|probe-set,.*|probe-set,$S/refs/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv|" \
      -e "s|,/mnt/pikachu/tenx_320k_scFFPE/320k_scFFPE_16-plex_GEM-X_FLEX_fastqs|,$S/fastqs|" \
      -e "s|,/home/lhhung/stage_320k_fastq|,$S/fastqs|" \
      $S/refs/benchmark_cr9_320k_config.csv > $R/cr9_config.csv
  cr_sampler=$(disk_sample_start cellranger)
  ( cd $R && cold && /usr/bin/time -v -o $R/cellranger.time \
      $S/cellranger-9.0.1/cellranger multi --id=cr9_320k --csv=$R/cr9_config.csv \
      --localcores=32 --localmem=200 > $R/cellranger.stdout 2> $R/cellranger.stderr )
  st=$?
  disk_sample_stop "$cr_sampler" cellranger
  du -sB1 $R/cr9_320k 2>/dev/null | cut -f1 > $R/cellranger.pipestance_bytes || true
  if [ $st -ne 0 ]; then
    record cellranger "exit=$st" "$R/cellranger.time"
    # Keep the pipestance intact: its _errors/_log are the only failure evidence.
    [ -f $R/cr9_320k/_errors ] && { log "--- Cell Ranger _errors (tail) ---"; tail -30 $R/cr9_320k/_errors; }
    fail cellranger "exited $st; pipestance retained at $R/cr9_320k for diagnosis"
  else
    record cellranger ok "$R/cellranger.time"; mark cellranger
    mkdir -p $R/cr9_logs
    for a in _log _errors _perf _versions _cmdline _finalstate; do
      [ -f $R/cr9_320k/$a ] && cp $R/cr9_320k/$a $R/cr9_logs/ 2>/dev/null
    done
    ( cd $R/cr9_320k && find outs -maxdepth 3 \
        \( -name 'metrics_summary.csv' -o -name 'web_summary.html' \) \
        -exec cp --parents {} $R/cr9_logs/ \; ) 2>/dev/null
    rm -rf $R/cr9_320k/SC_MULTI_CS   # free pipeline scratch ONLY after success; outs/ retained
  fi
fi

# CBQ arms depend on a complete, verified encode.
if runnable cbq_encode; then
  if [ ! -x $S/tools/cbq_ordered_encoder ]; then
    # Untimed producer tooling: not part of the installed package, so build it
    # from the released source tarball (a single translation unit).
    log "building CBQ encoder from released source"
    sudo apt-get install -y -qq build-essential zlib1g-dev || die "build prerequisites failed"
    curl -fL -o $S/tools/src.tar.gz \
      "https://github.com/morphic-bio/STAR-suite/releases/download/${SRC_TAG}/${SRC_TARBALL}" \
      || die "source tarball download failed"
    mkdir -p $S/tools/src && tar -xzf $S/tools/src.tar.gz -C $S/tools/src --strip-components=1 \
      || die "source extraction failed"
    ( cd $S/tools/src/core/legacy/source && make cbq-ordered-encoder ) \
      && cp $S/tools/src/core/legacy/source/cbq_ordered_encoder $S/tools/ \
      || die "CBQ encoder build failed"
  fi
  log "encoding CBQ (untimed producer step)"; pids=(); i=0
  for L in L001 L002 L003 L004; do
    $S/tools/cbq_ordered_encoder \
      --readFilesIn $S/fastqs/16-plex_GEM-X_FLEX_S1_${L}_R2_001.fastq.gz \
                    $S/fastqs/16-plex_GEM-X_FLEX_S1_${L}_R1_001.fastq.gz \
      --outFile $S/cbq/lane_00${i}.cbq > $R/encode_${L}.log 2>&1 & pids+=($!); i=$((i+1))
  done
  enc_ok=0; for p in "${pids[@]}"; do wait "$p" || enc_ok=1; done
  nfiles=$(ls $S/cbq/lane_*.cbq 2>/dev/null | wc -l)
  nsmall=$(find $S/cbq -name 'lane_*.cbq' -size -1M 2>/dev/null | wc -l)
  if [ $enc_ok -ne 0 ] || [ "$nfiles" -ne 4 ] || [ "$nsmall" -ne 0 ]; then
    fail cbq_encode "encoder status=$enc_ok files=$nfiles undersized=$nsmall; CBQ arms will be skipped"
    rm -f $S/cbq/lane_*.cbq   # never leave partial CBQs where a later arm could read them
  else
    mark cbq_encode
  fi
fi

if done_p cbq_encode; then
  CBQ=$(ls $S/cbq/lane_*.cbq | sort | paste -sd,)
  star_run flex_cbq_ram --soloBucketMode ram --flexNoAlign 1 \
    --readFilesType Binseq PE --readFilesCbqRangeMode range --readFilesIn "$CBQ"
  cyto_run cyto_cbq -g '[gex][:18][probe] | [barcode][umi:12]' $(ls $S/cbq/lane_*.cbq | sort)
else
  log "SKIP flex_cbq_ram and cyto_cbq: no verified CBQ inputs"
fi

if [ "$RUN_SPILL" = 1 ]; then
  star_run flex_bgzf_spill --soloBucketMode spill --soloBucketMemGB 32 \
    --soloBucketSpillDir $R/flex_bgzf_spill/bucket_spill --flexNoAlign 1 \
    --readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1 --readFilesIn "$R2" "$R1"
fi

# ---------- Output identity ----------
: > $R/identity.txt
for pair in "flex_bgzf_ram flex_cbq_ram" "flex_bgzf_ram flex_bgzf_spill"; do
  set -- $pair
  if ! { done_p "$1" && done_p "$2"; }; then
    echo "SKIPPED: $1 vs $2 (a run is absent)" >> $R/identity.txt; continue
  fi
  (cd $R/$1/per_sample && find . -type f -exec md5sum {} \; | sort -k2) > $R/$1.md5
  (cd $R/$2/per_sample && find . -type f -exec md5sum {} \; | sort -k2) > $R/$2.md5
  if diff -q $R/$1.md5 $R/$2.md5 >/dev/null; then
    echo "IDENTICAL: $1 vs $2" | tee -a $R/identity.txt
  else
    echo "DIFFER: $1 vs $2" | tee -a $R/identity.txt
    fail "identity:$1-vs-$2" "per-sample outputs are not byte-identical"
  fi
done

# ---------- Ship and report ----------
log "matrix finished; results table:"; column -t "$RESULTS"
for f in $R/*.peak_disk_gib; do
  [ -f "$f" ] && log "scratch high-water $(basename "$f" .peak_disk_gib): $(cat "$f") GiB"
done
PREFIX="$STAGING_BUCKET/results-$(date -u +%Y%m%dT%H%M%SZ)"
if aws s3 sync $R "$PREFIX" \
     --exclude "*/per_sample/*" --exclude "*/bucket_spill/*" \
     --exclude "cr9_320k/SC_MULTI_CS/*" --exclude "*/counts/*" --only-show-errors; then
  log "results shipped to $PREFIX (Cell Ranger logs included under cr9_logs/)"
else
  log "ERROR: results sync to $PREFIX FAILED — do NOT terminate the instance; retrieve $R by other means"
  FAILED+=("results_sync")
fi

if [ ${#FAILED[@]} -gt 0 ]; then
  uniq_failed=$(printf '%s\n' "${FAILED[@]}" | sort -u | paste -sd' ')
  log "MATRIX INCOMPLETE — failed steps: $uniq_failed"
  log "Investigate before terminating the instance; scratch storage is ephemeral."
  exit 1
fi
log "MATRIX COMPLETE — every arm verified. Remember to terminate the instance."
