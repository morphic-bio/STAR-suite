#!/bin/bash
# run_cloud_320k_matrix.sh — reproduce the STAR Suite 320k scFFPE benchmark matrix
# on a fresh AWS r6id.8xlarge (us-west-2), entirely from public data and the
# released STAR Suite binary. Companion narrative: docs/RUNBOOK_CLOUD_320K_SPILL_20260902.md
#
# Usage (on the instance):
#   export STAGING_BUCKET=s3://<your-staging-bucket>   # refs/tools staged per runbook
#   export RELEASE_TAG=vX.Y.Z BINARY_ASSET=<released asset: the _amd64.deb installs\n#     via a single apt command (Ubuntu AMI); a -glibc*.tar.gz extracts instead>
#   bash run_cloud_320k_matrix.sh            # add RUN_SPILL=1 for the optional spill run
#                                            # add SKIP_CELLRANGER=1 to omit run 5
# Resumable: each completed step drops a marker in /scratch/runs/.done/ and is
# skipped on re-invocation. Results accumulate in /scratch/runs/results.tsv.
set -uo pipefail

: "${STAGING_BUCKET:?set STAGING_BUCKET}"
INSTALL_MODE="${INSTALL_MODE:-asset}"   # asset (GitHub release) | ppa (Launchpad)
if [ "$INSTALL_MODE" = ppa ]; then : "${PPA:?set PPA=ppa:<owner>/star-suite}";
else : "${RELEASE_TAG:?set RELEASE_TAG}"; : "${BINARY_ASSET:?set BINARY_ASSET}"; fi
RUN_SPILL="${RUN_SPILL:-0}"; SKIP_CELLRANGER="${SKIP_CELLRANGER:-0}"
SRC_TAG="${SRC_TAG:-v1.8.1}"; SRC_TARBALL="${SRC_TARBALL:-star-suite_1.8.1.orig.tar.gz}"
S=/scratch; R=$S/runs; D=$R/.done; RESULTS=$R/results.tsv
TENX_TAR="s3://10x.files/samples/cell-exp/9.0.0/320k_scFFPE_16-plex_GEM-X_FLEX_Multiplex/320k_scFFPE_16-plex_GEM-X_FLEX_Multiplex_fastqs.tar"

log() { echo "[$(date -u +%FT%TZ)] $*"; }
mark() { touch "$D/$1"; }
done_p() { [ -e "$D/$1" ]; }
cold() { sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'; }
record() { # name status timefile
  local wall rss cpu
  wall=$(grep -oP 'Elapsed \(wall clock\).*: \K.*' "$3" 2>/dev/null || echo NA)
  rss=$(grep -oP 'Maximum resident set size \(kbytes\): \K\d+' "$3" 2>/dev/null || echo NA)
  cpu=$(grep -oP 'Percent of CPU this job got: \K.*' "$3" 2>/dev/null || echo NA)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$wall" "$rss" "$cpu" "$(date -u +%FT%TZ)" >> "$RESULTS"
  log "RESULT $1: status=$2 wall=$wall rss_kb=$rss cpu=$cpu"
}

# ---------- Step 0: filesystem + inputs ----------
if ! done_p setup; then
  if ! mountpoint -q $S; then sudo mkfs.xfs -f /dev/nvme1n1 && sudo mkdir -p $S && sudo mount /dev/nvme1n1 $S; fi
  sudo chown "$(id -u)" $S; mkdir -p $S/{fastqs,cbq,refs,tools,tmp} $R $D
  command -v aws >/dev/null || { sudo apt-get update -qq && sudo apt-get install -y -qq awscli; }
  printf 'run\tstatus\twall\tmax_rss_kb\tcpu\tutc\n' > "$RESULTS"
  if [ "$INSTALL_MODE" = ppa ]; then
    # The demonstration's opening line: the suite installs from a public apt repository.
    sudo apt-get update -qq && sudo apt-get install -y -qq software-properties-common
    sudo add-apt-repository -y "$PPA" && sudo apt-get install -y star-suite
    ln -sf "$(command -v STAR)" $S/tools/STAR
  else
  curl -L -o $S/tools/release_asset."${BINARY_ASSET##*.}" \
    "https://github.com/morphic-bio/STAR-suite/releases/download/${RELEASE_TAG}/${BINARY_ASSET}"
  case "$BINARY_ASSET" in
    *.deb)    sudo apt-get install -y "$S/tools/release_asset.deb" \
                && ln -sf "$(command -v STAR)" $S/tools/STAR ;;   # one command installs the suite
    *.tar.gz) mkdir -p $S/tools/release && tar -xzf $S/tools/release_asset.gz -C $S/tools/release \
                && cp "$(find $S/tools/release -type f -name STAR | head -1)" $S/tools/STAR && chmod +x $S/tools/STAR ;;
    *) cp "$S/tools/release_asset.${BINARY_ASSET##*.}" $S/tools/STAR && chmod +x $S/tools/STAR ;;
  esac
  fi
  $S/tools/STAR --version | tee $R/star_version.txt
  aws s3 sync "$STAGING_BUCKET/refs"  $S/refs  --only-show-errors
  aws s3 sync "$STAGING_BUCKET/tools" $S/tools --only-show-errors
  chmod +x $S/tools/cyto $S/tools/cbq_ordered_encoder 2>/dev/null || true
  if [ "$SKIP_CELLRANGER" != 1 ]; then
    aws s3 sync "$STAGING_BUCKET/cellranger-9.0.1" $S/cellranger-9.0.1 --only-show-errors
    chmod +x $S/cellranger-9.0.1/cellranger; chmod -R +x $S/cellranger-9.0.1/bin 2>/dev/null || true
  fi
  command -v uvx >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh
  { date -u; curl -s http://169.254.169.254/latest/meta-data/instance-type; uname -a; } > $R/preflight.txt
  mark setup
fi
export PATH="$HOME/.local/bin:$PATH"

if ! done_p fastqs; then
  log "streaming 320k FASTQs from 10x public bucket (index reads dropped)"
  aws s3 cp --no-sign-request "$TENX_TAR" - \
    | tar -x -C $S/fastqs --strip-components=1 --wildcards '*_R1_001.fastq.gz' '*_R2_001.fastq.gz'
  [ "$(ls $S/fastqs/*.fastq.gz | wc -l)" -eq 8 ] || { log "FATAL: expected 8 FASTQs"; exit 1; }
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
  local name=$1; shift
  done_p "$name" && { log "skip $name (done)"; return 0; }
  mkdir -p $R/$name; cold
  /usr/bin/time -v -o $R/$name.time $S/tools/STAR "${common[@]}" "$@" \
    --soloFlexOutputPrefix $R/$name/per_sample --outFileNamePrefix $R/$name/ \
    > $R/$name.stdout 2> $R/$name.stderr
  local st=$?; record "$name" "$st" "$R/$name.time"; [ $st -eq 0 ] && mark "$name"
}

cyto_run() { # name geometry-args... -- inputs...
  local name=$1; shift
  done_p "$name" && { log "skip $name (done)"; return 0; }
  mkdir -p $S/tmp; cold
  TMPDIR=$S/tmp /usr/bin/time -v -o $R/$name.time $S/tools/cyto workflow gex \
    --gex $S/refs/gene_probes_v1.1.0.tsv \
    -w $S/tools/dot-cyto/737K-fixed-rna-profiling.txt.gz \
    -p $S/tools/dot-cyto/probe-barcodes-fixed-rna-profiling-rna.txt \
    --probe-regex 'BC0(0[1-9]|1[0-6])' -T 32 -o $R/$name "$@" \
    > $R/$name.stdout 2> $R/$name.stderr
  local st=$?; record "$name" "$st" "$R/$name.time"; [ $st -eq 0 ] && mark "$name"
}

# ---------- The matrix ----------
star_run flex_bgzf_ram  --soloBucketMode ram --flexNoAlign 1 \
  --readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1 --readFilesIn "$R2" "$R1"
star_run flex_gzip_ram  --soloBucketMode ram --flexNoAlign 1 \
  --readFilesBgzfMode off --readFilesIn "$R2" "$R1"
star_run flex_align_ram --soloBucketMode ram --flexNoAlign 0 \
  --readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1 --readFilesIn "$R2" "$R1"
cyto_run cyto_fastq --preset gex-v1 $(ls $S/fastqs/*.fastq.gz | sort)

if [ "$SKIP_CELLRANGER" != 1 ] && ! done_p cellranger; then
  sed -e "s|reference,.*|reference,$S/refs/refdata-gex-GRCh38-2024-A|" \
      -e "s|probe-set,.*|probe-set,$S/refs/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv|" \
      -e "s|,/mnt/pikachu/tenx_320k_scFFPE/320k_scFFPE_16-plex_GEM-X_FLEX_fastqs|,$S/fastqs|" \
      $S/refs/benchmark_cr9_320k_config.csv > $R/cr9_config.csv
  ( cd $R && cold && /usr/bin/time -v -o $R/cellranger.time \
      $S/cellranger-9.0.1/cellranger multi --id=cr9_320k --csv=$R/cr9_config.csv \
      --localcores=32 --localmem=200 > $R/cellranger.stdout 2> $R/cellranger.stderr )
  st=$?; record cellranger "$st" "$R/cellranger.time"; [ $st -eq 0 ] && mark cellranger
  rm -rf $R/cr9_320k/SC_MULTI_CS   # free pipeline scratch; outs/ retained
fi

if ! done_p cbq_encode; then
  if [ ! -x $S/tools/cbq_ordered_encoder ]; then
    # Untimed producer tooling: the encoder is not part of the installed package,
    # so build it from the released source tarball (single translation unit).
    log "building CBQ encoder from released source"
    sudo apt-get install -y -qq build-essential zlib1g-dev
    curl -L -o $S/tools/src.tar.gz \
      "https://github.com/morphic-bio/STAR-suite/releases/download/${SRC_TAG}/${SRC_TARBALL}"
    mkdir -p $S/tools/src && tar -xzf $S/tools/src.tar.gz -C $S/tools/src --strip-components=1
    ( cd $S/tools/src/core/legacy/source && make cbq-ordered-encoder ) \
      && cp $S/tools/src/core/legacy/source/cbq_ordered_encoder $S/tools/ \
      || { log "FATAL: CBQ encoder build failed"; exit 1; }
  fi
  log "encoding CBQ (untimed producer step)"; pids=(); i=0
  for L in L001 L002 L003 L004; do
    $S/tools/cbq_ordered_encoder \
      --readFilesIn $S/fastqs/16-plex_GEM-X_FLEX_S1_${L}_R2_001.fastq.gz \
                    $S/fastqs/16-plex_GEM-X_FLEX_S1_${L}_R1_001.fastq.gz \
      --outFile $S/cbq/lane_00${i}.cbq > $R/encode_${L}.log 2>&1 & pids+=($!); i=$((i+1))
  done
  ok=0; for p in "${pids[@]}"; do wait "$p" || ok=1; done
  [ $ok -eq 0 ] && [ "$(ls $S/cbq/lane_*.cbq | wc -l)" -eq 4 ] && mark cbq_encode || log "WARN: CBQ encode incomplete"
fi

CBQ=$(ls $S/cbq/lane_*.cbq 2>/dev/null | sort | paste -sd,)
[ -n "$CBQ" ] && star_run flex_cbq_ram --soloBucketMode ram --flexNoAlign 1 \
  --readFilesType Binseq PE --readFilesCbqRangeMode range --readFilesIn "$CBQ"
[ -n "$CBQ" ] && cyto_run cyto_cbq -g '[gex][:18][probe] | [barcode][umi:12]' $(ls $S/cbq/lane_*.cbq | sort)

[ "$RUN_SPILL" = 1 ] && star_run flex_bgzf_spill --soloBucketMode spill --soloBucketMemGB 32 \
  --soloBucketSpillDir $R/flex_bgzf_spill/bucket_spill --flexNoAlign 1 \
  --readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1 --readFilesIn "$R2" "$R1"

# ---------- Output-identity checks + shipping ----------
for pair in "flex_bgzf_ram flex_cbq_ram" "flex_bgzf_ram flex_bgzf_spill"; do
  set -- $pair
  [ -d $R/$1/per_sample ] && [ -d $R/$2/per_sample ] || continue
  (cd $R/$1/per_sample && find . -type f -exec md5sum {} \; | sort -k2) > $R/$1.md5
  (cd $R/$2/per_sample && find . -type f -exec md5sum {} \; | sort -k2) > $R/$2.md5
  diff -q $R/$1.md5 $R/$2.md5 >/dev/null && echo "IDENTICAL: $1 vs $2" | tee -a $R/identity.txt \
    || echo "DIFFER: $1 vs $2" | tee -a $R/identity.txt
done

log "matrix complete"; column -t "$RESULTS"
aws s3 sync $R "$STAGING_BUCKET/results-$(date -u +%Y%m%dT%H%M%SZ)" \
  --exclude "*/per_sample/*" --exclude "*/cr9_320k/*" --exclude "*/bucket_spill/*" --only-show-errors
log "results shipped — REMEMBER: terminate the instance"
