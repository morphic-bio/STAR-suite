# Runbook: RNA mate reader threads for STAR Suite 1.11.0 (2026-09-29)

Status: **approved; implementation in progress.** M1 and M2 are committed
and G-R0 passes; M0 and M3-M5 (G-R1, G-R2, G-R3, G-S2, G-M1) are done
with no difference outside the B0 variance classes; stopped before G-S1
(see the handoff). Handoff (current state, supersedes this
file where they differ): `docs/handoffs/HANDOFF_RNA_MATE_THREADS_20260929.md`.

- Worktree `/mnt/pikachu/STAR-suite-rna-mate-threads-20260929`, branch
  `design/rna-mate-threads-20260929` from `origin/master` `9090fb4`
  (v1.10.0 plus the 1.9.5.b Launchpad changes). Under `core/`, `9090fb4`
  differs from the `v1.10.0` tag (`67211a2`) only by the `--build-features`
  flag (`Parameters.cpp:997-1001`), which does not touch mapping.
- Proposed external work root (large outputs, untracked):
  `W=/mnt/pikachu/star_rna_mate_threads_20260929`.
- Release: STAR Suite 1.11.0, a minor release (new option).

## 1. Goal and non-goals

### Finding that shapes the plan

The lane-1 RNA arm is reader-bound, but the cause is not that both mates are
decompressed on one thread. STAR already decompresses each mate in its own
process: one forked zlib helper per mate for plain `.gz` input, or one `zcat`
per mate with `--readFilesCommand zcat`. The serial part is the **parse**. Each
mapping thread in turn holds `g_threadChunks.mutexInRead` and parses both
mates, pair by pair, with iostream calls into its chunk buffer.

Evidence from the full-depth lane-1 run
(`/mnt/pikachu/perf_lane_throughput_20260929/full_baseline/out/lane01/`):

- `Log.out:161-171`: "using internal zlib streaming path" and
  "mate 1/2 using internal gzip FIFO helper", so one helper process per mate.
- `Log.final.out`, PIPELINE DIAGNOSTICS: 17,356 chunks; chunk read time
  2,476.35 thread-seconds, while RNA mapping took 2,498 s (`phases.json`).
  The input lock was held about 99% of the mapping phase.
- The average chunk read took 142.68 ms for about 55,030 pairs (2.59 µs per
  pair, a ceiling near 386k pairs/s; the observed mean was 382k pairs/s).
  The average chunk mapping took 1,943.50 ms, so at most about 13.6 chunks can
  be mapping at once. This matches the 12-16 permits seen in use.
- `zcat` alone takes 19.7 s (R2) and 10.8 s (R1) per 20M pairs. Because the
  mates are already decompressed in parallel, the limit this sets is about
  1.0M pairs/s (R2 on one core), not the 650k pairs/s of the sum. It becomes
  the next ceiling once the parse leaves the lock.

### Goal

Move FASTQ/FASTA parsing for the RNA path off the input lock. One reader
thread per mate parses that mate's stream, from the start, into record batches.
The two readers agree on the batch size and count records from the start of
each file. Inside the lock, the mapping thread only pairs batch *k* of each mate
and writes the chunk text. That text must be byte-identical to what v1.10.0
writes, so mapping, Solo and every output are unchanged. The model is the Flex
reader (`8219555`, `752c339`).

### Non-goals

- No change to decompression. The forked gzip helper, `--readFilesCommand`
  and the BGZF pipe group produce the same bytes into the same FIFOs.
  In-thread zlib is decision D4.
- No change to Flex (the fused pipeline and the Flex BGZF core reader),
  CBQ/Binseq, SAM input, the feature-arm (pf-multi) readers, SLAM or
  TranscriptVB. The new reader stands down for all of them (section 3.5).
- No change to chunk boundaries, mapping, Solo or output files.
- No benchmark claims. Measurements are correctness-first and informal.
- No reader design beyond the constraint below.

### Reader design constraint

- Each mate's file is read on its own thread, from the start of that file.
- The per-mate readers are coordinated only by counting records from the start of each file (ordinal pairing, agreed batch sizes).
- The model is STAR Suite's existing Flex reader, commit `8219555` "Read the two FASTQ mates on separate threads" (public since v1.9.0; see also `752c339`). Its mates agree on the batch size and count records from the file start; pairing truncates to the shorter side.
- No other reader design is in scope. Follow the exclusions in the maintainers' private notes.

The plan below follows this constraint. Each reader consumes its own stream
strictly in order, with no seeking and no positions exchanged between mates.

## 2. Current behaviour, per input case

Paths below are relative to `core/legacy/source/` unless stated.

### Common to every Fastx case

- `Parameters_readFilesInit.cpp:343-371`: the record-at-a-time
  `FastxInputModule` is validated but never activated (`fastxInputActive =
  false`). The comment explains that it serialized STARsolo input under the
  lock. Log line: "Fastx input path: direct STAR chunk reader".
- Read files are opened at the end of parameter parsing
  (`Parameters.cpp:2194`, reached from `STAR.cpp:597`). They are reopened by
  passes that rewind (`STAR.cpp:431-432, 1123-1125, 2313-2319`, two-pass
  `twoPassRunPass1.cpp:102-103`).
- The parse: `ReadAlignChunk_processChunks.cpp`.
  - Every mapping thread loops (942). It locks `mutexInRead` (952), fills its
    chunk and unlocks (1393).
  - The legacy fill loop is at 1105-1338:
    - Loop test (1106): mate-0 and mate-1 bytes below `chunkInSizeBytes`, and
      both streams `good()`. The loop also stops at `readMapNumber` (1108).
    - FASTQ (1195-1241): mate-0 ID token and header extras (1197-1206).
      Each other mate's ID token is read and discarded, keeping only its extras
      (1208-1212). Every mate's header line is built from **mate 0's** ID and
      Illumina filter flag, plus `iReadAll` and `readFilesIndex` (1216-1224).
      Then come the sequence, "+" and quality lines per mate (1227-1241).
    - FASTA and multi-line FASTA (1242-1275); end of stream (1276-1278);
      the `FILE n` lane marker, or the "wrong read ID line format" error
      (1279-1293); the new-lane handling and the "Starting to map file" log
      (1295-1336).
    - Helpers: `fastqHeaderExtraFromCurrentLine` and
      `illuminaFilterFlagFromHeaderExtra` (56-81); `fastqReadOneLine` and
      `removeStringEndControl` (1515-1533).
  - After the loop: the unequal-mate check (1339-1360), the chunk index
    (1371-1379) and the optional chunk trace (1391; `STAR_INPUT_CHUNK_TRACE`,
    149-210).
  - The MAP permit is acquired only **after** the fill (1417-1427). The serial
    parse has never been permit-accounted.
- Mapping threads re-read the chunk text with `readLoad`. For each mate it takes
  `iReadAll`, the filter and `readFilesIndex` from that mate's own header line,
  and the last mate wins (`readLoad.cpp:27`, `ReadAlign_oneRead.cpp:82`). So
  the mate-1 header text matters and must be reproduced exactly. `clipChunk`
  also runs on the chunk text (`ReadAlignChunk_mapChunk.cpp:80-90`).

### Case 1: plain `--readFilesIn R2.gz R1.gz`, no `--readFilesCommand`

This is the Multiomics lane runner's form
(`multiomics-suite/scripts/run_dogmaplex_lane.py:18-25`).

- `Parameters_readFilesInit.cpp`: all names end in `.gz` (237-250), so
  `readFilesUseInternalGzip = true` (254) and the command string is the
  sentinel `INTERNAL_GZIP` (261-264).
- `Parameters_openReadsFiles.cpp`:
  - It sniffs each file for BGZF (431-448). Plain gzip gives no pipe group.
  - For each mate it creates a FIFO (625-636) and forks a child
    (672-701). The child inflates that mate's lanes with zlib and writes
    `FILE i\n` plus text (`streamGzipMateToFifo`, 90-117). The parent opens an
    `ifstream` on the FIFO (694).
- Result: decompression runs in one process per mate. The parse is serial
  under the lock.
- Variants:
  - Uncompressed, single file: opened directly, with no FIFO and no marker
    (603-616).
  - Uncompressed, several lanes: `cat` script (`Parameters_readFilesInit.cpp:
    274-276` → `Parameters_openReadsFiles.cpp:703-761`).
  - A mix of `.gz` and plain files is not internal gzip.
  - `--readFilesLegacyZcat Yes`: a `zcat` script, as in case 2
    (`Parameters_readFilesInit.cpp:265-267`).

### Case 2: `--readFilesCommand zcat`

- The command string is built at `Parameters_readFilesInit.cpp:215, 256-260`.
- For each mate, `Parameters_openReadsFiles.cpp:703-761` writes a script:
  `exec > fifo`, then `echo FILE i` and `zcat "<lane i>"` for each lane. The
  script runs with `vfork`/`execlp`. The result is one `zcat` per mate, with
  lanes in sequence, and the same serial parse.
- Hosted runs differ. Multiomics drops `--readFilesCommand zcat` from the argv
  (`multiomics-suite/scripts/run_star_composed.py:77-80`, applied at
  `run_multiomics.py:61`). Recipes that pass `zcat` (PBMC, CAT-ATAC) therefore
  take case 1 when hosted, and case 2 when STAR runs standalone.

### Case 3: comma-separated multi-lane lists

- One list per mate, with equal counts enforced
  (`Parameters_readFilesInit.cpp:87-107`). Each mate's producer writes a
  `FILE i` marker before lane *i* (`Parameters_openReadsFiles.cpp:98-113`
  internal gzip; 719-727 command; `BgzfPipeGroup.h:77-79`).
- When mate 0 reaches its marker, the parse reads it and skips the whole marker
  line of every other mate without checking it (1329-1332). A chunk can span
  lanes.

### Case 4: BGZF inputs

- **Non-Flex:** a mate file that sniffs as BGZF (mode not `off`, internal
  gzip, no legacy zcat, not TranscriptVB) creates a `BgzfPipeGroup`
  (`Parameters_openReadsFiles.cpp:421-472`).
  - One producer thread per mate decodes that file's BGZF blocks with parallel
    inflate workers. The workers take MAP-domain BGZF decode permits. The
    producer writes the blocks in file order, from the start, with lane
    markers, into the mate's FIFO (660-671; `BgzfPipeGroup.h:67-111`). This
    already exists.
  - The same serial parse follows.
- **Flex:** the fused Flex pipeline closes STAR's streams and reads its own
  lanes, with a mate-reader thread per lane (`mapThreadsSpawn.cpp:296-299`;
  `FlexPipeline.cpp:791-1004`).
  - The Flex BGZF core reader (`Parameters_openReadsFiles.cpp:181-278,
    406-419`; `ReadAlignChunk_processChunks.cpp:960`) applies only to Flex with
    coordinate-sorted BAM.
  - Neither is used, changed or extended here.
- **Feature arms** (pf-multi `crAssign`) read feature FASTQs with their own
  reader (`core/features/process_features/src/pf_bgzf_input.cpp`). That reader
  is not part of the RNA path and is not changed.

### Case 5: STAR hosted by `multiomics`

- `star::host::runMain` runs the same code in-process (`StarHost.cpp`,
  `docs/HOST_API.md`). Multiomics passes the STAR argv and appends
  `--readFilesBgzfMode auto --crAssignBgzfMode auto --bgzfReaderThreads 0`
  (`run_star_composed.py:60-91`).
- Lane 1 is case 1. Evidence:
  - The exact argv is in the full-baseline `lane.stdout`,
    `out/lane01/run_record.json` ("argv") and `out/lane01/Log.out:8`
    (effective command line 157).
  - No `timed_command` record exists under `full_baseline/`, because that run
    was untimed.
- STAR opens the read files, forking the per-mate helpers, during parameter
  parsing (`STAR.cpp:597`). This is before the host's preflight and start hooks
  (`STAR.cpp:1349-1363`).
- The permit pool is `runThreadN` plus the host's extra threads
  (`mapThreadsSpawn.cpp:926-938`), 64 on lane 1. The helper processes and the
  serial parse are outside it.
- pikachu has 32 logical CPUs. While ATAC maps, RNA and Chromap share them.

## 3. Code-level plan

### 3.1 Design in one paragraph

A new module `input/FastxMateReaders.{h,cpp}` holds one reader thread per mate.

- Each thread owns its mate's existing input stream (`P.inOut->readIn[m]`,
  which is the FIFO or the file that STAR already opens). It runs the same
  per-stream calls that the legacy loop makes on that stream, and appends the
  parsed fields (ID token, header extras, sequence and quality lines) to
  fixed-size record batches.
- Under `mutexInRead`, the mapping thread takes records from mate 0's current
  batch and mate 1's current batch in step. It applies the legacy rules for
  `iReadAll`, `readFilesIndex`, lane markers and chunk boundaries, and writes
  the same bytes into `chunkIn`.
- Parsing, including the per-read `istringstream` for the filter flag, runs on
  the reader threads, one per mate and in parallel with mapping. Only
  `memcpy`-level composition stays in the lock.
- The legacy loop is left byte-for-byte as it is, as the `off` path.

### 3.2 Batch contract (how batches stay identical)

- Both mates use the same batch size `B`, a compile-time constant. The
  proposal is 2,048 records, the Flex value (`FlexPipeline.h:272`).
- Each mate reader counts records from the start of each lane file: after that
  lane's `FILE` marker, or from the start of the stream when there is no
  marker. Batch *k* of lane *L* holds records `[kB, kB + n)` of that lane, with
  `n = B` except for the lane's last batch. A batch never spans a lane.
- The filler pairs batches in queue order and checks that each pair carries
  the same lane and the same batch ordinal. Nothing else crosses between the
  readers: no offsets, no positions, no shared cursor. This mirrors Flex
  `processOneLane` (`FlexPipeline.cpp:815-866` reader, `877-929` pairing).
- There are `D` batches in flight per mate (proposed 32, about 64k records,
  roughly one legacy chunk ahead). Memory is about 16 MB for mate 0 and 6 MB
  for mate 1. A reader blocks when its free pool is empty.
- Chunk boundaries do not depend on `B`. The filler walks records one by one
  and stops exactly where the legacy loop stops (3.4).

### 3.3 Files and functions, step by step

1. **Chunk-trace digest (M1, debug only).** In `writeInputChunkTrace`
   (`ReadAlignChunk_processChunks.cpp:167-210`), add two columns when
   `STAR_INPUT_CHUNK_TRACE_DIGEST=1` is set: the zlib `crc32` of each mate's
   chunk bytes. Also report source `mate-threads` in `inputChunkTraceSource`
   (154-165). Without the variable, nothing changes.
2. **New module** `input/FastxMateReaders.h/.cpp`:
   - `MateRecord`: offsets and lengths into the batch arena for the ID
     (mate 0; for FASTA, every mate), extras, sequence bytes and quality bytes
     exactly as the legacy loop would copy them, plus the filter flag
     (mate 0) and the format (FASTQ/FASTA).
   - `MateBatch`: lane, batch ordinal, `laneStart`, `laneEnd` and `inputEnd`
     flags, records, arena, and an in-band error (message plus record ordinal).
   - `FastxMateReader`, one per mate: the thread loop, bounded ready/free
     queues, a stop flag, and counters.
   - `FastxMateReaderGroup`: `ensureStarted()`, `nextPair()` (cursor over the
     current batches, lane transitions, truncation and error propagation),
     `requestStop()`, `join()` and `summary()`.
   - The module carries its own copies of the four legacy helpers (3.4) and a
     small bounded queue. It does not include `FlexPipeline.h`, so Flex files
     stay untouched.
   - Each end-of-input batch records whether the stream had already failed
     before the record start. The filler uses this to write the legacy
     "end of input stream" Log.out line only where the legacy loop wrote it:
     the legacy loop tested mates 1 and 2 for `good()` before each pair, so a
     file without a final newline, or a FASTA file whose last peek reached
     EOF, ended silently.
3. **Reader loop, per mate.**
   - Mate 0: at a record start, `peek()`. `@` means a FASTQ record and `>` a
     FASTA record. A space, a newline or a stream that is not `good()` means
     end of input. Otherwise it reads a word: `FILE` means a lane marker
     (`>> lane`, then ignore the rest of the line); anything else becomes an
     in-band error carrying the legacy message and read number.
   - Mates 1 and 2: skip whitespace; EOF means end; a `FILE` token means a
     marker; anything else is a record. The legacy loop never validates these
     IDs.
   - The FASTQ and FASTA field parsing repeats the legacy calls in the same
     order on the same stream (3.4).
   - A batch is closed and pushed when it is full, at a lane end or at the
     input end. The stop flag is checked between records, and a queue close
     wakes a blocked push.
4. **Filler.** In `ReadAlignChunk::processChunks`, add a new branch before the
   legacy `else` at `ReadAlignChunk_processChunks.cpp:1105`:
   `else if (P.fastxMateReaders) { … }`.
   - It calls `ensureStarted()`.
   - It keeps the legacy loop test (mate-0 and mate-1 bytes below
     `chunkInSizeBytes`) and the `readMapNumber` stop.
   - It calls `nextPair()`:
     - end of input: log the legacy "end of input stream" line and break;
     - lane start: set `P.readFilesIndex` and write the legacy "Starting to map
       file # …" lines (1327-1334);
     - error: `exitWithError` with the legacy text and `EXIT_CODE_INPUT_FILES`;
     - a pair: `++P.iReadAll`, then write both mates' text directly into
       `chunkIn` (3.4).
   - The existing tail (1339-1393) runs unchanged. The cursor state lives in
     the group and is touched only under `mutexInRead`, as the CBQ path does
     with `cbqInputPendingBatch` (962-1045). Fully consumed batches go back to
     the free queues.
5. **Parameters.**
   - `Parameters.h` (input block near 157-189): `string readFilesMateThreads
     = "auto"`, `std::shared_ptr<star::input::FastxMateReaderGroup>
     fastxMateReaders` and `string fastxMateReadersReason`.
   - `Parameters.cpp`: register the option next to `readFilesBgzfMode`
     (255-256) and validate `auto|off|on` next to 1737-1752.
   - `parametersDefault`: text after `readFilesBgzfMode` (444-449). Regenerate
     `parametersDefault.xxd`; `tests/test_parameters_default_generation.py`
     must pass.
   - `shared_ptr` matters because two-pass copies `Parameters`
     (`twoPassRunPass1.cpp:17`); `P` and `P1` then share one group, as they
     already share `bgzfPipes`.
6. **Open.** At the end of the FIFO/direct branch of
   `Parameters::openReadsFiles` (after 762), evaluate the gate (3.5) and log
   either "Fastx mate readers: active (N mates, batch B, depth D, threads
   outside --runThreadN)" or "not active (<reason>)". If the gate accepts,
   create the group bound to `inOut->readIn[0..readNends-1]`, without starting
   it. Threads start lazily on the first fill, so a Flex run that takes over at
   `mapThreadsSpawn.cpp:1008` never starts them.
7. **Close.** At the top of `Parameters::closeReadsFiles`
   (`Parameters_closeReadsFiles.cpp:13`), when a reader group exists:
   1. end the helper children first (the existing kill loop, moved into a
      helper with unchanged behaviour), so a reader blocked in `read()` gets
      EOF;
   2. request stop, close the queues and join the reader threads;
   3. log `summary()` to `Log.out` if the readers started;
   4. close the streams;
   5. join the BGZF producers;
   6. run the helper-kill loop again at the end, as before (a no-op when step
      1 ran).

   Readers are also stopped and joined at the end of `mapThreadsSpawn`, after
   the mapping threads join and before MAP is marked complete. After a
   `--readMapNumber` stop or two-pass pass 1, readers may have read ahead;
   this way no reader holds or waits for a permit at the exit invariant. A
   reader that has been asked to stop does not wait for a new permit.
   `openReadsFiles` stops any previous group before it reopens. Reader
   threads therefore never exist across a later `fork()` in
   `openReadsFiles`.
8. **Build and test wiring.**
   - `Makefile`: add `FastxMateReaders.o` to `OBJECTS` (input objects near
     line 228), a rule like line 695, a `fastx-mate-reader-harness` target and
     `clean` entries (739). The host library picks the object up through
     `STAR_HOST_OBJECTS` (`OBJECTS` minus `STARmain.o`); verified with
     `ar t libstar_suite.a`.
   - New `input/fastx_mate_reader_harness.cpp` and
     `tests/run_fastx_mate_reader_harness_smoke.sh`. Add a manifest row
     `input-fastx  fastx-mate-readers  contract` to
     `tests/production_module_regression_manifest.tsv`.
9. **Docs (M6).** `parametersDefault` text, `docs/RELEASE_NOTES_v1.11.0.md`
   (draft; a local 1.10.1 notes file already exists for unreleased Launchpad
   content), and the README option list if it lists input options. Keep the
   runbook and handoff current.

### 3.4 Chunk text that must be reproduced (the byte rules)

**FASTQ pair** (mate-0 record starts with `@`):

- `iReadAll += 1`.
- The ID is mate 0's first token, read with `operator>>` semantics.
  - If its last byte, compared as signed `char`, is below 33, that byte is
    dropped (`removeStringEndControl`).
  - With `--outSAMreadID Number` the ID is `@<iReadAll>`.
- `extra_m` is the rest of mate *m*'s header line. Leading spaces and tabs are
  removed, and trailing bytes below 33 (compared as unsigned) are removed.
- The filter flag is `Y` when the first whitespace-delimited token of
  `extra_0` has length ≥ 4 with `[1]==':'`, `[2]=='Y'` and `[3]==':'`;
  otherwise it is `N`.
- Every mate *m* gets the header line
  `ID SP iReadAll SP filter SP readFilesIndex [SP extra_m] \n`.
- The sequence and quality lines follow `fastqReadOneLine` semantics:
  - `getline` with limit `DEF_readNameSeqLengthMax+1`;
  - one trailing byte below 33 (signed) is dropped;
  - a `\n` is appended.
- The `+` line is skipped with `ignore(DEF_readNameSeqLengthMax, '\n')` and
  written as `+\n`.

**FASTA pair** (`>`):

- Each mate writes its own first token, or `>` plus `iReadAll` when
  `outSAMreadID == "Number"`. The rest of the header is ignored.
- Then ` iReadAll N readFilesIndex \n` (note the trailing space).
- Sequence lines are joined until a line that starts with `@`, `>`, a space or
  a newline, or EOF. Each line loses its terminator and one trailing byte below
  33 (signed). A `\n` ends the record.

**Boundaries:**

- Before each pair the filler stops if mate-0 or mate-1 bytes are at least
  `chunkInSizeBytes`, or if `iReadAll == readMapNumber`. Mate 2 of a
  three-file input is not part of this test, as in legacy.
- End of input is a mate-0 record start at a stream that is not good, or a
  line that starts with a space and is not blank. Blank lines where a read
  header is expected follow the rule in 3.7 in every mate (a deliberate
  change: v1.10.0 ended the whole input at a blank line in mate 1 and skipped
  blank lines silently in mates 2 and 3).
- The chunk end (`\n` per mate) and `iChunkIn` are unchanged.

**Proof of identity:** the harness oracle is a verbatim copy of lines
1106-1337 in the harness file, run on the same inputs. The harness also pins
the new module's copies of the four helpers. On real data, the M1 trace digests
prove identity chunk by chunk (4.4).

### 3.5 When the new reader is active (gate)

Active only when all of these hold:

- `readFilesTypeN == 1` (Fastx);
- no Flex (`pSolo.flexMode` or `flexModeStr`), which keeps Flex unchanged;
- no Flex BGZF core reader (`bgzfCoreActive`);
- no SLAM (per-file skip and stop and auto-trim passes read the streams directly,
  `ReadAlignChunk_processChunks.cpp:851-937`);
- no TranscriptVB, which is order-sensitive online learning. This follows the
  existing BGZF precedent at `Parameters_openReadsFiles.cpp:422-430`;
- no batch mode;
- no spatial raw-R1 tap;
- 1-3 mates.

Some passes set their SLAM or TranscriptVB flags after the files were opened.
Before the readers start (on the first chunk fill, under `mutexInRead`),
`fastxMateReadersUsable` checks these flags again. If one is set, it logs
"not active (<reason>)", stops the unstarted group and falls back to the
legacy loop. The streams are untouched until the readers start, so the
fallback is exact. With `on`, this late rejection is fatal as well.

The reader does work with `--readFilesCommand`, `--readFilesManifest`,
comma-separated lanes, BGZF via the pipe group, two-pass (reopen),
`--outFilterType BySJout` (stage 1 reads input normally; stage 2 reads
per-thread files), `--readMapNumber`, `--outSAMreadID Number`, Y/noY FASTQ
emission and any `--runThreadN`.

Modes:

- `auto`: gate as above; a rejection logs the reason and uses the legacy loop.
- `on`: a rejection is fatal and names the reason.
- `off`: the legacy loop, exactly v1.10.0.

### 3.6 Threads and the permit allocator

- Reader threads take **MAP-domain decode permits** (`PermitWork::BGZF`)
  through the same hooks as the BGZF inflate workers
  (`mapDecodePermitHooks` in `Parameters_openReadsFiles.cpp`), and report
  their queue state with `mapPermitObserveDecode`. Releases go to decode
  accounting, not to completed pairs. Decision D3 (as changed by the author)
  covers this.
- A reader holds a permit **only while it parses bytes already in memory**.
  Its stream buffer gives the permit back before any blocking read and takes
  it again afterwards. It never holds one while waiting for queue space or
  input. This avoids a deadlock with a small pool, where the BGZF inflaters
  need permits to feed the same FIFO. Mapping threads hold their permit only
  inside `mapChunk`, never while waiting for `mutexInRead`.
- The pool size, the domain floors, the External domain and
  `docs/HOST_API.md` are unchanged. There is no allocator or host API
  change.
- A run has `runThreadN` mapping threads, plus one reader thread per mate,
  plus the existing decompressor processes. `Log.out` states this.
- **Expected effect:** MAP-domain use rises from about 13.6 toward
  `runThreadN`. In the hosted lane, RNA then competes harder with Chromap for
  the 32 CPUs while ATAC maps, and the allocator's floors govern that split.
  M4 reports ATAC mapping seconds next to RNA, informally.
- **Summary** (in `Log.out` only, never `Log.final.out`, which other tools
  parse), per mate:
  - records and batches;
  - time blocked reading the stream (decompressor-bound);
  - time blocked on a full queue (back-pressure);
  - time the filler waited for a batch.

  These show whether the ceiling has moved to the decompressor, the reader or
  mapping.

### 3.7 Errors and early EOF

| Situation | v1.10.0 | With mate threads |
|---|---|---|
| Well-formed input | chunks as described | identical chunks |
| Mate-0 record start neither `@`, `>` nor `FILE` | fatal "wrong read ID line format", read *N* | same text and *N*; raised by the filler when it reaches that record |
| Single blank line where a read header is expected, followed by a header in the file's format | mate 1: ends all input; mates 2-3: skipped silently (FASTA with `--outSAMreadID Number`: read as the header, mates misalign) | **deliberate change, every mate:** skipped; WARNING at the first one per file; count per file in `Log.out` |
| Blank line where a read header is expected, followed by anything else | mate 1: ends all input; mates 2-3: skipped, the next line read as a header | **deliberate change, every mate:** fatal, naming the file and line |
| Two or more blank lines in a row before the end of a file | mate 1: ends all input; mates 2-3: skipped silently | **deliberate change, every mate:** fatal at the second blank line |
| Blank lines at the end of a file (then end of input or the next file) | mate 1: ends all input; mates 2-3: skipped silently | **deliberate change, every mate:** skipped, any number; WARNING and count |
| Line starting with a space (not blank) at a mate-1 record start | ends all input | same |
| Mate 0 ends first | silent truncation | truncation plus a WARNING with both counts (D2) |
| Mate 1 ends first | the loop checks mate 1 only before each pair, so a malformed final pair is appended | clean truncation plus a WARNING (D2) |
| Per-lane counts differ | lanes silently misalign | truncate that lane to the shorter mate, WARNING, realign at the next marker (D2) |
| Different lane-marker structure | undefined | fatal, naming lane and mate |
| Mates in different formats (FASTA vs FASTQ) | mate ≥1 parsed in mate 0's format | fatal |
| Decompressor helper fails mid-file | stream ends, treated as end of input (exit status not checked, `Parameters_closeReadsFiles.cpp:51-76`) | unchanged; pre-existing, see D7 |
| `readMapNumber` reached | filler stops; helpers killed at close | same; readers stopped and joined at close |
| Reader-thread exception (allocation) | not applicable | in-band error, fatal `EXIT_CODE_INPUT_FILES` |

The unequal-bytes check (1339-1360) cannot fire, because both mates always
contribute the same record count per chunk.

**Blank lines where a read header is expected (the author's decisions of 29
and 30 Sep; a deliberate change from v1.10.0 on malformed input).** A blank
line must never be taken as a read header, because that is the safe choice
(29 Sep). The same rule applies to every mate (1, 2 and 3), for FASTQ and
FASTA:

- A blank line is a line of only spaces, tabs or carriage returns. STAR never
  reads one as a read header.
- **A single blank line followed immediately by a read header** in that
  file's format (`@` for FASTQ, `>` for FASTA; either if the file has no read
  yet) is skipped, with a WARNING and a count (29 Sep).
- **A blank line followed by anything else** (another kind of line, a header
  in the other format, or a line that starts with whitespace) is fatal: STAR
  stops with an error that names the file and the line and says the input is
  malformed (29 Sep).
- **A second blank line in a row, anywhere before the end of the file**, is
  fatal, reported at the second blank line, whatever follows it (30 Sep).
- **Blank lines at the end of a file**, followed by end of input or by the
  next input file, are skipped with a WARNING and counted, whatever their
  number (30 Sep).
- WARNING and counts: at the first skipped blank line in each input file,
  `Log.out` (and stderr) gets a WARNING that names the file, the mate and the
  line number. Later blank lines in the file are counted, not reported. When
  the readers are joined, the reader summary in `Log.out` gives each file's
  count.
- Line numbers count from 1 at the start of each input file.
- Where v1.10.0 already consumed such a line in another way, nothing changes.
  This covers a blank line inside a FASTQ record, and a CRLF blank line after
  FASTA sequence lines, which v1.10.0 and the readers both read as an empty
  sequence line.
- Harness: `fastx_mate_reader_harness` part 3 covers mates 1-3 in FASTQ and
  FASTA:
  - one blank line then a header: skipped, one WARNING with file and line,
    output equal to the input without it;
  - a blank line then garbage, an indented header or the other format: fatal,
    naming file and line;
  - two blank lines in a row (also with spaces and tabs) before a header or
    garbage: fatal at the second blank line;
  - several single blank lines plus three at the end of the file: one
    WARNING, exact count, output equal to the input without them;
  - two input files with blank lines at the end of each (before the next
    file's marker, and before end of input): a WARNING and a count per file.

  The identity part (inputs without blank lines) still matches the copy of
  the v1.10.0 loop, including its `Log.out` lines, with no new WARNING.
- Pre-existing and unchanged (left documented): with FASTA, a `FILE` marker
  line right after a read's sequence lines is read as sequence (the FASTA
  loop stops only at `@`, `>`, a space or a newline), in v1.10.0 and here
  alike. The harness's two-file FASTA cases therefore end file 0 with one
  blank line in every mate.

### 3.8 Multi-lane lists

Each mate's producer still writes one `FILE i` marker per lane. Each reader
records lane starts, including empty lanes, as zero-record `laneStart`
batches, so the filler writes the same "Starting to map file # i" lines at the
same point as legacy. A lane end closes the current batch. Chunks keep spanning
lanes, as before.

### 3.9 `--readFilesCommand` and other producers

- Nothing changes on the producer side. The reader consumes whatever writes the
  mate's FIFO: the internal gzip helper, a user command (`zcat`, `cat`, a
  parallel decompressor per mate as measured in `752c339`), the `cat` script
  for plain multi-lane input, or the BGZF pipe group.
- A single plain file is read directly, as before.
- `--readFilesLegacyZcat` behaves as before.

### 3.10 Option name and default

`--readFilesMateThreads auto|off|on`, lower-case like `--readFilesBgzfMode`.
The default is `auto` (D1).

Proposed `parametersDefault` text:

```
readFilesMateThreads        auto
    string: auto|off|on - parse each FASTX mate on its own reader thread
                            auto ... use for supported Fastx runs (not Flex, SLAM, TranscriptVB, batch mode)
                            off  ... parse all mates in the mapping threads, one chunk at a time (v1.10.0)
                            on   ... require mate reader threads; unsupported runs are fatal
```

Multiomics needs no change; the option passes through its argv adapter.
Optionally, Multiomics could later own or pin it as it does the BGZF modes
(Multiomics follow-up, not part of this change).

## 4. Output-identity gates against v1.10.0

### 4.1 Builds

- **B0:** `9090fb4` (v1.10.0 code). Built from `0f9701a`, which has no `core/`
  changes after `9090fb4`.
- **B1:** B0 plus the M1 trace digest only.
- **C:** the candidate.
- All three are built the same way from fresh clones at a detached commit in
  `$W/src_<label>`, with no branches or worktrees. For M3 onward, use the gate
  container (`star-suite-gate-v1100:20260928`,
  `/mnt/pikachu/star_suite_v1100_gates_20260928/tools/build_in_container.sh`).
  Development builds run on the host: `nice -n 10`, at most `-j16`.
- B1 outputs must equal B0 outputs on F1 and F2 before B1 traces are used as
  references.

### 4.2 Comparison rules

- **Byte identity:** every file except logs. `.gz` files are compared after
  decompression. BAM is compared as `samtools view -h --no-PG` records without
  `@PG`/`@CO`: in order for coordinate-sorted BAM, as a multiset for unsorted
  BAM (the SLAM harness `ORDER_ONLY` rule). Run and tree roots are replaced by
  placeholders.
- **Tools:** the G-S1 comparator
  (`/mnt/pikachu/star_suite_v1100_gates_20260928/tools/compare_outputs.py`)
  and, for feature arms, the Multiomics variance lists.
- **Always byte-compared:**
  - Solo outputs: `Solo.out/{Gene,GeneFull}/{raw,filtered}/*`, `Summary.csv`,
    `Features.stats`, `Barcodes.stats` and every other `Solo.out` file;
  - `SJ.out.tab`, `ReadsPerGene.out.tab`, `outs/` when present;
  - `Log.final.out` count lines, excluding the time, speed, "ms" and PIPELINE
    DIAGNOSTICS lines.
- **Repeatability control:** run B0 twice per fixture at the fixture's thread
  count first. Only files that differ between the two B0 runs may differ for C,
  and only in the same way (order-only or the named content-varying fields).
  The expected set is the `cr_assign/` feature files already listed in
  `multiomics-suite/tests/reference/{order_only,content_varying}_*.txt`;
  `Solo.out` is expected to be exact.
- **Single-thread runs:** at `--runThreadN 1`, everything must be byte
  identical.
- **Input-layer identity:** C's `STAR_INPUT_CHUNK_TRACE` with digests equals
  B1's row for row. Compare `chunk_index`, read range and count, mate bytes,
  both CRCs, `read_files_index` and `no_reads_left`; ignore `thread` and
  `source`.
- **Log.out assertions:** "Fastx mate readers: active" where expected, and
  "not active (<reason>)" for Flex, SLAM, TranscriptVB and CBQ.

### 4.3 Fixtures (all on pikachu)

STAR-only argv means the source argv minus the host options (`--chromapAtac*`,
`--multiomeAtac*`); STAR 1.10+ rejects them without a host.

| ID | Set | Inputs | Cases covered | Argv source | Threads |
|---|---|---|---|---|---|
| F1 | PBMC 3k 100k GEX, 2 lanes | `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/gex/pbmc_unsorted_3k_S01_L00{3,4}_R{2,1}_001.fastq.gz`; genome `/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star`; whitelist per the smoke | 2+3 (`zcat`); 1+3 (no command); 4 (`bgzip -c` copies in `$W/fixtures`); plain multi-lane (`cat`); plain single lane (direct); `--readFilesManifest` once | `multiomics-suite/tests/run_star_chromap_macs3_lowmem_smoke_100k.sh:62-75` (sorted BAM with CB/UB) | 8 and 1 |
| F2 | DOGMA-plex lane 1, 2M | `/mnt/pikachu/single_binary_reference_20260928/dogmaplex_lane1_2m_sidecar/run1/input_fastqs/{gex,adt,hto,ct}/`, staged as `$W/fixtures/dp2m/lane_01/<lib>/` symlinks | 1, plus pf-multi feature arms | `run_dogmaplex_lane.py --dry-run --fastq-root $W/fixtures/dp2m …` renders the argv and config; drop host options | 16 |
| F3 | DOGMA-plex lane 1, 20M | `/mnt/pikachu/multiomics-suite-single-binary-20260928/artifacts/single_binary_20260929/g_m2/inputs/lane_01/{gex,adt,hto,ct}/` | 1 at scale | same, `--threads 32` | 32 |
| F4 | HIV DOGMA 100k | `/mnt/pikachu/hiv_dogma_gse239916/star_four_arm_downsample_inputs_100k/fastq/gex_R{2,1}.fastq.gz` (plus ADT) | 1, `--readMapNumber 100000` | `…/single_binary_reference_20260928/hiv_dogma_four_arm_100k/run1/star_run_20260928_165218/Log.out:8` | 8 |
| F5 | CAT-ATAC 100k | `/mnt/pikachu/catatac_gse288996/fastq/GEX/SRR32265752_{2,1}.fastq.gz` (full files) | 2 and 1, with `--readMapNumber 100000` on full files (early stop while readers run ahead) | `…/single_binary_reference_20260928/catatac_trimodal_100k/run1/RUN_STAR_TRIMODAL_SMOKE.sh` | 8 |
| F6 | Bulk PE, SRR4422207, 500k pairs | `/tmp/starsuite-public-fixture.ZuSe1s/SRR4422207_{1,2}.fastq.gz` (**in /tmp: copy to `$W/fixtures/bulk` first and record sha256**); genome `/storage/autoindex_110_44/bulk_index` | 1; sorted BAM plus `GeneCounts`; unsorted BAM; `--twoPassMode Basic`; `--outFilterType BySJout`; single-end (mate 1 only); plain FASTQ | new argv, no TranscriptVB | 8 and 1 |
| F7 | scRNA 100k regression, 2 lanes, `zcat` | `tests/run_scrna_gex_100k_regression.py` defaults | 2+3 | the script (a G-S1 row) | default |
| F8 | Flex | `tests/run_flex_half_probe_100k_smoke.sh`, `tests/run_flex_tiny_public_smoke.sh` | reader stands down; outputs unchanged | the scripts | default |
| F9 | Synthetic | harness inputs; `tests/test_scrna_gex_counts.py` fixture (G-S2) | edge cases; host | — | — |

### 4.4 Gates

- **G-R0, harness (M2).** The oracle and the new module must produce equal
  chunk bytes, boundaries, `iReadAll` and `readFilesIndex` sequences over:
  - SE, PE and three-mate inputs; FASTQ, FASTA and multi-line FASTA;
  - CRLF; tabs and empty extras; filter flags `Y` and `N`;
    `--outSAMreadID Number`;
  - lanes, including an empty lane; blank lines in mate 0 and in mate 1;
    no final newline; lines at the `DEF_*` limits; bytes ≥ 0x80;
  - chunk sizes of one record, odd sizes and the default; `readMapNumber`
    limits.

  Mismatched-mate and structural-error inputs check the new, intended
  behaviour (3.7) separately.
- **G-R1, input layer (M3, M4).** Trace digests: C equals B1 on F1 (all input
  variants), F2, F3, F4, F5, F6 (all variants) and F7.
- **G-R2, STAR standalone outputs (M3, M4).** C equals B0 under 4.2 on
  F1-F7, at the fixture thread counts and at 1 thread for F1 and F6.
- **G-R3, Flex unchanged (M3).** F8 outputs are byte-identical to B0, and
  `Log.out` shows the reader not active (Flex).
- **G-R4, host (M5).**
  - G-S2 (`tests/host_api/run_host_api_tests.sh`) on C: a host with no
    callbacks, and a host with a dummy External domain. The `.gz` fixture makes
    the reader active.
  - Multiomics G-M1, five fixtures: a Multiomics development build against C,
    compared with `/mnt/pikachu/single_binary_reference_20260928` using
    `scripts/compare_composition_outputs.py` and the committed variance lists,
    as `run_g_m1.py` does (D6).
- **G-R5, release regression (M6).** G-S1 (24 non-multiome manifest rows
  plus Tier A) on B0 vs C, with the new default in force
  (`/mnt/pikachu/star_suite_v1100_gates_20260928/gate/run_gs1.sh`,
  `tools/compare_gs1.sh`). Gated-off modes (SLAM, TranscriptVB, CBQ, Flex rows)
  must be unchanged.

### 4.5 Informal measurements (no claims)

- **M0 diagnostic** on F2 with B0: `.gz` input vs the same content as plain
  FASTQ. Compare "Avg chunk read time (ms)". If they are similar, the parse
  dominates, as the lane-1 numbers suggest.
- **M4** on F3 with B0 vs C: alternate the order, two runs each, under the
  lock. Record the PIPELINE DIAGNOSTICS lines, the reader summary, MAP permits
  in use (telemetry) and wall time. Label the results informal.
- A full-depth lane run only on the author's request, ideally after Chromap
  1.2.0 lands.

## 5. Milestones

Every milestone ends with **stop and report**: results, paths and the next
proposed step. Commit only after the author approves.

| M | Work | Stop-and-report contents |
|---|---|---|
| M0 | Stage `$W` (fixtures, sha256 of inputs, STAR-only argv files). Build B0. Run the B0 repeatability pairs for F1, F2 and F4-F7 and derive the variance classes. Run the M0 diagnostic. | variance classes; diagnostic numbers; any fixture problem |
| M1 | Trace digest (3.3 step 1). Build B1; check B1 = B0 on F1 and F2; collect B1 traces for all fixtures. | B1 = B0 evidence; trace set; diff for review |
| M2 | Module, filler, option, gate, open and close changes, harness (steps 2-8); G-R0. | diff for review; harness report |
| M3 | Build C. G-R1, G-R2 and G-R3 on the small sets (F1 in all variants, F2, F4-F8). | per-fixture verdicts |
| M4 | G-R1 and G-R2 on F3 (20M); the informal M4 measurements; propose whether M4b is needed. | 20M identity; reader summary; ceiling analysis |
| M4b | Optional (D5): buffered parser in the reader threads replacing the iostream calls, proven by an extended G-R0 plus a rerun of G-R1 and G-R2. | as M3 and M4 |
| M5 | G-R4: G-S2 on C; Multiomics development build against C, then G-M1 (after D6). | host verdicts |
| M6 | Default (D1), docs, G-R5 (G-S1 plus Tier A); final handoff. **No tag, push or release.** | release-readiness note |

**Stop conditions** (stop, report, wait):

- Any output or trace difference outside the B0 variance classes. Any output
  difference means stop and report; do not fix and continue silently.
- A harness mismatch that is not an intended 3.7 behaviour.
- M0 shows that plain-FASTQ input is as slow as `.gz` input *and* that the
  chunk read is not dominated by the parse (the design premise fails).
- Any need to change Flex files, the Flex BGZF core reader, CBQ, the feature
  readers or the permit allocator.
- Any idea that goes beyond the reader design constraint (section 1). Do not
  prototype it; stop and ask.
- Any need to read vendor code or excluded material.
- The host is not released, or another job holds the lock: wait, and do not
  start heavy runs.
- Usage near 90%: stop, update the runbook and handoff, clean up `$W` scratch.

## 6. Open author decisions

| # | Decision | Recommendation |
|---|---|---|
| D1 | Option name and default: `--readFilesMateThreads auto\|off\|on`, default `auto` or `off` | **Decided: `auto`.** Multiomics gains without argv changes; the identity gates prove equivalence; `off` is the exact v1.10.0 path. Use `auto` in the branch from M2 so every gate runs the release configuration. |
| D2 | Mate-count mismatch (whole file or per lane): truncate to the shorter mate with a WARNING, or fatal | **Decided: truncate plus WARNING**, per lane. This matches the Flex model and what v1.10.0 does when mate 0 is shorter; when mate 1 is shorter, v1.10.0 appended a malformed pair. Structural mismatches (lane markers, formats) are fatal. |
| D3 | Reader threads outside the permit pool, or reserved permits | **Decided (changed by the author): reader threads take MAP-domain decode permits, only while parsing bytes already in memory** (3.6). The design-time recommendation was "outside". Reader time is reported in `Log.out`. No allocator or host API change. |
| D4 | Keep the forked helper and FIFO (this plan), or read gzip in-thread per mate (Flex-style `gzopen`/`gzdopen`) | **Decided: keep for 1.11.0.** It covers all five cases with one reader and the least new code. Revisit only if M4 shows the FIFO path to be the ceiling. |
| D5 | M4b fast parser | **Decided (author, 30 Sep): implement M4b in 1.11.0.** It is simple and will help in some cases. FASTQ is not text, so DOS versus Linux line endings need no handling: `\r` is an ordinary byte, as on the current getline path. Everything else must match the current reader exactly (line-length limit, the 30 Sep blank-line rules, a missing final newline, lane markers, FASTA and FASTQ, chunk text, `iReadAll`, `readFilesIndex`, Log.out lines); a reader holds a permit only while parsing bytes already in memory. The harness keeps the iostream path as an oracle next to the copy of the old loop; G-R0 is extended (length limits, `\r` bytes, records split across buffer refills); then G-R1, G-R2 and G-R3 are rerun with the new C. |
| D6 | Host gate through Multiomics now or in the Multiomics cycle | **Decided: now (M5), with conditions:** run the Multiomics host gate in M5; a local, unpushed tag only in a scratch clone (never `/mnt/pikachu/STAR-suite` or another shared checkout) and only for the development build; a scratch manifest via `--manifest`; delete the scratch clone afterwards and never push the tag. Background: the motivation is the hosted lane. `build_multiomics.py` checks commit and tree (`scripts/build_multiomics.py:25-29`) and, for a development build, a local tag in the override checkout (43-61). That means a scratch manifest (`--manifest`), a committed C, and a local, unpushed tag in a scratch clone. The rules otherwise forbid tags. |
| D7 | Pre-existing: a helper decompressor failing mid-file ends input without an error (exit status not checked) | **Decided: separate change after 1.11.0.** Fixing it here changes exit codes and would mix behaviour changes into an identity-gated release. |
| D8 | Reader policy text in repository docs | **Resolved.** The constraint above is stated neutrally. Follow the exclusions in the maintainers' private notes. |

## 7. Rules

- **Git.** Work only on `design/rna-mate-threads-20260929`. Commit only when
  the author approves. Plain messages, with no AI or Claude attribution. No
  push, no tag, no release, no merge. Do not create other branches or
  worktrees, and do not touch other checkouts (build sources are fresh clones
  in `$W`).
- **Runbook and handoff.** Keep this runbook and the handoff current at every
  milestone. Near 90% of the usage limit: stop, write the handoff, clean up.
- **Clean room.** Never read, grep or summarize 10x Genomics code. Whitelists
  used as test inputs are data. Stop and ask if a question needs vendor source.
- **Excluded material.** Follow the exclusions in the maintainers' private
  notes.
- **Reader design.** Implement only the design in this runbook: each mate
  read from the start of its own file, coordinated by record counts.
- **Host sharing.**
  - A speed check was running on 29 Sep. Start no builds or runs until the
    author releases the host.
  - Heavy runs use `flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock`.
  - Builds use `nice -n 10` and at most `-j16`.
  - Check `uptime` and `pgrep -af 'STAR|multiomics'` before runs.
  - Never use `pkill -f` with a pattern from your own command line.
- **Writing.** Do not name external users in any doc, log or commit.

## 8. Effort estimate

The estimates are agent working time plus machine time under the lock, assuming
the host is available and gates pass first time.

| M | Working time | Machine time |
|---|---|---|
| M0 | 0.5 day | about 1.5 h (B0 build about 15 min; repeatability pairs; 2M diagnostic) |
| M1 | 0.5 day | about 1 h |
| M2 | 2.5-3 days | builds only |
| M3 | 1 day | 3-4 h |
| M4 | 0.5 day | 1-1.5 h |
| M4b (optional) | 1.5-2 days | about 4 h (re-gates) |
| M5 | 1 day | 2-3 h (Multiomics build plus G-M1) |
| M6 | 1 day | about 2 h (G-S1 took 40-50 min per tree on 28 Sep) |
| **Total** | **about 7-7.5 days** (9-9.5 with M4b) | **about 11-13 h** (15-17 h with M4b) |
