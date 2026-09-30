# STAR Suite v1.11.0 Release Notes

Unreleased. The candidate source is branch `design/rna-mate-threads-20260929`
(the validated commit is named under Validation); tagging, the version bump
and packaging are part of the separate release step. Upstream STAR remains
`2.7.11b`, genome-index compatibility remains `2.7.4a`, and legacy
compatibility remains `2.7.1a`. Existing indexes do not need rebuilding.

STAR Suite 1.11.0 moves FASTQ/FASTA parsing out of the mapping threads' input
lock. Each mate is parsed on its own reader thread, and the mapping threads
only pair the parsed reads, so mapping sees the same input as before. A new
option controls it. Malformed input is now handled the same way in every
mate.

## FASTQ/FASTA mates parsed on their own reader threads

- New option `--readFilesMateThreads auto|off|on`, default `auto`.
  - With `auto` or `on`, one reader thread per mate parses that mate's input
    from the start of each file, in batches of reads counted from the start
    of the file. The mapping threads pair read *i* of each mate and write the
    same chunk text as before, so for well-formed input mapping sees the same
    input as in v1.10.0.
  - The readers use a buffered parser that scans the raw bytes for line ends
    instead of iostream calls. Its results are those of the v1.10.0 parsing
    loop, including the line-length limits; `\r` is an ordinary byte, as
    before.
  - Decompression is unchanged: the per-mate gzip helper, a
    `--readFilesCommand` producer or the BGZF reader feeds each mate as
    before.
  - `off` keeps the v1.10.0 parsing loop.
- The readers stand down, and the v1.10.0 loop runs, for Flex, SLAM,
  TranscriptVB, batch mode, the spatial raw-R1 tap, SAM and CBQ input, and
  more than three mates. With `on`, those runs are fatal. `Log.out` states
  whether the readers are active and why not.
- Reader threads take MAP-domain decode permits (as the BGZF inflate workers
  do), only while parsing bytes already in memory; a reader never holds a
  permit while it waits for input or queue space. The host interface and
  the permit allocator are unchanged. The end-of-mapping permit telemetry
  line in `Log.out` counts acquisitions up to the end of GEX mapping only;
  with faster GEX mapping it includes less of any concurrent feature-arm
  work, which is otherwise unchanged.
- `Log.out` gets a per-mate reader summary when the input closes.
- Debug: with `STAR_INPUT_CHUNK_TRACE`, setting
  `STAR_INPUT_CHUNK_TRACE_DIGEST=1` adds the CRC32 of each mate's chunk text
  to the trace.

## Deliberate changes on malformed input

These inputs are malformed. v1.10.0 handled them silently or inconsistently
between mates. With `--readFilesMateThreads off` the v1.10.0 behaviour
remains.

- **Blank lines where a read header is expected** (the author's decisions
  of 29 and 30 Sep). The rule is the same for every mate, in FASTQ and FASTA.
  A blank line (only spaces, tabs or carriage returns) is never read as a
  read header.
  - A single blank line followed by a read header in that file's format (`@`
    or `>`) is skipped.
  - Blank lines at the end of a file, followed by end of input or by the next
    input file, are skipped, whatever their number.
  - Skipped lines get a WARNING in `Log.out`, naming the file, mate and line,
    at the first one in each file, and each file's count of skipped blank
    lines when the input closes.
  - A second blank line in a row before the end of the file, or a blank line
    followed by anything other than a read header, stops STAR with a fatal
    error that names the file and line and says the input is malformed.
  - v1.10.0 ended the whole input at such a line in mate 1, skipped it
    silently in mates 2 and 3, and with FASTA plus `--outSAMreadID Number`
    could read it as mate 2's header.
- **Mates with different numbers of reads in an input file.** STAR maps the
  reads all mates have, skips the rest of that file, and writes a WARNING
  with the counts. v1.10.0 truncated silently when mate 1 was shorter, and
  appended a malformed pair when mate 2 was.
- **Input files out of step between mates** (different lane-marker
  structure) and **mates in different formats** (FASTQ in one, FASTA in
  another) are now fatal.

## Known limitations

- Unchanged from v1.10.0: with several FASTA input files per mate, a lane
  marker right after a read's sequence lines is read as sequence (the FASTA
  loop stops only at `@`, `>`, a space or a newline). With the readers
  active, ending each FASTA file with a blank line avoids it.
- Unchanged from v1.10.0: a decompression helper that fails mid-file ends the
  input without an error (its exit status is not checked). A fix is planned
  as a separate change.

## Validation

The validated code is commit `81e03b6` on branch
`design/rna-mate-threads-20260929`; later commits change documentation only.
Comparisons are against the v1.10.0 code at `0f9701a` (under `core/` it
differs from the `v1.10.0` tag only by the `--build-features` flag, which
does not affect mapping). Runs used the shared host lock; no timing
results are claimed.

- **Reader harness (G-R0).** `fastx_mate_reader_harness` passes: 336
  identity runs against a copy of the v1.10.0 parsing loop (SE, PE and three
  mates; FASTQ, FASTA and multi-line FASTA; CRLF and other `\r` bytes; lane
  markers; line-length limits; chunk sizes; `--readMapNumber`;
  `--outSAMreadID Number`; with and without a one-permit pool), including
  the `Log.out` lines, with no new WARNING; the mismatched-mate and
  blank-line cases in each of three mates; input buffers of 1 to 4,096
  bytes, so every record is split across refills; and the buffered parser
  against the iostream parser on every case, including malformed inputs.
  Repeated runs and a ThreadSanitizer build are clean.
- **Input identity (G-R1).** On every fixture below, the chunk trace with
  per-mate CRC32s is equal row for row to that of the v1.10.0 code with only
  the trace digest added.
- **Output identity (G-R2).** Outputs equal those of the v1.10.0 code:
  byte for byte at one thread and for PBMC 3k, bulk paired-end and the
  scRNA 100K regression at eight threads; for fixtures with feature
  libraries, only in files that also vary between two runs of the v1.10.0
  code (feature-arm row order and match positions). Fixtures: PBMC 3k 100K
  GEX with `zcat`, internal gzip, BGZF, plain lanes, a single plain file, a
  manifest and one thread; DOGMA-plex lane 1 at 2M and 20M reads; HIV DOGMA
  100K; CAT-ATAC GEX with an early `--readMapNumber` stop, with and without
  `zcat`; bulk SRR4422207 500K pairs with sorted BAM and gene counts, one
  thread, unsorted BAM, two-pass, `--outFilterType BySJout`, single-end and
  plain FASTQ; the three scRNA 100K regression profiles.
- **Flex unchanged (G-R3).** The Flex half-probe 100K smoke reproduces all
  121 outputs of the preserved reference; the readers report "not active".
- **Host (G-S2, G-M1).** The host API tests pass, with the readers active
  in all 31 runs. A Multiomics Suite development build against the
  candidate passed G-M1 on all five fixtures (run on the candidate before
  the buffered parser; the parser change was then re-gated with G-R0 to
  G-R3).
- **Regression (G-S1).** Pending.
