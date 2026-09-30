# STAR Suite v1.11.0 Release Notes

Draft. Unreleased; work in progress on branch
`design/rna-mate-threads-20260929`. The release gates against v1.10.0 have
not run yet, and this file will be completed before release.

## FASTQ/FASTA mates parsed on their own reader threads

- New option `--readFilesMateThreads auto|off|on`, default `auto`.
  - With `auto` or `on`, one reader thread per mate parses that mate's input
    from the start of each file, in batches of reads counted from the start
    of the file.
  - The mapping threads pair read *i* of each mate and write the same chunk
    text as before, so for well-formed input mapping sees the same input as
    in v1.10.0.
  - Decompression is unchanged.
- The readers stand down for Flex, SLAM, TranscriptVB, batch mode, the
  spatial raw-R1 tap and non-FASTQ/FASTA input. With `on`, those runs are
  fatal; `off` keeps the v1.10.0 parsing loop.
- Reader threads take MAP-domain decode permits, only while parsing bytes
  already in memory.

## Deliberate changes on malformed input

These inputs are malformed. v1.10.0 handled them silently or inconsistently
between mates.

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
