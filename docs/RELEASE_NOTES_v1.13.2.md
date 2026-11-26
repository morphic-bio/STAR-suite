# STAR Suite v1.13.2

## Batch mode with sorted BAM output

With `--batchMode 1` and coordinate-sorted BAM output, each sample's sort
files stayed open until the whole batch finished. Open files grew by about the
number of threads times the number of sort bins for every sample, about 200 per
sample at 4 threads. A batch of six samples at 4 threads therefore needed more
than 1,024 open files, a common default limit, and failed. The disk space held
by those files was also not released until the batch ended.

Each sample now closes its sort files when it finishes. The same six-sample
batch stays at about 212 open files. Alignments are unchanged, and batches
that finished before were not affected.

## SLAM per-file trimming

`--trimScope per-file` had two problems when a run had more than one input
file. Runs with one input file and the default `--trimScope first` were not
affected.

- The detection pass for each later file reused the temporary file names of
  mapping thread 0. When thread 0 had already mapped an earlier file, which is
  always the case with `--runThreadN 1`, this discarded that output. The
  coordinate-sorted BAM could be unreadable or incomplete while STAR still
  reported that it finished successfully. Detection now uses its own temporary
  folder.
- Each file is mapped with the same thread workspaces. When a thread wrote less
  unmapped-read or chimeric output for a later file than for an earlier one,
  the end of the earlier output was written again. `Unmapped.out.mate1`,
  `Unmapped.out.mate2`, `Chimeric.out.junction` and `Chimeric.out.sam` could
  then hold repeated and partial records. These temporary files are now
  cleared before each later file.

Rerun analyses that used `--trimScope per-file` with more than one input file.
Without `--autoTrim`, per-file and first-file scopes now give the same
alignments, gene counts, junctions, unmapped reads and chimeric output.

## Regression coverage

`tests/batch_bam_handles_regression.py` compares outputs and peak open files
with a baseline build on single, merged and batch runs.
`tests/slam_per_file_tmp_regression.py` compares per-file SLAM output with a
first-file-scope run at one and four threads, with unmapped reads and
chimeric junction and SAM output.

Existing genome indexes can be reused.
