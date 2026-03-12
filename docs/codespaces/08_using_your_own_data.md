# Using your own data

You do not need to use the small public demo files.

You can replace them with your own files at any time.

## The short version

For most real datasets, you will usually need to change only these things:
- the input FASTQ files
- the STAR index
- the reference files used to build that index
- the whitelist and barcode settings for single-cell data
- the feature reference or probe set for perturb or Flex

## 1. If you use your own reference

If your data is not meant to run against the small demo reference, use your own full reference instead.

You will need:
- a genome FASTA file
- a matching GTF file
- a STAR index built from that FASTA and GTF

In practice:
- replace the demo `genome.fa` and `genes.gtf`
- build a new STAR index
- point the run script at your new `--genome-dir`

Important rule:
- the FASTA and GTF should match each other
- the STAR index should be built from that same FASTA and GTF pair

## 2. Bulk and SLAM

For bulk or SLAM runs, the main changes are simple:
- replace the public FASTQ files with your own FASTQ files
- point `--genome-dir` to the STAR index you want to use

Use the demo guides mainly as examples of the command shape.

## 3. Single-cell data

For single-cell data, there are two extra things to think about:
- whitelist
- barcode and UMI positions

### Whitelist

The whitelist must match the chemistry you are using.

If you switch to a different assay or chemistry, you may need a different whitelist file.

### Barcode and UMI settings

These settings tell STAR where the cell barcode and UMI are found in the read.

If your assay uses different positions or lengths, you need to change:
- `soloCBstart`
- `soloCBlen`
- `soloUMIstart`
- `soloUMIlen`

The small demo uses values that match its public source data.
Your real dataset may need different values.

## 4. Perturb

For real perturb data, you will usually need:
- your own GEX FASTQs
- your own guide FASTQs
- the right whitelist for the assay
- the right barcode and UMI settings
- the right feature reference

If you are using a Cell Ranger-style config, make sure it points to the real files you want to run.

## 5. Flex

For real Flex data, you will usually need:
- your own FASTQs
- the right whitelist for the assay
- the right barcode and UMI settings
- the right probe set
- a filtered reference and STAR index built from the reference and probe set you actually want to use

## 6. A safe way to switch from demo data to real data

A simple order is:
1. get the demo guide working first
2. replace the STAR index
3. replace the FASTQ files
4. for single-cell runs, update whitelist and barcode settings
5. for perturb or Flex, update the feature reference or probe set
6. rerun and check the written `RUN_COMMAND.sh`

## 7. When in doubt

If you are unsure what to change:
- start from the guide that is closest to your assay
- run the demo in `--dry-run` mode first
- read the `RUN_COMMAND.sh` file it writes
- replace one thing at a time
