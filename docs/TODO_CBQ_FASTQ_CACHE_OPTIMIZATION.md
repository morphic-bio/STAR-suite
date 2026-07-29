# Moved: parallel FASTQ sharding

The sharding method is not a CBQ optimization. Its design and work plan now
live in `TODO_PARALLEL_GZIPPED_FASTQ_SHARDING.md`.

CBQ remains an optional downstream sink that can retain the already-decoded
records in an existing indexed binary format. It is not the headlining method
or a prerequisite for using the parallel FASTQ reader.
