# Bridges-2 execution profile

This profile maps portable STAR Suite resource classes to PSC Bridges-2 without
embedding an account, allocation, login host, or private path. PSC uses Slurm
and requires production computation on compute nodes. At the time of this
profile snapshot, Regular Memory nodes expose 128 AMD EPYC cores, usually 256
GB RAM, and 3.84 TB local NVMe; RM-shared supports partial-node jobs of up to 64
cores, while RM and RM-512 are full-node resources.

For formal measurements, request a full node, launch the worker with `srun`,
bind CPUs to cores, and record the exact Slurm request and accounting output.
For throughput runs, use local NVMe for appropriate transient partition state
and move final outputs and provenance to durable project storage before job
cleanup.

Queues, limits, charging, and hardware availability can change. Treat
[the official PSC Bridges-2 user guide](https://www.psc.edu/resources/bridges-2/user-guide/)
as authoritative at submission time.

