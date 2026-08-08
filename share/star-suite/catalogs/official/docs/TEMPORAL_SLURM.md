# Temporal and Slurm execution

Temporal consumes a resolved bundle; it does not discover a catalog or select
between recipes. An upstream CLI or Workbench session must lock the recipe,
helper digests, execution profile, site profile, and partition manifest first.

Each activity retains its partition ordinal across retries. Slurm receives the
locked resource request, and `srun` binds the task to allocated cores. Gather is
a separate activity with dependencies on every required terminal record. A
retry never fetches a newer catalog branch.

Recommended worker environment:

```text
OMP_PROC_BIND=close
OMP_PLACES=cores
srun --cpus-per-task=<locked-value> --cpu-bind=cores ...
```

Site catalogs may specialize queue, account, container runtime, storage root,
and walltime. Those values are execution policy and do not alter recipe
identity.
