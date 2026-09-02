# Lustre scatter/gather recipe

The portable rule is to separate high-churn partition work from durable
results. A scheduler expands a locked partition manifest, gives each ordinal a
private work directory, and writes one compact terminal record per ordinal.
Gather starts only after all required records are successful.

For Lustre deployments:

1. Keep the manifest, locks, and compact logs on durable project storage.
2. Prefer node-local scratch for temporary sorting, decompression, and other
   high-churn intermediate state.
3. Use one output directory per partition. Avoid shared append-only logs.
4. Apply striping only to newly created large sequential outputs. There is no
   universal stripe count; validate the site-selected value for the file size,
   concurrency, and OST availability.
5. Gather on a resource sized for aggregate BAM throughput and final sorting.
6. Write the final artifact and checksum before deleting temporary partitions.

The profile supplies a generic `lfs setstripe` template but intentionally does
not hard-code a stripe count. Site policy owns that decision.

