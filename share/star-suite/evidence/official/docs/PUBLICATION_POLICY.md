# Public provenance policy

Publish evidence when it helps an independent user interpret a STAR Suite
release, recipe, performance result, or reproducibility claim and can be shared
without exposing operational state.

Include:

- stable record identity and UTC publication date;
- exact public repository revisions or release versions;
- hardware class and thread count when relevant;
- normalized metrics with explicit units;
- parity or validation outcome;
- caveats that bound the interpretation;
- logical artifact identities and checksums when available.

Exclude:

- absolute workstation, cluster, or storage paths;
- scheduler accounts, login hosts, transfer collections, and task IDs;
- credentials, tokens, keys, cookies, and private image registries;
- collaborator sample/delivery manifests and human identifiers;
- raw FASTQ/BAM/MEX/h5ad outputs and large logs;
- unpublished partition-provider commands, formats, architecture, or benchmark
  claims.

Project provenance remains the authority for an executed project run. This
repository publishes a selective, normalized view suitable for release and
paper citation.

