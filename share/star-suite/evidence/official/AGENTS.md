# Agent instructions

This repository is an append-only public evidence ledger.

- Add compact, machine-readable records; do not copy complete operational run
  directories or collaborator handoff packets.
- Do not publish raw paths, accounts, login hosts, endpoint or task IDs,
  credentials, private sample manifests, or unpublished implementation details.
- Prefer aggregate metrics and logical artifact URIs with checksums.
- Record the source commit, hardware class, method, limitations, and validation
  result needed to interpret each claim.
- Never rewrite a published result silently. Add a correction record and use
  `supersedes`.
- Do not execute biological or performance benchmarks unless the user grants a
  specific execution authorization. Static record validation is safe.

