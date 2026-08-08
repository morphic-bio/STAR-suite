# STAR Suite Provenance

This public repository contains compact, normalized reproducibility evidence
for STAR Suite releases, paper results, recipe validation, and portable HPC
guidance. Records are immutable after publication; corrections are new records
that identify what they supersede.

This is not a project run registry. Collaborator manifests, raw locations,
transfer tasks, account names, credentials, full logs, and large outputs stay in
their project provenance repositories or controlled storage.

Every record follows `starsuite.provenance_record/v1` and passes the publication
validator:

```bash
python3 scripts/validate_records.py
```

The schemas, validators, and records are open source under the MIT License. See
[docs/PUBLICATION_POLICY.md](docs/PUBLICATION_POLICY.md).

