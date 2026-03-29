# Handoff: Post-UCSF Run Cleanup And Test Work

Date: 2026-03-29  
Repo: `/mnt/pikachu/STAR-suite`  
Run-sensitive branch: `production/ucsf-run-20260329`

## Scope

This handoff is for follow-up work that must **not** change STAR behavior while the UCSF production run is active.

Do not modify:

- `core/legacy/source/*.cpp`
- `core/legacy/source/STAR`
- `core/legacy/source/STAR.release`
- `scripts/run_ucsf_perturb_yremove_batch.sh`
- `scripts/paper/run_ucsf_corrected_production_workflow.sh`
- downstream wrappers used by the live run

If work must continue before the UCSF run finishes, do it on a new branch from the current production branch and keep the changes limited to docs, tests, and non-runtime cleanup.

## Current Runtime Constraint

The active UCSF production run is using:

- branch: `production/ucsf-run-20260329`
- release binary: [`STAR.release`](/mnt/pikachu/STAR-suite/core/legacy/source/STAR.release)
- corrected UCSF wrapper: [`run_ucsf_corrected_production_workflow.sh`](/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_corrected_production_workflow.sh)
- batch runner: [`run_ucsf_perturb_yremove_batch.sh`](/mnt/pikachu/STAR-suite/scripts/run_ucsf_perturb_yremove_batch.sh)

The batch runner was updated so that:

- STAR runs complete across the batch first
- downstream starts only after the STAR loop finishes

That ordering is now implemented in [`run_ucsf_perturb_yremove_batch.sh`](/mnt/pikachu/STAR-suite/scripts/run_ucsf_perturb_yremove_batch.sh).

## Safe Work Items Before Run Completion

These are safe while the run is live:

1. Document the Velocyto/BAM `packedReadInfo` lifetime fix and bisect findings.
2. Clean up docs for prerelease planning and merge sequencing.
3. Build a public or allowable smoke-test fixture plan.
4. Prepare non-runtime smoke harnesses that do not touch UCSF-restricted fixtures.
5. Tidy local-only UCSF smoke documentation so restricted artifacts are clearly marked host-local.

## Priority TODOs To Continue

Start from [`docs/todos`](/mnt/pikachu/STAR-suite/docs/todos), especially the top items:

1. Merge the current Velocyto/BAM lifetime fix set to `master` before any UCSF follow-on production work.
2. Cut only a prerelease after UCSF production, not a stable release.
3. Create acceptable/public regression fixtures that can ship with the repo or be fetched from public sources.
4. Do not make a stable/public release until acceptable/public regression coverage exists.

## Recommended Next Steps After UCSF Run Finishes

1. Re-run the local UCSF 100K smoke on the frozen production branch.
2. Validate one full-sample downstream output set, including h5ad surfaces.
3. Validate BAM CB/UB tags again on a small smoke.
4. Create a clean merge branch from `master` or `perturb` and port only the intended fixes.
5. Prepare the `v0.91.0` prerelease only after the public-fixture test plan is in place.

## Notes

- The UCSF-specific smoke remains host-local because the fixture surface is not redistributable.
- If a new smoke is added before release, prefer MSK or a public A375-compatible surface over UCSF-restricted data.
