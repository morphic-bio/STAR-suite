# Fused FLEX CBQ read-limit defect

Discovered during RAM prototypes based on v1.9.3 (`8276814`). Keep this fix
separate from the RAM representation changes so their benchmark comparison
remains interpretable.

`core/legacy/source/FlexPipeline.cpp`, `flexPrepareCbqRangeTasks`, discards its
`Parameters& P` argument and partitions the sum of all lane records. It does
not honor `P.readMapNumber`. The general planner in
`Parameters_openReadsFiles.cpp` applies that limit, producing contradictory
messages in the same run.

Observed on the immutable v1.9.3 binary, SHA256
`1c6f049a561ea1aafddd29664e8643313ea86d65ae89edef0668d9c1b04752b3`:

- argv requested `--readMapNumber 100000000` on L004 CBQ.
- General reader logged 48 tasks and 100,000,000 records.
- Fully fused FLEX logged 48 ranges and 1,823,648,323 records.
- Both input-read and hash-evaluated totals were 1,823,648,323.
- The experiment wrapper rejected the run's claimed prototype scope; it is
  preserved as a failed validation, not accepted as a 100M experiment.

Evidence: `/scratch/flex_ram_20260912/runs/baseline_proto100m_cbq/` on
`i-037683d7964956c63`; corresponding archive under
`s3://star-suite-320k-benchmark-alt-171440768238-us-west-2-20260904/analysis-tools/flex_ram_20260912/`.
Local orchestration evidence:
`/mnt/pikachu/star_suite_paper/analysis/flex_ram_20260912/diagnose_proto_limit.response.json`.

Workaround used for this experiment: a physically shorter CBQ file containing
99,996,705 pairs in 27,180 complete blocks, with unchanged compressed payloads
and a rebuilt index. Both mate data remain in each original paired block.
`subset/manifest.json` records its boundaries and hash. These runs use scope
`proto100mfile` and omit the ineffective command-line limit.

Future fix should apply the cap before fused task partitioning and test limits
inside a block, across a lane boundary, at zero, and above the total. Compare
read ordinals/counts and counted output against physically restricted input.
Check the other fused input planners for the same parameter contract. Full-set
runs without a read limit are unaffected by this specific issue.
