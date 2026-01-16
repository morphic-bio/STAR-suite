## GEDI `*.snpdata` Coordinate Convention (Do Not Forget)

### Summary

GEDI SLAM SNP detection (`SlamDetectSnps`) writes `*.snpdata` lines like:

`Location\tCoverage\tMismatches\tP value`

Where **Location is 0-based**.

This is easy to misinterpret because it is printed as `chr:pos` without explicitly saying 0-based.

### Proof (instrumentation-backed)

We verified this by overriding GEDI’s `gedi.slam.javapipeline.SlamDetectSnps` via classpath and logging the internal `pos` used in:

`int pos = r2.map(mpos);`

That `pos` is exactly what GEDI prints as `Location` in `*.snpdata`.

At the debug locus:
- BAM (1-based): `chr1:26429152` has **147** A→G mismatches out of **209** coverage.
- GEDI `snpdata` reports: `1:26429151\t147.0\t209.0\t2.98e-33`

So GEDI is clearly emitting **0-based** positions.

### Correct conversion: GEDI `snpdata` → BED

To convert GEDI `Location` to BED:
- **BED start = Location**
- **BED end = Location + 1**

**Do NOT subtract 1** from GEDI’s `Location`.

Also note: GEDI `*.snpdata` may contain duplicate loci; sort+unique before overlap stats.

### Related scripts

- `tests/filter_gedi_snpdata.sh`: updated to treat GEDI `Location` as 0-based.
- `tests/gedi_debug/run_gedi_snp_debug.sh`: reproduces the instrumentation run.

