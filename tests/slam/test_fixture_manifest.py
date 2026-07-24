#!/usr/bin/env python3
"""Unit tests for the pinned SLAM fixture verifier."""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("verify_fixture_manifest.py")
SPEC = importlib.util.spec_from_file_location("verify_fixture_manifest", MODULE_PATH)
assert SPEC and SPEC.loader
VERIFIER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFIER)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class FixtureManifestTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fastq = self.root / "reads.fastq.gz"
        self.fastq_content = b"@r1\nACGT\n+\nIIII\n"
        with gzip.GzipFile(filename=str(self.fastq), mode="wb", mtime=0) as handle:
            handle.write(self.fastq_content)

        self.snp_bed = self.root / "snps.bed"
        self.snp_bed.write_bytes(b"chr2\t2\t3\nchr1\t1\t2\n")
        self.reference = self.root / "reference.tsv.gz"
        self.reference_content = b"Gene\tMAP\nG1\t0.1\n"
        with gzip.GzipFile(filename=str(self.reference), mode="wb", mtime=0) as handle:
            handle.write(self.reference_content)

        self.index = self.root / "index"
        self.index.mkdir()
        (self.index / "Genome").write_bytes(b"genome")
        (self.index / "SA").write_bytes(b"sa")
        (self.index / "SAindex").write_bytes(b"sai")
        (self.index / "geneInfo.tab").write_bytes(b"1\nG1\tG1\tgene\n")
        (self.index / "genomeParameters.txt").write_text(
            "versionGenome\t2.7.4a\n"
            "genomeSAindexNbases\t14\n"
            "genomeChrBinNbits\t18\n"
            "genomeSAsparseD\t1\n"
            "sjdbOverhang\t100\n",
            encoding="utf-8",
        )

        sorted_bed = b"".join(sorted(self.snp_bed.read_bytes().splitlines(keepends=True)))
        self.manifest = {
            "fixture": "unit-test",
            "alignment_contract": {
                "required_options": {
                    "--clip3pAdapterSeq": "AGATCGGAAGAG",
                    "--clip3pAdapterMMp": "0.1",
                    "--outSAMtype": "None",
                    "--readFilesCommand": "zcat",
                    "--slamGrandSlamOut": "1",
                    "--slamQuantMode": "1",
                },
                "required_flags": [
                    "--slamQcReport",
                    "--slamDumpBinary",
                    "--slamDumpWeights",
                    "--slamSnpMaskIn",
                ],
                "forbidden_options": ["--autoTrim", "--slamCompatTrim5p"],
                "forbidden_option_values": {"--slamQcReport": ["-"]},
            },
            "parity_contract": {
                "ref_corr_min": 0.99,
                "exact_min_pearson": 0.999999,
                "exact_max_abs_delta": 0.000001,
                "align_min_pearson": 0.98,
            },
            "artifacts": {
                "fastq": {
                    "hash_mode": "gzip_content",
                    "sha256": digest(self.fastq_content),
                    "fastq_records": 1,
                },
                "snp_bed": {
                    "hash_mode": "sorted_lines",
                    "sha256": digest(sorted_bed),
                    "line_count": 2,
                },
                "reference": {
                    "hash_mode": "gzip_content",
                    "sha256": digest(self.reference_content),
                },
            },
            "index": {
                "files": {
                    "Genome": {"size": 6},
                    "SA": {"size": 2},
                    "SAindex": {"size": 3},
                    "geneInfo.tab": {
                        "sha256": digest((self.index / "geneInfo.tab").read_bytes())
                    },
                },
                "parameters": {
                    "versionGenome": "2.7.4a",
                    "genomeSAindexNbases": "14",
                    "genomeChrBinNbits": "18",
                    "genomeSAsparseD": "1",
                    "sjdbOverhang": "100",
                },
            },
        }

    def tearDown(self) -> None:
        self.temp.cleanup()

    def command(self) -> list[str]:
        return [
            "STAR",
            "--genomeDir",
            str(self.index),
            "--readFilesIn",
            str(self.fastq),
            "--readFilesCommand",
            "zcat",
            "--outSAMtype",
            "None",
            "--clip3pAdapterSeq",
            "AGATCGGAAGAG",
            "--clip3pAdapterMMp",
            "0.1",
            "--slamQuantMode",
            "1",
            "--slamSnpMaskIn",
            str(self.snp_bed),
            "--slamGrandSlamOut",
            "1",
            "--slamQcReport",
            "qc/star_slam",
            "--slamDumpBinary",
            "dump.bin",
            "--slamDumpWeights",
            "weights.bin",
        ]

    def verify(
        self,
        command: list[str] | None = None,
        thresholds: dict[str, float] | None = None,
    ) -> list[str]:
        return VERIFIER.verify(
            self.manifest,
            self.fastq,
            self.snp_bed,
            self.reference,
            self.index,
            command or self.command(),
            thresholds
            or {
                "ref_corr_min": 0.99,
                "exact_min_pearson": 0.999999,
                "exact_max_abs_delta": 0.000001,
                "align_min_pearson": 0.98,
            },
        )

    def test_accepts_matching_fixture_and_command(self) -> None:
        self.assertEqual([], self.verify())

    def test_rejects_missing_adapter_contract(self) -> None:
        command = self.command()
        start = command.index("--clip3pAdapterSeq")
        del command[start : start + 2]
        errors = self.verify(command)
        self.assertTrue(any("--clip3pAdapterSeq" in error for error in errors))

    def test_rejects_forbidden_legacy_trim(self) -> None:
        errors = self.verify(self.command() + ["--autoTrim", "-"])
        self.assertTrue(any("forbidden fixture option" in error for error in errors))

    def test_rejects_default_qc_output_prefix(self) -> None:
        command = self.command()
        command[command.index("--slamQcReport") + 1] = "-"
        errors = self.verify(command)
        self.assertTrue(any("--slamQcReport has forbidden value" in error for error in errors))

    def test_rejects_index_drift(self) -> None:
        (self.index / "Genome").write_bytes(b"changed")
        errors = self.verify()
        self.assertTrue(any("Genome size mismatch" in error for error in errors))

    def test_rejects_relaxed_parity_threshold(self) -> None:
        errors = self.verify(
            thresholds={
                "ref_corr_min": 0.9,
                "exact_min_pearson": 0.999999,
                "exact_max_abs_delta": 0.000001,
                "align_min_pearson": 0.98,
            }
        )
        self.assertTrue(any("ref_corr_min cannot be relaxed" in error for error in errors))

    def test_gzip_metadata_does_not_change_content_fingerprint(self) -> None:
        second = self.root / "same-content.fastq.gz"
        with gzip.GzipFile(filename=str(second), mode="wb", mtime=1234) as handle:
            handle.write(self.fastq_content)
        self.assertEqual(
            VERIFIER.artifact_digest(self.fastq, "gzip_content"),
            VERIFIER.artifact_digest(second, "gzip_content"),
        )


if __name__ == "__main__":
    unittest.main()
