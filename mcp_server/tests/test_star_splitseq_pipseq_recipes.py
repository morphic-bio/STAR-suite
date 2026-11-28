"""Tests for the star_splitseq_rna and star_pipseq_rna recipes.

Contract: schema load -> validate -> render -> assert the argv.
End to end: the rendered command is run with a built STAR binary on a small
generated genome and generated reads that use the bundled allowlists, and the
count matrices are checked. The end-to-end tests are skipped when no STAR
binary is found (STAR_BIN, or core/legacy/source/STAR of this repository);
with STAR_RECIPE_E2E=required a missing binary is an error instead.
No external data is needed.
"""

import csv
import gzip
import os
import random
import subprocess
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.workflows import (
    render_workflow_command,
    validate_workflow_parameters,
)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SPLITSEQ = "star_splitseq_rna"
PIPSEQ = "star_pipseq_rna"
LISTS = REPO_ROOT / "share" / "star-suite" / "barcode_lists"
SPLITSEQ_LISTS = LISTS / "splitseq_v3" / "48_wells"
PIPSEQ_LISTS = LISTS / "pipseq_readderived"
STAR = Path(os.environ.get("STAR_BIN", REPO_ROOT / "core" / "legacy" / "source" / "STAR"))
if os.environ.get("STAR_RECIPE_E2E") == "required" and not STAR.is_file():
    raise RuntimeError(f"STAR_RECIPE_E2E=required, but there is no STAR binary at {STAR}")
needs_star = pytest.mark.skipif(not STAR.is_file(), reason="no STAR binary (set STAR_BIN or build core/legacy/source/STAR)")


@pytest.fixture
def recipes(tmp_path):
    workflows = []
    for workflow_id in (SPLITSEQ, PIPSEQ):
        yaml_path = REPO_ROOT / "mcp_server" / "workflows" / f"{workflow_id}.yaml"
        with open(yaml_path) as f:
            wf_data = yaml.safe_load(f)
        assert wf_data["id"] == workflow_id
        workflows.append(
            {
                "id": workflow_id,
                "title": wf_data["title"],
                "entry_script": wf_data["entry_script"],
                "kind": wf_data["kind"],
                "schema_file": str(yaml_path),
                "visibility": "public",
            }
        )
    config_dict = {
        "server": {"host": "127.0.0.1", "port": 9999, "auth_token": "", "transport": "http"},
        "paths": {
            "repo_root": str(REPO_ROOT),
            "artifact_log_root": str(tmp_path / "artifacts"),
            "temp_root": str(tmp_path / "tmp"),
        },
        "trusted_roots": [str(REPO_ROOT), str(tmp_path), "/tmp"],
        "datasets": [],
        "scripts": [],
        "test_suites": [],
        "required_binaries": [],
        "workflows": workflows,
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_dict, f)
    load_config(config_path)
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


def splitseq_params(genome_dir, cdna, barcode, out_prefix):
    return {
        "genome_dir": str(genome_dir),
        "read_files_cdna": str(cdna),
        "read_files_barcode": str(barcode),
        "read_files_command": "zcat",
        "cb_position_round1": "50_8",
        "cb_position_round2": "30_8",
        "cb_position_round3": "10_8",
        "umi_position": "0_10",
        "whitelist_round1": str(SPLITSEQ_LISTS / "round1.tsv"),
        "whitelist_round2": str(SPLITSEQ_LISTS / "round2.tsv"),
        "whitelist_round3": str(SPLITSEQ_LISTS / "round3.tsv"),
        "solo_strand": "Forward",
        "out_file_name_prefix": str(out_prefix),
        "run_thread_n": 2,
    }


def pipseq_params(genome_dir, cdna, barcode, out_prefix):
    return {
        "genome_dir": str(genome_dir),
        "read_files_cdna": str(cdna),
        "read_files_barcode": str(barcode),
        "read_files_command": "zcat",
        "whitelist_tier1": str(PIPSEQ_LISTS / "readderived_tier1.tsv"),
        "whitelist_tier2": str(PIPSEQ_LISTS / "readderived_tier2.tsv"),
        "whitelist_tier3": str(PIPSEQ_LISTS / "readderived_tier3.tsv"),
        "whitelist_tier4": str(PIPSEQ_LISTS / "readderived_tier4.tsv"),
        "solo_strand": "Forward",
        "out_file_name_prefix": str(out_prefix),
        "run_thread_n": 2,
    }


# ---------------------------------------------------------------- contract


def test_both_recipes_are_registered_as_public_star_commands():
    with open(REPO_ROOT / "mcp_server" / "config.yaml") as f:
        registered = {w["id"]: w for w in yaml.safe_load(f)["workflows"]}
    for workflow_id in (SPLITSEQ, PIPSEQ):
        entry = registered[workflow_id]
        assert entry["kind"] == "star_cli"
        assert entry["visibility"] == "public"
        assert entry["entry_script"] == "core/legacy/source/STAR"
        assert entry["schema_file"] == f"mcp_server/workflows/{workflow_id}.yaml"
        assert "xperimental" in entry["title"] and "STAR-perturb" in entry["summary"]


def test_splitseq_renders_one_star_command(recipes, tmp_path):
    genome = tmp_path / "genome"
    genome.mkdir()
    params = splitseq_params(genome, "cDNA.fastq.gz", "barcode.fastq.gz", tmp_path / "out/")
    assert validate_workflow_parameters(SPLITSEQ, params).valid
    rendered = render_workflow_command(SPLITSEQ, params)
    assert rendered.argv == [
        str(REPO_ROOT / "core" / "legacy" / "source" / "STAR"),
        "--runMode", "alignReads",
        "--genomeDir", str(genome),
        "--readFilesIn", "cDNA.fastq.gz", "barcode.fastq.gz",
        "--readFilesCommand", "zcat",
        "--soloSplitseqMode", "Fixed",
        "--soloSplitseqCBposition", "50_8", "30_8", "10_8",
        "--soloSplitseqUMIposition", "0_10",
        "--soloSplitseqWhitelist",
        str(SPLITSEQ_LISTS / "round1.tsv"),
        str(SPLITSEQ_LISTS / "round2.tsv"),
        str(SPLITSEQ_LISTS / "round3.tsv"),
        "--soloStrand", "Forward",
        "--soloFeatures", "GeneFull",
        "--soloCrMultimapRescue", "yes",
        "--soloCellFilter", "EmptyDrops_CR",
        "--outSAMtype", "None",
        "--outFileNamePrefix", str(tmp_path / "out/"),
        "--runThreadN", "2",
    ]
    assert rendered.output_root == str(tmp_path / "out/")
    edit2 = render_workflow_command(SPLITSEQ, {**params, "solo_splitseq_mode": "FixedEdit2",
                                               "solo_splitseq_min_posterior": 0.975,
                                               "solo_splitseq_frequency_backup": "no"})
    assert edit2.argv[edit2.argv.index("--soloSplitseqMode") + 1] == "FixedEdit2"
    assert edit2.argv[edit2.argv.index("--soloSplitseqMinPosterior") + 1] == "0.975"
    assert edit2.argv[edit2.argv.index("--soloSplitseqFrequencyBackup") + 1] == "no"
    assert not validate_workflow_parameters(SPLITSEQ, {**params, "solo_splitseq_mode": "FixedEdit2",
                                                       "solo_splitseq_min_posterior": 95}).valid


@pytest.mark.parametrize(
    "missing",
    ["cb_position_round1", "cb_position_round3", "umi_position", "whitelist_round2", "solo_strand"],
)
def test_splitseq_layout_lists_and_strand_have_no_default(recipes, tmp_path, missing):
    genome = tmp_path / "genome"
    genome.mkdir()
    params = splitseq_params(genome, "cDNA.fastq", "barcode.fastq", tmp_path / "out/")
    del params[missing]
    result = validate_workflow_parameters(SPLITSEQ, params)
    assert not result.valid
    assert any(missing in error for error in result.errors)


def test_splitseq_rejects_a_missing_list_and_a_bam_request(recipes, tmp_path):
    genome = tmp_path / "genome"
    genome.mkdir()
    params = splitseq_params(genome, "cDNA.fastq", "barcode.fastq", tmp_path / "out/")
    params["whitelist_round1"] = str(tmp_path / "no_such_list.tsv")
    assert not validate_workflow_parameters(SPLITSEQ, params, check_paths=True).valid
    params = splitseq_params(genome, "cDNA.fastq", "barcode.fastq", tmp_path / "out/")
    params["out_sam_type"] = "BAM"
    assert not validate_workflow_parameters(SPLITSEQ, params).valid


def test_pipseq_renders_one_star_command(recipes, tmp_path):
    genome = tmp_path / "genome"
    genome.mkdir()
    params = pipseq_params(genome, "R2.fastq.gz", "R1.fastq.gz", tmp_path / "out/")
    assert validate_workflow_parameters(PIPSEQ, params).valid
    rendered = render_workflow_command(PIPSEQ, params)
    assert rendered.argv == [
        str(REPO_ROOT / "core" / "legacy" / "source" / "STAR"),
        "--runMode", "alignReads",
        "--genomeDir", str(genome),
        "--readFilesIn", "R2.fastq.gz", "R1.fastq.gz",
        "--readFilesCommand", "zcat",
        "--soloPipseqMode", "Illumina",
        "--soloPipseqWhitelist",
        str(PIPSEQ_LISTS / "readderived_tier1.tsv"),
        str(PIPSEQ_LISTS / "readderived_tier2.tsv"),
        str(PIPSEQ_LISTS / "readderived_tier3.tsv"),
        str(PIPSEQ_LISTS / "readderived_tier4.tsv"),
        "--soloStrand", "Forward",
        "--soloFeatures", "GeneFull_Ex50pAS",
        "--soloCellFilter", "EmptyDrops_CR",
        "--outSAMtype", "None",
        "--outFileNamePrefix", str(tmp_path / "out/"),
        "--runThreadN", "2",
    ]
    assert rendered.output_root == str(tmp_path / "out/")


def test_pipseq_optional_settings_render_when_changed(recipes, tmp_path):
    genome = tmp_path / "genome"
    genome.mkdir()
    params = pipseq_params(genome, "R2.fastq", "R1.fastq", tmp_path / "out/")
    del params["read_files_command"]
    params.update(solo_pipseq_trim="None", solo_pipseq_ipm=1.35, solo_pipseq_multi_match="FirstListed",
                  solo_cell_filter="None")
    assert validate_workflow_parameters(PIPSEQ, params).valid
    argv = render_workflow_command(PIPSEQ, params).argv
    assert "--readFilesCommand" not in argv
    for flag, value in (("--soloPipseqTrim", "None"), ("--soloPipseqIPM", "1.35"),
                        ("--soloPipseqMultiMatch", "FirstListed"), ("--soloCellFilter", "None")):
        assert argv[argv.index(flag) + 1] == value
    params["solo_pipseq_ipm"] = -1
    assert not validate_workflow_parameters(PIPSEQ, params).valid


@pytest.mark.parametrize("missing", ["whitelist_tier1", "whitelist_tier4", "solo_strand"])
def test_pipseq_lists_and_strand_have_no_default(recipes, tmp_path, missing):
    genome = tmp_path / "genome"
    genome.mkdir()
    params = pipseq_params(genome, "R2.fastq", "R1.fastq", tmp_path / "out/")
    del params[missing]
    result = validate_workflow_parameters(PIPSEQ, params)
    assert not result.valid
    assert any(missing in error for error in result.errors)


def read_list(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if line.strip()]


def test_bundled_lists_have_the_shape_the_recipes_describe():
    for wells, round1_wells in (("48_wells", 48), ("96_wells", 96)):
        folder = LISTS / "splitseq_v3" / wells
        round1 = read_list(folder / "round1.tsv")
        assert all(len(row) == 2 and len(row[0]) == 8 and len(row[1]) == 8 for row in round1)
        assert len({row[1] for row in round1}) == round1_wells
        for name in ("round2.tsv", "round3.tsv"):
            rows = read_list(folder / name)
            assert len(rows) == 96 and all(len(row) == 1 and len(row[0]) == 8 for row in rows)
    for tier, length in ((1, 8), (2, 6), (3, 6), (4, 8)):
        rows = read_list(PIPSEQ_LISTS / f"readderived_tier{tier}.tsv")
        assert len(rows) == 96 and all(len(row) == 1 and len(row[0]) == length for row in rows)
        assert len({row[0] for row in rows}) == 96
    for folder in (LISTS / "splitseq_v3", PIPSEQ_LISTS):
        assert (folder / "PROVENANCE.md").is_file()


# -------------------------------------------------------------- end to end


def write_fastq_gz(path, sequences):
    with gzip.open(path, "wt") as handle:
        for i, sequence in enumerate(sequences):
            handle.write("@r%d\n%s\n+\n%s\n" % (i, sequence, "I" * len(sequence)))


def small_genome(folder, seed):
    """An 8 kb random genome with one gene on each strand, and its STAR index."""
    rng = random.Random(seed)
    genome = "".join(rng.choice("ACGT") for _ in range(8000))
    (folder / "genome.fa").write_text(">chrSynthetic\n" + genome + "\n")
    (folder / "genes.gtf").write_text(
        'chrSynthetic\tfixture\texon\t501\t3500\t.\t+\t.\tgene_id "PLUS"; transcript_id "PLUS.1"; gene_name "PLUS";\n'
        'chrSynthetic\tfixture\texon\t4501\t7500\t.\t-\t.\tgene_id "MINUS"; transcript_id "MINUS.1"; gene_name "MINUS";\n'
    )
    index = folder / "index"
    index.mkdir()
    subprocess.run(
        [str(STAR), "--runMode", "genomeGenerate", "--genomeDir", str(index),
         "--genomeFastaFiles", str(folder / "genome.fa"), "--sjdbGTFfile", str(folder / "genes.gtf"),
         "--sjdbOverhang", "49", "--genomeSAindexNbases", "3", "--genomeChrBinNbits", "12",
         "--runThreadN", "1", "--outFileNamePrefix", str(folder / "index_")],
        check=True, capture_output=True,
    )
    return genome, index


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def raw_counts(folder):
    """{(barcode, gene): count} of a raw MEX folder."""
    barcodes = (folder / "barcodes.tsv").read_text().split()
    genes = [line.split("\t")[0] for line in (folder / "features.tsv").read_text().splitlines()]
    lines = [line for line in (folder / "matrix.mtx").read_text().splitlines() if not line.startswith("%")]
    return {(barcodes[int(c) - 1], genes[int(g) - 1]): int(n) for g, c, n in (line.split() for line in lines[1:])}


@needs_star
def test_splitseq_rendered_command_counts_reads_with_the_bundled_lists(recipes, tmp_path):
    genome, index = small_genome(tmp_path, 2703)
    round1 = read_list(SPLITSEQ_LISTS / "round1.tsv")
    round2 = [row[0] for row in read_list(SPLITSEQ_LISTS / "round2.tsv")]
    round3 = [row[0] for row in read_list(SPLITSEQ_LISTS / "round3.tsv")]
    # two entries of one round-1 well (they share the replacement target) and an entry of another well
    target = round1[0][1]
    same_well = [row[0] for row in round1 if row[1] == target]
    assert len(same_well) >= 2
    other = next(row for row in round1 if row[1] != target)

    def barcode_read(first, second, third, umi):
        # v3 layout: UMI 0-9, round 3 10-17, linker, round 2 30-37, linker, round 1 50-57
        return umi + third + "ATGAGGGGTCAG" + second + "TCCAACCACCTC" + first

    plus, minus = genome[700:750], reverse_complement(genome[5000:5050])
    umis = ["ACGTACGTAC", "GCTAGCTAGC", "TTGACCATGA"]
    rows = []
    rows += [(plus, barcode_read(same_well[0], round2[0], round3[0], umis[0]))] * 3  # one molecule, three reads
    rows += [(plus, barcode_read(same_well[0], round2[0], round3[0], umis[1]))]
    rows += [(plus, barcode_read(same_well[1], round2[0], round3[0], umis[2]))]      # same well, other primer
    rows += [(minus, barcode_read(same_well[0], round2[0], round3[0], umis[0]))]
    rows += [(plus, barcode_read(other[0], round2[1], round3[2], umis[0]))]          # a second cell
    rows += [(plus, barcode_read("N" + same_well[0][1:], round2[0], round3[0], umis[1]))]  # one N: same molecule
    rows += [(plus, barcode_read(same_well[0], "N" * 8, round3[0], umis[0]))]        # no valid barcode
    cdna, barcode = tmp_path / "cDNA.fastq.gz", tmp_path / "barcode.fastq.gz"
    write_fastq_gz(cdna, [row[0] for row in rows])
    write_fastq_gz(barcode, [row[1] for row in rows])

    out = tmp_path / "run"
    out.mkdir()
    params = splitseq_params(index, cdna, barcode, str(out) + "/")
    assert validate_workflow_parameters(SPLITSEQ, params, check_paths=True).valid
    rendered = render_workflow_command(SPLITSEQ, params)
    argv = [str(STAR)] + rendered.argv[1:]
    result = subprocess.run(argv, capture_output=True, text=True, cwd=out)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]

    cell_a = target + round2[0] + round3[0]
    cell_b = other[1] + round2[1] + round3[2]
    solo = out / "Solo.out" / "GeneFull"
    assert raw_counts(solo / "raw") == {(cell_a, "PLUS"): 3, (cell_a, "MINUS"): 1, (cell_b, "PLUS"): 1}
    assert (solo / "filtered" / "matrix.mtx").is_file()
    metrics = dict(csv.reader((out / "Solo.out" / "Splitseq.metrics.tsv").open(), delimiter="\t"))
    assert metrics["reads_total"] == str(len(rows))
    assert metrics["barcode_exact"] == str(len(rows) - 2)
    assert metrics["barcode_corrected"] == "1"


@needs_star
def test_pipseq_rendered_command_counts_reads_with_the_bundled_lists(recipes, tmp_path):
    genome, index = small_genome(tmp_path, 7319)
    tiers = [read_list(PIPSEQ_LISTS / f"readderived_tier{i}.tsv")[0][0] for i in (1, 2, 3, 4)]
    body = tiers[0] + "ATG" + tiers[1] + "GAG" + tiers[2] + "TCGAG" + tiers[3]
    cell = "".join(tiers)

    def binning_index(n):
        return "".join("ACGT"[(n >> shift) & 3] for shift in (4, 2, 0))

    rows = []
    # PLUS: 11 bins, two fragment starts per bin, three reads per fragment, at changing barcode phase
    for b in range(11):
        for start in (700, 900):
            for copy in range(3):
                rows.append((genome[start:start + 50], "G" * ((b + copy) % 4) + body + binning_index(b) + "TTTT"))
    # MINUS: 10 bins, two fragment ends per bin
    for b in range(10):
        for end in (5100, 5300):
            rows.append((reverse_complement(genome[end - 50:end]), "G" * (b % 4) + body + binning_index(b) + "TTTT"))
    rows.append((genome[700:750], "N" * 39 + "AAA" + "TTTT"))  # no valid barcode
    cdna, barcode = tmp_path / "R2.fastq.gz", tmp_path / "R1.fastq.gz"
    write_fastq_gz(cdna, [row[0] for row in rows])
    write_fastq_gz(barcode, [row[1] for row in rows])

    out = tmp_path / "run"
    out.mkdir()
    params = pipseq_params(index, cdna, barcode, str(out) + "/")
    # these generated cDNA reads carry no kit sequence
    params["solo_pipseq_trim"] = "None"
    assert validate_workflow_parameters(PIPSEQ, params, check_paths=True).valid
    rendered = render_workflow_command(PIPSEQ, params)
    argv = [str(STAR)] + rendered.argv[1:]
    result = subprocess.run(argv, capture_output=True, text=True, cwd=out)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]

    solo = out / "Solo.out" / "GeneFull_Ex50pAS"
    ledger = list(csv.DictReader((solo / "Pipseq.groups.tsv").open(), delimiter="\t"))
    assert {row["gene_id"]: int(row["occupied_bins"]) for row in ledger} == {"PLUS": 11, "MINUS": 10}
    assert {row["gene_id"]: int(row["unique_bi_imi"]) for row in ledger} == {"PLUS": 22, "MINUS": 20}
    assert raw_counts(solo / "raw") == {(cell, "PLUS"): 11, (cell, "MINUS"): 10}
    assert (solo / "filtered" / "barcodes.tsv").read_text().split() == [cell]
    metrics = dict(csv.reader((solo / "Pipseq.metrics.tsv").open(), delimiter="\t"))
    assert metrics["reads_total"] == str(len(rows))
    assert metrics["barcode_exact"] == str(len(rows) - 1)
