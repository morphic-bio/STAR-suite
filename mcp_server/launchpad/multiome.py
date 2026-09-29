"""Input and installed-runtime checks for the native Multiome recipe."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import re
import subprocess
import shutil
import hashlib

WORKFLOW_IDS = {"morphic_multiome", "starsuite.official/multiome"}
FASTQ_GROUPS = [("gex_r1", "gex_r2"), ("atac_r1", "atac_barcode", "atac_r2")]
CBQ_INPUTS = ("gex_cbq", "atac_read_pair_cbq", "atac_barcode_cbq")
REFERENCES = ("genome_dir", "gex_whitelist", "chromap_ref", "chromap_index",
              "atac_whitelist", "atac_to_gex")


def prepare_schema(schema):
    """Apply the same form contract to the built-in alias and pinned catalog."""
    from ..schemas.workflow import WorkflowParameterDef, WorkflowParameterGroup
    schema.title = "10x Multiome v1 — RNA + ATAC"
    schema.summary = "Run paired RNA and ATAC libraries to produce RNA matrices, ATAC peaks and a peak-by-cell matrix. Omit I1/I2 sample index files."
    labels = {"gex_r1": "RNA R1 — cell barcode + UMI", "gex_r2": "RNA R2 — cDNA",
              "atac_r1": "ATAC R1 — genomic mate 1", "atac_barcode": "ATAC R2 — cell barcode",
              "atac_r2": "ATAC R3 — genomic mate 2", "out_dir": "New output directory",
              "genome_dir": "STAR genome index", "gex_whitelist": "ARC v1 RNA whitelist",
              "chromap_ref": "Reference genome FASTA", "chromap_index": "Chromap genome index",
              "atac_whitelist": "ARC v1 ATAC whitelist", "atac_to_gex": "ATAC-to-RNA barcode translation"}
    for param in schema.parameters:
        if param.name in labels:
            param.label = labels[param.name]
        if param.name in sum((list(g) for g in FASTQ_GROUPS), []):
            param.widget_hint = "textarea"
            param.help = "Absolute FASTQ paths, comma-separated or one per line, in matching lane order within each library."
        if param.name == "atac_barcode":
            param.description = "Raw 10x ATAC R2: extract bases 9–24 and reverse-complement."
        if param.name == "atac_r2":
            param.description = "The sequencing file named R3 contains ATAC genomic mate 2."
        if param.name == "dry_run":
            param.default = False
        if param.name in ("dry_run", "skip_build", "stop_after_local_mex", "force"):
            param.widget_hint = "hidden"
        if param.name in REFERENCES:
            param.default = None
            param.required = True
        if param.name in ("threads", "chromap_threads"):
            param.min_value = 1
    names = {p.name for p in schema.parameters}
    for name, env, label in (("star_bin", "STAR_BIN", "Chromap-enabled STAR executable"),
                             ("atac_peak_mex_bin", "BUILD_ATAC_MEX_NATIVE", "ATAC peak-matrix executable")):
        if name not in names:
            schema.parameters.append(WorkflowParameterDef(name=name, cli_flag=name, type="file", env_var=env,
                source="runtime",
                label=label, description="Optional override. Otherwise use this installation or PATH.",
                path_must_exist=True, must_be_executable=True))
            schema.rendering.flag_order.append(name)
    schema.parameter_groups.append(WorkflowParameterGroup(name="runtime", title="Installed runtime",
                                   parameters=["star_bin", "atac_peak_mex_bin"]))
    return schema


def runtime_environment(params, root):
    root = Path(root)
    env = {}
    for field, key, candidates in (
        ("star_bin", "STAR_BIN", ("core/legacy/source/STAR", "bin/STAR", "STAR")),
        ("atac_peak_mex_bin", "BUILD_ATAC_MEX_NATIVE", (
            "core/features/libchromap_contract/star_multiome_atac_peak_mex", "bin/star_multiome_atac_peak_mex", "star_multiome_atac_peak_mex")),
    ):
        explicit = params.get(field) or os.environ.get(key)
        if explicit:
            if "/" not in explicit and not explicit.startswith("~"):
                explicit = shutil.which(explicit) or explicit
            path = Path(explicit).expanduser()
            path = (root / path).resolve() if not path.is_absolute() else path.resolve()
            env[key] = str(path)
            continue
        for candidate in candidates:
            found = root / candidate if "/" in candidate else Path(shutil.which(candidate) or root / candidate)
            if found.is_file() and os.access(found, os.X_OK):
                env[key] = str(found.resolve())
                break
    return env


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def runtime_check(params, root):
    env = runtime_environment(params, root)
    for key, title in (("STAR_BIN", "Chromap-enabled STAR"), ("BUILD_ATAC_MEX_NATIVE", "ATAC peak-matrix helper")):
        if key not in env or not Path(env[key]).is_file() or not os.access(env[key], os.X_OK):
            raise ValueError(f"Select an installed {title} executable in the runtime fields, or put it on PATH.")
        if any(char in env[key] for char in ('"', '$', '`', '\\', '\r', '\n', '\0')):
            raise ValueError(f"{title} path contains characters unsupported by the recipe command file.")
    result = subprocess.run([env["STAR_BIN"], "--build-features"], capture_output=True, text=True, timeout=30)
    try:
        features = json.loads(result.stdout)
        chromap = features["chromap_atac"] is True
    except (ValueError, TypeError, KeyError):
        raise ValueError("Select STAR Suite 1.9.5.b or newer with --build-features support and Chromap enabled.") from None
    if result.returncode or not chromap:
        raise ValueError("This STAR binary was built without Chromap. Select a build made with WITH_CHROMAP=1; the portable STAR package cannot run Multiome.")
    return {"binary": env["STAR_BIN"], "binary_sha256": _sha256(env["STAR_BIN"]),
            "features": features, "atac_helper": env["BUILD_ATAC_MEX_NATIVE"],
            "atac_helper_sha256": _sha256(env["BUILD_ATAC_MEX_NATIVE"]), "env": env}


def validate_multiome(params, *, check_paths=False):
    errors = []
    fastq = params.get("input_format", "fastq") == "fastq"
    inputs = [name for group in FASTQ_GROUPS for name in group] if fastq else list(CBQ_INPUTS)
    inactive = list(CBQ_INPUTS) if fastq else [name for group in FASTQ_GROUPS for name in group]
    counts = {}
    for name in inputs:
        value = params.get(name, "")
        if not isinstance(value, str) or not value.strip():
            errors.append((name, f"Choose {name} for {params.get('input_format', 'fastq')} input."))
            continue
        paths = [part.strip() for part in re.split(r"[,\n]", value.strip())]
        if not all(paths):
            errors.append((name, f"{name} contains an empty lane path."))
        params[name] = ",".join(paths)
        counts[name] = len(paths)
        if len(set(paths)) != len(paths):
            errors.append((name, f"{name} repeats a lane path."))
    for name in inactive:
        if params.get(name):
            errors.append((name, f"Clear {name}; it belongs to the other input format."))
    if fastq:
        for group in FASTQ_GROUPS:
            if all(name in counts for name in group) and len({counts[name] for name in group}) != 1:
                errors.append((group[0], f"Matching lane counts are required for {', '.join(group)}."))
    for name in (*inputs, *REFERENCES, "out_dir"):
        value = params.get(name)
        if not value:
            continue  # required reference/output checks are in the shared schema
        if not isinstance(value, str):
            errors.append((name, f"{name} must be a path string."))
            continue
        paths = value.split(",") if name in inputs else [value]
        for path in paths:
            if not isinstance(path, str) or not Path(path).is_absolute():
                errors.append((name, f"Use an absolute server path for {name}."))
            elif any(char in path for char in ('"', '$', '`', '\\', '\r', '\n', '\0')):
                # The existing recipe writes a quoted Bash command file.
                errors.append((name, f"{name} contains characters unsupported by the recipe's command file."))
            elif any(char in str(Path(path).resolve()) for char in ('"', '$', '`', '\\', '\r', '\n', '\0')):
                errors.append((name, f"{name} resolves to a path unsupported by the recipe's command file."))
            elif check_paths and name in inputs and not Path(path).is_file():
                errors.append((name, f"Input file does not exist: {path}"))
            elif check_paths and name in REFERENCES:
                valid = Path(path).is_dir() if name == "genome_dir" else Path(path).is_file()
                if not valid:
                    errors.append((name, f"Reference has the wrong path type: {path}"))
    if not params.get("stop_after_local_mex", True):
        errors.append(("stop_after_local_mex", "Launchpad runs the local matrix/peak workflow; keep 'Stop after matrices' enabled."))
    if not params.get("skip_build", True):
        errors.append(("skip_build", "Build the STAR/Chromap runtime before starting Launchpad; keep skip_build enabled."))
    for name in ("threads", "chromap_threads"):
        value = params.get(name)
        if isinstance(value, int) and value < 1:
            errors.append((name, "Thread counts must be positive."))
    qvalue = params.get("chromap_macs3_frag_qvalue")
    if isinstance(qvalue, (int, float)) and not math.isfinite(qvalue):
        errors.append(("chromap_macs3_frag_qvalue", "Use a finite q-value between 0 and 1."))
    return errors
