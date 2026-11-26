"""Checks for the dev-branch TEST process in community/lrsomatic-dev/0.1.

The pipeline schema is read from the pinned commit in the fork checkout (git show <sha>:...), so the
tests do not depend on which branch the fork has checked out. Set LRSOMATIC_FORK to another checkout.
"""
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DEV = ROOT / "community" / "lrsomatic-dev" / "0.1"
V1 = ROOT / "community" / "lrsomatic" / "1.0"
FORK = Path(os.environ.get("LRSOMATIC_FORK", "/net/bmc-lab3/data/bcc/yannvrb/GBM_TeamLab/Cirro_LRsomatic/lrsomatic-fork"))
PARAM_PREFIX = "$.dataset.params."
CIRRO_ONLY_KEYS = {"extra_params_json", "monochrome_logs"}
NEW_BOOLEANS = ["skip_coral", "skip_ampliconclassifier", "modkit_phased", "skip_dmr",
                "skip_padfoot", "skip_reconplot", "skip_shatterseek"]


def load(name):
    return json.loads((DEV / name).read_text())


def pinned_sha():
    return load("process-definition.json")["code"]["version"]


def fork_show(path):
    result = subprocess.run(["git", "-C", str(FORK), "show", f"{pinned_sha()}:{path}"],
                            capture_output=True, text=True)
    if result.returncode != 0:
        pytest.skip(f"{pinned_sha()} is not in the fork checkout {FORK}: {result.stderr.strip()}")
    return result.stdout


def schema_properties():
    schema = json.loads(fork_show("nextflow_schema.json"))
    props = dict(schema.get("properties", {}))
    for group in schema.get("$defs", schema.get("definitions", {})).values():
        props.update(group.get("properties", {}))
    return props


def form_property(form, dotted):
    node = form["form"]
    for part in dotted.split("."):
        assert node.get("type") == "object" and part in node.get("properties", {}), f"{dotted}: no {part!r} in form"
        node = node["properties"][part]
    return node


def param_paths(mapping):
    return {k: v[len(PARAM_PREFIX):] for k, v in mapping.items()
            if isinstance(v, str) and v.startswith(PARAM_PREFIX)}


def version_tuple(text):
    return tuple(int(part) for part in text.split("."))


def test_definition_is_a_separate_test_process_pinned_to_the_fork():
    d = load("process-definition.json")
    assert d["id"] == "process-intgenomicslab-lrsomatic-dev-test"
    assert d["id"] != json.loads((V1 / "process-definition.json").read_text())["id"]
    assert d["code"]["uri"] == "yannvrb/lrsomatic"
    assert re.fullmatch(r"[0-9a-f]{40}", d["code"]["version"]), "pin the fork by full commit SHA"
    assert d["code"]["repository"] == "GITHUBPUBLIC"
    assert d["parentProcessIds"] == ["ingest_aligned_bam", "ingest_unaligned_bam"]
    assert d["allowMultipleSources"] is True and d["usesSampleSheet"] is True


def test_executor_version_satisfies_the_pipeline_manifest():
    # dev declares nextflowVersion '!>=26.04.0' (the 25.10 compiler rejects the LRSOMATIC workflow body)
    manifest = fork_show("nextflow.config")
    required = re.search(r"nextflowVersion\s*=\s*'!?>=([\d.]+)'", manifest)
    assert required, "no nextflowVersion in the pinned manifest"
    assert version_tuple(load("process-definition.json")["code"]["executorVersion"]) >= version_tuple(required.group(1))


def test_every_mapped_param_exists_in_the_pinned_schema():
    mapping = load("process-input.json")
    props = schema_properties()
    unknown = {k for k in mapping if k not in props and k not in CIRRO_ONLY_KEYS}
    assert not unknown, f"mapped params not in the pinned schema: {sorted(unknown)}"
    for flag in NEW_BOOLEANS:
        assert flag in mapping, f"{flag} not mapped"
        assert props[flag]["type"] == "boolean"


def test_every_dataset_param_path_resolves_and_types_match():
    mapping = load("process-input.json")
    form = load("process-form.json")
    props = schema_properties()
    for key, dotted in param_paths(mapping).items():
        form_type = form_property(form, dotted).get("type")
        if key in props:
            schema_type = props[key].get("type")
            assert form_type == schema_type or (form_type == "integer" and schema_type == "number"), \
                f"{key}: form type {form_type} vs schema type {schema_type}"
    for flag in NEW_BOOLEANS:
        assert form_property(form, param_paths(mapping)[flag]).get("default") is False


def test_staged_signatures_volume_is_used_not_downloaded():
    # dev defaults download_sigprofiler_genome to true; with a staged volume it must be off
    mapping = load("process-input.json")
    assert mapping["sigprofiler_genome_dir"] == "s3://btc-refdata/sigprofiler/volume/"
    assert mapping["download_sigprofiler_genome"] is False
    assert schema_properties()["download_sigprofiler_genome"]["default"] is True


def test_staged_references_and_paths_match_version_1():
    new, old = load("process-input.json"), json.loads((V1 / "process-input.json").read_text())
    for key in ("input", "outdir", "sigprofiler_genome_dir", "vep_revel", "vep_revel_tbi", "monochrome_logs", "genome"):
        assert new[key] == old[key], key


def test_preprocess_and_output_are_copies_of_version_1():
    # a copy, not a fork of the logic: a fix in one must be made in both
    assert (DEV / "preprocess.py").read_text() == (V1 / "preprocess.py").read_text()
    assert load("process-output.json") == {"commands": []}


def test_file_mapping_rules_cover_the_new_outputs():
    rules = {r["description"]: r for r in load("process-definition.json")["fileMappingRules"]}
    classifier = rules["Amplicon (ecDNA/BFB) classification"]["fileNamePatterns"][0]
    dmr = rules["Haplotype differential methylation regions"]["fileNamePatterns"][0]
    for pattern in (classifier, dmr):
        regex = pattern["sampleMatchingPattern"].replace("(?P<", "(?<")
        py = pattern["sampleMatchingPattern"]
        match = re.search(py, pattern["exampleName"])
        assert match and match.group("sampleName") == "sampleName", pattern["exampleName"]
        assert regex.count("(?<") == 1
    assert re.search(dmr["sampleMatchingPattern"], "s1/methylation/normal/dmr/s1.nearestGene.bed")
    assert not re.search(dmr["sampleMatchingPattern"], "s1/methylation/normal/dmr/s1.log")


def _config_code():
    text = (DEV / "process-compute.config").read_text()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return "\n".join(re.sub(r"//.*$", "", line) for line in text.splitlines())


def test_compute_config_is_minimal_and_safe():
    code = _config_code()
    assert not re.search(r"^\s*params\s*\{", code, flags=re.M)
    for word in ("slurm", "clusterOptions", "singularity", "queueSize"):
        assert word not in code
    assert re.search(r"resourceLimits\s*=\s*\[\s*cpus:\s*64,\s*memory:\s*488\.GB,\s*time:\s*240\.h\s*\]", code)
    strategy = re.findall(r"errorStrategy\s*=\s*(.+)", code)
    assert len(strategy) == 1 and "task.exitStatus" in strategy[0] and "Integer.MAX_VALUE" in strategy[0]
    assert re.findall(r"maxRetries\s*=\s*(\d+)", code) == ["2"]
    assert re.search(r"withName:\s*'\.\*'\s*\{\s*queue\s*=\s*env\('PW_ONDEMAND_JOB_QUEUE'\)", code)
    # dev's own sizing is kept: no label or per-process size is overridden here
    assert not re.search(r"\b(cpus|memory|time)\s*=", code.replace("resourceLimits", ""))


def test_compute_config_turns_copy_attributes_off():
    # dev sets workflow.output.copyAttributes = true; nf-amazon 3.9.2 rejects COPY_ATTRIBUTES on S3 to S3 copies
    # (run b94197a1: "options are not supported: [COPY_ATTRIBUTES]"), so Cirro's config must override it
    assert re.findall(r"^\s*workflow\.output\.copyAttributes\s*=\s*(\w+)", _config_code(), flags=re.M) == ["false"]


def test_compute_config_drops_clair3_records_with_a_bad_ref():
    # run 6a568a0c: Clair3 wrote 31 records with N in REF (chr3:16902879 TGTGN); bcftools norm exits on the first
    # mismatch by default. Same selector string as dev's conf/modules.config so this config (loaded last) replaces ext.args;
    # dev's own options stay.
    code = _config_code()
    match = re.search(r"withName:\s*'(\.\*:PHASING_HAPLOTYPING:REJOIN_GERMLINE)'\s*\{\s*ext\.args\s*=\s*\{\s*'([^']*)'\s*\}", code)
    assert match, "no REJOIN_GERMLINE ext.args override"
    args = match.group(2).split()
    assert args == ["-m", "+any", "--output-type", "z", "--write-index=tbi", "--check-ref", "x"]


def test_tumor_pileup_rule_exposes_the_files_to_the_timepoint_process():
    rules = {r["description"]: r for r in load("process-definition.json")["fileMappingRules"]}
    pattern = rules["Tumor methylation pileup (bedMethyl)"]["fileNamePatterns"][0]
    assert rules["Tumor methylation pileup (bedMethyl)"]["min"] == 0
    assert pattern["exampleName"] == "sampleName/methylation/tumor/modkit_pileup/sampleName.bed.gz"
    regex = pattern["sampleMatchingPattern"]
    assert re.search(regex, pattern["exampleName"]).group("sampleName") == "sampleName"
    assert not re.search(regex, "s1/methylation/normal/modkit_pileup/s1.bed.gz")
    # the stricter whole-sample selection happens in the timepoint preprocess, not here
    assert re.search(regex, "s1/methylation/tumor/modkit_pileup/s1_combined.bed.gz")
    assert re.search(regex, "s1/methylation/tumor/modkit_pileup/s1_hp1.bed.gz")
