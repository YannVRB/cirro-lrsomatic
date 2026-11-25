"""Checks for the timepoint DMR process folder community/timepoint-dmr/0.1."""
import importlib.util
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TP = ROOT / "community" / "timepoint-dmr" / "0.1"
DEV = ROOT / "community" / "lrsomatic-dev" / "0.1"
FORK = Path(os.environ.get("LRSOMATIC_FORK", "/net/bmc-lab3/data/bcc/yannvrb/GBM_TeamLab/Cirro_LRsomatic/lrsomatic-fork"))
PARAM_PREFIX = "$.dataset.params."
ZERO_SHA = "0" * 40


def load(name, folder=TP):
    return json.loads((folder / name).read_text())


def preprocess_module():
    spec = importlib.util.spec_from_file_location("timepoint_preprocess_for_process_tests", TP / "preprocess.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def form_property(form, dotted):
    node = form["form"]
    for part in dotted.split("."):
        assert node.get("type") == "object" and part in node.get("properties", {}), f"{dotted}: no {part!r} in form"
        node = node["properties"][part]
    return node


def version_tuple(text):
    return tuple(int(part) for part in text.split("."))


def config_code(folder):
    text = (folder / "process-compute.config").read_text()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    lines = (re.sub(r"//.*$", "", line).rstrip() for line in text.splitlines())
    return "\n".join(line for line in lines if line.strip())


def test_definition():
    d = load("process-definition.json")
    assert d["id"] == "process-intgenomicslab-timepoint-dmr-dev-test"
    assert d["parentProcessIds"] == ["process-intgenomicslab-lrsomatic-dev-test"]
    assert d["parentProcessIds"][0] == load("process-definition.json", DEV)["id"]
    assert d["code"]["uri"] == "yannvrb/lrsomatic"
    assert d["code"]["repository"] == "GITHUBPUBLIC"
    assert d["code"]["script"] == "timepoint_dmr.nf"
    assert re.fullmatch(r"[0-9a-f]{40}", d["code"]["version"])
    assert d["usesSampleSheet"] is False
    assert d["allowMultipleSources"] is True
    assert version_tuple(d["code"]["executorVersion"]) >= (26, 4, 0)
    sha = d["code"]["version"]
    if sha != ZERO_SHA:
        shown = subprocess.run(["git", "-C", str(FORK), "show", f"{sha}:nextflow.config"], capture_output=True, text=True)
        if shown.returncode != 0:
            pytest.skip(f"{sha} is not in the fork checkout {FORK}")
        required = re.search(r"nextflowVersion\s*=\s*'!?>=([\d.]+)'", shown.stdout)
        assert required, "no nextflowVersion in the pinned manifest"
        assert version_tuple(d["code"]["executorVersion"]) >= version_tuple(required.group(1))


def test_every_dataset_param_path_resolves_to_a_form_property():
    mapping, form = load("process-input.json"), load("process-form.json")
    paths = {k: v[len(PARAM_PREFIX):] for k, v in mapping.items() if isinstance(v, str) and v.startswith(PARAM_PREFIX)}
    assert set(paths) == {"baseline_sample", "extra_params_json"}
    for dotted in paths.values():
        form_property(form, dotted)


def test_baseline_is_required_text():
    form = load("process-form.json")
    assert "baseline_sample" in form["form"]["required"]
    assert form["form"]["properties"]["baseline_sample"]["type"] == "string"
    assert form["form"]["properties"]["baseline_sample"]["title"] == "Baseline sample"


def test_form_description_names_every_protected_param():
    description = load("process-form.json")["form"]["properties"]["extra_params_json"]["description"]
    missing = [k for k in sorted(preprocess_module()._PROTECTED_PARAMS)
               if not re.search(rf"(?<![\w]){re.escape(k)}(?![\w])", description)]
    assert not missing, f"extra_params_json description does not list: {missing}"


def test_input_mapping():
    mapping = load("process-input.json")
    assert set(mapping) == {"baseline_sample", "extra_params_json", "comparisons", "outdir", "genome", "monochrome_logs"}
    assert mapping["comparisons"] == "$.dataset.s3|/config/comparisons.csv"
    assert mapping["outdir"] == "$.dataset.dataPath"
    assert mapping["genome"] == "GRCh38"
    assert mapping["monochrome_logs"] is True


def test_output_rules_match_examples():
    rules = {r["description"]: r["fileNamePatterns"][0] for r in load("process-definition.json")["fileMappingRules"]}
    assert len(rules) == 2
    for pattern in rules.values():
        match = re.search(pattern["sampleMatchingPattern"], pattern["exampleName"])
        assert match and match.group("sampleName"), pattern["exampleName"]
    nearest = [p for p in rules.values() if "nearestGene" in p["sampleMatchingPattern"]][0]
    assert not re.search(nearest["sampleMatchingPattern"], "s1_vs_s0/dmr/s1_vs_s0.log")
    assert re.search(nearest["sampleMatchingPattern"], "s1_vs_s0/dmr/s1_vs_s0.nearestGene.bed")


def test_compute_config_matches_dev_test_process():
    assert config_code(TP) == config_code(DEV)


def test_preprocess_loads_and_exposes_the_contract():
    assert preprocess_module().COMPARISON_COLUMNS == [
        "comparison", "treatment_name", "treatment_bed", "control_name", "control_bed"]
