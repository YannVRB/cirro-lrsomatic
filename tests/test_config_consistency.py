import importlib.util
import json
import os
import re
from pathlib import Path

PROCESS_DIR = Path(__file__).resolve().parents[1] / "community" / "lrsomatic" / "1.0"
FORK = Path(os.environ.get("LRSOMATIC_FORK", "/net/bmc-lab3/data/bcc/yannvrb/GBM_TeamLab/Cirro_LRsomatic/lrsomatic-fork"))

# Keys that are mapped for Cirro but are not pipeline parameters:
#   extra_params_json - free-text JSON textarea, consumed by Cirro (as in nf-core/sarek)
#   monochrome_logs   - fixed boolean so Nextflow logs carry no ANSI codes in Cirro's log viewer
CIRRO_ONLY_KEYS = {"extra_params_json", "monochrome_logs"}

SKIP_FLAGS = [
    "skip_ascat", "skip_savana", "skip_wakhan", "skip_vep", "skip_vep_plugins",
    "skip_fiber", "skip_normalfiber", "skip_m6a", "skip_modcall", "skip_qc",
]
PARAM_PREFIX = "$.dataset.params."


def load(name):
    path = PROCESS_DIR / name
    assert path.exists(), f"missing {path}"
    return json.loads(path.read_text())


def fork_schema_properties():
    path = FORK / "nextflow_schema.json"
    assert path.exists(), f"fork schema missing: {path} (set LRSOMATIC_FORK)"
    schema = json.loads(path.read_text())
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


def test_every_mapped_param_exists_in_fork_schema():
    assert CIRRO_ONLY_KEYS == {"extra_params_json", "monochrome_logs"}
    # the fork is pinned by SHA in process-definition.json, so users do not choose a version
    assert "workflow_version" not in load("process-input.json")
    assert "workflow_version" not in load("process-form.json")["form"]["properties"]
    mapping = load("process-input.json")
    schema_props = fork_schema_properties()
    unknown = {k for k in mapping if k not in schema_props and k not in CIRRO_ONLY_KEYS}
    assert not unknown, f"mapped params not in fork schema: {sorted(unknown)}"
    for flag in SKIP_FLAGS + ["skip_signatures"]:
        assert flag in mapping, f"{flag} not mapped"


def test_every_dataset_param_path_resolves_to_a_form_property():
    mapping = load("process-input.json")
    form = load("process-form.json")
    paths = param_paths(mapping)
    assert paths, "no $.dataset.params paths mapped"
    for key, dotted in paths.items():
        form_property(form, dotted)


def test_boolean_fields_are_boolean_in_form_and_schema():
    mapping = load("process-input.json")
    form = load("process-form.json")
    schema_props = fork_schema_properties()
    for key, dotted in param_paths(mapping).items():
        form_type = form_property(form, dotted).get("type")
        if key in schema_props:
            schema_type = schema_props[key].get("type")
            assert form_type == schema_type or (form_type == "integer" and schema_type == "number"), \
                f"{key}: form type {form_type} vs schema type {schema_type}"
    for flag in SKIP_FLAGS + ["skip_signatures"]:
        assert form_property(form, param_paths(mapping)[flag]).get("type") == "boolean"
        assert schema_props[flag]["type"] == "boolean"
        assert form_property(form, param_paths(mapping)[flag]).get("default") is False
    assert mapping["monochrome_logs"] is True
    assert schema_props["monochrome_logs"]["type"] == "boolean"


def test_signatures_default_on():
    mapping = load("process-input.json")
    prop = form_property(load("process-form.json"), param_paths(mapping)["skip_signatures"])
    assert prop["default"] is False
    assert prop["title"] == "Skip mutational signatures"


def test_genome_enum_is_grch38_only():
    mapping = load("process-input.json")
    prop = form_property(load("process-form.json"), param_paths(mapping)["genome"])
    assert prop["enum"] == ["GRCh38"]
    assert prop["default"] == "GRCh38"


def test_input_path_is_config_pipeline_samplesheet():
    mapping = load("process-input.json")
    assert mapping["input"] == "$.dataset.s3|/config/pipeline_samplesheet.csv"
    assert mapping["outdir"] == "$.dataset.dataPath"


def test_public_references_are_left_to_fork_defaults():
    # iGenomes and the VEP cache are read from public S3 buckets, so they are not mapped. The fork's
    # defaults are pinned here so a fork change cannot silently point them somewhere else.
    mapping = load("process-input.json")
    assert "igenomes_base" not in mapping
    assert "vep_cache" not in mapping
    props = fork_schema_properties()
    assert props["igenomes_base"]["default"] == "s3://ngi-igenomes/igenomes/"
    assert props["vep_cache"]["default"] == "s3://annotation-cache/vep_cache/"
    assert props["vep_cache_version"]["default"] == 113


def test_staged_references_come_from_references_bucket():
    mapping = load("process-input.json")
    assert mapping["sigprofiler_genome_dir"].startswith("$.referencesBucketUri|/sigprofiler/")
    assert mapping["vep_revel"].startswith("$.referencesBucketUri|/vep_revel/")
    assert mapping["vep_revel_tbi"].startswith("$.referencesBucketUri|/vep_revel/")


def test_staged_reference_paths_match_the_upload_request():
    # These are the S3 paths requested from the references-bucket admin; change both together.
    mapping = load("process-input.json")
    assert mapping["sigprofiler_genome_dir"] == "$.referencesBucketUri|/sigprofiler/volume/"
    assert mapping["vep_revel"] == "$.referencesBucketUri|/vep_revel/revel_grch38.tsv.gz"
    assert mapping["vep_revel_tbi"] == "$.referencesBucketUri|/vep_revel/revel_grch38.tsv.gz.tbi"
    assert "TODO" not in json.dumps(mapping)


def _protected_params():
    spec = importlib.util.spec_from_file_location("preprocess", PROCESS_DIR / "preprocess.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return set(module._PROTECTED_PARAMS)


def test_every_reference_mapped_param_is_protected():
    mapping = load("process-input.json")
    expected = {k for k, v in mapping.items()
                if isinstance(v, str) and v.startswith("$.referencesBucketUri")} | {"outdir", "input"}
    missing = expected - _protected_params()
    assert not missing, f"mapped but not protected from extra_params_json: {sorted(missing)}"


def test_form_description_names_every_protected_param():
    description = load("process-form.json")["form"]["properties"]["extra_params_json"]["description"]
    # whole-word match, so e.g. `vep_revel` cannot pass because `vep_revel_tbi` is named
    missing = [k for k in sorted(_protected_params())
               if not re.search(rf"(?<![\w]){re.escape(k)}(?![\w])", description)]
    assert not missing, f"extra_params_json description does not list: {missing}"
