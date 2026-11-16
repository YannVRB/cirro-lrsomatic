import json
from pathlib import Path

PROCESS_DIR = Path(__file__).resolve().parents[1] / "community" / "lrsomatic" / "1.0"
FORK_SHA = "42458112ee927a37c1879a3f23e0fb82ab7b8fdb"


def load(name):
    path = PROCESS_DIR / name
    assert path.exists(), f"missing {path}"
    return json.loads(path.read_text())


def readme():
    path = PROCESS_DIR / "README.md"
    assert path.exists(), f"missing {path}"
    return path.read_text()


def version_tuple(text):
    return tuple(int(part) for part in text.split("."))


def test_code_pins_fork():
    code = load("process-definition.json")["code"]
    assert code["uri"] == "yannvrb/lrsomatic"
    assert code["version"] == FORK_SHA
    assert code["repository"] == "GITHUBPUBLIC"
    assert code["script"] == "main.nf"
    assert version_tuple(code["executorVersion"]) >= (25, 10, 4)


def test_allow_multiple_sources_true():
    assert load("process-definition.json")["allowMultipleSources"] is True


def test_uses_sample_sheet_true():
    assert load("process-definition.json")["usesSampleSheet"] is True


def test_parents_are_the_two_bam_intake_types():
    # custom_dataset is excluded on purpose: it carries no sample association, so pairing cannot work
    parents = load("process-definition.json")["parentProcessIds"]
    assert parents == ["ingest_aligned_bam", "ingest_unaligned_bam"]


def test_output_json_is_empty_commands():
    assert load("process-output.json") == {"commands": []}


def test_readme_lists_every_staged_reference_path():
    text = readme()
    paths = [v for v in load("process-input.json").values()
             if isinstance(v, str) and v.startswith("s3://btc-refdata/")]
    assert paths, "process-input.json maps no btc-refdata path"
    for path in paths:
        assert path in text, f"README does not list {path}"


def test_readme_documents_metadata_columns():
    text = readme()
    for column in ["patient", "status", "sex", "platform", "fiber"]:
        assert column in text, f"README does not mention {column}"


def test_vep_revel_is_prepared_gz_with_matching_tbi():
    mapping = load("process-input.json")
    revel = mapping["vep_revel"]
    assert revel.lower().endswith(".tsv.gz")
    assert mapping.get("vep_revel_tbi") == revel + ".tbi"
