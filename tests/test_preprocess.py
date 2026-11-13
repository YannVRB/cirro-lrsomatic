import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PREPROCESS = Path(__file__).resolve().parents[1] / "community" / "lrsomatic" / "1.0" / "preprocess.py"


def load_module():
    spec = importlib.util.spec_from_file_location("preprocess", PREPROCESS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module



@pytest.fixture(scope="module")
def pre():
    assert PREPROCESS.exists(), f"missing {PREPROCESS}"
    return load_module()


COLUMNS = ["sample", "bam_tumor", "bam_normal", "platform", "sex", "fiber"]


def uri(name, ds="ds1"):
    return f"s3://proj/{ds}/{name}.bam"


def make(rows, columns=("sample", "patient", "status", "sex", "platform", "fiber"), datasets=None):
    """rows: (sample, patient, status, sex, platform, fiber) -> (files, samples)."""
    samples = pd.DataFrame(rows, columns=list(columns))
    key = samples.columns[0]
    files = pd.DataFrame({"sample": samples[key], "file": [uri(s) for s in samples[key]]})
    if datasets is not None:
        files["dataset"] = datasets
    return files, samples


P3 = [
    ("PT3.S1.L2", "P3", "tumor", "female", "ont", "n"),
    ("PT3.S6.L1", "P3", "tumor", "female", "ont", "n"),
    ("P3_PBMC", "P3", "normal", "female", "ont", "n"),
]
P4 = [
    ("PT4.S1.L1", "P4", "tumor", "male", "ont", "n"),
    ("PT4.S4.L2", "P4", "tumor", "male", "ont", "n"),
    ("PT4.S6.L6", "P4", "tumor", "male", "ont", "n"),
    ("P4_PBMC", "P4", "normal", "male", "ont", "n"),
]


def test_p3_layout_reproduces_samplesheet(pre):
    files, samples = make(P3)
    out = pre.build_pipeline_samplesheet(files, samples)
    assert list(out.columns) == COLUMNS
    assert len(out) == 2
    assert out["sample"].tolist() == ["PT3.S1.L2", "PT3.S6.L1"]
    assert (out["bam_normal"] == uri("P3_PBMC")).all()
    assert out["bam_tumor"].tolist() == [uri("PT3.S1.L2"), uri("PT3.S6.L1")]
    assert (out["sex"] == "female").all()
    assert (out["platform"] == "ont").all()
    assert (out["fiber"] == "n").all()


def test_p4_layout_three_tumors_one_normal(pre):
    files, samples = make(P4)
    out = pre.build_pipeline_samplesheet(files, samples)
    assert len(out) == 3
    assert (out["bam_normal"] == uri("P4_PBMC")).all()
    assert (out["sex"] == "male").all()


def test_two_patients_never_cross_paired(pre):
    files, samples = make(P4 + P3)
    out = pre.build_pipeline_samplesheet(files, samples)
    assert out["sample"].tolist() == sorted(out["sample"].tolist())
    normals = dict(zip(out["sample"], out["bam_normal"]))
    for s in ("PT3.S1.L2", "PT3.S6.L1"):
        assert normals[s] == uri("P3_PBMC")
    for s in ("PT4.S1.L1", "PT4.S4.L2", "PT4.S6.L6"):
        assert normals[s] == uri("P4_PBMC")


def test_tumor_and_normal_in_different_datasets(pre):
    files, samples = make(P3)
    files["dataset"] = ["ds1", "ds1", "ds2"]
    files["file"] = [uri("PT3.S1.L2", "ds1"), uri("PT3.S6.L1", "ds1"), uri("P3_PBMC", "ds2")]
    out = pre.build_pipeline_samplesheet(files, samples)
    assert len(out) == 2
    assert (out["bam_normal"] == uri("P3_PBMC", "ds2")).all()


def test_patient_without_normal_raises(pre):
    files, samples = make(P3[:2])
    with pytest.raises(ValueError, match="P3"):
        pre.build_pipeline_samplesheet(files, samples)


def test_patient_with_two_normals_raises(pre):
    files, samples = make(P3 + [("P3_PBMC2", "P3", "normal", "female", "ont", "n")])
    with pytest.raises(ValueError, match="P3"):
        pre.build_pipeline_samplesheet(files, samples)


def test_sample_with_two_bams_raises(pre):
    files, samples = make(P3)
    extra = pd.DataFrame({"sample": ["PT3.S1.L2"], "file": [uri("PT3.S1.L2.again")]})
    files = pd.concat([files, extra], ignore_index=True)
    with pytest.raises(ValueError, match="PT3.S1.L2"):
        pre.build_pipeline_samplesheet(files, samples)


def test_metadata_sample_without_bam_raises(pre):
    files, samples = make(P3)
    files = files[files["sample"] != "PT3.S6.L1"]
    with pytest.raises(ValueError, match="PT3.S6.L1"):
        pre.build_pipeline_samplesheet(files, samples)


def test_bam_without_metadata_raises(pre):
    files, samples = make(P3)
    samples = samples[samples["sample"] != "PT3.S6.L1"]
    with pytest.raises(ValueError, match="PT3.S6.L1"):
        pre.build_pipeline_samplesheet(files, samples)


def test_files_without_sample_column_raises(pre):
    files, samples = make(P3)
    files = files.rename(columns={"sample": "name"})
    with pytest.raises(ValueError, match="sample"):
        pre.build_pipeline_samplesheet(files, samples)


def test_missing_metadata_column_raises(pre):
    files, samples = make(P3)
    samples = samples.drop(columns=["fiber"])
    with pytest.raises(ValueError, match="fiber"):
        pre.build_pipeline_samplesheet(files, samples)


def test_bad_sex_names_sample(pre):
    rows = [("PT3.S1.L2", "P3", "tumor", "unknown", "ont", "n")] + P3[1:]
    files, samples = make(rows)
    with pytest.raises(ValueError, match="PT3.S1.L2") as exc:
        pre.build_pipeline_samplesheet(files, samples)
    assert "sex" in str(exc.value)


def test_bad_platform_raises(pre):
    rows = [("PT3.S1.L2", "P3", "tumor", "female", "illumina", "n")] + P3[1:]
    files, samples = make(rows)
    with pytest.raises(ValueError, match="platform"):
        pre.build_pipeline_samplesheet(files, samples)


def test_platform_mismatch_within_patient_raises(pre):
    rows = [("PT3.S1.L2", "P3", "tumor", "female", "pb", "n")] + P3[1:]
    files, samples = make(rows)
    with pytest.raises(ValueError, match="P3") as exc:
        pre.build_pipeline_samplesheet(files, samples)
    assert "platform" in str(exc.value)


def test_sex_mismatch_within_patient_raises(pre):
    rows = [P3[0], P3[1], ("P3_PBMC", "P3", "normal", "male", "ont", "n")]
    files, samples = make(rows)
    with pytest.raises(ValueError, match="P3") as exc:
        pre.build_pipeline_samplesheet(files, samples)
    message = str(exc.value)
    assert "sex" in message
    for token in ("PT3.S1.L2", "PT3.S6.L1", "P3_PBMC", "female", "male"):
        assert token in message, f"{token!r} not in: {message}"


def test_sex_mismatch_is_reported_with_a_platform_mismatch(pre):
    rows = [("PT3.S1.L2", "P3", "tumor", "female", "pb", "n"), P3[1],
            ("P3_PBMC", "P3", "normal", "male", "ont", "n")]
    files, samples = make(rows)
    with pytest.raises(ValueError) as exc:
        pre.build_pipeline_samplesheet(files, samples)
    message = str(exc.value)
    assert "platform mismatch" in message and "sex mismatch" in message


def test_fiber_difference_is_not_an_error_and_tumor_fiber_is_used(pre):
    rows = [("PT3.S1.L2", "P3", "tumor", "female", "ont", "y"),
            ("P3_PBMC", "P3", "normal", "female", "ont", "n")]
    files, samples = make(rows)
    out = pre.build_pipeline_samplesheet(files, samples)
    assert out["fiber"].tolist() == ["y"]
    rows = [("PT3.S1.L2", "P3", "tumor", "female", "ont", "n"),
            ("P3_PBMC", "P3", "normal", "female", "ont", "y")]
    files, samples = make(rows)
    assert pre.build_pipeline_samplesheet(files, samples)["fiber"].tolist() == ["n"]


def test_normal_only_patients_lists_patients_without_a_tumor(pre):
    rows = P3 + [("P5_PBMC", "P5", "normal", "male", "ont", "n"),
                 ("P2_PBMC", "P2", "pbmc", "male", "ont", "n")]
    files, samples = make(rows)
    assert pre.normal_only_patients(files, samples) == ["P2", "P5"]


def test_normal_only_patients_is_empty_when_every_patient_has_a_tumor(pre):
    files, samples = make(P3 + P4)
    assert pre.normal_only_patients(files, samples) == []


def test_normal_only_patients_skips_patients_without_a_valid_status(pre):
    rows = P3 + [("P5_PBMC", "P5", "unknown", "male", "ont", "n"),
                 ("P6_X", "P6", "", "male", "ont", "n")]
    files, samples = make(rows)
    assert pre.normal_only_patients(files, samples) == []


def test_normal_only_patients_ignores_samples_without_a_bam(pre):
    rows = P3 + [("P5_PBMC", "P5", "normal", "male", "ont", "n")]
    files, samples = make(rows)
    files = files[files["sample"] != "P5_PBMC"]
    assert pre.normal_only_patients(files, samples) == []


def test_all_problems_reported_together(pre):
    rows = [
        ("PT3.S1.L2", "P3", "tumor", "unknown", "ont", "n"),
        ("PT3.S6.L1", "P3", "tumor", "female", "illumina", "n"),
        ("P3_PBMC", "P3", "normal", "female", "ont", "n"),
    ]
    files, samples = make(rows)
    with pytest.raises(ValueError) as exc:
        pre.build_pipeline_samplesheet(files, samples)
    msg = str(exc.value)
    assert "PT3.S1.L2" in msg and "PT3.S6.L1" in msg


def test_status_spellings(pre):
    assert pre.normalize("status", "Tumour", "s") == "tumor"
    assert pre.normalize("status", " tumor ", "s") == "tumor"
    assert pre.normalize("status", "Tumor", "s") == "tumor"
    assert pre.normalize("status", "Normal", "s") == "normal"
    assert pre.normalize("status", "PBMC", "s") == "normal"
    with pytest.raises(ValueError, match="s1.*status|status.*s1"):
        pre.normalize("status", "blood", "s1")


def test_sex_spellings(pre):
    assert [pre.normalize("sex", v, "s") for v in ("M", "f", "Male", "FEMALE")] == [
        "male", "female", "male", "female"]


def test_pacbio_spelling(pre):
    assert pre.normalize("platform", "Pacbio", "s") == "pb"
    assert pre.normalize("platform", "pb", "s") == "pb"
    assert pre.normalize("platform", "ONT", "s") == "ont"


def test_fiber_yes_no(pre):
    assert pre.normalize("fiber", "Yes", "s") == "y"
    assert pre.normalize("fiber", "y", "s") == "y"
    assert pre.normalize("fiber", "No", "s") == "n"
    assert pre.normalize("fiber", " N ", "s") == "n"


@pytest.mark.parametrize("value", [None, float("nan"), "", "   "])
def test_blank_and_nan_values_invalid(pre, value):
    with pytest.raises(ValueError, match="sex"):
        pre.normalize("sex", value, "s1")


def test_column_names_case_insensitive(pre):
    files, samples = make(P3, columns=("Sample", "Patient", "STATUS", "Sex", "PLATFORM", "Fiber"))
    out = pre.build_pipeline_samplesheet(files, samples)
    assert len(out) == 2


def test_spellings_normalised_end_to_end(pre):
    rows = [
        ("PT3.S1.L2", "P3", " Tumour ", "F", "Pacbio", "Yes"),
        ("P3_PBMC", "P3", "PBMC", "F", "pb", "no"),
    ]
    files, samples = make(rows)
    out = pre.build_pipeline_samplesheet(files, samples)
    assert out.iloc[0].tolist() == [
        "PT3.S1.L2", uri("PT3.S1.L2"), uri("P3_PBMC"), "pb", "female", "y"]


def test_blank_patient_raises(pre):
    rows = [("PT3.S1.L2", "", "tumor", "female", "ont", "n")] + P3[1:]
    files, samples = make(rows)
    with pytest.raises(ValueError, match="patient") as exc:
        pre.build_pipeline_samplesheet(files, samples)
    assert "PT3.S1.L2" in str(exc.value)


def test_nan_patient_raises(pre):
    rows = [("PT3.S1.L2", None, "tumor", "female", "ont", "n")] + P3[1:]
    files, samples = make(rows)
    with pytest.raises(ValueError, match="patient"):
        pre.build_pipeline_samplesheet(files, samples)


def test_sample_name_with_space_raises(pre):
    rows = [("PT3 S1", "P3", "tumor", "female", "ont", "n")] + P3[1:]
    files, samples = make(rows)
    with pytest.raises(ValueError, match="PT3 S1"):
        pre.build_pipeline_samplesheet(files, samples)


def test_sample_whose_only_file_is_not_a_bam_raises(pre):
    # a non-.bam row is ignored, so this sample ends up with no BAM and the error names it
    files, samples = make(P3)
    files.loc[0, "file"] = "s3://proj/ds1/PT3.S1.L2.bam.bai"
    with pytest.raises(ValueError, match="PT3.S1.L2") as exc:
        pre.build_pipeline_samplesheet(files, samples)
    assert "no BAM file" in str(exc.value)


def with_bai(files):
    """Add a `<bam>.bai` row beside every `.bam` row."""
    bai = files.assign(file=files["file"] + ".bai")
    return pd.concat([files, bai], ignore_index=True)


def test_bai_rows_beside_each_bam_are_ignored(pre):
    files, samples = make(P3 + P4)
    expected = pre.build_pipeline_samplesheet(files, samples)
    out = pre.build_pipeline_samplesheet(with_bai(files), samples)
    pd.testing.assert_frame_equal(out, expected)


def test_ignored_files_lists_the_non_bam_rows(pre):
    files, _ = make(P3)
    mixed = with_bai(files)
    assert pre.ignored_files(mixed) == [f"{p}.bai" for p in files["file"]]
    assert pre.ignored_files(files) == []


def test_bam_extension_is_case_sensitive(pre):
    files, samples = make(P3)
    files.loc[0, "file"] = "s3://proj/ds1/PT3.S1.L2.BAM"
    assert pre.ignored_files(files) == ["s3://proj/ds1/PT3.S1.L2.BAM"]
    with pytest.raises(ValueError, match="PT3.S1.L2"):
        pre.build_pipeline_samplesheet(files, samples)


def test_no_bam_files_at_all_raises(pre):
    files, samples = make(P3)
    files["file"] = files["file"] + ".bai"
    with pytest.raises(ValueError, match="no .bam files"):
        pre.build_pipeline_samplesheet(files, samples)


def test_no_tumor_rows_raises(pre):
    files, samples = make([P3[2]])
    with pytest.raises(ValueError, match="no tumor|No tumor|tumor"):
        pre.build_pipeline_samplesheet(files, samples)


def test_module_imports_without_cirro():
    code = (
        "import importlib.util, sys; "
        f"s = importlib.util.spec_from_file_location('preprocess', {str(PREPROCESS)!r}); "
        "m = importlib.util.module_from_spec(s); s.loader.exec_module(m); "
        "assert 'cirro' not in sys.modules and 'boto3' not in sys.modules"
    )
    result = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_duplicate_metadata_columns_differing_in_case_raise(pre):
    files, samples = make(P3)
    samples = pd.concat([samples, samples[["sex"]].rename(columns={"sex": "Sex"})], axis=1)
    with pytest.raises(ValueError, match="sex") as exc:
        pre.build_pipeline_samplesheet(files, samples)
    assert "metadata" in str(exc.value)


def test_duplicate_files_columns_differing_in_case_raise(pre):
    files, samples = make(P3)
    files = pd.concat([files, files[["sample"]].rename(columns={"sample": "Sample"})], axis=1)
    with pytest.raises(ValueError, match="sample") as exc:
        pre.build_pipeline_samplesheet(files, samples)
    assert "files" in str(exc.value)


def test_extra_params_parses_json_object(pre):
    raw = '{"skip_cramino": true, "vep_cache_version": 113}'
    assert pre.extra_params(raw) == {"skip_cramino": True, "vep_cache_version": 113}


@pytest.mark.parametrize("raw", [None, "", "   ", float("nan")])
def test_extra_params_empty_is_empty_dict(pre, raw):
    assert pre.extra_params(raw) == {}


def test_extra_params_invalid_json_raises(pre):
    with pytest.raises(ValueError, match="extra_params_json"):
        pre.extra_params("{not json")


@pytest.mark.parametrize("raw", ["[1, 2]", '"x"'])
def test_extra_params_non_object_raises(pre, raw):
    with pytest.raises(ValueError, match="extra_params_json"):
        pre.extra_params(raw)


PROTECTED = [
    "input", "outdir", "igenomes_base", "vep_cache", "sigprofiler_genome_dir", "vep_revel", "vep_revel_tbi",
    "monochrome_logs",
    "custom_config_base", "custom_config_version", "validate_params", "igenomes_ignore", "genome",
    "download_vep_cache", "download_sigprofiler_genome",
]


@pytest.mark.parametrize("key", PROTECTED)
def test_extra_params_protected_keys_raise(pre, key):
    with pytest.raises(ValueError, match="extra_params_json") as exc:
        pre.extra_params(f'{{"{key}": "x"}}')
    assert key in str(exc.value)


def test_extra_params_lists_all_protected_keys_in_one_error(pre):
    with pytest.raises(ValueError, match="extra_params_json") as exc:
        pre.extra_params('{"outdir": "s3://elsewhere", "skip_vep": true, "input": "x"}')
    assert "outdir" in str(exc.value) and "input" in str(exc.value)
    assert "skip_vep" not in str(exc.value)


def _object_of_length(n):
    prefix, suffix = '{"skip_vep": "', '"}'
    return prefix + "a" * (n - len(prefix) - len(suffix)) + suffix


def test_extra_params_oversize_raises(pre):
    raw = _object_of_length(10001)
    assert len(raw) == 10001
    with pytest.raises(ValueError, match="10000|10,000") as exc:
        pre.extra_params(raw)
    assert "extra_params_json" in str(exc.value)


def test_extra_params_at_size_limit_is_accepted(pre):
    raw = _object_of_length(10000)
    assert len(raw) == 10000
    assert pre.extra_params(raw) == {"skip_vep": "a" * (10000 - len('{"skip_vep": "') - 2)}
