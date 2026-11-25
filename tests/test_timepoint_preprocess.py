import importlib.util
import re
from pathlib import Path

import pandas as pd
import pytest

PREPROCESS = Path(__file__).resolve().parents[1] / "community" / "timepoint-dmr" / "0.1" / "preprocess.py"


@pytest.fixture(scope="module")
def pre():
    assert PREPROCESS.exists(), f"missing {PREPROCESS}"
    spec = importlib.util.spec_from_file_location("timepoint_preprocess", PREPROCESS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frame(*paths):
    return pd.DataFrame({"file": list(paths)})


def pile(sample, suffix=".bed.gz", kind="tumor", prefix="s3://b/ds1/data"):
    return f"{prefix}/{sample}/methylation/{kind}/modkit_pileup/{sample}{suffix}"


def test_three_samples_give_two_comparisons(pre):
    files = frame(pile("S1"), pile("S4"), pile("S6"))
    result = pre.build_comparisons(files, "S1")
    assert list(result.columns) == pre.COMPARISON_COLUMNS
    assert result["comparison"].tolist() == ["S4_vs_S1", "S6_vs_S1"]
    assert result["treatment_name"].tolist() == ["S4", "S6"]
    assert result["control_name"].tolist() == ["S1", "S1"]
    assert result["control_bed"].tolist() == [pile("S1"), pile("S1")]
    assert result["treatment_bed"].tolist() == [pile("S4"), pile("S6")]


def test_combined_file_is_accepted(pre):
    files = frame(pile("S1", "_combined.bed.gz"), pile("S4", "_combined.bed.gz"))
    result = pre.build_comparisons(files, "S1")
    assert result["comparison"].tolist() == ["S4_vs_S1"]
    assert result["treatment_bed"].tolist() == [pile("S4", "_combined.bed.gz")]


def test_baseline_is_an_exact_match(pre):
    files = frame(pile("S1"), pile("S10"), pile("S4"))
    result = pre.build_comparisons(files, "S1")
    assert result["comparison"].tolist() == ["S10_vs_S1", "S4_vs_S1"]
    with pytest.raises(ValueError) as err:
        pre.build_comparisons(files, "S")
    message = str(err.value)
    assert "not found" in message
    assert all(name in message for name in ("S1", "S10", "S4"))


def test_same_sample_in_two_datasets_is_an_error(pre):
    first, second = pile("S4", prefix="s3://b/ds1/data"), pile("S4", prefix="s3://b/ds2/data")
    with pytest.raises(ValueError) as err:
        pre.build_comparisons(frame(pile("S1"), first, second), "S1")
    assert first in str(err.value) and second in str(err.value)


def test_hp_and_normal_files_are_never_used(pre):
    extras = [pile("S4", "_hp1.bed.gz"), pile("S4", "_hp2.bed.gz"), pile("S4", kind="normal")]
    files = frame(pile("S1"), pile("S4"), *extras)
    result = pre.build_comparisons(files, "S1")
    assert result["treatment_bed"].tolist() == [pile("S4")]
    pileups, ignored = pre.tumor_pileups(files)
    assert pileups == {"S1": [pile("S1")], "S4": [pile("S4")]}
    for path in extras:
        assert path in ignored


def test_plain_and_combined_for_one_sample_is_an_error(pre):
    files = frame(pile("S1"), pile("S4"), pile("S4", "_combined.bed.gz"))
    with pytest.raises(ValueError) as err:
        pre.build_comparisons(files, "S1")
    assert pile("S4") in str(err.value) and pile("S4", "_combined.bed.gz") in str(err.value)


def test_missing_baseline_blank_baseline_no_files_no_other_sample(pre):
    with pytest.raises(ValueError, match="not found"):
        pre.build_comparisons(frame(pile("S1"), pile("S4")), "S9")
    with pytest.raises(ValueError, match="blank"):
        pre.build_comparisons(frame(pile("S1"), pile("S4")), "  ")
    with pytest.raises(ValueError, match=r"no .*pileup"):
        pre.build_comparisons(frame("s3://b/ds1/data/S1/variants/x.vcf.gz"), "S1")
    with pytest.raises(ValueError, match="no other tumor"):
        pre.build_comparisons(frame(pile("S1")), "S1")


def test_all_problems_are_listed_together(pre):
    with pytest.raises(ValueError) as err:
        pre.build_comparisons(frame(), "")
    message = str(err.value)
    assert message.startswith("Invalid timepoint DMR input:\n- ")
    assert "blank" in message and re.search(r"no .*pileup", message)


def test_relative_paths_work(pre):
    rel = lambda s: f"data/{s}/methylation/tumor/modkit_pileup/{s}.bed.gz"
    result = pre.build_comparisons(frame(rel("S1"), rel("S4")), "S1")
    assert result["comparison"].tolist() == ["S4_vs_S1"]
    assert result["control_bed"].tolist() == [rel("S1")]
