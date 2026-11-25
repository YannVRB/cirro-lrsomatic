#!/usr/bin/env python3
"""Cirro preprocess step for the timepoint DMR process.

Turns the files of the selected lrsomatic-dev result datasets into the comparison table
comparison,treatment_name,treatment_bed,control_name,control_bed: every tumor sample except the named
baseline is compared with the baseline (treatment = that sample, control = baseline).

`cirro` is imported only inside main() so that the pure functions can be tested without Cirro installed.
"""

import re

import pandas as pd

COMPARISON_COLUMNS = ["comparison", "treatment_name", "treatment_bed", "control_name", "control_bed"]

# The tumor's bedMethyl in an lrsomatic-dev result: data/<sample>/methylation/tumor/modkit_pileup/<stem>.bed.gz.
# Only a stem equal to <sample> (no phasing) or <sample>_combined (with "Phased methylation pileup") is a
# whole-sample file; <sample>_hp1 and <sample>_hp2 are single haplotypes. The normal's files (.../normal/...)
# are named after the tumor and never match.
_PILEUP = re.compile(r"(?:^|/)data/(?P<sample>[^/]+)/methylation/tumor/modkit_pileup/(?P<stem>[^/]+)\.bed\.gz$")


def _text(value: object) -> str:
    """Stripped string, or '' for None/NaN."""
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def tumor_pileups(files: pd.DataFrame) -> tuple[dict[str, list[str]], list[str]]:
    """Whole-sample tumor pileups by sample name, and every other path (ignored)."""
    files = files.rename(columns={c: str(c).strip().lower() for c in files.columns})
    if "file" not in files.columns:
        raise ValueError(f"files table is missing column 'file'; found {list(files.columns)}")
    pileups: dict[str, list[str]] = {}
    ignored: list[str] = []
    for path in (_text(p) for p in files["file"]):
        match = _PILEUP.search(path)
        if match and match["stem"] in (match["sample"], match["sample"] + "_combined"):
            pileups.setdefault(match["sample"], []).append(path)
        else:
            ignored.append(path)
    return pileups, ignored


def build_comparisons(files: pd.DataFrame, baseline: str) -> pd.DataFrame:
    """One row per tumor sample other than `baseline`, sorted by sample name.

    Raises one ValueError listing every problem found.
    """
    baseline = _text(baseline)
    pileups, _ = tumor_pileups(files)
    problems = []

    if baseline == "":
        problems.append("the baseline sample name is blank")
    if not pileups:
        problems.append(
            "no tumor methylation pileup files in the input; the datasets are probably not lrsomatic-dev "
            "results, or were run without 'Phased methylation pileup'"
        )
    for sample, paths in sorted(pileups.items()):
        if len(paths) > 1:
            problems.append(f"sample {sample!r} has {len(paths)} candidate pileup files, expected 1: {sorted(paths)}")
    if baseline and pileups and baseline not in pileups:
        problems.append(f"baseline sample {baseline!r} not found; samples found: {sorted(pileups)}")
    others = sorted(sample for sample in pileups if sample != baseline)
    if baseline in pileups and not others:
        problems.append(f"no other tumor sample besides the baseline {baseline!r}; nothing to compare")

    if problems:
        raise ValueError("Invalid timepoint DMR input:\n- " + "\n- ".join(problems))

    control = pileups[baseline][0]
    rows = [
        {
            "comparison": f"{sample}_vs_{baseline}",
            "treatment_name": sample,
            "treatment_bed": pileups[sample][0],
            "control_name": baseline,
            "control_bed": control,
        }
        for sample in others
    ]
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)
