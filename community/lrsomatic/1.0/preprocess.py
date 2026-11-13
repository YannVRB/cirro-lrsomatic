#!/usr/bin/env python3
"""Cirro preprocess step for LRsomatic.

Turns a Cirro dataset (BAM file list + per-sample metadata) into the pipeline
samplesheet: sample,bam_tumor,bam_normal,platform,sex,fiber (one row per tumor).

`cirro` is imported only inside main() so that the pure function
build_pipeline_samplesheet can be tested without Cirro installed.
"""

import json
import re

import pandas as pd

# Accepted spellings per metadata field -> canonical value. Keys are matched
# after strip() + lower(); anything else is an error, never a default.
SPELLINGS = {
    "status": {
        "tumor": "tumor", "tumour": "tumor",
        "normal": "normal", "pbmc": "normal",
    },
    "sex": {
        "male": "male", "m": "male",
        "female": "female", "f": "female",
    },
    "platform": {
        "ont": "ont",
        "pb": "pb", "pacbio": "pb",
    },
    "fiber": {
        "y": "y", "yes": "y",
        "n": "n", "no": "n",
    },
}

OUTPUT_COLUMNS = ["sample", "bam_tumor", "bam_normal", "platform", "sex", "fiber"]
METADATA_COLUMNS = ["patient", "status", "sex", "platform", "fiber"]


def _text(value: object) -> str:
    """Stripped string, or '' for None/NaN."""
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def normalize(field: str, value: object, sample: str) -> str:
    """Canonical value of `field` for `value`, or ValueError naming sample and field."""
    key = _text(value).lower()
    try:
        return SPELLINGS[field][key]
    except KeyError:
        accepted = ", ".join(sorted(SPELLINGS[field]))
        raise ValueError(
            f"sample {sample!r}: invalid {field} {_text(value)!r} (accepted: {accepted})"
        ) from None


def _lower_columns(df: pd.DataFrame, table: str) -> pd.DataFrame:
    df = df.rename(columns={c: str(c).strip().lower() for c in df.columns})
    duplicated = sorted(set(df.columns[df.columns.duplicated()]))
    if duplicated:
        raise ValueError(
            f"{table} table has duplicate column(s) {duplicated} "
            f"(column names are matched case-insensitively)"
        )
    return df


def _is_bam(files: pd.DataFrame) -> pd.Series:
    """Boolean mask of the rows whose `file` ends in `.bam` (case-sensitive, as the pipeline expects)."""
    return files["file"].map(_text).str.endswith(".bam")


def ignored_files(files: pd.DataFrame) -> list[str]:
    """The `file` values that build_pipeline_samplesheet ignores because they do not end in `.bam`."""
    files = _lower_columns(files, "files")
    if "file" not in files.columns:
        return []
    return [_text(p) for p in files.loc[~_is_bam(files), "file"]]


def normal_only_patients(files: pd.DataFrame, samples: pd.DataFrame) -> list[str]:
    """Sorted patients whose samples (with a BAM and a valid status) are all normal.

    build_pipeline_samplesheet gives such a patient no rows; main() warns about each one.
    """
    files = _lower_columns(files, "files")
    samples = _lower_columns(samples, "metadata")
    with_bam = {_text(s) for s in files.loc[_is_bam(files), "sample"]}
    statuses = {}
    for row in samples.to_dict("records"):
        patient = _text(row["patient"])
        if patient == "" or _text(row["sample"]) not in with_bam:
            continue
        try:
            statuses.setdefault(patient, set()).add(normalize("status", row["status"], ""))
        except ValueError:
            continue
    return sorted(p for p, found in statuses.items() if found == {"normal"})


def build_pipeline_samplesheet(files: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Pair each tumor BAM with its patient's normal BAM.

    files:   columns `sample`, `file` (optionally `dataset`, ignored).
    samples: column `sample` plus patient, status, sex, platform, fiber (any case).
    Raises one ValueError listing every problem found.
    """
    files = _lower_columns(files, "files")
    samples = _lower_columns(samples, "metadata")

    missing = [c for c in ("sample", "file") if c not in files.columns]
    if missing:
        raise ValueError(
            f"files table is missing column(s) {missing}; found {list(files.columns)}"
        )
    missing = [c for c in ["sample"] + METADATA_COLUMNS if c not in samples.columns]
    if missing:
        raise ValueError(
            f"sample metadata is missing column(s) {missing}; found {list(samples.columns)}"
        )

    problems = []

    # Only `.bam` rows are BAMs; other rows (e.g. `<sample>.bam.bai`) are ignored, see ignored_files()
    files = files[_is_bam(files)]
    if files.empty:
        problems.append("no .bam files in the input")

    # One BAM per sample, every BAM has metadata and vice versa
    bams = {}
    for name, group in files.groupby(files["sample"].map(_text), sort=True):
        paths = [_text(p) for p in group["file"]]
        if name == "":
            problems.append(f"{len(paths)} BAM file(s) have a blank sample name: {paths}")
        elif len(paths) > 1:
            problems.append(f"sample {name!r} has {len(paths)} BAM files, expected 1: {paths}")
        else:
            bams[name] = paths[0]

    meta_names = [_text(s) for s in samples["sample"]]
    for name in sorted({n for n in meta_names if meta_names.count(n) > 1}):
        problems.append(f"sample {name!r} appears more than once in the metadata")
    for name in sorted(set(meta_names) - set(bams) - {""}):
        if name not in set(files["sample"].map(_text)):
            problems.append(f"sample {name!r} is in the metadata but has no BAM file")
    for name in sorted(set(bams) - set(meta_names)):
        problems.append(f"sample {name!r} has a BAM file but no metadata row")
    if "" in meta_names:
        problems.append("a metadata row has a blank sample name")

    # Per-sample validation
    records = []
    for row in samples.to_dict("records"):
        name = _text(row["sample"])
        if name == "":
            continue
        rec = {"sample": name, "bam": bams.get(name), "ok": True}
        if not re.fullmatch(r"\S+", name):
            problems.append(f"sample name {name!r} must not contain whitespace")
        patient = _text(row["patient"])
        if patient == "":
            problems.append(f"sample {name!r}: patient is blank")
            rec["ok"] = False
        rec["patient"] = patient
        for field in ("status", "sex", "platform", "fiber"):
            try:
                rec[field] = normalize(field, row[field], name)
            except ValueError as err:
                problems.append(str(err))
                rec["ok"] = False
        records.append(rec)

    # Pair tumors with their patient's normal
    rows = []
    valid = pd.DataFrame([r for r in records if r["ok"]])
    if not valid.empty:
        for patient, group in valid.groupby("patient"):
            tumors = group[group["status"] == "tumor"].to_dict("records")
            normals = group[group["status"] == "normal"].to_dict("records")
            if len(normals) != 1:
                names = [n["sample"] for n in normals]
                problems.append(
                    f"patient {patient!r} has {len(normals)} normal sample(s) {names}, expected exactly 1"
                )
                continue
            normal = normals[0]
            for tumor in tumors:
                # platform and sex must agree within a patient (a sex mismatch is the one free sign of
                # a sample filed under the wrong patient); fiber may differ (skip_normalfiber), the
                # tumor's value is used
                mismatch = False
                for field in ("platform", "sex"):
                    if tumor[field] != normal[field]:
                        problems.append(
                            f"patient {patient!r}: {field} mismatch, tumor {tumor['sample']!r} is "
                            f"{tumor[field]} but normal {normal['sample']!r} is {normal[field]}"
                        )
                        mismatch = True
                if mismatch:
                    continue
                rows.append({
                    "sample": tumor["sample"],
                    "bam_tumor": tumor["bam"],
                    "bam_normal": normal["bam"],
                    "platform": tumor["platform"],
                    "sex": tumor["sex"],
                    "fiber": tumor["fiber"],
                })
        if not rows and not problems:
            problems.append("no tumor samples found; nothing to run")
    elif not problems:
        problems.append("no tumor samples found; nothing to run")

    if problems:
        raise ValueError("Invalid samplesheet input:\n- " + "\n- ".join(problems))

    return (
        pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
        .sort_values("sample", kind="stable")
        .reset_index(drop=True)
    )


# Params that extra_params_json must never override: those set by process-input.json (the
# samplesheet/output locations, the staged reference paths, monochrome_logs), the form's genome
# (GRCh38 only in v1; the form field is the only way to set it), the remote config include
# (custom_config_base, custom_config_version), parameter validation, igenomes_ignore, and the
# downloads of references that are staged instead.
_PROTECTED_PARAMS = frozenset({
    "input", "outdir", "igenomes_base", "vep_cache",
    "sigprofiler_genome_dir", "vep_revel", "vep_revel_tbi", "monochrome_logs",
    "custom_config_base", "custom_config_version", "validate_params", "igenomes_ignore", "genome",
    "download_vep_cache", "download_sigprofiler_genome",
})
# The form's maxLength is client-side only, so it is enforced here too.
_MAX_EXTRA_PARAMS_CHARS = 10_000


def extra_params(raw: object) -> dict[str, object]:
    """Parse the extra_params_json textarea: a JSON object string -> dict, blank -> {}."""
    text = _text(raw)
    if not text:
        return {}
    if len(text) > _MAX_EXTRA_PARAMS_CHARS:
        raise ValueError(
            f"extra_params_json is {len(text)} characters, over the {_MAX_EXTRA_PARAMS_CHARS}-character limit"
        )
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"extra_params_json is not valid JSON: {exc}") from exc
    if not isinstance(parsed, dict):
        raise ValueError(f"extra_params_json must be a JSON object, got {type(parsed).__name__}")
    protected = sorted(_PROTECTED_PARAMS & parsed.keys())
    if protected:
        raise ValueError(f"extra_params_json cannot set protected parameter(s): {', '.join(protected)}")
    return parsed


def main() -> None:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    sheet = build_pipeline_samplesheet(ds.files, ds.samplesheet)
    ignored = ignored_files(ds.files)
    if ignored:
        ds.logger.warning(f"Ignored {len(ignored)} file(s) that do not end in .bam: {ignored}")
    ds.logger.info(f"Pipeline samplesheet: {sheet.shape[0]} tumor sample(s)\n{sheet.to_csv(index=False)}")
    for patient in normal_only_patients(ds.files, ds.samplesheet):
        ds.logger.warning(f"Patient {patient!r} has only normal sample(s); skipped (no tumor to analyse)")
    # Must go to the path mapped in process-input.json (an S3 URI), not the working directory
    sheet.to_csv(ds.params["input"], index=None)

    # User-supplied extra params override earlier values (as in nf-core/sarek 3.1)
    for key, value in extra_params(ds.params.get("extra_params_json")).items():
        ds.add_param(key, value, overwrite=True)
    # Cirro-only field, never a Nextflow parameter
    ds.remove_param("extra_params_json", force=True)


if __name__ == "__main__":
    main()
