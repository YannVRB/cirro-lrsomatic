"""Check 2: rebuild the real P3 and P4 samplesheets with build_pipeline_samplesheet.

Not a pytest test. Run: python -I docs/verification/verify_real_samplesheets.py

Reads the two real samplesheets (read-only), turns them into a Cirro-style `files` table and a
metadata table (tumors and the shared normal in different datasets), calls
build_pipeline_samplesheet once for both patients, and compares the result row by row.
This proves the pairing logic on the real P3/P4 layout. It does not prove the real shape of
Cirro's ds.files or of the Cirro metadata.
"""
import importlib.util
import posixpath
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PREPROCESS = ROOT / "community" / "lrsomatic" / "1.0" / "preprocess.py"
BASE = Path("/net/bmc-lab3/data/bcc/yannvrb/GBM_TeamLab")
REAL = {"P3": BASE / "P3.LRsomatic" / "samplesheet.csv", "P4": BASE / "P4.LRsomatic" / "samplesheet.csv"}
FIELDS = ["platform", "sex", "fiber"]


def uri(dataset, local_path):
    return f"s3://proj/{dataset}/{posixpath.basename(local_path)}"


def main():
    spec = importlib.util.spec_from_file_location("preprocess", PREPROCESS)
    pre = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pre)

    real = {}
    files_rows, meta_rows = [], []
    for patient, path in REAL.items():
        sheet = pd.read_csv(path, dtype=str)
        real[patient] = sheet
        normals = set(sheet["bam_normal"])
        if len(normals) != 1:
            sys.exit(f"{patient}: expected one distinct bam_normal in {path}, found {len(normals)}")
        normal_path = normals.pop()
        shared = sheet.iloc[0][FIELDS].to_dict()
        for _, r in sheet.iterrows():
            files_rows.append((r["sample"], uri(f"{patient}_tumors", r["bam_tumor"]), f"{patient}_tumors"))
            meta_rows.append((r["sample"], patient, "tumor", r["sex"], r["platform"], r["fiber"]))
        files_rows.append((f"{patient}_PBMC", uri(f"{patient}_normal", normal_path), f"{patient}_normal"))
        meta_rows.append((f"{patient}_PBMC", patient, "normal", shared["sex"], shared["platform"], shared["fiber"]))

    files = pd.DataFrame(files_rows, columns=["sample", "file", "dataset"]).drop_duplicates()
    samples = pd.DataFrame(meta_rows, columns=["sample", "patient", "status", "sex", "platform", "fiber"])
    out = pre.build_pipeline_samplesheet(files, samples).set_index("sample")

    failures = 0
    expected_samples = set()
    for patient, sheet in real.items():
        for _, r in sheet.iterrows():
            name = r["sample"]
            expected_samples.add(name)
            problems = []
            if name not in out.index:
                problems.append("missing from the built samplesheet")
            else:
                got = out.loc[name]
                for f in FIELDS:
                    if got[f] != r[f]:
                        problems.append(f"{f}: built {got[f]!r}, real {r[f]!r}")
                want_normal = uri(f"{patient}_normal", r["bam_normal"])
                if got["bam_normal"] != want_normal:
                    problems.append(f"bam_normal: built {got['bam_normal']!r}, want {want_normal!r}")
                want_tumor = uri(f"{patient}_tumors", r["bam_tumor"])
                if got["bam_tumor"] != want_tumor:
                    problems.append(f"bam_tumor: built {got['bam_tumor']!r}, want {want_tumor!r}")
            print(f"{'PASS' if not problems else 'FAIL'} {patient} {name}" + ("" if not problems else ": " + "; ".join(problems)))
            failures += bool(problems)
    extra = sorted(set(out.index) - expected_samples)
    if extra:
        print(f"FAIL built samplesheet has samples not in the real ones: {extra}")
        failures += 1
    total = len(expected_samples)
    print(f"TOTAL {total - min(failures, total)}/{total} PASS, {len(out)} rows built, {failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
