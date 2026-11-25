"""Register the LRsomatic process in ONE Cirro project (not tenant-wide).

Without --create this only builds the request from <folder>/process-definition.json (default community/lrsomatic/1.0) and prints it (no
login, no network). With --create it logs in to Cirro, creates the custom process, then syncs its configuration from
the config repo (form, input mapping, preprocess, compute config).

    python docs/registration/register_lrsomatic.py --project PROJECT_ID            # dry run
    python docs/registration/register_lrsomatic.py --project PROJECT_ID --create   # register and sync

Needs the `cirro` Python package (pip install cirro). The first login prints a Cirro link to open in a browser; the
Cirro instance comes from ~/.cirro/config.ini (run `cirro configure`) or the CIRRO_BASE_URL variable. On an old Linux
(glibc 2.17) pip may try to compile pandas; install pandas>=2.2,<3 from conda first, then pip install cirro.

Run with `python -u` so the login link is not held back by output buffering.

Verified once, against the dev project on 2026-11-16: the create call succeeded and the sync reported SUCCESSFUL. A
launch of the pipeline from Cirro has not been tested.
"""
import argparse
import json
from pathlib import Path

from cirro_api_client.v1.models import (
    CustomPipelineSettings, CustomProcessInput, Executor, FileMappingRule, PipelineCode, RepositoryType,
)

CONFIG_REPO = "yannvrb/cirro-lrsomatic"
CONFIG_BRANCH = "main"
DEFAULT_FOLDER = "community/lrsomatic/1.0"  # the dev test process is community/lrsomatic-dev/0.1
ROOT = Path(__file__).resolve().parents[2]


def build_request(project_id: str, folder: str = DEFAULT_FOLDER) -> CustomProcessInput:
    d = json.loads((ROOT / folder / "process-definition.json").read_text())
    code = d["code"]
    return CustomProcessInput(
        id=d["id"],
        name=d["name"],
        description=d["desc"],
        executor=Executor.NEXTFLOW,
        category=d["category"],
        data_type=d["dataType"],
        child_process_ids=d["childProcessIds"],
        parent_process_ids=d["parentProcessIds"],
        documentation_url=f"https://github.com/YannVRB/cirro-lrsomatic/tree/{CONFIG_BRANCH}/{folder}",
        linked_project_ids=[project_id],
        is_tenant_wide=False,
        allow_multiple_sources=d["allowMultipleSources"],
        uses_sample_sheet=d["usesSampleSheet"],
        pipeline_code=PipelineCode(
            repository_path=code["uri"],
            version=code["version"],
            repository_type=RepositoryType.GITHUB_PUBLIC,
            entry_point=code["script"],
            executor_version=code["executorVersion"],
        ),
        custom_settings=CustomPipelineSettings(
            repository=CONFIG_REPO,
            branch=CONFIG_BRANCH,
            folder=folder,
            repository_type=RepositoryType.GITHUB_PUBLIC,
        ),
        # The create API rejects Python-style named groups "(?P<name>...)" (the form the catalogue files use);
        # it wants "(?<name>...)". Only the request is converted, not process-definition.json.
        file_mapping_rules=[FileMappingRule.from_dict(json.loads(json.dumps(rule).replace("(?P<", "(?<")))
                            for rule in d["fileMappingRules"]],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, help="Cirro project ID to link the process to")
    parser.add_argument("--folder", default=DEFAULT_FOLDER, help="config folder in the repository (default: %(default)s)")
    parser.add_argument("--create", action="store_true", help="really register (default: print only)")
    args = parser.parse_args()

    request = build_request(args.project, args.folder)
    if not args.create:
        print(json.dumps(request.to_dict(), indent=2))
        print("\nDRY RUN: nothing was sent. Add --create to register.")
        raise SystemExit(0)

    from cirro.cirro_client import CirroApi

    cirro = CirroApi()
    created = cirro.processes.create_custom_process(request)
    print("created:", created)
    settings = cirro.processes.sync_custom_process(request.id)
    print("sync:", settings)
