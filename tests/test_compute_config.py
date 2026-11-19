import re
from pathlib import Path

CONFIG = Path(__file__).resolve().parents[1] / "community" / "lrsomatic" / "1.0" / "process-compute.config"


def config_code():
    """Config text with comments removed.

    Comments are stripped before matching so that a comment which mentions a
    forbidden word (e.g. "no params block", "no singularity settings") does not
    trip the checks below.
    """
    assert CONFIG.exists(), f"missing {CONFIG}"
    text = CONFIG.read_text()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return "\n".join(re.sub(r"//.*$", "", line) for line in text.splitlines())


def test_no_params_block():
    # The pipeline requires parameters on the CLI, never via -c.
    assert not re.search(r"^\s*params\s*\{", config_code(), flags=re.M)


def test_no_slurm_settings():
    code = config_code()
    for word in ("slurm", "clusterOptions", "singularity", "queueSize"):
        assert word not in code, f"{word} must not appear in the Cirro compute config"


def test_has_resource_limits():
    assert "resourceLimits" in config_code()


def test_error_strategy_is_conditional_and_retries_are_capped():
    code = config_code()
    strategies = re.findall(r"errorStrategy\s*=\s*(.+)", code)
    assert strategies, "errorStrategy must be set"
    for value in strategies:
        # no unconditional 'retry': the decision must depend on the exit status
        assert "task.exitStatus" in value, f"unconditional errorStrategy: {value}"
    retries = re.findall(r"maxRetries\s*=\s*(\S+)", code)
    assert retries, "maxRetries must be set"
    for value in retries:
        assert re.fullmatch(r"\d+", value) and 1 <= int(value) <= 3, f"maxRetries out of range: {value}"


def test_error_strategy_also_retries_unknown_exit_status():
    # a spot reclaim may surface without a container exit code (Nextflow reports Integer.MAX_VALUE)
    strategies = re.findall(r"errorStrategy\s*=\s*(.+)", config_code())
    assert strategies
    for value in strategies:
        assert "task.exitStatus == Integer.MAX_VALUE" in value, f"unknown exit status not retried: {value}"


def test_every_task_uses_the_on_demand_queue():
    # Spot hosts were lost 18 times in the first full-size run. A plain process-level `queue` was ignored by
    # Cirro (run c6012ed5), so the queue is set in a withName: '.*' block that covers every task.
    code = config_code()
    assert re.search(r"withName:\s*'\.\*'\s*\{\s*queue\s*=\s*env\('PW_ONDEMAND_JOB_QUEUE'\)\s*\}", code)
    assert not re.search(r"^    queue\s*=", code, flags=re.M)
    assert not re.search(r"withName:[^{]*MINIMAP2_ALIGN", code)
