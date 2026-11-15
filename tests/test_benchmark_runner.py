"""Failed benchmark cases must remain visible and return a failing status."""

import json
import sys
from types import SimpleNamespace

import pytest

from benchmarks import run_benchmarks as runner


@pytest.mark.parametrize("keep_going", [False, True])
def test_failed_case_is_recorded_without_fabricating_timings(
    monkeypatch, tmp_path, keep_going
):
    output = tmp_path / "measurements.json"
    calls = []
    monkeypatch.setattr(
        runner, "CASES", [SimpleNamespace(id="bad"), SimpleNamespace(id="good")]
    )

    def measure(case, repeat):
        calls.append(case.id)
        if case.id == "bad":
            raise RuntimeError("Invalid density")
        return dict(
            max_density_fit_relative_error=0.0, seconds={"solve": {"median": 1.0}}
        )

    monkeypatch.setattr(runner, "measure", measure)
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_benchmarks.py", "--output", str(output)]
        + (["--keep-going"] if keep_going else []),
    )
    if keep_going:
        with pytest.raises(SystemExit) as error:
            runner.main()
        assert error.value.code == 1
    else:
        with pytest.raises(RuntimeError, match="Invalid density"):
            runner.main()
    report = json.loads(output.read_text())
    assert report["failures"]["bad"]["message"] == "Invalid density"
    assert "bad" not in report["cases"]
    assert ("good" in report["cases"]) is keep_going
    assert calls == (["bad", "good"] if keep_going else ["bad"])
