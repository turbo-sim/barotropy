"""The report preserves measured data and safely embeds user labels."""

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "render_report", ROOT / "benchmarks/render_report.py"
)
renderer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(renderer)


def test_multiple_runs_are_embedded_without_external_dependencies(tmp_path):
    first = tmp_path / "baseline.json"
    second = tmp_path / "release.json"
    reports = [
        dict(schema_version=1, cases={"co2": {"samples": [1, 2]}}),
        dict(schema_version=1, cases={"co2": {"samples": [3, 4]}}),
    ]
    for path, data in zip((first, second), reports):
        path.write_text(json.dumps(data))
    output = renderer.render_report(
        [("Baseline", first), ("Release", second)], tmp_path / "report.html"
    )
    html = output.read_text(encoding="utf-8")
    embedded = html.split('<script id="benchmark-data" type="application/json">', 1)[
        1
    ].split("</script>", 1)[0]
    data = json.loads(embedded)
    assert [run["label"] for run in data] == ["Baseline", "Release"]
    assert [run["report"] for run in data] == reports
    assert "<script src=" not in html
    assert "__REPORT_" not in html


def test_label_cannot_close_embedded_data_script(tmp_path):
    source = tmp_path / "input.json"
    source.write_text(json.dumps(dict(schema_version=1, cases={"case": {}})))
    label = '</script><script>alert("injected")</script>'
    output = renderer.render_report([(label, source)], tmp_path / "report.html")
    html = output.read_text(encoding="utf-8")
    assert label not in html
    embedded = html.split('<script id="benchmark-data" type="application/json">', 1)[
        1
    ].split("</script>", 1)[0]
    assert json.loads(embedded)[0]["label"] == label


@pytest.mark.parametrize(
    "report",
    [{}, {"schema_version": 2, "cases": {"a": {}}}, {"schema_version": 1, "cases": {}}],
)
def test_unsupported_report_is_rejected(tmp_path, report):
    source = tmp_path / "invalid.json"
    source.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="schema-version-1"):
        renderer.render_report([("Invalid", source)], tmp_path / "report.html")
