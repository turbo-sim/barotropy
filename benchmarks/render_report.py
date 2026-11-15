"""Build an offline HTML comparison from one or more benchmark JSON files."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def render_report(runs, output):
    """Embed measured inputs; no dependencies, server, or network are required."""
    payload = []
    for label, path in runs:
        report = json.loads(Path(path).read_text(encoding="utf-8"))
        if report.get("schema_version") != 1 or not report.get("cases"):
            raise ValueError(f"{path}: expected a nonempty schema-version-1 benchmark")
        payload.append(dict(label=label, filename=Path(path).name, report=report))
    if not payload:
        raise ValueError("At least one measured run is required")
    if len({item["label"] for item in payload}) != len(payload):
        raise ValueError("Run labels must be unique")
    # JSON is data in an HTML script element, so escape HTML-significant chars.
    data = (
        json.dumps(payload, allow_nan=False)
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
    )
    page = (HERE / "report.html").read_text(encoding="utf-8")
    page = page.replace(
        "__REPORT_JS__", (HERE / "report.js").read_text(encoding="utf-8")
    )
    page = page.replace(
        "__GENERATED__", datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    )
    page = page.replace("__REPORT_DATA__", data)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding="utf-8")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        action="append",
        required=True,
        metavar="LABEL=JSON",
        help="Repeat for each version, in baseline-first order",
    )
    parser.add_argument("--output", type=Path, default=HERE / "results/report.html")
    args = parser.parse_args()
    runs = []
    for item in args.run:
        label, separator, path = item.partition("=")
        if not separator or not label.strip() or not path:
            parser.error("--run requires LABEL=path/to/report.json")
        runs.append((label.strip(), path))
    try:
        print(render_report(runs, args.output).resolve())
    except (ValueError, OSError) as error:
        parser.exit(1, f"Cannot build benchmark report: {error}\n")


if __name__ == "__main__":
    main()
