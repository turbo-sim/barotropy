"""Repeatable downstream benchmarks; run from any working directory."""

import argparse
from datetime import datetime, timezone
from dataclasses import asdict
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import barotropy as bp
import jaxprop as jp
from validation.process_cases import Case, FLUIDS
from benchmarks.case_diagnostics import collect_diagnostics

CASES = [Case(f, "flashing") for f in FLUIDS] + [
    Case("CO2", "gas_compression", efficiency=0.8),
    Case("Water", "wet_liquid", efficiency=0.8),
    Case("CO2", "meta_liquid", "metastable"),
    Case("R134a", "meta_vapor", "metastable"),
    Case("CO2", "blend_liquid", "blending"),
    Case("Novec649", "blend_vapor", "blending"),
]


def git_info(path):
    def git(*args):
        return subprocess.check_output(
            ["git", "-C", str(path), *args], text=True, stderr=subprocess.DEVNULL
        ).strip()

    try:
        return dict(
            revision=git("rev-parse", "HEAD"), dirty=bool(git("status", "--porcelain"))
        )
    except (OSError, subprocess.CalledProcessError):
        return None


def measure(case, repeat):
    # Parameter/oracle preparation is excluded from model timings.
    params = case.parameters()
    variables = ["density", "speed_of_sound", "void_fraction", "vapor_quality"]
    if case.fluid != "Novec649":
        variables.append("viscosity")
    params["polynomial_variables"] = variables
    if case.family == "flashing":
        params.update(
            polynomial_degree=10,
            ODE_max_step=abs(params["p_out"] - params["p_in"]) / 500,
        )
    samples = []
    for _ in range(repeat + 1):
        timings = {}

        def timed(name, call):
            start = time.perf_counter()
            value = call()
            timings[name] = time.perf_counter() - start
            return value

        model = timed("construct", lambda: bp.BarotropicModel(**params))
        timed("solve", model.solve)
        if not all(sol.success for sol in model.ode_solution):
            raise RuntimeError(f"Integration failed: {case.id}")
        for name in ("density", "temperature", "speed_of_sound"):
            if not np.all(np.isfinite(model.states[name])) or not np.all(
                model.states[name] > 0
            ):
                raise RuntimeError(f"Invalid {name}: {case.id}")
        np.testing.assert_allclose(
            model.states["pressure"][[0, -1]], [model.p_max, model.p_min], rtol=2e-6
        )
        timed("fit", model.fit_polynomials)
        grid = np.linspace(model.p_min, model.p_max, 501)
        prediction = timed(
            "evaluate_501",
            lambda: model.poly_fitter.evaluate_polynomial(grid, "density"),
        )
        if not np.all(np.isfinite(prediction)) or not np.all(prediction > 0):
            raise RuntimeError(f"Invalid density fit: {case.id}")
        with tempfile.TemporaryDirectory(prefix="barotropy-benchmark-") as directory:
            timed(
                "export",
                lambda: (
                    model.export_expressions_fluent(directory),
                    model.export_expressions_cfx(directory),
                ),
            )
        samples.append(timings)
    summary = {
        name: dict(
            first=samples[0][name],
            median=statistics.median(s[name] for s in samples[1:]),
            minimum=min(s[name] for s in samples[1:]),
            maximum=max(s[name] for s in samples[1:]),
        )
        for name in samples[0]
    }
    return dict(
        case=asdict(case),
        parameters=params,
        seconds=summary,
        samples=samples,
        raw_endpoints={
            key: model.states[key][[0, -1]].tolist()
            for key in ("pressure", "temperature", "density", "enthalpy")
        },
        max_density_fit_relative_error=float(
            np.max(
                np.abs(
                    model.poly_fitter.evaluate_polynomial(
                        model.states["pressure"], "density"
                    )
                    / model.states["density"]
                    - 1
                )
            )
        ),
        nfev=[sol.nfev for sol in model.ode_solution],
        pressure=grid.tolist(),
        density=prediction.tolist(),
        diagnostics=collect_diagnostics(
            model,
            grid,
            dict(
                source_version=jp.__version__,
                import_path=jp.__file__,
                source_revision=os.environ.get("JAXPROP_SOURCE_REVISION"),
                timing_scope="Captured after all timing samples",
            ),
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument(
        "--case", action="append", help="Exact case ID; repeat to select several"
    )
    parser.add_argument("--list", action="store_true")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "benchmarks/results/latest.json"
    )
    parser.add_argument(
        "--accuracy-rtol",
        type=float,
        default=1e-4,
        help="Raw endpoint comparison tolerance",
    )
    parser.add_argument(
        "--compare", type=Path, help="Previous JSON on the same machine"
    )
    parser.add_argument(
        "--diagnostics-only",
        type=Path,
        metavar="JSON",
        help="Replay a historical report with the same jaxprop version; preserve timings",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Record failed cases and continue; still exit nonzero on failures",
    )
    args = parser.parse_args()
    if args.list:
        print("\n".join(c.id for c in CASES))
        return
    if args.repeat < 1:
        parser.error("--repeat must be at least 1")
    if args.diagnostics_only:
        report = json.loads(args.diagnostics_only.read_text(encoding="utf-8"))
        if report["jaxprop"]["source_version"] != jp.__version__:
            parser.error(
                "Diagnostics replay requires the original jaxprop source version"
            )
        for case_id, result in report["cases"].items():
            model = bp.BarotropicModel(**result["parameters"])
            model.solve()
            for key, expected in result["raw_endpoints"].items():
                np.testing.assert_allclose(
                    model.states[key][[0, -1]],
                    expected,
                    rtol=1e-6,
                    atol=1e-6,
                    err_msg=f"Historical replay changed {case_id}: {key}",
                )
            model.fit_polynomials()
            result["diagnostics"] = collect_diagnostics(
                model,
                np.asarray(result["pressure"]),
                dict(
                    source_version=jp.__version__,
                    import_path=jp.__file__,
                    source_revision=os.environ.get("JAXPROP_SOURCE_REVISION"),
                    timing_scope="Replayed separately; original timing samples preserved",
                ),
            )
            print(f"Captured diagnostics: {case_id}", flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        return
    selected = [c for c in CASES if args.case is None or c.id in args.case]
    if args.case and set(args.case) - {c.id for c in CASES}:
        parser.error("Unknown --case; use --list")
    report = dict(
        schema_version=1,
        measured_at=datetime.now(timezone.utc).isoformat(),
        cohort=os.environ.get("BAROTROPY_BENCHMARK_COHORT"),
        platform=platform.platform(),
        processor=platform.processor(),
        python=sys.version,
        executable=sys.executable,
        dependencies={
            p: importlib.metadata.version(p)
            for p in ("jaxprop", "CoolProp", "numpy", "scipy", "jax", "jaxlib", "pysolver_view")
        },
        execution={name: os.environ.get(name) for name in ("JAX_ENABLE_X64", "JAX_PLATFORM_NAME", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
        jaxprop=dict(
            source_revision=os.environ.get("JAXPROP_SOURCE_REVISION"),
            wheel_sha256=os.environ.get("JAXPROP_WHEEL_SHA256"),
            source_version=jp.__version__,
            path=jp.__file__,
            git=git_info(Path(jp.__file__).parents[1]),
        ),
        barotropy=dict(
            source_version=bp.__version__, path=bp.__file__, git=git_info(ROOT)
        ),
        cases={},
        failures={},
    )
    old = json.loads(args.compare.read_text()) if args.compare else None
    for case in selected:
        try:
            result = measure(case, args.repeat)
            if result["max_density_fit_relative_error"] > 0.01:
                raise RuntimeError(f"Density fit exceeds one percent: {case.id}")
            if old:
                previous = old["cases"][case.id]
                if previous["parameters"] != result["parameters"]:
                    raise ValueError(f"Cannot compare changed parameters: {case.id}")
                np.testing.assert_allclose(
                    result["pressure"], previous["pressure"], rtol=1e-10
                )
                for key, values in result["raw_endpoints"].items():
                    np.testing.assert_allclose(
                        values,
                        previous["raw_endpoints"][key],
                        rtol=args.accuracy_rtol,
                        atol=1e-5,
                        err_msg=f"{case.id}: {key} regression",
                    )
                result["comparison"] = dict(
                    median_ratios={
                        key: value["median"] / previous["seconds"][key]["median"]
                        for key, value in result["seconds"].items()
                    },
                    max_density_relative_change=float(
                        np.max(
                            np.abs(
                                np.array(result["density"]) / previous["density"] - 1
                            )
                        )
                    ),
                )
            report["cases"][case.id] = result
        except Exception as error:
            report["failures"][case.id] = dict(
                type=type(error).__name__, message=str(error)
            )
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
            )
            print(f"FAILED {case.id}: {type(error).__name__}: {error}", flush=True)
            if not args.keep_going:
                raise
            continue
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(
            f"{case.id}: solve {result['seconds']['solve']['median']:.4f} s", flush=True
        )

    if report["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
