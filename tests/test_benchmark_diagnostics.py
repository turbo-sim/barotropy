"""Scientific diagnostic export: direction, units, and zero-reference errors."""

from types import SimpleNamespace

import numpy as np

from benchmarks.case_diagnostics import collect_diagnostics, numbers


def test_compression_trajectory_and_zero_reference_fit_errors():
    class Fitter:
        variables = ["density", "vapor_quality"]
        poly_breakpoints = [1.0, 0.5]
        p_in = 2e5

        def evaluate_polynomial(self, pressure, variable):
            if variable == "density":
                return np.asarray(pressure) / 1e5 + 0.1
            return np.full_like(pressure, 0.01, dtype=float)

    model = SimpleNamespace(
        process_type="compression",
        fluid_name="CO2",
        poly_fitter=Fitter(),
        states=dict(
            pressure=np.array([2e5, 1e5]),
            temperature=np.array([350.0, 300.0]),
            entropy=np.array([1200.0, 1100.0]),
            density=np.array([2.0, 1.0]),
            vapor_quality=np.array([0.0, 0.0]),
        ),
    )
    data = collect_diagnostics(
        model, np.linspace(1e5, 2e5, 5), {"source_version": "test"}
    )
    assert data["trajectory"]["pressure"] == [1e5, 2e5]
    assert data["trajectory"]["temperature"] == [300.0, 350.0]
    density = data["fits"]["density"]
    np.testing.assert_allclose(density["absolute_error"], [0.1, 0.1])
    np.testing.assert_allclose(density["relative_error_percent"], [5.0, 10.0])
    assert data["fits"]["vapor_quality"]["relative_error_percent"] == [None, None]
    np.testing.assert_allclose(
        data["fits"]["vapor_quality"]["absolute_error"], [0.01, 0.01]
    )
    assert data["breakpoints"] == [2e5, 1e5]
    assert len(data["saturation"]["temperature"]) == 180
    assert np.all(
        np.asarray(data["saturation"]["entropy_vapor"])
        > data["saturation"]["entropy_liquid"]
    )


def test_unavailable_values_are_json_null_gaps():
    assert numbers([1.0, np.nan, np.inf, -np.inf]) == [1.0, None, None, None]
