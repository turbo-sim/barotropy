"""Downstream fits, exports, absent transitions, and pressure assembly."""

import numpy as np
import pytest
import barotropy as bp
from validation.process_cases import Case, FLUIDS


@pytest.mark.parametrize("fluid", FLUIDS)
@pytest.mark.parametrize("mode", ["equilibrium", "metastable", "blending"])
def test_public_fit_and_export(fluid, mode, tmp_path):
    case = Case(fluid, "gas_compression", mode)
    args = case.parameters()
    args["polynomial_variables"] = [
        "density",
        "speed_of_sound",
        "vapor_quality",
        "void_fraction",
    ]
    model = bp.BarotropicModel(**args)
    model.solve()
    original = {key: value.copy() for key, value in model.states.items()}
    model.fit_polynomials()
    # Midpoints use an independent entropy flash, rather than training nodes.
    import CoolProp.CoolProp as CP

    p = np.linspace(model.p_min, model.p_max, 41)
    entropy = model.state_in["entropy"]
    expected = CP.PropsSI("D", "P", p, "S", entropy, fluid)
    actual = model.poly_fitter.evaluate_polynomial(p, "density")
    np.testing.assert_allclose(actual, expected, rtol=3e-4)
    assert (
        max(len(poly.coef) for poly in model.poly_fitter.poly_handles["density"]) <= 5
    )
    for key, value in original.items():
        np.testing.assert_array_equal(model.states[key], value)
    model.export_expressions_fluent(str(tmp_path))
    model.export_expressions_cfx(str(tmp_path))
    files = list(tmp_path.glob("*.txt"))
    assert len(files) == 2
    for path in files:
        text = path.read_text().lower()
        assert "density" in text
        assert "nan" not in text and "inf" not in text


@pytest.mark.parametrize("family", ["blend_liquid", "blend_vapor"])
def test_complete_blending_fit(family):
    model = bp.BarotropicModel(**Case("CO2", family, "blending").parameters())
    model.solve()
    model.fit_polynomials()
    fitter = model.poly_fitter
    for boundary in fitter.poly_breakpoints[1:-1]:
        index = fitter.poly_breakpoints.index(boundary)
        np.testing.assert_allclose(
            fitter.poly_handles["density"][index - 1](boundary),
            fitter.poly_handles["density"][index](boundary),
            rtol=1e-8,
        )
    predicted = fitter.evaluate_polynomial(model.states["pressure"], "density")
    np.testing.assert_allclose(predicted, model.states["density"], rtol=2e-3)


def test_two_component_smoke():
    model = bp.BarotropicModel(
        ["Water", "Nitrogen"],
        p_in=2e5,
        T_in=600.0,
        p_out=1e5,
        mixture_ratio=10.0,
        polynomial_variables=["density"],
    )
    model.solve()
    model.fit_polynomials()
    assert all(sol.success for sol in model.ode_solution)
    assert np.all(model.states["density"] > 0)


@pytest.mark.parametrize("region", ["before", "early", "late", "after"])
def test_partial_blending_fit(region):
    model = bp.BarotropicModel(**Case("CO2", "blend_liquid", "blending").parameters())
    model.solve()
    x = model.states["x"]
    masks = {"before": x >= 1, "early": x > 0.4, "late": x < 0.6, "after": x <= 0}
    states = {key: value[masks[region]] for key, value in model.states.items()}
    fitter = bp.PolynomialFitter(
        states, model.state_in, ["density"], 4, "blending", "one-component"
    )
    fitter.fit_polynomials()
    fitted = fitter.evaluate_polynomial(states["pressure"], "density")
    np.testing.assert_allclose(fitted, states["density"], rtol=2e-3)


@pytest.mark.parametrize("fluid", FLUIDS)
def test_flashing_fit_at_independent_pressures(fluid):
    import CoolProp.CoolProp as CP

    args = Case(fluid, "flashing").parameters()
    args.update(
        polynomial_degree=10,
        ODE_max_step=abs(args["p_out"] - args["p_in"]) / 500,
        polynomial_variables=["density"],
    )
    model = bp.BarotropicModel(**args)
    model.solve()
    model.fit_polynomials()
    pressure = np.linspace(model.p_min, model.p_max, 401)
    expected = CP.PropsSI("D", "P", pressure, "S", model.state_in["entropy"], fluid)
    actual = model.poly_fitter.evaluate_polynomial(pressure, "density")
    np.testing.assert_allclose(actual, expected, rtol=0.01)
