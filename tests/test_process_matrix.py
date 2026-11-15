"""Independent numerical checks of the public barotropic process API."""

import CoolProp.CoolProp as CP
import numpy as np
import pytest
from scipy.integrate import solve_ivp
import barotropy as bp
from validation.process_cases import CASES


@pytest.mark.integration
@pytest.mark.parametrize("case", CASES, ids=lambda c: c.id)
def test_process(case):
    args = case.parameters()
    model = bp.BarotropicModel(**args)
    model.solve()
    states = model.states
    p, h, rho, temperature = (
        states[k] for k in ("pressure", "enthalpy", "density", "temperature")
    )
    assert all(sol.success for sol in model.ode_solution)
    np.testing.assert_allclose(p[[0, -1]], [model.p_max, model.p_min], rtol=2e-8)
    assert np.all(np.diff(p) < 0)
    for name in ("density", "temperature", "speed_of_sound"):
        assert np.all(np.isfinite(states[name])), name
        assert np.all(states[name] > 0), name
    if case.fluid == "Novec649":
        # CoolProp HEOS has no viscosity correlation for this fluid.
        assert np.all(np.isnan(states["viscosity"]))
    else:
        assert np.all(np.isfinite(states["viscosity"]))
        assert np.all(states["viscosity"] > 0)
    # Check actual inlet/outlet in physical process order, including compression.
    physical = slice(None) if args["p_out"] < args["p_in"] else slice(None, None, -1)
    phase = states["is_two_phase"][physical]
    if case.mode == "equilibrium":
        if case.family in ("flashing", "vapor_entry"):
            assert phase[0] == 0 and phase[-1] == 1
        elif case.family in ("wet_liquid", "wet_vapor"):
            assert phase[0] == 1 and phase[-1] == 0
        elif case.family == "wet_inside":
            assert np.all(phase == 1)
        else:
            assert np.all(phase == 0)
        oracle = CP.AbstractState("HEOS", case.fluid)
        oracle.update(CP.DmassP_INPUTS, args["rho_in"], args["p_in"])
        h0, s0 = oracle.hmass(), oracle.smass()
        if case.efficiency == 1:
            expected_h, expected_rho = [], []
            for pressure in p:
                oracle.update(CP.PSmass_INPUTS, pressure, s0)
                expected_h.append(oracle.hmass())
                expected_rho.append(oracle.rhomass())
        else:
            factor = (
                case.efficiency if args["p_out"] < args["p_in"] else 1 / case.efficiency
            )

            def rhs(pressure, enthalpy):
                oracle.update(CP.HmassP_INPUTS, enthalpy[0], pressure)
                return [factor / oracle.rhomass()]

            ref = solve_ivp(
                rhs,
                (args["p_in"], args["p_out"]),
                [h0],
                method="DOP853",
                rtol=1e-10,
                atol=1e-6,
                max_step=abs(args["p_in"] - args["p_out"]) / 40,
                dense_output=True,
            )
            assert ref.success
            expected_h = ref.sol(p)[0]
            expected_rho = []
            for pressure, enthalpy in zip(p, expected_h):
                oracle.update(CP.HmassP_INPUTS, enthalpy, pressure)
                expected_rho.append(oracle.rhomass())
        np.testing.assert_allclose(h, expected_h, rtol=2e-5, atol=0.02)
        np.testing.assert_allclose(rho, expected_rho, rtol=1e-4)
    if case.mode == "metastable":
        # Evaluate the homogeneous EOS directly, independently of the h-p inverse.
        oracle = CP.AbstractState("HEOS", case.fluid)
        oracle.specify_phase(
            CP.iphase_gas
            if "vapor" in case.family or "gas" in case.family
            else CP.iphase_liquid
        )
        for pressure, enthalpy, density, temp in zip(p, h, rho, temperature):
            oracle.update(CP.DmassT_INPUTS, density, temp)
            np.testing.assert_allclose(
                [pressure, enthalpy], [oracle.p(), oracle.hmass()], rtol=2e-6, atol=0.02
            )
    factor = case.efficiency if args["p_out"] < args["p_in"] else 1 / case.efficiency
    if case.mode != "equilibrium":
        np.testing.assert_allclose(
            h[-1] - h[0], factor * np.trapezoid(1 / rho, p), rtol=0.01, atol=0.05
        )
    if case.mode == "metastable" and case.family.startswith("meta_"):
        # The equilibrium projection crosses the dome while the homogeneous
        # branch remains on its chosen side of critical density.
        eq_phase = CP.PropsSI("Q", "P", p, "H", h, case.fluid)[physical]
        assert eq_phase[0] < 0 and 0 < eq_phase[-1] < 1
        critical_density = CP.PropsSI("rhocrit", case.fluid)
        assert (
            np.all(rho > critical_density)
            if case.family == "meta_liquid"
            else np.all(rho < critical_density)
        )
    if case.mode != "blending":
        entropy = states["entropy"][physical]
        if case.efficiency == 1:
            np.testing.assert_allclose(entropy, entropy[0], rtol=0, atol=0.02)
        else:
            assert entropy[-1] > entropy[0]
    else:
        assert np.all(np.isfinite(states["x"]))
        if case.family.startswith("blend_"):
            assert np.any(states["x"] >= 1)
            assert np.any((states["x"] > 0) & (states["x"] < 1))
            assert np.any(states["x"] <= 0)


@pytest.mark.integration
@pytest.mark.parametrize("fluid", ["CO2", "Water", "Nitrogen", "R134a", "Novec649"])
@pytest.mark.parametrize("family", ["blend_liquid", "blend_vapor"])
def test_blending_constitutive_density(fluid, family):
    import jaxprop as jp
    from validation.process_cases import Case

    args = Case(fluid, family, "blending").parameters()
    model = bp.BarotropicModel(**args)
    model.solve()
    x = model.states["x"]
    interior = np.flatnonzero((x > 0) & (x < 1))
    assert len(interior)
    source = jp.Fluid(fluid)
    for index in (0, int(interior[len(interior) // 2]), len(x) - 1):
        p, h = model.states["pressure"][index], model.states["enthalpy"][index]
        eq_rho = CP.PropsSI("D", "P", p, "H", h, fluid)
        limited = np.clip(x[index], 0, 1)
        weight = limited**2 * (3 - 2 * limited)
        if weight:
            q = 0 if family == "blend_liquid" else 1
            guess = (
                CP.PropsSI("D", "P", p, "Q", q, fluid),
                CP.PropsSI("T", "P", p, "Q", q, fluid),
            )
            meta = source.get_state_metastable(
                "enthalpy", h, "pressure", p, rhoT_guess=guess, solver_tolerance=1e-9
            )
            expected = (1 - weight) * eq_rho + weight * meta.rho
        else:
            expected = eq_rho
        np.testing.assert_allclose(model.states["density"][index], expected, rtol=2e-6)
