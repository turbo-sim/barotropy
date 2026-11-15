"""Capture untimed scientific curves for the offline benchmark inspector."""

import numpy as np
import CoolProp.CoolProp as CP

UNITS = {
    "density": "kg/m³",
    "speed_of_sound": "m/s",
    "viscosity": "Pa·s",
    "void_fraction": "—",
    "vapor_quality": "—",
}


def numbers(values):
    """JSON has no NaN: preserve unavailable diagnostics as null gaps."""
    return [float(x) if np.isfinite(x) else None for x in np.asarray(values)]


def collect_diagnostics(model, grid, provenance):
    """Called after all timing samples; never contributes to benchmark costs."""
    states = model.states
    order = (
        slice(None, None, -1) if model.process_type == "compression" else slice(None)
    )
    trajectory = {
        key: numbers(states[key][order])
        for key in ("pressure", "temperature", "entropy", "density")
    }
    source = CP.AbstractState("HEOS", model.fluid_name)
    temperatures = np.linspace(
        source.Ttriple() + 0.05, source.T_critical() * (1 - 1e-6), 180
    )
    liquid, vapor = [], []
    for temperature in temperatures:
        source.update(CP.QT_INPUTS, 0, temperature)
        liquid.append(source.smass())
        source.update(CP.QT_INPUTS, 1, temperature)
        vapor.append(source.smass())
    fits = {}
    for variable in model.poly_fitter.variables:
        reference = np.asarray(states[variable])
        at_nodes = model.poly_fitter.evaluate_polynomial(states["pressure"], variable)
        absolute = at_nodes - reference
        scale = max(float(np.max(np.abs(reference))), np.finfo(float).tiny)
        relative = np.full_like(absolute, np.nan)
        np.divide(
            100 * absolute,
            reference,
            out=relative,
            where=np.abs(reference) > 1e-12 * scale,
        )
        fits[variable] = dict(
            unit=UNITS.get(variable, "SI"),
            pressure=numbers(states["pressure"]),
            reference=numbers(reference),
            fitted_at_nodes=numbers(at_nodes),
            absolute_error=numbers(absolute),
            relative_error_percent=numbers(relative),
            curve_pressure=numbers(grid),
            curve_fitted=numbers(model.poly_fitter.evaluate_polynomial(grid, variable)),
        )
    return dict(
        schema_version=1,
        provenance=provenance,
        trajectory=trajectory,
        saturation=dict(
            temperature=numbers(temperatures),
            entropy_liquid=numbers(liquid),
            entropy_vapor=numbers(vapor),
        ),
        fits=fits,
        breakpoints=numbers(
            np.asarray(model.poly_fitter.poly_breakpoints) * model.poly_fitter.p_in
        ),
        error_basis="Polynomial minus barotropic ODE state at integration nodes; not an independent EOS validation.",
    )
