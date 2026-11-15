# Changelog

## 0.4.5 - 2026-11-15

Barotropy now uses the published JaxProp 0.6 API. Installations require
`jaxprop>=0.6.0,<0.7.0`; earlier JaxProp versions are no longer supported.
The Barotropy model construction, solve, fit and export workflow is retained.

- Fix compression endpoint assembly and preserve exact integration pressures.
- Correct blending branch selection for dry fluids such as Novec649 and retain
  the transition coordinate in Barotropy's public results.
- Respect requested polynomial degrees, avoid modifying source arrays during
  fitting, include equilibrium outlets and handle partial or absent transitions.
- Add optional `ODE_max_step` control for resolving narrow transitions.
- Reuse a configured JaxProp binary-mixture model during integration.
- Add regression coverage for CO2, Water, Nitrogen, R134a and Novec649 across
  equilibrium, metastable and blending processes, expansion and compression,
  reversible and irreversible paths, and saturation-boundary crossings.
- Add reproducible stage timings and an offline interactive HTML report with
  clickable timing bars, T-s trajectories and polynomial-fit diagnostics.
- Gate publication on the full numerical and demo suite on Linux and Windows.

Local validation passed all 271 tests against the published JaxProp 0.6.0 wheel.
On the same machine and Barotropy measurement code, ten comparable cases showed
observed solve speedups of 6.0-9.8 times over JaxProp 0.5.8. These timings used
three repeats per case and the prepared 0.6.0 wheel, whose Python modules match
the published wheel; they are not a fresh timing of this final patch or a
performance guarantee. R134a metastable timings are excluded from speed ratios
because the numerical endpoints differ.

See [validation results](benchmarks/VALIDATION_0_6_0.md) and
[benchmark/report instructions](benchmarks/README.md).
