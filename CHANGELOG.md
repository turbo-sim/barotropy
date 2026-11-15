# Changelog

## Unreleased

- Exclude the standalone `barotropic_simulations_2D.py` Fluent simulation from
  generated API documentation. Importing it during a documentation build tried
  to read local case files; the corrected documentation build and deployment
  passed. This fix follows the 0.4.5 release and is available on `main`.

## 0.4.5 - 2026-11-15

Barotropy now uses the published JaxProp 0.6 API. Installations require
`jaxprop>=0.6.0,<0.7.0`; earlier JaxProp versions are no longer supported.
The Barotropy model construction, solve, fit and export workflow is retained.
The dependency range for `pysolver_view` is now `>=0.7.1,<1.0.0`.

- Fix compression endpoint assembly and preserve exact integration pressures.
- Correct blending branch selection for dry fluids such as Novec649 and retain
  the transition coordinate in Barotropy's public results.
- Respect requested polynomial degrees, avoid modifying source arrays during
  fitting, include equilibrium outlets and handle partial or absent transitions.
- Add optional `ODE_max_step` control for resolving narrow transitions.
- Reuse a configured JaxProp binary-mixture model during integration.
- Reject nonfinite or nonpositive binary-mixture ratios and preserve Barotropy's
  public state aliases when reading JaxProp's canonical fields.
- Add regression coverage for CO2, Water, Nitrogen, R134a and Novec649 across
  equilibrium, metastable and blending processes, expansion and compression,
  reversible and irreversible paths, and saturation-boundary crossings.
- Add reproducible stage timings and an offline interactive HTML report with
  clickable timing bars, T-s trajectories and polynomial-fit diagnostics. Saved
  runs can be added for comparison; failed or numerically incompatible cases
  do not receive misleading timing ratios.
- Add a manual Linux/Windows workflow to build and validate a candidate JaxProp
  wheel, run demo smoke tests, and upload benchmark JSON and HTML artifacts.
- Separate demo smoke tests from the default numerical suite and ignore
  generated test, benchmark, package-build and local-environment artifacts.
- Gate publication on the full numerical and demo suite on Linux and Windows.
- Migrate narrative documentation to MyST Markdown, refresh the developer guide
  and Fluent/CFX tutorials, and revise the documentation builder to use the
  active environment, rebuild bibliography entries and fail on build warnings.
- Standardize package descriptions and add `CITATION.cff`; version bumps now
  update the citation metadata alongside package metadata.

All **271 tests passed locally and on both Linux and Windows** against the
published JaxProp 0.6.0 dependency in the
[0.4.5 release workflow](https://github.com/turbo-sim/barotropy/actions/runs/43720439666).
The published Barotropy wheel also passed CO2 compression and Water/Nitrogen
solve, polynomial-fit and Fluent/CFX-export smoke checks outside the source
checkout; dependency checking passed.

On the same machine and Barotropy measurement code, ten comparable cases showed
observed solve speedups of 6.0-9.8 times over JaxProp 0.5.8. These timings used
three repeats per case and the prepared 0.6.0 wheel, whose Python modules match
the published wheel; they are not a fresh timing of this final patch or a
performance guarantee. R134a metastable timings are excluded from speed ratios
because the numerical endpoints differ.

See [validation results](benchmarks/VALIDATION_0_6_0.md) and
[benchmark/report instructions](benchmarks/README.md).

Published artifacts: [PyPI 0.4.5](https://pypi.org/project/barotropy/0.4.5/) and
[GitHub release](https://github.com/turbo-sim/barotropy/releases/tag/v0.4.5).
