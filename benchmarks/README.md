# Release validation and performance

Run commands from the barotropy repository in an environment with its dependencies and pytest installed.

```sh
python -m pytest
python benchmarks/run_benchmarks.py --list
python benchmarks/run_benchmarks.py --repeat 5 --output benchmarks/results/baseline.json
```

The numerical suite covers five fluids, equilibrium/metastable/blending, expansion/compression, efficiencies 1 and 0.8, both saturation boundaries, wet inlets, and wet-to-dry paths. Cases live in `validation/process_cases.py` and assert their actual phase topology. Novec649 is dry; R134a irreversible condensation needs different conditions. Saturation-based recipes exclude triple-point and spinodal limits. Metastable wet-mixture inlets are deliberately undefined: a homogeneous liquid or vapor branch must be chosen first. Wet-inlet blending cases exercise an inlet already on the equilibrium branch.

Equilibrium isentropes use direct CoolProp P-S flashes. Irreversible equilibrium paths use an independent DOP853 P-h integration. Homogeneous branches are checked with direct density-temperature EOS calls, entropy and energy balances. Blending checks include branch recovery and cubic constitutive weights. Fits are checked at additional pressures, at transition boundaries, and with partial/absent transitions; exports must contain finite expressions. HEOS has no Novec649 viscosity correlation: its missing viscosity is asserted explicitly and excluded from benchmark fits.

The default suite excludes demo smoke tests, which generate plots and exports. Run them separately:

```sh
python -m pytest -m examples
```

## Test unreleased jaxprop

Use an isolated environment for each dependency version. An editable checkout is useful during development:

```sh
python -m pip install -e . pytest
python -m pip install -e ../jaxprop
python -m pytest
```

For the actual release gate, build a candidate wheel, install it into a clean environment together with barotropy and pytest, then run the suite. Do not rely on `poetry install` alone: the lockfile selects the released dependency. Verify `jaxprop.__file__` and `jaxprop.__version__`; editable distribution metadata can be stale. The manual **Validate jaxprop candidate** GitHub workflow builds the requested Git ref into a wheel, installs it with barotropy, runs tests on Linux and Windows, and uploads timings. It does not publish either package.

## Compare performance

After recording a baseline, install the candidate in a separate environment and run on the same machine:

```sh
python -m pytest
python benchmarks/run_benchmarks.py --repeat 5 --compare benchmarks/results/baseline.json --output benchmarks/results/candidate.json
```

The eleven representative cases share definitions with the tests. Reports contain parameters, raw endpoint properties, fitted density on 501 common pressures, fit error, solver evaluation counts, raw timing samples, first invocation, repeated median/minimum/maximum, versions, import paths, Git revisions and dirty state. Construction, solve, fit, vector evaluation and export are timed separately. Case preparation and correctness checks are excluded from timings. `first` means first invocation of that case in this process, not a fresh-interpreter cold start. Run each `--case ID` in a separate process when investigating cold behavior.

Comparison requires matching parameters and pressure grids, rejects raw endpoint changes beyond `--accuracy-rtol` (default 0.0001), and reports fitted-density changes and median timing ratios (candidate/baseline). It does not impose wall-clock CI thresholds. Always run correctness tests first; timing success alone is not scientific validation. Keep hardware, power settings, thread settings, dependencies and solver tolerances consistent. Local reports under `benchmarks/results/` are ignored by Git.

Flashing fit tests and benchmarks use degree 10 and a maximum pressure step of 1/500 of the pressure span. Their density approximation must stay within 1% on independent isentropic pressure samples; the benchmark also rejects errors above 1% at its integration nodes. Coarse low-degree fits are not assumed accurate merely because the ODE converged. Adjust sampling and polynomial degree when designing other operating ranges.

## Visual HTML reports

Generate an offline report from the saved timings. No extra Python packages, npm packages, web server, or internet connection are needed:

```sh
python benchmarks/render_report.py --run "Local candidate=benchmarks/results/local-candidate.json" --output benchmarks/results/report.html
```

Open `benchmarks/results/report.html` in a browser. The grouped bar chart lets you select solve, construction, fit, vector evaluation, export, or their total. Choose repeated medians (with observed min/max whiskers) or first invocations. A second chart shows the stage breakdown for a selected process. The ledger keeps timing, solver work and fit accuracy together, and expandable provenance retains the original environment and source information.

When the released jaxprop has been measured on the same machine, add another series:

```sh
python benchmarks/render_report.py --run "Local candidate=benchmarks/results/local-candidate.json" --run "PyPI release=benchmarks/results/release.json" --output benchmarks/results/report.html
```

Each `--run` adds a colour and a bar per available case. The first run is the baseline for the time-ratio view; values below 1 mean faster. Alternatively, use **Add benchmark JSON** inside the report to add measured runs, edit their labels, and **Save HTML snapshot** to retain the comparison in a new standalone file. Browser-loaded files stay local. Missing cases have no bars; changed parameters have no baseline ratio. Differences in machine metadata and other dependency versions are flagged. Source parameters are preserved even when you edit a display label. No synthetic release results are generated.

For totals, the report sums the timed stages within each repetition, then takes the median. The headline suite cost sums per-case median totals; it is not wall-clock time for the complete benchmark command. Stacked stage proportions use mean stage contributions scaled to the displayed median total. Whiskers show observed repeated min/max, not confidence intervals. HTML reports remain local ignored artifacts; the candidate CI workflow uploads both its JSON and HTML.

Report-specific checks:

```sh
python -m pytest tests/test_benchmark_report.py
node --test benchmarks/test_report.cjs
```

The JavaScript tests need Node only when running those developer checks; viewing or generating the HTML does not need Node.

## Click a bar to inspect the process

Each timing bar opens an offline case-inspector window within the report. It shows a temperature–entropy diagram with the saturation dome, inlet/outlet markers and process direction; the selected polynomial against the original ODE states; and relative or absolute fit errors with segment boundaries. Choose another saved run in the inspector to examine the same process after a dependency update. Density, sound speed, quality, void fraction and available viscosity fits can be selected. Relative errors are undefined at zero reference values; use absolute errors there. Error plots compare fits with their ODE training states and do not replace the independent numerical tests.

New benchmark runs capture these curves after all timed samples. Existing JSON without diagnostic curves still displays timing bars and explains why the case plots are unavailable. You can enrich a historical report without changing its timing samples, using its original frozen jaxprop installation:

```sh
python benchmarks/run_benchmarks.py --diagnostics-only benchmarks/results/local-candidate.json --output benchmarks/results/candidate-with-diagnostics.json
python benchmarks/render_report.py --run "Local candidate=benchmarks/results/candidate-with-diagnostics.json" --output benchmarks/results/report.html
```

Replay checks the jaxprop source version and original endpoint properties before attaching curves. Use the exact original wheel when available: matching version strings alone do not prove identical source code. The attached diagnostic provenance distinguishes replayed curves from the historical timings. The current HTML report has been enriched this way for all eleven cases, preserving the original timing samples and other recorded results.

## Jaxprop 0.6.0 validation

The [local verification record](VALIDATION_0_6_0.md) contains the tested wheel hashes, test outcomes, measured speedups and comparison limitations.

Barotropy accepts jaxprop 0.5.8 through 0.6.x. The comprehensive process matrix is validated against the prepared 0.6.0 wheel; the published 0.5.8 package is retained as a historical performance baseline and has known limitations under these stronger checks. Until 0.6.0 is published, the Poetry lockfile deliberately resolves the published 0.5.8 dependency. Install the prepared wheel explicitly, or use the candidate-wheel CI workflow, for the complete 0.6.0 test run. After publication, run `poetry update jaxprop` to move the default locked environment to 0.6.0.

For historical versions with failed benchmark cases, `--keep-going` records errors and continues measuring the remaining cases, while still returning a failing exit status. The HTML lists failed cases and omits their bars. Ratios also require matching endpoint pressure, temperature, density and enthalpy within the comparison budget; numerical changes are flagged instead of presented as speedups. Reports retain a measurement timestamp, optional cohort label, wheel hash, and CPU/x64/thread settings.
