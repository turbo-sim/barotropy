# Jaxprop 0.6.0 downstream verification

Local Windows verification completed on 15 November 2026 against the prepared jaxprop 0.6.0 wheel. Barotropy measurement code is committed as `3302d5d`; jaxprop source is `3192cade32b299059a19cd80a5f6a0294b217279`.

The full numerical, fitting, exporter, diagnostic-report and demo suite passed **269 tests**. Two subsequently added benchmark failure-recording tests and five JavaScript calculation checks also passed. The barotropy wheel builds, and the project and lockfile checks pass. Remote Linux/Windows CI execution is not claimed.

## Matched timing comparison

Both wheels used the same Python 3.13.2 interpreter and machine, CoolProp 7.2.0, NumPy 2.4.4, SciPy 1.17.1 and JAX 0.6.2. CPU execution and x64 were enabled; `OMP_NUM_THREADS` and `MKL_NUM_THREADS` were set to 1. The published baseline was measured first, followed by the candidate, with three repeated invocations per case and a separately recorded first invocation. There was no counterbalanced repeat or control over other system activity; these are observed timings, not confidence bounds or performance guarantees.

These benchmarks measure eager HEOS/SciPy barotropic processes. Flashing cases use degree-10 polynomial fits and a maximum pressure step of 1/500 of the pressure span. They do not measure compiled bicubic or perfect-gas kernels. Solve costs include the process integration and its postprocessing; case preparation, fit/export and diagnostic-curve collection are outside that stage.

| Case | PyPI 0.5.8 solve (s) | Local 0.6.0 solve (s) | Solve speedup |
| --- | ---: | ---: | ---: |
| CO2 equilibrium flashing | 29.3694 | 3.0066 | 9.77× |
| Water equilibrium flashing | 28.9902 | 3.1965 | 9.07× |
| Nitrogen equilibrium flashing | 25.4311 | 2.9367 | 8.66× |
| R134a equilibrium flashing | 25.5025 | 3.0250 | 8.43× |
| Novec649 equilibrium flashing | 30.7000 | 3.2186 | 9.54× |
| CO2 irreversible gas compression | 0.4585 | 0.0763 | 6.01× |
| Water wet-to-liquid compression | 1.0723 | 0.1656 | 6.48× |
| CO2 metastable liquid continuation | 0.9454 | 0.1008 | 9.38× |
| R134a metastable vapor continuation | 0.8147 | 0.1385 | Excluded: different endpoints |
| CO2 liquid-side blending | 14.5955 | 1.9924 | 7.33× |
| Novec649 vapor-side blending | 15.0284 | 2.0483 | 7.34× |

All eleven benchmarks completed for both versions and satisfied the benchmark density-fit check. Ten cases had matching parameters, pressure grids and endpoint properties within the comparison budget. Their total measured-stage speedups ranged from **5.16× to 9.73×**. R134a metastable endpoint temperature and density differed beyond the budget, so its cost remains visible but its speed ratio is omitted. The published version also fails independent R134a metastable density checks; timings alone do not certify the older implementation's accuracy. No scientific assertions were relaxed for this comparison.

## Reproducibility

Wheel SHA-256 values:

- Published 0.5.8: `3b6db9bfa5bd9ee37a4c3441f92f86b772e0c88ce8609bf76d32c322fbc795a3`.
- Prepared local 0.6.0: `5031e6bd7e00dc745fd5998e836e55666d0e243d2cfe8dc228b84a278dd374f0`.

Local ignored artifacts are in `benchmarks/results/release-comparison/`: raw run reports, normalized comparison reports, manifests, the candidate test XML, and `comparison-summary.json` / `.csv`. `benchmarks/results/report.html` embeds both versions and their T–s/fit diagnostic curves. Baseline measurement began before the code commit; timed model and measurement functions were unchanged by committing. Raw source/dirty-state metadata is retained.

Regenerate the offline comparison from those local artifacts:

```sh
python benchmarks/render_report.py --run "PyPI 0.5.8=benchmarks/results/release-comparison/pypi-0.5.8.json" --run "Local 0.6.0=benchmarks/results/release-comparison/jaxprop-0.6.0.json" --output benchmarks/results/report.html
```

This record describes the prepared wheel and the measurement code at `3302d5d`. Barotropy now requires jaxprop 0.6.0 or newer within the 0.6 series; earlier versions are unsupported. Historical 0.5.8 measurements require the earlier code checkout. See [benchmark instructions](README.md) for current installation and candidate-wheel CI.

## Published-release follow-up

Jaxprop 0.6.0 was published on 15 November 2026 after its hosted Linux and Windows HEOS release gates passed. Its PyPI wheel SHA-256 is `0a8dfb68ca4d21cf97b7483e1da776b946cec108a94d50ebb770bb9848127d3c`; all 43 Python modules match the prepared downstream-validated code after line-ending normalization. Barotropy now requires `jaxprop>=0.6.0,<0.7.0`, and its lockfile pins the published 0.6.0 artifacts.

Old-version mixture/blending adapters were removed. Binary processes now reuse a configured model through the public `Fluid(method="binary_mixture")` API, with an independent mixture volume-balance assertion. The complete **271-test** suite passed in a fresh environment with the published PyPI wheel installed directly, NumPy 2.5.4, SciPy 1.17.1 and pysolver_view 0.8.2; `pip check` also passed. Local XML evidence is `benchmarks/results/release-comparison/published-0.6.0-tests.xml`.
