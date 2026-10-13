# tests/local

Tests that sample with a real backend: parameter recovery, `posterior_predict()`
and `conditional_effects()` on real fits, refits through `update()`, and
comparisons against reference estimates. They are too slow for `tests/testthat`
and never run in R CMD check (the folder is build-ignored, see ADR 0004 in
`.dev/decisions/`).

Run them from the package root with CmdStan installed:

```sh
Rscript tests/local/run.R
```

`run.R` loads the package and turns the parallel workers off (the package
sets `Config/testthat/parallel` for `tests/testthat`). The weekly `local-tests`
workflow runs the same script on Ubuntu; it is not a required check. `setup.R` is sourced first and provides `expect_range()`, `SW()`
and `SM()`.

Files:

- `run.R`: the runner.
- `test-*.R`: the tests. Keep each fit short (one or two chains, a few hundred
  iterations) and assert with `expect_range()` or a tolerance, never with exact
  equality to a stored fit.
- `make-fixtures.R`: regenerates the small fitted models in
  `tests/testthat/assets/` (at most 300 KB and 200 draws each). Run it by hand
  when a fixture no longer loads.
- `conditional_effects_comprehensive.R`: a script that fits every model and
  prints `conditional_effects()` output; not yet a testthat file.

A test that only compiles a Stan program (parity or gradient checks) belongs in
`tests/testthat` behind `skip_unless_stan_tests()`, not here.
