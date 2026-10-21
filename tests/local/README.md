# tests/local

Tests that sample with a real backend: parameter recovery, `posterior_predict()`
and `conditional_effects()` on real fits, refits through `update()`, and
comparisons against reference estimates. They are too slow for `tests/testthat`
and never run in R CMD check: the folder is build-ignored (listed in
`.Rbuildignore`, so not in the built package). See ADR 0004 in
`.dev/decisions/`.

Run them from the package root with CmdStan installed
(`cmdstanr::install_cmdstan()`; `bmm_setup()` checks the installation):

```sh
Rscript tests/local/run.R
```

`run.R` loads the package and turns the parallel workers off (the package
sets `Config/testthat/parallel` for `tests/testthat`). The weekly `local-tests`
workflow runs the same script on Ubuntu; it is not a required check.

Files:

- `run.R`: the runner.
- `setup.R`: sourced before the tests. It sets the backend options and defines
  three helpers (`expect_range()` and `SM()` exist only in `tests/local`; the
  package itself has an internal `SW()`):
  - `expect_range(object, lower, upper)`: passes when every value of `object`
    lies inside `[lower, upper]`.
  - `SW()`: short for `suppressWarnings()`.
  - `SM()`: short for `suppressMessages()`.
- `helper-mpt.R`: sources `tests/testthat/helper-mpt.R`, so the MPT fits use
  the same tree and data builders as the unit tests.
- `helper-reference-fits.R`: `reference_cases()`, one small case per model in
  `model_names()` (model, formula, data, seed; every formula gives a parameter
  a two-level predictor), and `reference_fit(name)`, which fits a case (two
  chains, 500 warmup and 500 sampling iterations) once per R session, so the
  test files share the fits.
- `test-*.R`: the tests, among them `test-kfold-update.R` (refits through
  `update()` on a slice of the fit's data, grouped `kfold()`),
  `test-update-fits.R` (`brms::combine_models()` on a fit and its refit, and
  the init an rstan refit without a recompile stores; that test needs rstan,
  not CmdStan),
  `test-conditional_effects-fits.R` (`conditional_effects()` on `sdt_rating()`
  and `m3()` fits), `test-conditional_effects-models.R` (`conditional_effects()`
  and its plots for every parameter with a predictor in each model's reference
  case, on both scales, with native estimates inside the range of the
  parameter's link), `test-references.R` (each reference case refit and
  compared with `references.rds`), `test-sdt-fits.R` (`sdt_ranking()` recovery,
  `sdt_rating()` `posterior_predict()` and a `log_distance` fit at six
  categories),
  `test-mpt-fits.R` (`mpt()` fits read through `conditional_effects()`, and
  the prior-only fit behind the default simplex prior) and
  `test-ppp-fits.R` (`ppp()` on well-specified, misspecified and
  impossible-category `mpt()` fits). Keep each fit
  short (one or two chains, a few hundred iterations) and assert with
  `expect_range()` or a tolerance, never with exact equality to a stored fit.
- `make-fixtures.R`: regenerates every fixture in `tests/testthat/assets/`
  (one chain, at most 300 KB and 200 draws each, ADR 0004 R5), all of them or
  the ones named on the command line
  (`Rscript tests/local/make-fixtures.R bmmfit_m3_ppcheck.rds`). Each fixture
  carries its seed and the bmm and brms versions in `attr(fit, "fixture")`.
  A new fixture needs its own block in this script, saved with
  `save_fixture()`, and a `fixture` case in its model's provider.
- `references.rds` and `make-references.R`: the reference summaries of ADR
  0004 (issue #123), a few KB: per model, the posterior mean and SD of every
  fixed effect of its reference case, with the bmm, brms and CmdStan versions
  and the date in the attributes of each model's entry. `test-references.R`
  passes when a refit's estimate lies within 0.5 reference SDs of the
  reference mean and its SD within 0.7 to 1.4 times the reference SD. After an
  intended change to a model's likelihood, priors or case, regenerate that
  model's entry (`Rscript tests/local/make-references.R sdm`) or all of them
  (`Rscript tests/local/make-references.R`); a new model needs a case in
  `helper-reference-fits.R` and an entry.

A test that only compiles a Stan program (parity or gradient checks) belongs in
`tests/testthat` behind `skip_unless_stan_tests()`, registered under rule R1
of ADR 0004, not here.
