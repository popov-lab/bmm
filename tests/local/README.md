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
- `test-*.R`: the tests, among them `test-kfold-update.R` (refits through
  `update()` on a slice of the fit's data, grouped `kfold()`),
  `test-update-fits.R` (`brms::combine_models()` on a fit and its refit),
  `test-conditional_effects-fits.R` (`conditional_effects()` on `sdt_rating()`
  and `m3()` fits), `test-sdt-fits.R` (`sdt_ranking()` recovery,
  `sdt_rating()` `posterior_predict()` and a `log_distance` fit at six
  categories),
  `test-mpt-fits.R` (`mpt()` fits read through `conditional_effects()`, and
  the prior-only fit behind the default simplex prior) and
  `test-ppp-fits.R` (`ppp()` on well-specified, misspecified and
  impossible-category `mpt()` fits). Keep each fit
  short (one or two chains, a few hundred iterations) and assert with
  `expect_range()` or a tolerance, never with exact equality to a stored fit.
- `make-fixtures.R`: regenerates four fixtures in `tests/testthat/assets/`:
  `bmmfit_ddm_ppcheck.rds`, `bmmfit_cswald_ppcheck.rds`,
  `bmmfit_ezdm3_ppcheck.rds` and `bmmfit_ezdm4_ppcheck.rds` (one chain, at
  most 300 KB and 200 draws each). Run it by hand when one of them no longer
  loads. A new fixture needs its own block in this script: fit the model and
  save it with `save_fixture()`, defined at the top of the script.
- `conditional_effects_comprehensive.R`: an older script kept for reference.
  It fits eight of the 13 models (no SDT model) and prints
  `conditional_effects()` output. It is out of date (its `ezdm()` call uses
  arguments `ezdm()` no longer has), and `run.R` does not run it.

A test that only compiles a Stan program (parity or gradient checks) belongs in
`tests/testthat` behind `skip_unless_stan_tests()`, registered under rule R1
of ADR 0004, not here.
