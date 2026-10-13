# 0004: Test tiers, a time budget, and rules for the test suite

- **Status:** accepted
- **Date:** 2026-10-13
- **Issues:** #422, #123

## Context

R CMD check on GitHub Actions has a 30-minute limit on the check step. On
the develop run of 2026-10-12 (commit 8e772d9c6) the tests alone took
26 min on Windows, 15 on Ubuntu and 15 on macOS; the rest of the check took
about 2 min. The MPT pull requests #486–#489 timed out on Windows. The
suite grew from 2,284 lines in 16 files (v1.2.0) to 24,617 lines in 42
files on develop (8e772d9c6), and the open pull requests for 1.4.0 add
thousands more.

Measured on 2026-10-13, on develop (8e772d9c6, or 42fb3e92 after #544)
unless stated:

- On a developer machine the suite ran sequentially in 458 s on 8e772d9c6
  and in 450 s on 42fb3e92. With four testthat workers it ran in 181 s on
  42fb3e92, which sets no start order, and in 126 to 129 s (three runs) on
  this record's pull request, which adds `Config/testthat/start-first` to
  start the long files first. No run had a failure.
- Twenty-one tests need CmdStan on every run (#545 added two more on
  2026-10-13): eighteen tests at the
  fourteen `cmdstanr::cmdstan_model()` sites (R-against-Stan parity and
  gradient checks), two `bmm_setup()` tests (one of them the smoke test),
  and one sampled fit with its refit in `test-update.R`. They skip on CRAN
  and where CmdStan is not installed. R CMD check installed CmdStan on all
  three OS for every event, so they ran everywhere; a compile is cheap on
  Linux and slow on the Windows runner.
- A mock-backend fit (`backend = "mock", mock_fit = 1`) costs 0.2 to 1.1 s
  on 50 rows and 0.2 to 1.3 s on all 15,200 rows of `oberauer_lin_2017`,
  depending on the model (`sdm()` fastest, `imm()` slowest of the four
  timed; measured on this record's pull request, which does not change the
  fit path). The data size adds at most about 0.3 s. Several files fit the
  full dataset.
- The fixture-based tests (`tests/testthat/assets/`) never run on CI: the
  folder is build-ignored and R CMD check tests the built tarball, so 141 of
  the 153 skips of the 2026-10-12 develop run are those tests.
- `tests/internal/` holds 34 MB of tracked reference fits from February
  2024 (bmm 0.2.1–0.4.8), compared by eye in an Rmd; ADR 0003 left its fate
  to #422. No test in `tests/testthat` reads them (the `restructure()` test
  that would read one always skips); `tests/local/test-summary-fits.R`
  reads one.

Issue #123 asks for tests that run like testthat but only locally, in the
way brms keeps `tests/local/`. Issue #422 schedules a consolidation of the
suite after the 1.4.0 feature freeze, because it touches the files the open
pull requests edit.

## Decision

Tests are placed in four tiers, each with a mechanism that decides where it
runs:

| Tier | Content | Runs | Mechanism |
| --- | --- | --- | --- |
| Unit and mock pipeline | pure R, mock-backend fits on small data, `stancode()` and `standata()` content, priors, inits | every pull request and push, three OS, inside R CMD check, four workers | `Config/testthat/parallel: true`, `TESTTHAT_CPUS=4` |
| Fixture post-processing | tests on the stored fits in `tests/testthat/assets/` | every pull request and push, Ubuntu, from the source tree | the `tests-source` workflow |
| Stan compile | parity and gradient checks, the two `bmm_setup()` tests that need CmdStan | every pull request on Ubuntu; every push on three OS | `skip_unless_stan_tests()` reads `BMM_TEST_STAN`; the Windows and macOS pull-request jobs set it to `false` and install no CmdStan |
| Sampling | recovery, `posterior_predict()` and `conditional_effects()` on real fits, refits, reference comparisons | a developer machine; the weekly `local-tests` workflow, not required | `tests/local/`, build-ignored |

The Stan compile tier stays automated because its tests are the only ones
that execute `inst/stan_chunks/*.stan`; running them on Windows and macOS
for every pull request bought only toolchain-specific breakage detection,
which the push-time run keeps. The gate is opt-out so a local
`devtools::test()` still runs them.

A time budget is enforced by the `tests-source` job
(`.github/scripts/test-budget.R`): the job fails when the suite, a file or a
test exceeds its cap. The caps are 1.5 times the maxima of the first
four-worker run on the runner (this record's pull request with develop at
dfde3d25 merged, 2026-10-13), rounded:

| Applies to | Measured maximum | Cap |
| --- | --- | --- |
| a test | 44.1 s | 66 s |
| a file | 450 s (`test-update.R`) | 680 s |
| the suite | 532 s | 800 s |

Two known offenders in `test-update.R` carry their own cap until the #422
consolidation removes them: "update() takes a slice of the model frame as
newdata (kfold(), reloo())" (87.6 s, cap 131 s) and "only a slice of a fit's
model frame counts as stored data" (55.5 s, cap 83 s).

History. The pull request's first CI run had one worker:
`r-lib/actions/setup-r` sets `options(Ncpus = 1)` by default, and testthat
reads `Ncpus` before `TESTTHAT_CPUS`; the script now clears `Ncpus`. Caps
taken from that run (35 s per test, 330 s per file, 600 s for the suite)
were breached by the first four-worker run: four workers on the 4-vCPU
runner cut the wall time from 978 s to 532 s but doubled each file's and
test's own time (the test times summed to 971 s on one worker and 2,032 s
on four). Every later change is recorded here.

The check step's `timeout-minutes` stays at 30 until the first push run on
develop is measured. Push runs compile Stan on Windows; on the develop run
of 2026-10-12 the check step took 28.7 min on one worker. Once a push run
is measured, the timeout becomes 1.5 times the larger of the pull-request
and push measurements.

Rules. `tests/testthat/test-conventions.R` checks R1, R2, R7, R8 and part
of R5 (size, draws, and that `make-fixtures.R` lists the fixture); it pins
counts per file as literals, so a new occurrence fails loudly and names the
rule. The budget job checks R6, and R3 only indirectly. R4's checker comes
with #422. R9 and the rest of R5 are checked in review.

- **R1 Tier placement.** `cmdstan_model()`, `backend = "cmdstanr"` and
  `$sample()` appear in `tests/testthat` only at the registered sites, each
  opening with `skip_unless_stan_tests()` (new sites now, the existing ones
  in #422); nothing in `tests/testthat` samples (that belongs in
  `tests/local`). Checker: three counts per file, of `cmdstan_model(`, of
  `backend = "cmdstanr"` in any spacing and either quote, and of
  `$sample(`. A `bmm()` call that samples with the default backend passes
  the checker; only the budget catches it, through its time. So do
  `backend = c("cmdstanr")`, an argument split over two lines, and
  `$sample (` with a space; review catches those.
- **R2 Gates are named.** `skip_unless_stan_tests()` and
  `skip_unless_fixtures()` are the tier gates, `skip_on_cran()` is the CRAN
  gate; `skip_on_ci()` is not a gate.
- **R3 Small data.** A mock fit uses at most two participants or 200 rows
  unless the test is about data size. Small data keeps the tests
  readable. It saves little time (see Context), so the file caps catch a
  breach only indirectly.
- **R4 Contract tests loop over the registry.** Cross-model behaviours
  (default priors, inits, stored frame and `update()`, `pp_check()`,
  `pp_observables()`, `conditional_effects()`, report extraction, links)
  read their cases from the model-cases registry introduced with the rdm
  model (`model_test_cases()`, one `helper-model-<name>.R` provider per
  model); no hand-typed model vector. Checker: every model in `bmm_models()`
  has a provider. Lands with #422; until then R4 is checked in review.
- **R5 Fixtures.** One chain, at most 300 KB and 200 draws, generated by
  `tests/local/make-fixtures.R` with the seed and the bmm and brms versions
  in the object. A new test loads a fixture through `load_fixture_fit()`;
  `skip_unless_fixtures()` is for a test that needs the assets folder
  without loading a fit. Most existing tests read their fixtures directly
  with `readRDS()`; seven files already use `load_fixture_fit()`.
- **R6 Budget.** The caps above.
- **R7 One behaviour, one home.** A `test-model_<name>.R` file tests the
  model's own contract: constructor validation, `check_data()`,
  `configure_model()`, `stancode()` content (never that the code is
  non-empty), its d/r functions, Stan parity under the gate. Cross-model
  behaviours live in their contract file only.
- **R8 Seeds.** A test seeds through `withr::local_seed()` only. The
  checker counts lines that start with `set.seed(`, so `base::set.seed()`
  or a call after other code on the line passes it. It scans
  `tests/testthat`; a script such as
  `tests/local/make-fixtures.R` may use `set.seed()`.
- **R9 New model checklist.** A model pull request adds
  `test-model_<name>.R`, a fixture only when post-processing needs one,
  and stays inside the budget; once the model-cases registry exists (R4),
  also a `helper-model-<name>.R` provider.

`tests/internal/`: the reference fits and the Rmd are to be deleted and
replaced by reference summaries (posterior mean and SD per parameter per
model, a few KB, written by `tests/local/make-references.R` and compared
with `expect_range()` in `tests/local/`). The deletion waits for the
maintainer's confirmation and first replaces the two references to the
fits, in `tests/local/test-summary-fits.R` and in the skipped
`restructure()` test; the record notes the last commit holding the blobs
when it happens. The manual test files, including
`test-kfold-update.R`, moved to `tests/local/` on 2026-10-13.

## Consequences

- Pull-request checks on Windows and macOS no longer compile Stan; a
  toolchain-specific Stan problem surfaces on the push to develop or a
  `dev-*` branch, after the merge.
- Ubuntu runs the suite three times per pull request: R CMD check,
  `tests-source`, and `test-coverage` (covr, without CmdStan). R CMD check
  and `tests-source` run concurrently.
- The fixture-based tests run on CI for the first time, in `tests-source`.
- The literal counts in `test-conventions.R` register the open pull
  requests' files ahead (MPT, rdm, circmix); any other new occurrence of a
  gated or forbidden pattern fails the suite until it is registered or
  removed. Registering ahead works only for a file that does not exist on
  develop yet, because a present file must match its count.
  `test-model_mixture2p.R` is not on develop and is registered ahead (one
  compile site, one `$sample()`); the circmix sites in the existing
  `test-model_imm.R` and `test-model_mixture3p.R` are registered by the
  circmix merge itself.
- A push to `master`, `develop` or a `dev-*` branch is no longer cancelled
  by a later push: `cancel-in-progress` applies to pull requests only, so
  every merge gets its R CMD check and `tests-source` results.
- The #422 session inserts `skip_unless_stan_tests()` at the registered
  sites, moves the five `skip_on_ci()` sampling tests and the sampled
  `update()` test to `tests/local/`, lifts `stored_frame_cases()` into the
  model-cases registry, regenerates the two oversized fixtures, shrinks the
  full-dataset mock fits, consolidates the duplicated per-model templates,
  and lowers the ratchets and caps.
- Open: whether `create_initfun.default` should reuse the configured
  standata and stancode instead of rebuilding them (40 % of a mock fit); a
  package change, decided separately.
