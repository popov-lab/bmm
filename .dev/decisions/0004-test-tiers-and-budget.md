# 0004: Test tiers, a time budget, and rules for the test suite

- **Status:** accepted
- **Date:** 2026-10-13
- **Issues:** #422, #123

## Context

R CMD check on GitHub Actions has a 30-minute limit on the check step. On
the last develop run before this record (2026-10-12) the tests alone took
26 min on Windows, 17 on Ubuntu and 12 on macOS; the rest of the check is
2–3 min. The MPT pull requests #486–#489 timed out on Windows. The suite
grew from 2,284 lines in 16 files (v1.2.0) to 24,617 lines in 42 files on
develop, and the open pull requests for 1.4.0 add thousands more.

Measured on develop (8e772d9c6) on 2026-10-13:

- The suite runs in 458 s sequentially on a developer machine and in 131 s
  with four testthat workers, with the same result.
- Nineteen tests compile a Stan program on every run: fourteen
  `cmdstanr::cmdstan_model()` sites (R-against-Stan parity and gradient
  checks), the `bmm_setup()` smoke test, and one sampled fit with its refit
  in `test-update.R`. They are guarded by `skip_on_cran()` only, so they run
  on all three OS; a compile is cheap on Linux and slow on the Windows
  runner.
- A mock-backend fit costs 0.3 s on 50 rows and 1.1 s on the full
  `oberauer_lin_2017`; several files fit the full dataset.
- The fixture-based tests (`tests/testthat/assets/`) never run on CI: the
  folder is build-ignored and R CMD check tests the built tarball, so 141 of
  the 153 CI skips are those tests.
- `tests/internal/` holds 27 MB of tracked reference fits from February
  2024 (bmm 0.2.1–0.4.8) that nothing in the suite reads, compared by eye in
  an Rmd; ADR 0003 left its fate to #422.

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
| Stan compile | parity and gradient checks, the `bmm_setup()` smoke test | every pull request on Ubuntu; every push on three OS | `skip_unless_stan_tests()` reads `BMM_TEST_STAN`; the Windows and macOS pull-request jobs set it to `false` and install no CmdStan |
| Sampling | recovery, `posterior_predict()` and `conditional_effects()` on real fits, refits, reference comparisons | a developer machine; the weekly `local-tests` workflow, not required | `tests/local/`, build-ignored |

The Stan compile tier stays automated because its tests are the only ones
that execute `inst/stan_chunks/*.stan`; running them on Windows and macOS
for every pull request bought only toolchain-specific breakage detection,
which the push-time run keeps. The gate is opt-out so a local
`devtools::test()` still runs them.

A time budget is enforced by the `tests-source` job
(`.github/scripts/test-budget.R`): the job fails when the suite, a file or a
test exceeds its cap. The caps in the script are provisional and are set to
1.5 times the maxima of the first measured run on the runner; every later
change is recorded here. Known offenders carry their own cap in the script
until the #422 consolidation removes them. After the first measured run the
check step's `timeout-minutes` drops from 30 to 20.

Rules, each with its checker (`tests/testthat/test-conventions.R` pins
counts per file as literals, so a new occurrence fails loudly and names the
rule):

- **R1 Tier placement.** `cmdstan_model()`, `backend = "cmdstanr"` and
  `$sample()` appear in `tests/testthat` only at the registered sites, each
  opening with `skip_unless_stan_tests()`; nothing in `tests/testthat`
  samples (that belongs in `tests/local`).
- **R2 Gates are named.** `skip_unless_stan_tests()` and
  `skip_unless_fixtures()` are the tier gates, `skip_on_cran()` is the CRAN
  gate; `skip_on_ci()` is not a gate.
- **R3 Small data.** A mock fit uses at most two participants or 200 rows
  unless the test is about data size. Enforced by the file caps.
- **R4 Contract tests loop over the registry.** Cross-model behaviours
  (default priors, inits, stored frame and `update()`, `pp_check()`,
  `pp_observables()`, `conditional_effects()`, report extraction, links)
  read their cases from the model-cases registry introduced with the rdm
  model (`model_test_cases()`, one `helper-model-<name>.R` provider per
  model); no hand-typed model vector. Checker: every model in `bmm_models()`
  has a provider. Lands with #422.
- **R5 Fixtures.** At most 300 KB and 200 draws, one chain, generated by
  `tests/local/make-fixtures.R` with the seed and the bmm and brms versions
  in the object, loaded only through `load_fixture_fit()`.
- **R6 Budget.** The caps above.
- **R7 One behaviour, one home.** A `test-model_<name>.R` file tests the
  model's own contract: constructor validation, `check_data()`,
  `configure_model()`, `stancode()` content (never that the code is
  non-empty), its d/r functions, Stan parity under the gate. Cross-model
  behaviours live in their contract file only.
- **R8 Seeds.** `withr::local_seed()` only.
- **R9 New model checklist.** A model pull request adds
  `test-model_<name>.R`, a `helper-model-<name>.R` provider, a fixture only
  when post-processing needs one, and stays inside the budget.

`tests/internal/`: the reference fits and the Rmd are to be deleted and
replaced by reference summaries (posterior mean and SD per parameter per
model, a few KB, written by `tests/local/make-references.R` and compared
with `expect_range()` in `tests/local/`). The deletion waits for the
maintainer's confirmation; the record notes the last commit holding the
blobs when it happens. The manual test files moved to `tests/local/` on
2026-10-13.

## Consequences

- Pull-request checks on Windows and macOS no longer compile Stan; a
  toolchain-specific Stan problem surfaces on the push to develop or a
  `dev-*` branch, after the merge.
- Ubuntu runs the suite twice per pull request (R CMD check and
  `tests-source`); the two jobs run concurrently.
- The fixture-based tests run on CI for the first time, in `tests-source`.
- The literal counts in `test-conventions.R` register the open pull
  requests' files ahead (MPT, rdm, circmix); any other new occurrence of a
  gated or forbidden pattern fails the suite until it is registered or
  removed.
- The #422 session inserts `skip_unless_stan_tests()` at the registered
  sites, moves the five `skip_on_ci()` sampling tests and the sampled
  `update()` test to `tests/local/`, lifts `stored_frame_cases()` into the
  model-cases registry, regenerates the two oversized fixtures, shrinks the
  full-dataset mock fits, consolidates the duplicated per-model templates,
  and lowers the ratchets and caps.
- Open: whether `create_initfun.default` should reuse the configured
  standata and stancode instead of rebuilding them (40 % of a mock fit); a
  package change, decided separately.
