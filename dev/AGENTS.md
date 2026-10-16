# bmm - Bayesian Measurement Models Development Guide

## Project Overview

bmm is an R package that translates domain-specific cognitive
measurement models into brms/Stan syntax for Bayesian hierarchical
modeling. Users specify models through S3 objects (e.g.,
[`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md),
[`imm()`](https://venpopov.com/bmm/dev/reference/imm.md),
[`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md)) which are
internally transformed into brms distributional families and Stan code.

**Key Architecture**: Model specification → Formula construction → brms
translation → Stan compilation

## Critical Workflows

### Development Cycle

``` r

devtools::load_all()      # Load package functions (NOT library())
devtools::document()      # Generate man/ files from roxygen
devtools::check()         # R CMD CHECK - must pass before PR
devtools::test()          # Run all tests
testthat::test_file()     # Run single test file
```

**ALWAYS** use
[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
instead of [`library(bmm)`](https://github.com/popov-lab/bmm) during
development.

### Package Requirements

- Must pass R CMD CHECK for CRAN (enforced by CI on every commit/PR)
- CI runs on every push and PR to `master`, `develop` and `dev-*`: R CMD
  check on three OS, and on Ubuntu `tests-source` (the suite from the
  source tree, with the time budget) and `test-coverage` (covr). On PRs,
  `file-size` also fails when any commit in the PR adds a file over 1
  MB, even if a later commit removes it.
- All code changes require corresponding tests in `tests/testthat/`
- Documentation generated from roxygen comments (never edit `man/` or
  `docs/` directly)

## Code Conventions

### Error Handling & Messaging

Use custom utility functions instead of base R equivalents:

``` r

stop2("Invalid parameter {param_name}")     # NOT stop()
warning2("Set size mismatch: {details}")    # NOT warning()
message2("Processing {n} observations")     # NOT message()
stopif(condition, "Error message")          # Conditional stop
warnif(condition, "Warning message")        # Conditional warning
```

These functions: - Automatically suppress call stacks for cleaner error
messages - Support glue syntax without explicit `glue()` calls - Accept
`env.frame` parameter for environment resolution

### Function Style

- **Implicit returns**: Last expression is returned (no explicit
  [`return()`](https://rdrr.io/r/base/function.html))
- **Vectorized operations**: Prefer vector operations over loops
- **Functional patterns**: Use functional programming where possible
- **Explicit namespacing**: Always use `package::function()` (e.g.,
  [`brms::lf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html),
  [`glue::glue()`](https://glue.tidyverse.org/reference/glue.html))
- Only use functions from packages in `DESCRIPTION` “Imports”
- **Pure functions**: Always write functions that preserve the global
  state
- NEVER use a `seed` argument and `set.seed` inside of functions
- **Self-documenting code**: Write clear, descriptive names for
  functions and variables; refactor unclear code rather than adding
  explanatory comments
- **Single-responsibility functions**: Extract complex logic into
  well-named helper functions
- **Early returns**: Use early returns to simplify conditional flow and
  reduce nesting

``` r

# Bad - never do this
myfun <- function(n, seed) {
   set.seed(seed)
   rnorm(n)
}
```

### String Formatting

Use [`glue::glue()`](https://glue.tidyverse.org/reference/glue.html)
instead of [`sprintf()`](https://rdrr.io/r/base/sprintf.html):

``` r

# Good
message2("Processing {nrow(data)} rows with set_size={ss}")

# Bad
sprintf("Processing %d rows with set_size=%d", nrow(data), ss)
```

Use the internal function `collapse_comma()` to create a comma separated
string from a character vector

``` r

message2("The following variables are missing from the data: {collapse_comma(missing_vars)}")

# Bad
message2("The following variables are missing from the data: {paste(missing_vars, collapse = ', ')}")
```

### Documentation Requirements

Exported functions must: 1. Have complete roxygen documentation 2. Be
listed in `_pkgdown.yml` reference section OR have one of these
keywords: - `@keywords extract_info` - for info extraction utilities -
`@keywords bmmodel` - for model specification functions -
`@keywords distribution` - for distribution functions -
`@keywords transform` - for transformation utilities -
`@keywords dataset` - for included datasets - `@keywords developer` -
for internal/developer functions

## Model Development Pattern

### Model Structure (S3 System)

Models are S3 objects created by internal `.model_*()` functions with
user-facing aliases:

**Internal function** (`R/model_mixture3p.R`):

``` r

.model_mixture3p <- function(resp_error, nt_features, set_size, ...) {
  structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(nt_features, set_size),
      parameters = list(kappa = "...", thetat = "..."),
      links = list(kappa = "log", thetat = "identity"),
      fixed_parameters = list(mu1 = 0),
      default_priors = list(...)
    ),
    class = c("bmmodel", "circular", "non_targets", "mixture3p")
  )
}
```

**User-facing alias**:

``` r

#' @keywords bmmodel
#' @export
mixture3p <- function(...) { .model_mixture3p(...) }
```

### One Constructor per Response Type

One model constructor handles exactly one response type/format. If a
conceptually related model applies to tasks with different response
structures (e.g., binary decisions vs. confidence ratings), create
separate constructors with independent distribution functions. - Good:
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md),
[`sdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating.md),
[`sdt_mafc()`](https://venpopov.com/bmm/dev/reference/sdt_mafc.md) -
each with its own `check_data`, `configure_model` and `d*/r*`
functions - Bad: a single `sdt()` with a `type` argument that changes
the response format and branches throughout the pipeline

### Model Configuration Pipeline

Each model needs a `configure_model.*` S3 method (one per version class
when the versions are configured differently, as for `imm`, `ezdm` and
`cswald`), in the model’s own file (`R/model_<name>.R`); the generic is
in `R/helpers-model.R`:

1.  **Preprocess data** (extract attributes from data checks)
2.  **Build formula** using
    [`bmf2bf()`](https://venpopov.com/bmm/dev/reference/bmf2bf.md) to
    convert bmmformula → brmsformula
3.  **Define family** and attach it to the formula
    (`formula$family <- ...`; a brms family such as
    [`brms::mixture()`](https://paulbuerkner.com/brms/reference/mixture.html)
    or
    [`brms::multinomial()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
    or a
    [`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html))
4.  **Return** `nlist(formula, data)`, plus `stanvars` when the model
    adds Stan code (a custom family, or Stan functions called from the
    formula, as in `sdt_rating` and `mpt`):
    `nlist(formula, data, stanvars)`. The list is passed on to
    [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html);
    default priors do not belong here,
    [`configure_prior()`](https://venpopov.com/bmm/dev/reference/configure_prior.md)
    builds them from the model’s `default_priors`. Priors beyond the
    defaults go in the model’s own `configure_prior.<model>()` method,
    as in `R/model_mixture3p.R` and `R/model_mpt.R`.

[`use_model_template()`](https://venpopov.com/bmm/dev/reference/use_model_template.md)
generates this skeleton.

Find more details about the typical model code structure in the
[developer notes](https://venpopov.com/bmm/dev-notes/).

### Stan Code Integration

Stan code snippets go in `inst/stan_chunks/*.stan` - these are injected
into brms-generated models via stanvars, for a custom family or for Stan
functions called from the formula.

## Testing Requirements

**Always add/update tests** for code changes - even if not explicitly
requested.

Test tiers (ADR 0004, the architecture decision record
`.dev/decisions/0004-test-tiers-and-budget.md`): - `tests/testthat/` -
unit tests, mock-backend pipeline tests, fixture tests and Stan compile
tests. R CMD check runs them on every PR and push with four parallel
workers (`Config/testthat/parallel` in DESCRIPTION). It tests the built
package, which lacks the fixtures: `tests/testthat/assets/` is
build-ignored (listed in `.Rbuildignore`, so not in the built package),
and the fixture tests skip there. The `tests-source` workflow runs the
suite from the source tree on Ubuntu, fixtures included, and enforces
the time budget. - Tests run on parallel workers (four on CI) in no
fixed order. A test must not depend on another file having run, or on
global state another test leaves behind (options, loaded packages,
seeds). - Stan compile tests (R-against-Stan parity and gradient checks,
the two
[`bmm_setup()`](https://venpopov.com/bmm/dev/reference/bmm_setup.md)
tests that need CmdStan) run on Ubuntu for PRs and on all three OS for
pushes. Each opens with `skip_unless_stan_tests()`
(`tests/testthat/helper-tiers.R`; `BMM_TEST_STAN=false` skips them
locally) and is registered under R1. - A fixture test loads its fit
through `load_fixture_fit()` (`helper-epred.R`), which skips when the
file is absent. A test that needs the assets folder without loading a
fit opens with `skip_unless_fixtures()` (`helper-tiers.R`). -
`tests/local/` - tests that sample with a real backend; build-ignored,
run by hand (`Rscript tests/local/run.R`) and by the weekly
`local-tests` workflow. See `tests/local/README.md`. -
`tests/internal/` - reference fits from 2024; nothing in
`tests/testthat` reads them since \#422, the reference summaries in
`tests/local/references.rds` replace them, and the folder’s deletion
awaits the maintainer (ADR 0004) - `bmm-dev/feature_tests/` -
exploratory/development scripts (synced via SWITCHdrive, not tracked by
git) - `local/` - machine-specific scratch scripts (not tracked by git)

Rules, with the ADR’s rule ids. `tests/testthat/test-conventions.R`
checks R1, R2, R4, R7, R8, the file-level part of R10 and part of R5
(size, draws, listed in `make-fixtures.R`) and runs with the rest of the
suite; a violation fails with the file name and the rule. The
`tests-source` budget job checks R6, and R3 only indirectly; its CRAN
budget step (`.github/scripts/cran-budget.R`) checks the cap of R10. R9
and the rest of R5 are checked in review. - R1: No sampling in
`tests/testthat/`. `cmdstan_model()`, `$sample()` and
`backend = "cmdstanr"` appear only at registered sites: files whose
count of each pattern is written as a literal in `test-conventions.R`. A
new site opens with `skip_unless_stan_tests()` and updates that count in
the same PR. A model’s parity test lives in its `test-model_<name>.R`
(R7). - R2: The tier gates are `skip_unless_stan_tests()` and
`skip_unless_fixtures()`; `skip_on_cran()` is the CRAN gate.
`skip_on_ci()` is not a gate. - R3: A mock fit
(`backend = "mock", mock_fit = 1`) uses at most two participants or 200
rows unless the test is about data size. - R4: Cross-model behaviours
(default priors, inits, stored frame and
[`update()`](https://rdrr.io/r/stats/update.html),
[`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md),
[`pp_observables()`](https://venpopov.com/bmm/dev/reference/pp_observables.md),
[`conditional_effects()`](https://venpopov.com/bmm/dev/reference/conditional_effects.bmmfit.md),
report extraction, links) loop over the model-cases registry, never over
a hand-typed list of models: `model_test_cases(kind)` in
`helper-model-cases.R` reads one provider per model,
`model_cases_<name>()` in `helper-model-<name>.R`. Every provider has a
`stored_frame`, `sd_default`, `init` and `links` case (checked);
`check_data`, `observables` and `fixture` cases where they apply. The
header of `helper-model-cases.R` lists what each kind holds. - R5: A
fixture has one chain, at most 300 KB and 200 draws, and is generated by
`tests/local/make-fixtures.R` (seed and versions stored in the object).
A new test loads it through `load_fixture_fit()`. - R6: The budget job
fails when a test takes over 66 s, a file over 680 s, or the suite over
800 s, all measured on the 4-vCPU Ubuntu runner with four workers; no
test has a cap of its own. A new test that takes more than a few seconds
says in a comment why it needs the time. Run
`devtools::test(reporter = "slow")` before opening a PR. - R7: A
`test-model_<name>.R` file tests that model’s contract only: constructor
validation,
[`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md),
[`configure_model()`](https://venpopov.com/bmm/dev/reference/configure_model.md),
[`stancode()`](https://venpopov.com/bmm/dev/reference/stancode.bmmformula.md)
content (never `nchar(code) > 0`), its d/r functions, Stan parity. - R8:
A test file seeds through
[`withr::local_seed()`](https://withr.r-lib.org/reference/with_seed.html),
never [`set.seed()`](https://rdrr.io/r/base/Random.html); a script such
as `tests/local/make-fixtures.R` may. - R9: A new model adds
`test-model_<name>.R`, a `helper-model-<name>.R` provider (R4), a
fixture only when post-processing needs one (then also a `fixture` case
and a block in `make-fixtures.R`), and stays inside the budget (R6). -
R10: CRAN tier. The three cross-model contract files
(`test-conditional_effects.R`, `test-default-priors.R`,
`test-helpers-init.R`) open with a file-level `skip_on_cran()`, the only
file-level ones (count pinned in `test-conventions.R`). A test that
takes more than 0.4 s in CRAN mode on a developer machine (the installed
tarball’s tests, `NOT_CRAN` unset, one worker) opens with
`skip_on_cran()`. `cran-budget.R` runs the tests CRAN runs, fails on a
failing test or when the wall time exceeds its cap, and writes
`cran-test-timings.csv` to find slow tests.

Test patterns from codebase:

``` r

test_that("check_data() produces expected errors", {
  expect_error(
    check_data(model, data, formula),
    "expected error message regex"
  )
  expect_warning(..., "expected warning")
  expect_silent(...) # for valid input
})
```

## Data and Formula System

### bmmformula vs brmsformula

bmm uses explicit parameter prediction (NOT implicit mu):

``` r

# brms: response on LHS implies mu parameter
brmsformula(rt | dec(response) ~ condition)

# bmm: explicit parameter names
bmmformula(drift ~ condition, bs ~ 1, ndt ~ 1)
# Shorthand: bmf()
```

### Fixed Parameters

Set parameters to constants using assignment syntax:

``` r

bmf(kappa ~ set_size, c ~ 1, a = 0.5)  # 'a' fixed to 0.5
```

## Directory Structure

    R/                    # All package code
      bmm.R              # Main fitting function
      bmmformula.R       # Formula specification
      model_*.R          # Model definitions (one file per model)
      helpers-*.R        # Shared utilities organized by function
      distributions.R    # Custom distribution functions
      utils.R            # stop2/warning2/message2, stopif/warnif

    inst/stan_chunks/    # Stan code snippets for models
    data/                # Package datasets (.rda files)
    tests/testthat/      # Unit, mock-pipeline, fixture and Stan compile tests (test-*.R pattern)
    tests/local/         # Sampling tests, run by hand and weekly (build-ignored)
    tests/internal/      # Reference fits from 2024, to be deleted (ADR 0004)
    vignettes/articles/  # Long-form tutorials (Rmd)
    man/                 # AUTO-GENERATED - do not edit
    docs/                # AUTO-GENERATED pkgdown site - do not edit
    _pkgdown.yml         # pkgdown configuration and function organization

    # Git-ignored development folders (not part of the package)
    bmm-dev/             # Developer-specific folder for exploratory/development scripts (git-ignored)
    local/               # Machine-local scratch files (not synced, not shared)

## Local Configuration

If a `LOCAL-AGENTS.md` file exists in the repository root, read it for
machine-specific configuration. This file is git-ignored, so each
developer maintains their own copy.

## Common Pitfalls

1.  **Never edit** `man/` or `docs/` directly - regenerated from
    roxygen/pkgdown
2.  **Namespace all imports** - don’t assume packages are loaded
3.  **Test coverage** - add tests for all new code paths
4.  **Implicit returns** - avoid explicit
    [`return()`](https://rdrr.io/r/base/function.html) statements
5.  **devtools::load_all()** - never use
    [`library(bmm)`](https://github.com/popov-lab/bmm) during
    development
6.  **Git branching** - feature branches → PR to `develop`, or to the
    release’s `dev-*` integration branch for features planned for a
    later release (see
    `.dev/decisions/0002-release-integration-branches.md`); never commit
    directly to `develop`, `dev-*` or `master`
7.  **Local `.gitignore` changes** - to ignore files locally without
    touching the tracked `.gitignore`, add patterns to
    `.git/info/exclude` instead (see
    [venpopov.com/posts/2024/git-local-ignore](https://venpopov.com/posts/2024/git-local-ignore/))
