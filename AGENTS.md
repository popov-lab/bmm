# bmm - Bayesian Measurement Models Development Guide

## Project Overview

bmm is an R package that translates domain-specific cognitive measurement models into brms/Stan syntax for Bayesian hierarchical modeling. Users specify models through S3 objects (e.g., `mixture3p()`, `imm()`, `sdm()`) which are internally transformed into brms distributional families and Stan code.

**Key Architecture**: Model specification → Formula construction → brms translation → Stan compilation

## Critical Workflows

### Development Cycle
```r
devtools::load_all()      # Load package functions (NOT library())
devtools::document()      # Generate man/ files from roxygen
devtools::check()         # R CMD CHECK - must pass before PR
devtools::test()          # Run all tests
testthat::test_file()     # Run single test file
```

**ALWAYS** use `devtools::load_all()` instead of `library(bmm)` during development.

### Package Requirements
- Must pass R CMD CHECK for CRAN (enforced by CI on every commit/PR)
- CI runs on every push and PR to `master`, `develop` and `dev-*`: R CMD check on three OS, and on Ubuntu `tests-source` (the suite from the source tree, with the time budget) and `test-coverage` (covr). On PRs, `file-size` also fails when any commit in the PR adds a file over 1 MB, even if a later commit removes it.
- All code changes require corresponding tests in `tests/testthat/`
- Documentation generated from roxygen comments (never edit `man/` or `docs/` directly)

## Code Conventions

### Error Handling & Messaging
Use custom utility functions instead of base R equivalents:
```r
stop2("Invalid parameter {param_name}")     # NOT stop()
warning2("Set size mismatch: {details}")    # NOT warning()
message2("Processing {n} observations")     # NOT message()
stopif(condition, "Error message")          # Conditional stop
warnif(condition, "Warning message")        # Conditional warning
```

These functions:
- Automatically suppress call stacks for cleaner error messages
- Support glue syntax without explicit `glue()` calls
- Accept `env.frame` parameter for environment resolution

### Function Style
- **Implicit returns**: Last expression is returned (no explicit `return()`)
- **Vectorized operations**: Prefer vector operations over loops
- **Functional patterns**: Use functional programming where possible
- **Explicit namespacing**: Always use `package::function()` (e.g., `brms::lf()`, `glue::glue()`)
- Only use functions from packages in `DESCRIPTION` "Imports"
- **Pure functions**: Always write functions that preserve the global state
- NEVER use a `seed` argument and `set.seed` inside of functions
- **Self-documenting code**: Write clear, descriptive names for functions and variables; refactor unclear code rather than adding explanatory comments
- **Single-responsibility functions**: Extract complex logic into well-named helper functions
- **Early returns**: Use early returns to simplify conditional flow and reduce nesting

```r
# Bad - never do this
myfun <- function(n, seed) {
   set.seed(seed)
   rnorm(n)
}
```

### String Formatting
Use `glue::glue()` instead of `sprintf()`:
```r
# Good
message2("Processing {nrow(data)} rows with set_size={ss}")

# Bad
sprintf("Processing %d rows with set_size=%d", nrow(data), ss)
```

Use the internal function `collapse_comma()` to create a comma separated string from a character vector
```r
message2("The following variables are missing from the data: {collapse_comma(missing_vars)}")

# Bad
message2("The following variables are missing from the data: {paste(missing_vars, collapse = ', ')}")
```


### Documentation Requirements
Exported functions must:
1. Have complete roxygen documentation
2. Be listed in `_pkgdown.yml` reference section OR have one of these keywords:
   - `@keywords extract_info` - for info extraction utilities
   - `@keywords bmmodel` - for model specification functions
   - `@keywords distribution` - for distribution functions
   - `@keywords transform` - for transformation utilities
   - `@keywords dataset` - for included datasets
   - `@keywords developer` - for internal/developer functions

## Model Development Pattern

### Model Structure (S3 System)
Models are S3 objects created by internal `.model_*()` functions with user-facing aliases:

**Internal function** (`R/model_mixture3p.R`):
```r
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
    class = c("bmmodel", "circular", "mixture3p")
  )
}
```

**User-facing alias**:
```r
#' @keywords bmmodel
#' @export
mixture3p <- function(...) { .model_mixture3p(...) }
```

### Model Configuration Pipeline
Each model class needs a `configure_model.*` S3 method in `R/helpers-model.R`:

1. **Preprocess data** (extract attributes from data checks)
2. **Build formula** using `bmf2bf()` to convert bmmformula → brmsformula
3. **Define family** (often brms mixture distributions)
4. **Return** `nlist(formula, data, family, prior, stanvars)`

Find more details about the typical model code structure in `code-structure.md`.

### Stan Code Integration
Stan code snippets go in `inst/stan_chunks/*.stan` - these are injected into brms-generated models via stanvars for custom families.

## Testing Requirements

**Always add/update tests** for code changes - even if not explicitly requested.

Test tiers (ADR 0004, the architecture decision record `.dev/decisions/0004-test-tiers-and-budget.md`):
- `tests/testthat/` - unit tests, mock-backend pipeline tests, fixture tests and Stan compile tests. R CMD check runs them on every PR and push with four parallel workers (`Config/testthat/parallel` in DESCRIPTION). It tests the built package, which lacks the fixtures: `tests/testthat/assets/` is build-ignored (listed in `.Rbuildignore`, so not in the built package), and the fixture tests skip there. The `tests-source` workflow runs the suite from the source tree on Ubuntu, fixtures included, and enforces the time budget.
- Tests run on parallel workers (four on CI) in no fixed order. A test must not depend on another file having run, or on global state another test leaves behind (options, loaded packages, seeds).
- Stan compile tests (R-against-Stan parity and gradient checks, the two `bmm_setup()` tests that need CmdStan) run on Ubuntu for PRs and on all three OS for pushes. Each opens with `skip_unless_stan_tests()` (`tests/testthat/helper-tiers.R`; `BMM_TEST_STAN=false` skips them locally) and is registered under R1.
- A fixture test loads its fit through `load_fixture_fit()` (`helper-epred.R`), which skips when the file is absent. A test that needs the assets folder without loading a fit opens with `skip_unless_fixtures()` (`helper-tiers.R`).
- `tests/local/` - tests that sample with a real backend; build-ignored, run by hand (`Rscript tests/local/run.R`) and by the weekly `local-tests` workflow. See `tests/local/README.md`.
- `tests/internal/` - reference fits from 2024, pending replacement by reference summaries (ADR 0004)
- `bmm-dev/feature_tests/` - exploratory/development scripts (synced via SWITCHdrive, not tracked by git)
- `local/` - machine-specific scratch scripts (not tracked by git)

Rules, with the ADR's rule ids. `tests/testthat/test-conventions.R` checks R1, R2, R4, R7, R8 and part of R5 (size, draws, listed in `make-fixtures.R`) and runs with the rest of the suite; a violation fails with the file name and the rule. The `tests-source` budget job checks R6, and R3 only indirectly. R9 and the rest of R5 are checked in review.
- R1: No sampling in `tests/testthat/`. `cmdstan_model()`, `$sample()` and `backend = "cmdstanr"` appear only at registered sites: files whose count of each pattern is written as a literal in `test-conventions.R`. A new site opens with `skip_unless_stan_tests()` and updates that count in the same PR. A model's parity test lives in its `test-model_<name>.R` (R7).
- R2: The tier gates are `skip_unless_stan_tests()` and `skip_unless_fixtures()`; `skip_on_cran()` is the CRAN gate. `skip_on_ci()` is not a gate.
- R3: A mock fit (`backend = "mock", mock_fit = 1`) uses at most two participants or 200 rows unless the test is about data size.
- R4: Cross-model behaviours (default priors, inits, stored frame and `update()`, `pp_check()`, `pp_observables()`, `conditional_effects()`, report extraction, links) loop over the model-cases registry, never over a hand-typed list of models: `model_test_cases(kind)` in `helper-model-cases.R` reads one provider per model, `model_cases_<name>()` in `helper-model-<name>.R`. Every provider has a `stored_frame`, `sd_default`, `init` and `links` case (checked); `check_data`, `observables` and `fixture` cases where they apply. The header of `helper-model-cases.R` lists what each kind holds.
- R5: A fixture has one chain, at most 300 KB and 200 draws, and is generated by `tests/local/make-fixtures.R` (seed and versions stored in the object). A new test loads it through `load_fixture_fit()`.
- R6: The budget job fails when a test takes over 66 s, a file over 680 s, or the suite over 800 s, all measured on the 4-vCPU Ubuntu runner with four workers; no test has a cap of its own. A new test that takes more than a few seconds says in a comment why it needs the time. Run `devtools::test(reporter = "slow")` before opening a PR.
- R7: A `test-model_<name>.R` file tests that model's contract only: constructor validation, `check_data()`, `configure_model()`, `stancode()` content (never `nchar(code) > 0`), its d/r functions, Stan parity.
- R8: A test file seeds through `withr::local_seed()`, never `set.seed()`; a script such as `tests/local/make-fixtures.R` may.
- R9: A new model adds `test-model_<name>.R`, a `helper-model-<name>.R` provider (R4), a fixture only when post-processing needs one (then also a `fixture` case and a block in `make-fixtures.R`), and stays inside the budget (R6).

Test patterns from codebase:
```r
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
```r
# brms: response on LHS implies mu parameter
brmsformula(rt | dec(response) ~ condition)

# bmm: explicit parameter names
bmmformula(drift ~ condition, bs ~ 1, ndt ~ 1)
# Shorthand: bmf()
```

### Fixed Parameters
Set parameters to constants using assignment syntax:
```r
bmf(kappa ~ set_size, c ~ 1, a = 0.5)  # 'a' fixed to 0.5
```

## Directory Structure

```
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
tests/internal/      # Reference fits from 2024, pending replacement (ADR 0004)
vignettes/articles/  # Long-form tutorials (Rmd)
man/                 # AUTO-GENERATED - do not edit
docs/                # AUTO-GENERATED pkgdown site - do not edit
_pkgdown.yml         # pkgdown configuration and function organization

# Git-ignored development folders (not part of the package)
bmm-dev/             # Developer-specific folder for exploratory/development scripts (git-ignored)
local/               # Machine-local scratch files (not synced, not shared)
```

## Local Configuration

If a `LOCAL-AGENTS.md` file exists in the repository root, read it for machine-specific
configuration. This file is git-ignored, so each developer maintains their own copy.

## Common Pitfalls

1. **Never edit** `man/` or `docs/` directly - regenerated from roxygen/pkgdown
2. **Namespace all imports** - don't assume packages are loaded
3. **Test coverage** - add tests for all new code paths
4. **Implicit returns** - avoid explicit `return()` statements
5. **devtools::load_all()** - never use `library(bmm)` during development
6. **Git branching** - feature branches → PR to `develop`, or to the release's `dev-*` integration branch for features planned for a later release (see `.dev/decisions/0002-release-integration-branches.md`); never commit directly to `develop`, `dev-*` or `master`
7. **Local `.gitignore` changes** - to ignore files locally without touching the tracked `.gitignore`, add patterns to `.git/info/exclude` instead (see [venpopov.com/posts/2024/git-local-ignore](https://venpopov.com/posts/2024/git-local-ignore/))
