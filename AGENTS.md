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
- A behaviour change (new or changed output or prior, or a change in what the emitted Stan code computes, which deleting an unused Stan function is not; a new error, warning or message, or a change in when one is raised or in its condition class) comes with a test in `tests/testthat/`; a reworded message or a refactor that keeps behaviour adds none, and a test that called a helper the refactor removes is rewritten against the function that now holds the code
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
- **Implicit returns**: Last expression is returned (no `return()` at the end)
- **Vectorized operations**: Prefer vector operations over loops
- **Functional patterns**: Use functional programming where possible
- **Explicit namespacing**: Always use `package::function()` (e.g., `brms::lf()`, `glue::glue()`)
- Only use functions from packages in `DESCRIPTION` "Imports"
- **Pure functions**: Always write functions that preserve the global state
- NEVER use a `seed` argument and `set.seed` inside of functions
- **Self-documenting code**: Write clear, descriptive names for functions and variables; refactor unclear code (better names, a flatter structure; extraction only under the helper rule below) rather than adding explanatory comments
- **Comments say why, not what**: a comment that restates the next line, or narrates steps (`# Compute log-likelihood`), is deleted. No banner blocks (`####...!` lines around a title): a file that needs sections is split, and a short `# ---- title ----` line is the most a section gets
- **Helpers are extracted only with two callers or a dispatch target** (S3 method, registered callback, a name looked up at run time with `do.call()` or `get()` such as the `.model_*()` constructors, or a probe of the machine or session that tests replace with `local_mocked_bindings()`, such as the `probe_*()` functions in `R/bmm_setup.R`). A single-use helper is inlined at its call site; a long single-use block stays inline, and a comment above it, if any, says why the block is there, not what it computes. A helper of four lines or fewer that is not a dispatch target exists only when it has three or more callers. Callers are call sites in `R/` (two calls from one function count as two), tests excluded, and a helper's length is its non-blank, non-comment lines from `function` to the closing brace; a numerical or parsing unit that R11 lets tests call directly (a density, a bound, an expression parser) may have one caller
- **Early returns**: Use early returns to simplify conditional flow and reduce nesting; `return()` appears only in an early return, never as the last expression
- **Scope**: the comment, helper and `return()` rules above and the rule on internal functions under Documentation Requirements bind the lines a pull request adds or changes; editing lines inside an existing function does not oblige the pull request to strip that function's roxygen block or to inline it, and a block it keeps stays in step with the function (a new argument gets its `@param`). Banners, roxygen blocks on internal functions and single-use helpers elsewhere stay until the work tracked in #699 removes them after 1.5.0; removing them by hand on an open branch conflicts with every other branch that edits the same file. A pull request opened on or before 11 November 2026 is not reworked to these rules or to the test rules added that day (R11, the second half of R3), including commits it gains later; once it has merged, its lines are left to that same work

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

Internal functions (not exported, not S3 methods) get no roxygen block: no `@param`, `@return`, `@noRd` or `@keywords internal`. What a reader needs is a plain `#` comment of one to three lines above the definition saying why the function exists and what it assumes; the argument names and the body say the rest. The exception is a generic with a page in the developer section of the reference, exported or not (`check_model()`, `check_formula()`, `pp_observables()`, `native_transform()`): it keeps its full roxygen block tagged `@keywords developer` or `@keywords internal developer`. An internal generic without such a page, even one with methods in the model files (`settable_links()`, `built_data_columns()`), follows the rule above.

### NEWS
A NEWS entry is written for a user: the symptom they saw or the feature they get, what changed in their call or output, and what to do, in under 100 words. Of bmm's functions it names only exported ones (listed in `NAMESPACE`), or ones the entry records as deprecated or removed; model parameters, links, arguments and data columns are named as the user writes them. Any other dot-prefixed or unexported function, a Stan chunk, a stanvar or how the fix works does not appear, except inside an error message quoted as the user saw it. A change that only a model author or contributor notices goes under `### Developer-facing changes`, one line, with the PR number; that line may also name a generic with a page in the developer section of the reference, such as `check_model()`.

## Model Development Pattern

### Model Structure (S3 System)
Models are S3 objects created by internal `.model_*()` functions with user-facing aliases:

**Internal function** (`R/model_mixture3p.R`):
```r
.model_mixture3p <- function(resp_error = NULL, nt_features = NULL, set_size = NULL,
                             version = "simple", call = NULL, ...) {
  structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(nt_features, set_size),
      parameters = list(mu = "...", kappa = "...", thetat = "...", thetant = "..."),
      links = list(mu = "tan_half", kappa = "log", thetat = "softmax", thetant = "softmax"),
      fixed_parameters = list(mu = 0),
      default_priors = list(...),
      version = version
    ),
    class = c("bmmodel", "circular", "non_targets", "mixture3p", paste0("mixture3p_", version)),
    call = call
  )
}
```

**User-facing alias**:
```r
#' @keywords bmmodel
#' @export
mixture3p <- function(resp_error, nt_features, set_size,
                      version = c("simple", "slot", "slot_averaging"), ...) {
  call <- match.call()
  .check_constructor_dots(list(...))
  stop_missing_args()
  version <- match.arg(version)
  .model_mixture3p(resp_error = resp_error, nt_features = nt_features,
                   set_size = set_size, version = version, call = call, ...)
}
```

### One Constructor per Response Type
One model constructor handles exactly one response type/format. If a conceptually related model applies to tasks with different response structures (e.g., binary decisions vs. confidence ratings), create separate constructors with independent distribution functions.
- Good: `sdt_yn()`, `sdt_rating()`, `sdt_mafc()` - each with its own `check_data`, `configure_model` and `d*/r*` functions
- Bad: a single `sdt()` with a `type` argument that changes the response format and branches throughout the pipeline

### Model Configuration Pipeline
Each model needs a `configure_model.*` S3 method (one per version class when the versions are configured differently, as for `ezdm` and `cswald`), in the model's own file (`R/model_<name>.R`); the generic is in `R/helpers-model.R`:

1. **Preprocess data** (extract attributes from data checks)
2. **Build formula** using `bmf2bf()` to convert bmmformula → brmsformula
3. **Define family** and attach it to the formula (`formula$family <- ...`; a brms family such as `brms::mixture()` or `brms::multinomial()`, or a `brms::custom_family()`)
4. **Return** `nlist(formula, data)`, plus `stanvars` when the model adds Stan code (a custom family, or Stan functions called from the formula, as in `sdt_rating` and `mpt`): `nlist(formula, data, stanvars)`. The list is passed on to `brms::brm()`; default priors do not belong here, `configure_prior()` builds them from the model's `default_priors`. Priors beyond the defaults go in the model's own `configure_prior.<model>()` method, as in `R/model_mpt.R`.

`use_model_template()` generates this skeleton; until the comment pass updates it, delete the banner blocks it writes.

Find more details about the typical model code structure in the [developer notes](https://venpopov.com/bmm/dev-notes/).

### Stan Code Integration
Stan code snippets go in `inst/stan_chunks/*.stan` - these are injected into brms-generated models via stanvars, for a custom family or for Stan functions called from the formula.

## Testing Requirements

**Always add/update tests** for a behaviour change (as defined under Package Requirements) - even if not explicitly requested. A refactor, a rename of an internal function, or an extracted or inlined helper adds no test: the tests that covered the behaviour before still cover it, and a test that called a removed helper is rewritten against the function that now holds the code. A reworded message adds no test either; a test that matched the old wording is updated to the new one.

Test tiers (ADR 0004, the architecture decision record `.dev/decisions/0004-test-tiers-and-budget.md`):
- `tests/testthat/` - unit tests, mock-backend pipeline tests, fixture tests and Stan compile tests. R CMD check runs them on every PR and push with four parallel workers (`Config/testthat/parallel` in DESCRIPTION). It tests the built package, which lacks the fixtures: `tests/testthat/assets/` is build-ignored (listed in `.Rbuildignore`, so not in the built package), and the fixture tests skip there. The `tests-source` workflow runs the suite from the source tree on Ubuntu, fixtures included, and enforces the time budget.
- Tests run on parallel workers (four on CI) in no fixed order. A test must not depend on another file having run, or on global state another test leaves behind (options, loaded packages, seeds).
- Stan compile tests (R-against-Stan parity and gradient checks, the two `bmm_setup()` tests that need CmdStan) run on Ubuntu for PRs and on all three OS for pushes. Each opens with `skip_unless_stan_tests()` (`tests/testthat/helper-tiers.R`; `BMM_TEST_STAN=false` skips them locally) and is registered under R1.
- A fixture test loads its fit through `load_fixture_fit()` (`helper-epred.R`), which skips when the file is absent. A test that needs the assets folder without loading a fit opens with `skip_unless_fixtures()` (`helper-tiers.R`).
- `tests/local/` - tests that sample with a real backend; build-ignored, run by hand (`Rscript tests/local/run.R`) and by the weekly `local-tests` workflow. See `tests/local/README.md`.
- `tests/internal/` - reference fits from 2024; nothing in `tests/testthat` reads them since #422, the reference summaries in `tests/local/references.rds` replace them, and the folder's deletion awaits the maintainer (ADR 0004)
- `bmm-dev/feature_tests/` - exploratory/development scripts (synced via SWITCHdrive, not tracked by git)
- `local/` - machine-specific scratch scripts (not tracked by git)

Rules, with the ADR's rule ids. `tests/testthat/test-conventions.R` checks R1, R2, R4, R7, R8, the file-level part of R10 and part of R5 (size, draws, listed in `make-fixtures.R`) and runs with the rest of the suite; a violation fails with the file name and the rule. The `tests-source` budget job checks R6, and R3 only indirectly; its CRAN budget step (`.github/scripts/cran-budget.R`) checks the cap of R10. R9, R11, the second half of R3 and the rest of R5 and R7 are checked in review. R11, the second half of R3 and the `bmm()` block in R7 were added on 2026-11-11 (#699) and are not in ADR 0004.
- R1: No sampling in `tests/testthat/`. `cmdstan_model()`, `$sample()` and `backend = "cmdstanr"` appear only at registered sites: files whose count of each pattern is written as a literal in `test-conventions.R`. A new site opens with `skip_unless_stan_tests()` and updates that count in the same PR. A model's parity test lives in its `test-model_<name>.R` (R7).
- R2: The tier gates are `skip_unless_stan_tests()` and `skip_unless_fixtures()`; `skip_on_cran()` is the CRAN gate. `skip_on_ci()` is not a gate.
- R3: A mock fit (`backend = "mock", mock_fit = 1`, spelled out: `mock = 1` reaches brms only by partial matching) uses at most two participants or 200 rows unless the test is about data size. A new `test-model_<name>.R` file has one block that runs every version through `bmm()`; a mock fit whose only expectation is `expect_silent()`/`expect_no_error()` is not added for a single version.
- R4: Cross-model behaviours (default priors, inits, stored frame and `update()`, `pp_check()`, `pp_observables()`, `conditional_effects()`, report extraction, links, `simulate_data()`) loop over the model-cases registry, never over a hand-typed list of models: `model_test_cases(kind)` in `helper-model-cases.R` reads one provider per model, `model_cases_<name>()` in `helper-model-<name>.R`. Every provider has a `stored_frame`, `sd_default`, `init` and `links` case, and a `simulate` case per version for the model's `simulate_data()` method (checked); `check_data`, `observables` and `fixture` cases where they apply. The header of `helper-model-cases.R` lists what each kind holds.
- R5: A fixture has one chain, at most 300 KB and 200 draws, and is generated by `tests/local/make-fixtures.R` (seed and versions stored in the object). A new test loads it through `load_fixture_fit()`.
- R6: The budget job fails when a test takes over 0.033 or a file over 0.23 of the summed time of all tests, or the suite over 1000 s of wall time, all measured on the 4-vCPU Ubuntu runner with four workers; no test has a cap of its own. A new test that takes more than a few seconds says in a comment why it needs the time. Run `devtools::test(reporter = "slow")` before opening a PR. The shares are read from the CI run; a local run on another machine gives different shares (`test-conditional_effects.R` 0.082 locally, 0.150 on CI), so the local run spots slow tests but does not check the caps.
- R7: A `test-model_<name>.R` file tests that model's contract only: constructor validation, `check_data()`, `configure_model()`, the every-version `bmm()` block of R3, `stancode()` content (never `nchar(code) > 0`), its d/r functions, Stan parity.
- R8: A test file seeds through `withr::local_seed()`, never `set.seed()`; a script such as `tests/local/make-fixtures.R` may.
- R9: A new model adds a recovery plan `tests/recovery/<model>.R` (bmmtools `validation_plan()`; `.github/CONTRIBUTING.md`), `test-model_<name>.R`, a `helper-model-<name>.R` provider (R4), a fixture only when post-processing needs one (then also a `fixture` case and a block in `make-fixtures.R`), and stays inside the budget (R6).
- R10: CRAN tier. The three cross-model contract files (`test-conditional_effects.R`, `test-default-priors.R`, `test-helpers-init.R`) open with a file-level `skip_on_cran()`, the only file-level ones (count pinned in `test-conventions.R`). A test that takes more than 0.4 s in CRAN mode on a developer machine (the installed tarball's tests, `NOT_CRAN` unset, one worker) opens with `skip_on_cran()`. `cran-budget.R` runs the tests CRAN runs, fails on a failing test or when the wall time exceeds its cap, and writes `cran-test-timings.csv` to find slow tests.
- R11: A test reaches code through an exported function, an S3 generic or the mock pipeline (`backend = "mock", mock_fit = 1`). An internal function is called directly only when it is a numerical or parsing unit that the exported surface cannot drive in a few lines (a density, a bound, an expression parser), and then the test names the behaviour, not the function (`test_that("each link sets the natural-scale bounds", ...)`, not `test_that(".circmix_bounds() returns ...", ...)`). A new test file is named after the behaviour or model it covers, not after the R file it reads, and no test exists only because a helper exists. R11 binds new test files only.

Test patterns from codebase:
```r
test_that("a response column missing from the data stops the fit", {
  model <- mixture2p(resp_error = "y")
  formula <- bmf(thetat ~ 1, kappa ~ 1)
  expect_error(
    check_data(model, data.frame(x = 1), formula),
    "not present in the data"
  )
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
tests/internal/      # Reference fits from 2024, to be deleted (ADR 0004)
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
3. **Test coverage** - add tests for new behaviour, not for every helper; a new test file follows R11
4. **Implicit returns** - `return()` only in an early return
5. **devtools::load_all()** - never use `library(bmm)` during development
6. **Git branching** - feature branches → PR to `develop`, or to the release's `dev-*` integration branch for features planned for a later release (see `.dev/decisions/0002-release-integration-branches.md`); never commit directly to `develop`, `dev-*` or `master`
7. **Local `.gitignore` changes** - to ignore files locally without touching the tracked `.gitignore`, add patterns to `.git/info/exclude` instead (see [venpopov.com/posts/2024/git-local-ignore](https://venpopov.com/posts/2024/git-local-ignore/))
