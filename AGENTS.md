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
- CI runs on push/PR to `master` and `develop` branches
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

Test organization:
- `tests/testthat/` - automated unit tests (run by CI)
- `tests/internal/` - long-running manual tests
- `bmm-dev/feature_tests/` - exploratory/development scripts (synced via SWITCHdrive, not tracked by git)
- `local/` - machine-specific scratch scripts (not tracked by git)

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
tests/testthat/      # Unit tests (test-*.R pattern)
tests/internal/      # Manual long-running tests
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
6. **Git branching** - feature branches → PR to `develop` (never commit directly to `develop` or `master`)
7. **Local `.gitignore` changes** - to ignore files locally without touching the tracked `.gitignore`, add patterns to `.git/info/exclude` instead (see [venpopov.com/posts/2024/git-local-ignore](https://venpopov.com/posts/2024/git-local-ignore/))
