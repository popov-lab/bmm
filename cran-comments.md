## Changes

See NEWS.md. A feature release: signal detection models for five response formats,
multinomial processing tree models, new helper functions for priors, parameters and data checks, and
bug fixes. New imports: graphics, mvtnorm and posterior. New suggested packages: loo, MPTinR,
pkgbuild, roxygen2, rstan and vctrs.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Suggests or Enhances not in mainstream repositories: cmdstanr
  Availability using Additional_repositories specification: cmdstanr yes https://stan-dev.r-universe.dev

  Comment: cmdstanr is an optional dependency not required for the package's main functionality.

## Test environments

* local macOS (arm64), R 4.6.1: `R CMD check --as-cran` with `_R_CHECK_CRAN_INCOMING_=true`, and a
  second run with `_R_CHECK_DEPENDS_ONLY_=true` (Suggests unavailable): same result, 0 errors,
  0 warnings, 1 note
* GitHub Actions (macOS, Windows, Ubuntu; R release): passing on the release commit
* win-builder, R-devel (2026-10-14 r90673, x86_64-w64-mingw32, cmdstanr not installed): 0 errors,
  0 warnings, 1 note (the same note)
* win-builder, R-release (R 4.6.1, x86_64-w64-mingw32, cmdstanr not installed): 0 errors,
  0 warnings, 1 note (the same note)
* R-hub could not be used: its dependency resolver does not honour `Additional_repositories`, so the
  optional dependency cmdstanr blocks the build before the check starts

## Reverse dependencies

One reverse dependency on CRAN, rtprep, which suggests bmm. R CMD check of rtprep 0.1.0 with this
version: Status OK.
