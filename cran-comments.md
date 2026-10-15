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

TO BE FILLED IN: confirm the counts and the note above against the final checks.

## Test environments

* local macOS (arm64), R x.y.z: `R CMD check --as-cran` with `_R_CHECK_CRAN_INCOMING_=true`, and a
  second run with `_R_CHECK_DEPENDS_ONLY_=true` (Suggests unavailable): TO BE FILLED IN
* GitHub Actions (macOS, Windows, Ubuntu; R release): TO BE FILLED IN
* win-builder, R-devel: TO BE FILLED IN
* win-builder, R-release: TO BE FILLED IN
* R-hub could not be used: its dependency resolver does not honour `Additional_repositories`, so the
  optional dependency cmdstanr blocks the build before the check starts

## Reverse dependencies

One reverse dependency on CRAN, rtprep, which suggests bmm. R CMD check of rtprep with this version:
TO BE FILLED IN.
