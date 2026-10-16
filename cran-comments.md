## Resubmission

This is 1.4.0 resubmitted as 1.4.1, with the version increased as the CRAN policy prefers. The
submission of 1.4.0 (2026-10-15) was returned by the incoming pretest with a NOTE on
r-devel-windows-x86_64: overall check time 16 min, of which the tests took 12 min (Debian: OK,
tests 302 s). The tests that run every model through the full pipeline now skip on CRAN and run on
the package's continuous integration instead. On win-builder R-devel the check of 1.4.1 took 462 s
with the tests at 172 s (1.4.0: 16 min and 720 s), Status OK. The package is otherwise unchanged.

## Changes

See NEWS.md. A feature release: signal detection models for five response formats,
multinomial processing tree models, new helper functions for priors, parameters and data checks, and
bug fixes. New imports: graphics, mvtnorm and posterior. New suggested packages: loo, MPTinR,
pkgbuild, roxygen2, rstan and vctrs.

## R CMD check results

0 errors | 0 warnings | 0 notes

The CRAN incoming feasibility check lists cmdstanr as a suggested package not in a mainstream
repository, available through `Additional_repositories` (https://stan-dev.r-universe.dev).
win-builder reports this as INFO; the local check of 1.4.0 with `_R_CHECK_CRAN_INCOMING_=true`
reported it as the one NOTE. cmdstanr is an optional dependency not required for the package's main
functionality.

## Test environments

* win-builder, R-devel (2026-10-14 r90673, x86_64-w64-mingw32, cmdstanr not installed): Status OK,
  check time 462 s, tests 172 s
* win-builder, R-release (R 4.6.1, x86_64-w64-mingw32, cmdstanr not installed): Status OK,
  check time 476 s, tests 175 s
* local macOS (arm64), R 4.6.1: `R CMD check --as-cran` on the 1.4.1 tarball: 0 errors, 0 warnings,
  0 notes, tests 34 s. On 1.4.0 the same check with `_R_CHECK_CRAN_INCOMING_=true`, and a second run
  with `_R_CHECK_DEPENDS_ONLY_=true` (Suggests unavailable), gave 0 errors, 0 warnings and the note
  above
* GitHub Actions (macOS, Windows, Ubuntu; R release): passing on the 1.4.1 commit
* The macOS builder and R-hub could not be used: the macOS builder stopped at package dependencies
  (brms and other imports not available there), and R-hub's dependency resolver does not honour
  `Additional_repositories`, so the optional dependency cmdstanr blocks the build before the check
  starts

## Reverse dependencies

One reverse dependency on CRAN, rtprep, which suggests bmm. R CMD check of rtprep 0.1.0 with
bmm 1.4.0: Status OK. 1.4.1 changes only which tests run on CRAN.
