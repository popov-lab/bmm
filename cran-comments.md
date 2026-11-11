## Changes

A patch release with bug fixes only, see NEWS.md. It comes three weeks after 1.4.1 (published
2026-10-16), sooner than the one to two months the policy suggests. We did not want to wait because
several bugs in 1.4.1 give users wrong results without any warning:

* `sdm()` fits overestimated the likelihood when `c` and `kappa` were both large, which biased the
  estimates upwards.
* `sdm()` fits with two or more predictors could evaluate a design cell with another cell's
  parameters, and `ppp()` pooled such cells.
* `psdm()` returned wrong probabilities, up to 1 off.
* A prior on a parameter fixed in the formula silently freed that parameter.
* `update()` kept initial values, sampler settings and default priors that belonged to the fit it
  started from.

Users fitting these models with 1.4.1 can get wrong estimates now, so we would rather not hold the
fixes back for the 1.5.0 feature release. No new exports, no new or changed dependencies.

## R CMD check results

0 errors | 0 warnings | 0 notes

The CRAN incoming feasibility check lists cmdstanr as a suggested package not in a mainstream
repository, available through `Additional_repositories` (https://stan-dev.r-universe.dev).
win-builder reports this as INFO. cmdstanr is an optional dependency not required for the package's
main functionality.

## Test environments

* win-builder, R-devel (2026-11-03 r90835, x86_64-w64-mingw32, cmdstanr not installed): Status OK,
  check time 536 s, tests 216 s (an earlier run the same morning: 574 s, tests 244 s)
* win-builder, R-release (R 4.6.2, x86_64-w64-mingw32, cmdstanr not installed): Status OK,
  check time 520 s, tests 204 s
* local macOS (arm64), R 4.6.1: `R CMD check --as-cran` on the 1.4.2 tarball: 0 errors, 0 warnings,
  0 notes, tests 46 s. The same check with cmdstanr made unloadable: 0 errors, 0 warnings, 0 notes,
  tests 46 s. The tests CRAN runs take 41.0 s here against 35.9 s for 1.4.1 (mean of two
  interleaved runs each)
* GitHub Actions (macOS, Windows, Ubuntu; R release): R CMD check passing; the tests CRAN runs take
  128 s on the hosted Ubuntu runner
* The macOS builder and R-hub could not be used, as for 1.4.1: the macOS builder stops at package
  dependencies (brms and other imports not available there), and R-hub's dependency resolver does
  not honour `Additional_repositories`, so the optional dependency cmdstanr blocks the build before
  the check starts

## Reverse dependencies

One reverse dependency on CRAN, rtprep, which suggests bmm. R CMD check of rtprep 0.1.0 with
bmm 1.4.2: Status OK, as with bmm 1.4.1.
