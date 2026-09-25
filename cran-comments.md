## Maintainer change

The maintainer of bmm changes from Vencislav Popov <vencislav.popov@gmail.com> to
Gidon T. Frischkorn <gfrischkorn@icloud.com>. Vencislav Popov confirms the change by a separate e-mail
to CRAN-submissions@R-project.org from the previous maintainer address. He remains an author of the
package. The repository URLs in DESCRIPTION now point to the GitHub organisation that hosts the package
(https://github.com/popov-lab/bmm); the previous URLs redirect there.

## Changes

See NEWS.md. Bug fixes, within-chain threading for one model, an opt-in link function, recalibrated
default priors for one model, and the URL and maintainer updates above. No new model and no change to
the exported API.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New maintainer: Gidon T. Frischkorn <gfrischkorn@icloud.com>
  Old maintainer(s): Vencislav Popov <vencislav.popov@gmail.com>

  Suggests or Enhances not in mainstream repositories: cmdstanr
  Availability using Additional_repositories specification: cmdstanr yes https://stan-dev.r-universe.dev

  Comment: cmdstanr is an optional dependency not required for the package's main functionality.

## Test environments

* local macOS (arm64), R 4.6.1: `R CMD check --as-cran` with `_R_CHECK_CRAN_INCOMING_=true`, and a
  second run with `_R_CHECK_DEPENDS_ONLY_=true` (Suggests unavailable): same result, 0 errors,
  0 warnings, 1 note
* GitHub Actions (macOS, Windows, Ubuntu; R release): passing on the develop content
* win-builder, R-devel (2026-09-15 r90540, x86_64-w64-mingw32, cmdstanr not installed): 0 errors,
  0 warnings, 1 note (the same note)
* win-builder, R-release (R 4.6.1, x86_64-w64-mingw32, cmdstanr not installed): 0 errors,
  0 warnings, 1 note (the same note)
* R-hub could not be used: its dependency resolver does not honour `Additional_repositories`, so the
  optional dependency cmdstanr blocks the build before the check starts

## Reverse dependencies

There are no reverse dependencies on CRAN.
