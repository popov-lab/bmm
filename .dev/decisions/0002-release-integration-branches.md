# 0002: Integration branches for the next two releases

- **Status:** accepted
- **Date:** 2026-10-04

## Context

`develop` holds the next CRAN release (1.4.0). Features for 1.5.0 and 1.6.0
are ready for review now. Kept on separate feature branches until their
release opens, each would drift from `develop` again and need another port.

## Decision

- `develop` stays the line for the next CRAN release.
- Two integration branches collect reviewed features for later releases:
  `dev-1.5.0` and `dev-1.6.0`.
- A feature pull request targets the branch of its release (see 0001). Bug
  fixes target `develop`.
- Changes flow forward only: `develop` → `dev-1.5.0` → `dev-1.6.0`, by merge
  commits, weekly and after any fix to shared code. Integration branches are
  never rebased or force-pushed.
- `DESCRIPTION` version lines change only on `develop`. Each integration
  branch keeps its NEWS bullets under its own heading
  (`# bmm 1.5.0 (development)`), above the heading it inherits from `develop`.
- When a release is on CRAN, the next integration branch is merged into
  `develop` and deleted: `dev-1.5.0` after 1.4.0, `dev-1.6.0` after 1.5.0. At
  most two integration branches exist at a time.
- R CMD check and test coverage run on pushes and pull requests to `dev-*`.
  The pkgdown site keeps building from `develop` only.

## Consequences

- Reviewed features merge when they are ready, not when their release opens.
- The cost is one forward merge per branch per week, and a NEWS conflict at the
  top of the file on each forward merge.
- Contributors need to know which branch to target; CONTRIBUTING says so.
- A feature that moves to another release (e.g. MPT if it misses the 1.4.0 cut
  line) is retargeted to that release's branch before it merges.
