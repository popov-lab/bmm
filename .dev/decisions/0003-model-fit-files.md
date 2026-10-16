# 0003: Large model-fit files: keep the history, stop the growth

- **Status:** accepted
- **Date:** 2026-10-04
- **Issue:** #438

## Context

The website articles load cached fits from `vignettes/articles/assets/`. On
2026-10-04, `develop` tracks 32 files there with 150.6 MiB, 13 of them
larger than 5 MiB, and the packed repository is 620.8 MiB. The fits added by the
SDT articles (#371, #372) are already part of the history. The CRAN package is
not affected, because the folder is excluded from the build.

Removing them from history would mean rewriting `develop` and force-pushing it.
Every open pull request and every existing clone would have to be rebuilt.

## Decision

- The history stays as it is.
- New and refitted article fits are hosted as assets of a GitHub release and
  downloaded when the website is built. They are not committed.
- A CI check fails a pull request that adds a tracked file above a size limit,
  or grows one. Files already tracked stay as they are until their article is
  next refitted.
- The reference fits in `tests/internal/` are handled in #422.

## Consequences

- The circular-mixture rebuild (#434) refits the mixture, IMM and
  `vwm_crt` article fits. These refits go to release assets; overwriting the
  tracked `.rds` files would add another copy of each to the history.
- The same applies to the articles for MPT, TCC, the multivariate models and
  the racing models.
- The size limit is 1 MiB per added file or new version of a file, in any
  commit of the pull request. The `file-size` workflow enforces it and is a
  required check on `develop` and the `dev-*` branches (#517). The fits live
  in the `article-fits` release.
