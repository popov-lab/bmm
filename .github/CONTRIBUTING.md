# Contributing to `bmm`

Thank you for helping with `bmm`. This page tells you where to branch, what
your pull request (PR) must contain, who reviews it and when it merges.

## Ways to contribute

- **Report a bug.** Open an [issue](https://github.com/popov-lab/bmm/issues)
  with the bug report template.
- **Fix a bug or improve the docs.** Open a PR. For anything larger than a
  few lines, open an issue first so we can agree on the approach.
- **Add a model.** Open an issue with the `new_model` template first. We
  decide together whether the model fits the package, and in which release.
- **Ask a question or float an idea** in the
  [Discussions](https://github.com/popov-lab/bmm/discussions).

We credit every contribution in `NEWS.md`. Contributors of a new model or a
larger feature are also listed with the role `ctb` in `DESCRIPTION`.

## Which branch to target

1. Fork the repo and create your branch from `develop`. If your change is a
   feature planned for a later release, branch from that release's integration
   branch (`dev-1.5.0`, `dev-1.6.0`) instead; the milestone of the issue tells
   you which release it is planned for. Bug fixes always go to `develop`.
2. Open the PR against the branch you started from.

The reasoning is in [ADR 0002](../.dev/decisions/0002-release-integration-branches.md).
If you are unsure, ask in the issue.

## Before you open a PR

Run these from the package root and fix what they report:

```r
devtools::document()   # regenerate man/ and NAMESPACE; never edit them by hand
devtools::test()       # all tests pass
devtools::check()      # 0 errors, 0 warnings; explain any note in the PR
```

Your PR also needs:

- **Tests** for new behaviour, in `tests/testthat/`. Test what your code does,
  not how R subsets a data frame.
- **A NEWS bullet** under `# bmm (development version)` for anything users can
  see. Keep it short: what changes for the user and what they should do. Do not
  describe how the fix works.
- **Documentation** for new or changed exported functions (roxygen).
- **For a new model:** read the
  [developer notes](https://venpopov.com/bmm/dev/dev-notes/) and follow the
  one-constructor-per-response-type rule in [AGENTS.md](../AGENTS.md). Show
  parameter recovery in the PR description (see below).

### Recovery evidence for new models

Until we have a `tests/recovery/` folder, the PR description shows the result
of a parameter recovery for the new model:

- true against recovered values, as a table or a plot;
- a hierarchical recovery with several subjects whose parameters vary;
- data simulated with an independent generator, not the model's own `r*`
  function, so that a bug in the density cannot cancel out.

Attach the script as a gist or in a collapsed block in the PR.

### Model fits for the website articles

The articles load cached fits through `bmm(..., file = )`. Do not commit new
or refitted fits: a check fails every PR that adds a file larger than 1 MiB,
or a new version of one
([ADR 0003](../.dev/decisions/0003-model-fit-files.md)).

- Point `file =` to `fits/<name>`. The folder `vignettes/articles/fits/` is
  gitignored. Fits that are already tracked in `assets/` stay there until their
  article is refitted.
- Say in the PR where we can get the fit. A maintainer uploads it to the
  `article-fits` release, and the website build downloads it from there.
- Give a refitted fit a new file name instead of replacing the old asset. The
  website builds from `develop` and reads whatever the release holds at that
  moment.

To build the articles locally, download the fits first:

```sh
gh release download article-fits --dir vignettes/articles/fits --skip-existing
```

## Review and merging

How many reviews a PR needs depends on what it touches.

| Tier | What it covers | Reviews |
|---|---|---|
| **Low-risk** | Typos, docs only, CI config, tests only, NEWS only | None. May merge on green CI by anyone with write access. |
| **Standard** | Everything else, including bug fixes and new models | One approving review from someone other than the author. |
| **Pipeline or exported API** | The shared pipeline files (`R/bmm.R`, `R/bmmformula.R`, `R/helpers-model.R`, `R/helpers-data.R`, `R/helpers-prior.R`, `R/helpers-postprocess.R`, `R/helpers-inits.R`, `R/update.R`), or the signature or return value of an exported function | Two approving reviews. |

A **new model** is a standard PR with one addition: the reviewer works through
the reviewer block in the PR template. They read the likelihood against the
cited source, run or inspect the recovery, and check the default priors on the
natural scale.

Three rules apply to all tiers:

- All three `R-CMD-check` jobs (`ubuntu-latest`, `macos-latest`,
  `windows-latest`, each on `release`) must pass. CI runs on pushes and PRs to
  `master`, `develop` and `dev-*`.
- Approvals stay valid after new pushes. If you review a PR and the author
  pushes afterwards, you can still re-request changes.
- We may decline a PR that conflicts with the principles of the package. Talk
  to us early, especially if your change touches many files.

### While there is one active maintainer

At the moment one maintainer does most of the reviewing. So that PRs are not
stuck, the maintainer may merge their own PR without a GitHub approval when
all of the following hold:

- CI is green.
- An independent review was done, and the PR says what was checked and what
  changed after it.
- For the pipeline or exported API tier, the PR stays open for at least 72
  hours with that summary posted, so others can object.

PRs from other contributors always need an approval from a maintainer. This
fallback ends when a second person takes on regular reviews.

## Coding style

- Name variables and functions in `snake_case`, without upper case letters.
- Use `stopif()`, `warnif()`, `stop2()` and `message2()`, not `stopifnot()` or
  the base R equivalents.
- Use implicit returns. Do not write `return()` at the end of a function.
- Write comments that explain why, not what. If code needs a comment to say
  what it does, rename or restructure the code.
- Validate arguments in exported functions only. Internal helpers trust their
  callers.
- Call functions as `package::function()`.

[AGENTS.md](../AGENTS.md) has the rest, and the
[developer notes](https://venpopov.com/bmm/dev/dev-notes/) explain which file
holds which step of the fitting pipeline. If you are unsure where code
belongs, ask.

## AI-assisted contributions

You may use AI tools. You are responsible for every line you submit and must be
able to explain it in review. Tell us in the PR if a large share was generated.

## Licence

By submitting a change you agree that it is released under the
[GPL-2](https://choosealicense.com/licenses/gpl-2.0/) licence that covers the
package. Contact us if that is a problem.

## Bug reports

Use the bug report template. A good report has:

- a short summary;
- steps to reproduce, ideally with a minimal code example;
- what you expected and what happened instead;
- notes on what you tried, or why you think it happens.

We appreciate thorough bug reports a lot.
