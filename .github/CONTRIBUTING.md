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
   branch (e.g. `dev-1.6.0`) instead; the milestone of the issue tells
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
  not how R subsets a data frame. A new model adds `test-model_<name>.R` and a
  provider `model_cases_<name>()` in `tests/testthat/helper-model-<name>.R`
  with a `stored_frame`, `sd_default`, `init` and `links` case, which the
  cross-model tests read. Gate tests by tier: `skip_unless_stan_tests()` for a
  test that compiles Stan (and update that file's count of `cmdstan_model()`
  and `$sample()` in `test-conventions.R`), `skip_on_cran()` for one over 0.4 s in CRAN mode
  (rules R1, R4, R9 and R10 in
  [ADR 0004](../.dev/decisions/0004-test-tiers-and-budget.md); short form in
  [AGENTS.md](../AGENTS.md), "Testing Requirements").
- **A NEWS bullet** for anything users can see, under the development heading of
  the branch you target: `# bmm (development version)` on `develop`,
  `# bmm 1.6.0 (development)` on `dev-1.6.0`. Keep it short: what changes for the
  user and what they should do. Do not describe how the fix works.
- **Documentation** for new or changed exported functions (roxygen).
- **For a new model:** read the
  [developer notes](https://venpopov.com/bmm/dev-notes/) and follow the
  one-constructor-per-response-type rule in [AGENTS.md](../AGENTS.md). Add a
  recovery plan, `tests/recovery/<model>.R` (see below).

### Recovery evidence for new models

Recovery evidence is a validation plan for the
[bmmtools](https://github.com/GidonFrischkorn/bmmtools) package. Commit the plan
as `tests/recovery/<model>.R`, a script that builds it with
`bmmtools::validation_plan()` (see `tests/recovery/README.md`). The plan must
cover:

- true against recovered values;
- a hierarchical recovery with several subjects whose parameters vary;
- data simulated with an independent generator, not the model's own `r*`
  function, so that a bug in the density cannot cancel out.

You run the plan locally (see `tests/recovery/README.md`); it takes hours of
sampling. Until a card site exists, attach
the rendered card.html to the PR.

bmmtools 0.3.0 ships defaults for `mixture2p` and `ezdm` only, and those simulate
with bmm's own `r*` function; other models with a bmmtools adapter (for example
`mixture3p`) fall back to it too, so always pass `generator =`. For any other
model, pass the model object and `pars`, not its name:
`bmmtools::validation_plan(bmm::<model>(...), pars = c(...), generator = ..., seed = 2026)`.
The plan is not run by `R CMD check` or CI.

### Model fits for the website articles

The articles load cached fits through `bmm(..., file = )`. Do not commit new
or refitted fits. The `file-size` check fails a PR when any of its commits adds
a file larger than 1 MiB, or a new version of one, even if a later commit
removes the file again, because a merge commit keeps every commit of the PR in
the history. A squash merge would not, but the check cannot know which one will
be used ([ADR 0003](../.dev/decisions/0003-model-fit-files.md)). A file whose
content is already in the history of `develop` or of the PR's base does not
count, wherever it sits. A changed version of a tracked large file counts as
new. If the check names a commit of yours, rewrite the branch without that
file, for example with an interactive rebase.

- Load a fit with `bmm(..., file = "fits/<name>")` or
  `readRDS(file = "fits/<name>.rds")`, with the path written as a string
  literal. The folder `vignettes/articles/fits/` is gitignored. Fits that are
  already tracked in `assets/` stay there until their article is refitted.
  Release fits work only for articles directly in `vignettes/articles/`. The
  build downloads them to `vignettes/articles/fits/`, and an article in a
  subfolder looks for `fits/` next to itself.
- Say in the PR where we can get the fit. A maintainer uploads it to the
  `article-fits` release, and the website build downloads it from there.
- The website build checks every relative `file = "..."` that an evaluated R
  chunk passes as a string literal to `bmm()`, `brm()`, `readRDS()` or
  `read_rds()`. A `fits/<name>` file must be in the release. Any other file,
  such as `assets/<name>`, must be tracked. Both must be readable. The build
  stops, naming the article and the file, when one is not.
- A string that names a `fits/` path and is not such a literal fails the
  build, because the check cannot read it. This includes
  `readRDS("fits/<name>")` without `file =`, a path built with `paste0()` or
  `file.path()`, and inline code. Paths passed to a writer, such as
  `saveRDS()`, `write.csv()`, `ggsave()` or `png()`, are ignored. A path that
  reaches a reader only through a variable or a computed value is not seen
  unless the string itself names `fits/`. Sourced code and child documents
  are not checked.
- Articles that pkgdown does not build (names starting with `_`,
  `vignettes/tutorials/`) are not checked.
- Chunks are found as knitr finds them. A chunk opens at any line like
  `` ```{r} ``, also inside a longer fenced block that only shows it, because
  knitr runs it there. `#| ` options count only at the top of the chunk.
- Upload the fit before the PR is merged into a branch that builds the site.
- Give a refitted fit a new file name instead of replacing the old asset. The
  website builds from `develop` and reads whatever the release holds at that
  moment.

To build the articles locally, download the fits first:

```sh
gh release download article-fits --dir vignettes/articles/fits --skip-existing
```

## Deprecating and removing features

Users should know how long their scripts keep working. From bmm 1.5.0 on,
three rules apply:

1. **Deprecate before you remove.** A deprecated function or argument keeps
   working, with a warning, for at least one minor release. Something
   deprecated in 1.5.0 still works in every 1.5.x release and is removed in
   1.6.0 at the earliest.
2. **List it in NEWS.** Every deprecation and every removal gets a bullet under
   the heading `### Deprecated functions and arguments`.
3. **Announce breaking changes.** A change is breaking when a script that ran
   before now fails, has to be edited, or gives different results for the same
   data and model. Removing a deprecated form counts. A bug fix that corrects
   results is not; it goes in NEWS under the bug fixes. We announce a breaking
   change in a GitHub issue before the release that contains it, with the old
   and the new code side by side (use the "Breaking change" issue template).

The rules cover what users type or read: exported functions and their
arguments, model parameter names, documented defaults, and the data layout a
model expects. Internal functions are not covered. Fits saved by an earlier
release are covered: they keep working with `update()` and the extraction
functions, through `restructure()` where they need converting, or the NEWS
bullet says they must be refitted.

A changed default has no old form to warn about. Announce it as a breaking
change (rule 3) and say in NEWS how to get the old default back. For a changed
data layout, accept both layouts for the window and give the old and the new
column names as `what` and `with`.

If your PR deprecates something:

- Deprecate in a minor release (1.5.0, 1.6.0), not in a patch release, so that
  users meet the warning before the release that may remove the form.
- Keep the old form working by passing it on to the new one. If the new form
  gives a different result, say so in `details`. Keep a deprecated function
  exported and mark its help page as deprecated (title `Deprecated: use new()`,
  `@keywords deprecated`). If a call gives both the old and the new argument,
  stop and name both.
- Warn with `warn_deprecated(what, with, since = "1.5.0")` where the old form
  is detected. `what` and `with` are the old and the new code as users would
  type it; pass `with = NULL` when nothing replaces the old form. The message
  is fixed: it names the release that deprecated the form, the earliest release
  that may remove it (the next minor release; set `removed_in` only to name a
  later one) and the replacement, if there is one. Use `details` only for
  something users would otherwise miss, such as a changed meaning.
- Call `warn_deprecated()` once per user call. If one call can meet many old
  forms, such as prior rows, collect them first and pass them as vectors, so
  users get one warning and not one per row. The message lists the first 20
  forms and counts the rest.
- Add the NEWS bullet under `### Deprecated functions and arguments` in the
  section of the release that deprecates the form. On a `dev-*` branch that is
  the branch's own section ([ADR 0002](../.dev/decisions/0002-release-integration-branches.md)).
- Test that the old form still works and warns.

When the window has passed, you may remove the old form. Remove its tests with
it, list the removal under the same heading in the release that removes it, and
announce it in an issue before that release (rule 3). Forms deprecated before
1.5.0 may be removed in any later release once the removal is announced the
same way.

## Review and merging

How many reviews a PR needs depends on what it touches.

| Tier | What it covers | Reviews |
|---|---|---|
| **Low-risk** | Typos, docs only, CI config, tests only, NEWS only | None. May merge on green CI by anyone with write access. |
| **Standard** | Everything else, including bug fixes and new models | One approving review from someone other than the author. |
| **Pipeline or exported API** | The shared pipeline files (`R/bmm.R`, `R/bmmformula.R`, `R/helpers-model.R`, `R/helpers-data.R`, `R/helpers-prior.R`, `R/helpers-postprocess.R`, `R/helpers-inits.R`, `R/update.R`), the signature or return value of an exported function, or deprecating or removing an exported function or argument | Two approving reviews. |

A **new model** is a standard PR with one addition: the reviewer works through
the reviewer block in the PR template. They read the likelihood against the
cited source, read the recovery plan and inspect its card, and check the default priors on the
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
[developer notes](https://venpopov.com/bmm/dev-notes/) explain which file
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
