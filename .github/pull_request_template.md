## What and why
<!-- Closes #… . What changes for users, in two or three sentences. -->

## Checks
- [ ] `devtools::test()` passes locally
- [ ] `devtools::check()`: 0 errors, 0 warnings (notes explained below)
- [ ] NEWS bullet added (or: not user-visible)
- [ ] Documentation regenerated with `devtools::document()`

## Scope
- [ ] Pipeline or exported API: touches the shared pipeline or an exported function's signature/return value (two reviews)
- [ ] Low-risk: docs, tests, CI, NEWS only (may merge on green CI)
- [ ] Targets the right branch: `develop` for fixes and the next release, `dev-x.y.0` for later features

## New model only
- [ ] Recovery plan in `tests/recovery/<model>.R` (true vs recovered, hierarchical, independent generator so that a bug in the density cannot cancel out); card.html attached to the PR, or linked once a card site exists
- [ ] One constructor per response type

<details><summary>Reviewer (new models)</summary>

- [ ] Likelihood checked against the cited source
- [ ] Recovery plan read and card inspected (or plan run)
- [ ] Default priors checked on the natural scale
</details>

## Review notes
<!-- What a reviewer should read first; open questions; what was not run. -->
<!-- If the maintainer merges without a GitHub approval: summarise the independent review (what was checked, what changed after it). -->
