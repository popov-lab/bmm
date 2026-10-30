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
- [ ] Recovery plan in `tests/recovery/<model>.R` (true vs recovered, hierarchical, independent generator); card link below when available
- [ ] One constructor per response type

<details><summary>Reviewer (new models)</summary>

- [ ] Likelihood checked against the cited source
- [ ] Plan file runs (or card inspected)
- [ ] Default priors checked on the natural scale
</details>

## Review notes
<!-- What a reviewer should read first; open questions; what was not run. -->
<!-- If the maintainer merges without a GitHub approval: summarise the independent review (what was checked, what changed after it). -->
