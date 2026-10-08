# 0001: Release scope for 1.4.0, 1.5.0 and 1.6.0

- **Status:** accepted
- **Date:** 2026-10-04
- **Replaces:** the release plan of 2026-09-26 (1.5.0 = MPT + circular mixtures)

## Context

Three dated events set the release deadlines. Workshops teach against the
package as installed from CRAN, so the CRAN date is the deadline, not the merge
date.

| Event | Date | Release |
|---|---|---|
| Psychonomics workshop (SDT, M3, DDM) | November 2026 | 1.4.0, CRAN submission 26–29 Oct 2026 |
| TeaP pre-conference workshop | late March 2027 | 1.5.0, feature freeze 15 Jan 2027, CRAN mid-February |
| MathPsych developer talk on the extension API | summer 2027 | 1.6.0, on CRAN at least four weeks before the talk |

As of 2026-10-04, all feature branches have been brought up to date with
`develop` (at most 10 commits behind). The work that remains is review and
validation, not porting. The MPT stack (#484–#489) is open for review.

## Decision

### 1.4.0: signal detection, MPT, usability

- Already on `develop`: the SDT suite, `bmm_setup()`, the website move, the
  renames to `prior_info()`, `parameter_info()` and `bmm_models()`, the
  extraction helpers and `posterior_epred()`.
- MPT models (#255). The package code (#484–#488) must be merged by
  **20 Oct 2026**. The article (#489) is website-only and can follow after the
  CRAN submission. If #484–#487 miss 20 Oct, MPT moves to 1.5.0 and 1.4.0
  ships without it.

### 1.5.0: visual working memory and perception

- Circular mixture models rebuilt as custom families (#434), with variable
  precision and capacity versions. The lifecycle policy (#441) is merged first,
  because #434 moves user priors from `nlpar =` to `dpar =`.
- Change detection (#136), built on #434.
- Target confusability competition (TCC), a new model. Committed for 1.5.0.
- Psychometric functions, a new model.
- Smaller items: index-safe random generation (#444), `report_methods()`
  (#432), `restructure()` for old fits (#458), easystats stages 1 and 2
  (#477, #478).

### 1.6.0: joint models and extensions

- Extension API (#442), a simulation generic (new issue) and parameter-name
  clashes (#443).
- Multivariate models (#394).
- Racing models: LBA (#349), LNR (#352), RDM (#353), and cswald
  non-decision-time variability (#387, PR #406).
- easystats stage 3 (#479).

### Future (issues filed, not scheduled)

Valuation models (discounting, cumulative prospect theory), a check of a
Gaussian choice rule for m3, general recognition theory, latent factor
structure on model parameters, and reliability of person-level parameters.

The m3 Gaussian choice rule is filed as a check, not a feature. Adding it needs
a custom-family code path for m3 and a re-aggregation layer for all
post-processing, and the quadrature likelihood costs about 90 times more per
gradient. The check asks first whether `c` and `a` estimates change at all
under the alternative rule in realistic designs.

### Not planned

Models without a tractable likelihood (e.g. diffusion models for conflict
tasks), models with recursion over trial sequences (reinforcement learning,
trial-by-trial category learning), free-recall architectures, and
simulation-based likelihoods.

## Rationale

- MPT is finished and validated; shipping it in 1.4.0 puts it on CRAN before
  TeaP, so the TeaP workshop can teach it from an installed release.
- 1.5.0 groups the visual working memory models. TCC is the main rival account
  to the mixture models, so 1.5.0 lets users compare all major accounts of
  continuous reproduction on the same data. Psychometric functions extend the
  perception side next to the SDT suite.
- Racing models move with the multivariate models to 1.6.0. Together they
  support joint models of memory and response-time tasks.

## Consequences

- 1.5.0 carries about 15–17 pull requests plus TCC as new work, in roughly nine
  working weeks between the 1.4.0 submission and the freeze.
- The circular-mixture validation (pointwise log-likelihood equivalence,
  parameter recovery, speed comparison with `brms::mixture()`) is on the
  critical path for both #434 and #136.
- TCC needs an early measurement of likelihood cost per gradient. Its product
  over all response options at every quadrature node is the main risk.
- Open: what yields first if 1.5.0 runs short. Proposed order of priority:
  #434, TCC, #136, psychometric functions, then the smaller items.
