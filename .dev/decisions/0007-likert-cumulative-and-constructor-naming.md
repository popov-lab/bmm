# 0007: Cumulative response models for Likert items, and a naming rule for constructors

- **Status:** accepted 2026-11-05
- **Date:** 2026-11-05
- **Amends:** [0001](0001-release-scope-1.4-to-1.6.md)
- **Tracking issue:** #662

## Context

bmm fits measurement models for experimental tasks. Questionnaire data from Likert-scale items are the most
common data in psychology and are fitted, when fitted properly, as cumulative ordinal models (brms
`cumulative()`); what those fits lack is the response-process components the psychometric literature has
established: careless and insufficient-effort responding as a person-by-item mixture with a content-free
category distribution (Ulitzsch et al., 2022), acquiescence as a shift or as a mixture over affirmative
categories (Plieninger & Heck, 2018), and extreme or midpoint response style as person-specific threshold
dispersion (Johnson, 2003; Lubbe & Schuster, 2020). bmm already has the pieces: `sdt_rating()`'s threshold
builders are generic in (criterion, spacing, deltas, K), its noise registry gives normal, logistic and Gumbel
distributions, and the psychometric function model (#507) is the K = 2 case with the latent location given by
a stimulus intensity.

brms's native `cumulative()` cannot host bmm's threshold structures or the mixtures (it owns the ordered
intercepts and refuses `0 +`, and `mixture(cumulative(), categorical())` errors in brms 2.23), and the
multinomial per-category route costs 7.8–9.2 µs per row (0006). A trial-level custom family is the only route
that evaluates one category per row and takes the mixtures as a `log_sum_exp`.

The constructor on #507 was called `psychometric()`. "Psychometric function" is the literature's term, but the
bare word names a field that, with this record, bmm starts to serve: a model called `psychometric` would sit
beside questionnaire models, reliability and factor structure (#394).

## Decision

Decided by GF on 2026-11-03 (names, releases) and 2026-11-05 (the rest):

- **`likert_cumulative()`** is added in 1.6.0 on `dev-1.6.0`: a custom family on trial-level rows (one row per
  person × item, `loop = TRUE`; brms `cumulative()` stays the parity reference), K − 1 thresholds from the shared
  builders centred on 0 with `sdt_rating()`'s `threshold_type` set and no `criterion` (default `log_distance`,
  which matches brms's K − 1 free thresholds), the latent location in one parameter `activation` that carries
  person and item terms, a noise distribution from the shared registry. Data are long rows: `response` as
  integer 1..K with `n_categories` or as an ordered factor, an `item` column, an optional `keying` column (±1,
  constant within item); raw categories are never recoded; K ≤ 2 is refused (binary rows are
  `psychometric_function()` or brms `bernoulli()`).
- **No item discrimination** in the first version. `spacing ~ 0 + item` is item-specific threshold spread
  (a location-scale model on the thresholds), not Samejima's discrimination, which multiplies the trait; the
  documentation says so. An explicit discrimination waits for the factor structure of #394.
- **Default priors are data-free**: the latent scale is the noise SD. `activation` main normal(0, 1);
  `spacing` as `sdt_rating()`; `log_distance` deltas main normal(log(2·log(K − 1)/(K − 2)), 0.5), the first
  K-dependent default in bmm, chosen by prior-predictive category probabilities at K ∈ {4, 5, 7} so that
  `log_distance` and the two-parameter types imply the same marginal categories and few degenerate items.
- **Process components** are family variants added in a second PR, after an identification study:
  `careless` as a logit-linked rate with a formula and a uniform category draw (an estimated category-preference
  simplex follows in a third PR); `acquiescence` as a logit-linked rate whose affirmative response is drawn
  from the person's own cumulative distribution truncated to the agree categories (the only form that keeps
  Plieninger & Heck's property that both agree categories gain at every trait level on a cumulative kernel);
  the mixtures nest careless outermost. `acquiescence` is refused without both keyings; `careless` warns below
  the identification conditions the study measures. The acquiescence *shift* is a formula, not a parameter:
  `activation ~ 0 + item + (0 + keying | id) + (1 | id)`.
- **Shared helpers leave the SDT files**: threshold builders, the ordinal category kernel and the noise
  registry get neutral names (`R/helpers-thresholds.R`, `threshold_functions.stan`, `ordinal_functions.stan`,
  `.noise_dists`), with the `.sdt_*` names kept as aliases while any branch reads them.
- **`psychometric()` becomes `psychometric_function()`** before it reaches CRAN (1.5.0); `psychometric_threshold()`
  keeps its name. The two constructors stay separate: at K = 2 the ordinal kernel's gumbel_min is the
  psychometric gumbel_max, and the guess/lapse mixture is asymmetric.
- **Vocabulary**: domain "Personality and individual differences", task "Questionnaire and rating-scale
  responses"; `psychometric_function()` stays under "Perception". The article uses `psychTools::bfi` from
  Suggests.
- **Naming rule**: a constructor names a model or a response format, never a field. The field appears in the
  model vocabulary and in the pkgdown groups ("Questionnaires and individual differences"). The prefix
  `likert_` is a family: a later tree model is `likert_tree()`.

## Consequences

- 1.6.0 gains one model and a helper lift that touches `sdt_rating`, `sdt_ranking`, `sdt_cdp` and
  `psychometric_function` through aliases; parity tests pin the lift (log-likelihoods unchanged to 1e-12).
- The careless and acquiescence components ship with warnings and one refusal derived from a measured
  identification study, not from the literature's rules of thumb alone.
- The K-dependent delta prior is a precedent: a default prior may depend on a data property the constructor
  already requires (`n_categories`), never on the observed responses.
- Open: explicit item discrimination (GRM), IRTree models, factor structure on person parameters (#394), a
  person scaling of item-specific free thresholds (Lubbe & Schuster's full model), an `acquiescence = "centre"`
  variant with a trait-free extremity, a `loop = FALSE` kernel under 0006's 2× rule.
