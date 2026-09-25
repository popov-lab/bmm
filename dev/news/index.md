# Changelog

## bmm (development version)

#### New models

- Add **Yes/No Signal Detection Theory** (`sdt_yn`) for detection and
  old/new recognition tasks with aggregated response counts. It
  estimates sensitivity (`d`) and response bias (`criterion`), and
  optionally the unequal-variance ratio (`sdratio`). Also adds
  [`dsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md),
  [`rsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md),
  [`sdt_d()`](https://venpopov.com/bmm/dev/reference/SDTdist.md) and
  [`sdt_criterion()`](https://venpopov.com/bmm/dev/reference/SDTdist.md).
  See [`?sdt_yn`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) for
  the parameters, links, default priors and the designs that identify
  `sdratio`. Thanks to
  [@GidonFrischkorn](https://github.com/GidonFrischkorn)

#### New datasets

- Add **`broeder_schuetz_2009_e3`**, binary old/new recognition data
  from Broeder & Schuetz (2009, Exp. 3), with five base-rate conditions
  from 40 subjects. See
  [`?broeder_schuetz_2009_e3`](https://venpopov.com/bmm/dev/reference/broeder_schuetz_2009_e3.md).

#### New features

- `bmm(file_refit = "on_change")` is now implemented and no longer warns
  and falls back to `"never"`. The cached fit saved under `file` is
  returned only while the Stan code, the Stan data, the factor levels of
  the model variables and the algorithm are unchanged; any change
  refits. The comparison happens where `brms` makes it — after the bmm
  configuration pipeline has produced the Stan code and data, before
  compilation — so a cache hit costs one run of the pipeline plus
  [`standata()`](https://venpopov.com/bmm/dev/reference/standata.bmmformula.md)
  and
  [`stancode()`](https://venpopov.com/bmm/dev/reference/stancode.bmmformula.md):
  about 0.4 s rather than 0.02 s for an **sdm** model of
  `oberauer_lin_2017`, far less than compiling and sampling but not
  free. As in `brms`, only those four things are compared, so sampler
  settings do not force a refit — `control = list(adapt_delta = )` in
  particular, and also `iter`, `warmup`, `chains`, `seed`, `init` and
  `save_pars`: rerunning with a higher `adapt_delta` after divergent
  transitions, or with `save_pars(all = TRUE)` for `loo()`, returns the
  cached fit unchanged
  ([\#411](https://github.com/popov-lab/bmm/issues/411)).
- `bmm_options(file_refit = )` accepts the same values as
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md). It
  previously required a logical, so the string forms could be set only
  through `options(bmm.file_refit = )`
  ([\#411](https://github.com/popov-lab/bmm/issues/411)).
- [`update()`](https://rdrr.io/r/stats/update.html) gained the `file`
  and `file_compress` arguments. `file` writes the updated fit in the
  same order [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)
  does, i.e. after the bmm postprocessing
  ([\#411](https://github.com/popov-lab/bmm/issues/411)).
- New function
  [`report_priors()`](https://venpopov.com/bmm/dev/reference/report_priors.md)
  reports, for each parameter of a fitted model, its link function, the
  prior actually used on the sampling scale, and whether it was a bmm
  default, a brms default, or user-specified; parameters left with
  improper flat priors are flagged. `format = "text"` produces
  methods-section-ready sentences. Parameters that exist only because
  the family machinery requires them are omitted, so the report lists
  the model’s own vocabulary: the `mu` brms forces on the response-time
  custom families, the `mu2`/`kappa2` second component of a
  [`brms::mixture()`](https://paulbuerkner.com/brms/reference/mixture.html)
  family, and the `theta2` mixture-weight reference the sampler holds at
  zero ([\#391](https://github.com/popov-lab/bmm/issues/391)).
- New function **native_parameters()** returns posterior draws of the
  model parameters on their native scale, evaluated over a grid of
  predictor values. Because draws are returned rather than summaries,
  contrasts between conditions are ordinary arithmetic on the draws. The
  inverse link transformation is applied to the draws before any summary
  is computed, so quantiles and credible intervals are exact. Mixture
  weights are transformed jointly through their softmax
  ([\#388](https://github.com/popov-lab/bmm/issues/388)).
- New developer generic **native_transform()** defines how a model maps
  its parameters from the sampling scale to the native scale. The
  default method covers every transformation that can be expressed
  through an elementwise `links` declaration in a `.model_*()`
  constructor, or through a single softmax group that is active on every
  row of the data, so most new models are supported without writing a
  method.
  [`native_transform.non_targets()`](https://venpopov.com/bmm/dev/reference/native_transform.md)
  ships as the worked example of a design-dependent method
  ([\#388](https://github.com/popov-lab/bmm/issues/388)).
- New function
  [`bmm_data_check()`](https://venpopov.com/bmm/dev/reference/bmm_data_check.md)
  prints a human-readable pre-fit data report. It runs the same
  validation pipeline as
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) without
  compiling the model, captures all errors, warnings and messages, and
  summarizes the response variables (with the coding the model expects),
  the data columns each parameter formula uses (including factor
  coding), the number of observations per design cell, and
  model-specific diagnostics for common data mistakes — e.g. circular
  responses in degrees or coded on \[0, 2\*pi), and misplaced `NA`
  values in `nt_features`/`nt_distances` for set size varying designs.
  Developers can extend the model-specific diagnostics by adding
  [`data_check_findings()`](https://venpopov.com/bmm/dev/reference/data_check_findings.md)
  methods, building each finding with the new
  [`data_check_finding()`](https://venpopov.com/bmm/dev/reference/data_check_finding.md)
  constructor ([\#389](https://github.com/popov-lab/bmm/issues/389)).
- Random-effects standard deviations now get domain-informed default
  priors instead of the `student_t(3, 0, 2.5)` default of `brms`, which
  on a log link lets individual parameters vary by a factor of 12 around
  the group value. Every model parameter declares an `sd` entry in its
  `default_priors` (next to `main` and `effects`), applied as a blanket
  `class = "sd"` prior whenever the parameter has random effects. The
  priors are `exponential()` on the link scale: rate 1 for
  memory-strength, mixing-weight and identity-linked drift parameters,
  rate 2 for boundary, non-decision time, start point, contaminant and
  log-linked drift parameters, rate 4 for the circular bias `mu`/`mu1`.
  Override by addressing the parameter with `dpar` or `nlpar`,
  e.g. `set_prior("exponential(2)", class = "sd", nlpar = "kappa")`
  ([\#342](https://github.com/popov-lab/bmm/issues/342)).
- Correlations among random effects now get an `lkj(2)` default prior
  instead of the uniform `lkj(1)` of `brms`. For two correlated effects
  this lowers the prior probability of a correlation beyond ±.9 from 10%
  to 1.45%, so estimated correlations shrink slightly towards zero. The
  default applies only to models that estimate a correlation matrix, not
  to `(1 | ID)` or `(x || ID)`. To return to the previous behaviour,
  pass `prior = set_prior("lkj(1)", class = "cor")`.
  [`report_priors()`](https://venpopov.com/bmm/dev/reference/report_priors.md)
  shows the prior as `lkj(2)`, class `cor`, as written in `set_prior()`
  ([\#417](https://github.com/popov-lab/bmm/issues/417)).
- [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) now starts
  every model from tight initial values for its random effects, and
  **mixture2p**, **mixture3p** and **imm** start their population-level
  parameters inside the central 50% of their default priors, wherever
  the predictors can hold the parameter at one value there — a design
  that cannot, such as `~ 0 + poly(x, 2)`, starts as close as it can;
  these models and **m3** previously used `init = 1`. This reduces the
  `lkj_corr_cholesky_lpdf: Random variable[k] is 0`,
  `von_mises_lpdf: Scale parameter is inf` and `Rejecting initial value`
  messages at the start of warmup. Posteriors are unaffected. A
  user-supplied `init` still replaces the bmm default.
- [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) and
  [`update()`](https://rdrr.io/r/stats/update.html) now start Stan’s
  step-size search at 0.01 instead of 1, which further reduces the
  `lkj_corr_cholesky_lpdf: Random variable[k] is 0` and
  `von_mises_lpdf: Scale parameter is inf` messages at the start of
  warmup. The adapted step size and the posterior do not change. A
  `step_size` (or `stepsize`) in your own `control` list wins, and
  `bmm_options(step_size = FALSE)` restores Stan’s default for new fits;
  [`update()`](https://rdrr.io/r/stats/update.html) keeps the step size
  a fit was run with.

#### Bug fixes

- A name in a model’s `links` argument that names no parameter of that
  model is now an error instead of being added:
  `sdm(resp_error = "y", links = list(kapa = "identity"))` reported a
  parameter `kapa` while `kappa` kept its default link. A name one edit
  from a parameter is read as that parameter, with a warning; a link
  allowing values the default excludes warns once. A link the model
  cannot apply is refused too, naming those it can:
  `links = list(bound = "loglog")` used to fail with
  `argument is of length zero`. **sdm**, **mixture2p**, **mixture3p**
  and **imm** take none
  ([\#420](https://github.com/popov-lab/bmm/issues/420)).
- **sdm** no longer fails before sampling when `mu` is predicted without
  an intercept
  (`Initial values for vectors are only specified for b-coefficients, sd and z parameters`).
  **sdm**, **ddm**, **ezdm** and **cswald** no longer fail with the
  rstan backend when a parameter has exactly one slope
  (`no more scalars to read`), nor when their random effects use
  `gr(..., by = )`.
- [`report_priors()`](https://venpopov.com/bmm/dev/reference/report_priors.md)
  no longer fails for a fit whose formula needs `data2`,
  e.g. `(1 | gr(ID, cov = A))` (`Object 'A' was not found in 'data2'`).
- [`report_priors()`](https://venpopov.com/bmm/dev/reference/report_priors.md)
  now names every kind of correlation prior the way `set_prior()`
  documents it, not only the one on group-level effects. A fit with two
  or more `me()` terms was reported as class `Lme` with
  `lkj_corr_cholesky(1)`, brms’s internal spelling, instead of `corme`
  with `lkj(1)`. The same applied to `rescor`, `lncor` and `cortime`.
- [`update()`](https://rdrr.io/r/stats/update.html) with `newdata` or a
  new `formula.` no longer fails before sampling for a model that starts
  from bmm’s initial values (`no more scalars to read` on rstan,
  `Fitting failed` on cmdstanr). The initial values are now built for
  the data and formula the update fits; an `init` passed to
  [`update()`](https://rdrr.io/r/stats/update.html) still wins
  ([\#415](https://github.com/popov-lab/bmm/issues/415)).
- Formulas with `mo()` or `s()` terms no longer start with a partial
  init list, which made cmdstanr print
  `Init values were only set for a subset of parameters`.
  [`extract_parameter_dimensions()`](https://venpopov.com/bmm/dev/reference/extract_parameter_dimensions.md)
  now reads a declaration sized by an element of a data array, such as
  `simplex[Jmo_c[1]]`, and simplex parameters get an initial value.
- `file_refit` no longer accepts a logical that is not a single
  `TRUE`/`FALSE`. `NA`, `logical(0)` and multi-element logicals were
  coerced to `"never"`, so `file_refit = cfg$refit` with a missing or
  `NA` config key silently returned the cached fit instead of erroring
  as it did before `"on_change"` was added
  ([\#411](https://github.com/popov-lab/bmm/issues/411)).
- Under `file_refit = "on_change"`, a cached fit that
  [`restructure()`](https://venpopov.com/bmm/dev/reference/restructure.bmmfit.md)
  cannot bring forward to the current bmm version is now refitted rather
  than aborting [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)
  with *“Unable to restructure the object… Please refit”* — that is the
  one thing `"on_change"` exists to do. A cached fit with no `algorithm`
  field also no longer aborts on `brms`’s bare
  [`stopifnot()`](https://rdrr.io/r/base/stopifnot.html); the comparison
  falls back to the Stan code, Stan data and factor levels
  ([\#411](https://github.com/popov-lab/bmm/issues/411)).
- Fix `file_refit` silently falling back to `"never"` for any
  capitalisation other than all-lowercase. `"Always"` passed validation,
  which lowercases, but the coercion that followed it did not, so
  `bmm(..., file = , file_refit = "Always")` returned the cached fit
  instead of refitting
  ([\#411](https://github.com/popov-lab/bmm/issues/411)).
- [`update()`](https://rdrr.io/r/stats/update.html) no longer returns a
  fit whose `file` field has been dropped, and no longer lets `brms`
  handle the file at all.
  [`brms::update.brmsfit()`](https://paulbuerkner.com/brms/reference/update.brmsfit.html)
  clears `file` and
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  writes it before any bmm postprocessing has run, so
  `update(fit, file = )` stored a plain `brmsfit` that
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) then refused
  to read back (`Object loaded via 'file' is not of class 'bmmfit'`).
  Worse, when the file already existed,
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) read
  it and returned its contents instead of fitting, so
  `update(fit, newdata = , file = )` — rerunning a script that caches
  its fits — silently discarded the update and returned the old fit with
  no message. [`update()`](https://rdrr.io/r/stats/update.html) now
  writes the file itself, after the postprocessing, as
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) does
  ([\#411](https://github.com/popov-lab/bmm/issues/411)).
- Fix `options(bmm.default_priors = FALSE)` making
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) fail for
  every model instead of fitting with flat priors. `set_default_prior()`
  returned `NULL`, which `combine_prior()` indexed unconditionally
  (`second argument must be a list`); `set_default_prior()` now returns
  an empty prior, and `combine_prior()` passes a `NULL` argument through
  ([\#391](https://github.com/popov-lab/bmm/issues/391)).
- [`update()`](https://rdrr.io/r/stats/update.html) now configures the
  likelihood for the threading spec that will actually be used.
  [`brms::update.brmsfit`](https://paulbuerkner.com/brms/reference/update.brmsfit.html)
  falls back to the original fit’s `threads` when the argument is not
  passed, but bmm only inspected the new request, so updating a threaded
  fit without repeating `threads` emitted the serial likelihood chunk
  into threaded Stan code — for the **sdm** model a compile error
  (`Identifier 'COSN' not in scope`), and for any `loop = FALSE` custom
  family a silently mis-sliced likelihood. The fallback follows brms in
  distinguishing an absent `threads` argument from an explicit
  `threads = NULL`, which turns threading off.
- [`update()`](https://rdrr.io/r/stats/update.html) no longer lets a
  global `options(brms.threads = )` reach the likelihood configuration.
  Updating an unthreaded fit under a session-wide threading option
  emitted the sliced likelihood chunk while brms generated serial Stan
  code (`Identifier 'start' not in scope`); the spec of the fit being
  updated now always wins, as it does in brms.
- [`update()`](https://rdrr.io/r/stats/update.html) now re-resolves
  every parameter whose constant the new formula changes.
  [`update()`](https://rdrr.io/r/stats/update.html) never called
  [`check_model()`](https://venpopov.com/bmm/dev/reference/check_model.md),
  so the constant was never resolved again and the original fit’s prior
  overrode the freshly configured one: `update(fit, bmf(..., mu ~ 1))`
  returned a model in which `mu` was still pinned,
  `update(fit, bmf(..., kappa = 5))` one that reported `kappa = 5` while
  Stan estimated `kappa` freely, and `update(fit, bmf(..., mu = 0.5))`
  one that reported `mu = 0.5` while Stan kept `mu` at the original
  value — all three with no error or warning.
- The **ezdm** likelihood no longer assumes that reaction times are
  normally distributed. Because response times are right-skewed, it
  treated each cell as more informative than it is, and posteriors,
  especially for `bound`, came out too narrow. Existing **ezdm** fits
  are not reproduced: in our simulations, credible intervals for
  subject-level `bound` estimates widen by 20% to 75%, and fitting takes
  longer (see
  [`?ezdm_dist`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md)).
  A fit cached with `bmm(file = )` is reused unless
  `file_refit = "on_change"`.
  [`rezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md) no
  longer truncates `mean_rt` at `ndt`, so with very few trials it can
  return `mean_rt <= 0`, which
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) rejects
  ([\#407](https://github.com/popov-lab/bmm/issues/407)).
- Fix numerical failures in **ezdm**. With large drift rates, the
  4-parameter model returned `NaN` in
  [`dezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md),
  [`rezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md),
  `log_lik()` and
  [`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md),
  and both models could return a log-likelihood of `-Inf` when the
  predicted accuracy was very close to 1 but a cell contained errors.
  Near zero drift, the likelihood jumped where the code switched between
  formulas, and the response counts said nothing about the direction of
  drift.
  [`dezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md) now
  rejects counts that are not whole numbers
  ([\#407](https://github.com/popov-lab/bmm/issues/407)).
- The mixing weights `thetat` (**mixture2p**) and `thetat`/`thetant`
  (**mixture3p**) had no `effects` prior, so any regression coefficient
  on them was flat. They now get `normal(0, 0.5)` on the logit/softmax
  scale ([\#305](https://github.com/popov-lab/bmm/issues/305)).
- The circular bias parameters `mu` (**sdm**) and `mu1` (**mixture2p**,
  **mixture3p**, **imm**) had a `student_t(1, 0, 1)` intercept prior,
  uniform over the circle under the `tan_half` link, and no `effects`
  prior, which gave non-reference levels a bimodal prior on the native
  scale. They now get `normal(0, 0.5)` on the intercept (95% of the
  prior bias within 89 degrees of the target), `normal(0, 0.25)` on
  effects (a condition difference with SD 24 degrees) and
  `exponential(4)` on random-effects SDs. `mu`/`mu1` stay fixed at 0 by
  default, so this only affects models that free them
  ([\#306](https://github.com/popov-lab/bmm/issues/306)).
- The **cswald** model no longer rejects extreme response times at high
  drift rates. Its likelihood returned `NaN` there, which the sampler
  reported as *“Log probability evaluates to log(0)”* or as divergent
  transitions, and could bias the posterior for data with long tails.
  Fits of such data change; refit to benefit
  ([\#387](https://github.com/popov-lab/bmm/issues/387)).
- Passing `threads = NULL` while `options(brms.threads = )` is set no
  longer produces Stan code that fails to compile
  (`Identifier 'start' not in scope`) for the **sdm** and **cswald**
  models. `threads = NULL` now turns parallelization off, as it does in
  `brms`.

#### Other changes

- The **cswald** likelihood now evaluates all observations in one call
  instead of one at a time. This makes fitting faster, improves the
  accuracy of the gradients the sampler uses, and adds support for
  within-chain parallelization: `bmm(..., threads = 2)` now works for
  **cswald** as it does for **sdm**. The posterior is unchanged
  ([\#387](https://github.com/popov-lab/bmm/issues/387)).

## bmm 1.3.2

CRAN release: 2026-09-16

#### Changes to default priors

- Recalibrated the default priors for the **m3** activation parameters
  (`a`, `c`) so that the `simple` and `softmax` choice rules imply a
  comparable, broad prior-predictive range of average performance, and
  so that the `softmax` defaults place equal prior means on general
  (`a`) and context (`c`) activation — centering the implied `c - a`
  prior at zero for fair comparisons under cell-means coding. The
  mis-scaled `normal(0, 2)` effect prior on `c` is replaced by the
  shared `normal(0, 0.5)`. Because the `simple` rule requires `c > a` to
  predict accurate recall, its defaults remain asymmetric; direct
  comparisons of context and general activation should use the `softmax`
  choice rule ([\#364](https://github.com/popov-lab/bmm/issues/364)).

#### New features

- The **sdm** model now supports within-chain parallelization via the
  `threads` argument (e.g. `bmm(..., threads = 2)`), reducing fitting
  time by up to ~45% in benchmarks
  ([\#374](https://github.com/popov-lab/bmm/issues/374)).
- Add `softplus` as an opt-in link function for positively-bounded
  parameters, as an alternative to the default `log` link.
  `softplus(x) = log(1 + exp(x))` keeps parameters positive while
  growing linearly for large values, avoiding the numerical blow-up of
  [`exp()`](https://rdrr.io/r/base/Log.html) and giving predictor
  effects an additive (rather than multiplicative) interpretation on the
  natural scale. Enable it per parameter via the model’s `links` list,
  e.g. `m3(...)$links <- list(c = "softplus", a = "softplus")` or
  `ddm(rt, response, links = list(bound = "softplus"))`
  ([\#363](https://github.com/popov-lab/bmm/issues/363)).
- [`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
  can now check every observable of a model’s likelihood, not just the
  primary response. brms’s
  [`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
  only plots the brms `Y` variable, so the **ddm** and **cswald**
  responses and the **ezdm** RT variance and accuracy went unchecked.
  The new `resp_var` argument selects the observable
  (e.g. `pp_check(fit, resp_var = "response")`), including derived ones
  (`"signed_rt"` for the RT models; `"mean_pc"`, the proportion of
  upper-boundary responses, for `ezdm`), and `resp_var = "all"` returns
  a panel of all checks drawn from one shared joint simulation.
  `pp_check_vars(fit)` lists the available checks, the plot type each
  uses by default, and the brms standata slots it reads; model authors
  declare observables via the
  [`pp_observables()`](https://venpopov.com/bmm/dev/reference/pp_observables.md)/[`pp_simulate()`](https://venpopov.com/bmm/dev/reference/pp_observables.md)
  S3 generics. The **ezdm** checks default to `type = "intervals"`,
  since each observation is one design cell and a density overlay of a
  handful of summary statistics is uninformative. Where an observable is
  undefined for some cells (an **ezdm** boundary reached by fewer than
  two responses), observations are dropped only when the *observed*
  value is undefined; undefined *simulated* values are absorbed by
  dropping those posterior draws, so the number of observations checked
  does not depend on `ndraws`
  ([\#401](https://github.com/popov-lab/bmm/issues/401)).

#### Bug fixes

- Fix
  [`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
  for the RT models silently producing a misleading plot with
  `negative_rt = TRUE`: brms forwarded the argument to
  `posterior_predict()` (signed predicted RTs) while the observed
  response times stayed unsigned, so the plot looked like severe misfit.
  [`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
  now checks the `"signed_rt"` observable (with a message) so both
  halves are signed, and errors for models without signed RTs
  ([\#401](https://github.com/popov-lab/bmm/issues/401)).
- Remove the unreachable `dv` argument of the internal **ezdm**
  `posterior_predict` functions and its documentation.
  `posterior_predict(fit, dv = "var_rt")` silently returned `mean_rt`:
  brms forwards `posterior_predict()` dots to `prepare_predictions()`
  only, never to the family’s prediction function, so `dv` was dropped
  without a warning and has never worked. Use
  `pp_check(fit, resp_var = ...)` instead
  ([\#401](https://github.com/popov-lab/bmm/issues/401)).
- Fix `rezdm(version = "4par")` crashing (recycling errors) when some
  but not all simulated cells produced fewer than 2 responses at the
  upper boundary — the upper-boundary branch did not subset its moments
  and non-decision time by the affected cells, while the lower-boundary
  branch did. The bug was invisible whenever every cell shared the same
  parameters (the unindexed vectors are then constant), and is reachable
  through
  [`rezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md)
  directly; fixing it is a prerequisite for the multi-observable
  [`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md),
  which calls
  [`rezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md) with
  per-draw posterior parameters
  ([\#401](https://github.com/popov-lab/bmm/issues/401)).
- Fix the grouped plot type auto-selected by `pp_check(group = )`
  (e.g. `dens_overlay_grouped`) being silently dropped when `type` was
  not supplied, so that brms fell back to the ungrouped type and
  bayesplot warned about an unrecognized `group` argument.
- Fix initial values being set in two places, where the `init` returned
  by
  [`configure_model()`](https://venpopov.com/bmm/dev/reference/configure_model.md)
  was silently overwritten by
  [`create_initfun()`](https://venpopov.com/bmm/dev/reference/create_initfun.md).
  This caused the **m3** model’s intended `init = 0` (needed for stable
  sampling with the `simple` choice rule and an `identity` link) to be
  lost, and left dead `init` code in the **sdm** model.
  [`create_initfun()`](https://venpopov.com/bmm/dev/reference/create_initfun.md)
  is now the single source of truth for initial values, with
  model-specific behaviour expressed through S3 methods
  ([\#375](https://github.com/popov-lab/bmm/issues/375)).
- Fix `.pwald()` returning `NaN`/`-Inf` in the upper tail of the
  shifted-Wald survival function, which propagated to
  [`dcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md)
  (and therefore `log_lik`/`posterior_predict`) for the **cswald** model
  at extreme reaction times. The R-side survival now uses the stable
  `log_diff_exp` form already used by the Stan likelihood
  (`swald_lccdf`)
  ([\#376](https://github.com/popov-lab/bmm/issues/376)).
- Fix [`print()`](https://rdrr.io/r/base/print.html) for model summaries
  selecting regression-coefficient rows by an unanchored substring
  match, so a parameter such as `a` could pull in rows of another
  parameter like `kappa` (e.g. in the **imm** model). Rows are now
  matched on the exact parameter prefix. This also fixes a crash when
  only a single coefficient row is shown
  ([\#379](https://github.com/popov-lab/bmm/issues/379),
  [\#369](https://github.com/popov-lab/bmm/issues/369)).
- [`create_initfun()`](https://venpopov.com/bmm/dev/reference/create_initfun.md)
  now matches Stan parameters to model parameters with a word-boundary
  regex (`(^|_)param(_|$)`) instead of a substring match, preventing
  collisions in models with short parameter names (e.g. `s`, `c`, `a`)
  that are substrings of longer ones (`sim`, `correct`, `activation`);
  the longest (most specific) match is selected when several apply
  ([\#354](https://github.com/popov-lab/bmm/issues/354),
  [\#355](https://github.com/popov-lab/bmm/issues/355)).
- [`create_initfun()`](https://venpopov.com/bmm/dev/reference/create_initfun.md)
  now resolves initialization terms from `nlpars` when a model parameter
  is not a distributional parameter, so models built as non-linear brms
  formulas (e.g. native-multinomial models whose parameters live in
  `bterms$nlpars`) no longer error with
  `no applicable method for 'has_intercept' applied to an object of class "NULL"`
  ([\#362](https://github.com/popov-lab/bmm/issues/362)).
- The **sdm** model now emits its serial likelihood under
  `threading(n, force = TRUE)`, matching the unsliced Stan code brms
  generates in that case. Previously the threaded chunk was emitted and
  the model failed to compile.
- [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) now warns
  when a predictor in the formula shares its name with both a predicted
  parameter and a column in the data. Such a predictor was silently
  treated as a non-linear term (emitted via `nlf()` instead of `lf()`),
  changing the likelihood without any error. Short parameter names (`c`,
  `a`, `s`, `b`) collide naturally with condition codes or columns like
  `accuracy`/`stimulus`
  ([\#378](https://github.com/popov-lab/bmm/issues/378)).

#### Other changes

- Added an internal consistency check in
  [`configure_prior()`](https://venpopov.com/bmm/dev/reference/configure_prior.md):
  if a model’s `fixed_parameters` includes a parameter that its
  [`configure_model()`](https://venpopov.com/bmm/dev/reference/configure_model.md)
  never wires into the formula (neither a dpar nor an nlpar), bmm now
  fails with a clear model-definition error instead of letting a
  malformed `b_Intercept ~ constant()` prior reach `brm()`. This is a
  safety net for model development; it cannot be reached through the
  normal [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)
  interface, where an unrecognized parameter is already caught earlier
  by
  [`check_formula()`](https://venpopov.com/bmm/dev/reference/check_formula.md)
  ([\#377](https://github.com/popov-lab/bmm/issues/377)).
- `print.bmmodel()` now also displays the required response variables
  (one per line, annotated with the expected coding — e.g. radians in
  \[-pi, pi\] for the circular models, seconds and 0/1 responses for the
  trial-wise RT models) and the default parameter links, answering “what
  does my data frame need to look like?” directly at the console
  ([\#392](https://github.com/popov-lab/bmm/issues/392)).
- The package maintainer contact is now Gidon T. Frischkorn, and the
  repository URLs point to `https://github.com/popov-lab/bmm`.

#### Developer-facing changes

- **[`use_model_template()`](https://venpopov.com/bmm/dev/reference/use_model_template.md)
  now scaffolds the current model-specification patterns.** It generates
  a flat `.{model}_defaults` block for unversioned models (like `ddm`)
  or, with the new `versions` argument, a `.{model}_version_table` block
  for versioned models (like `cswald`). The generated constructor spells
  out every field of the model object inline (referencing the
  defaults/version table for `parameters`, `links`, `fixed_parameters`,
  `default_priors`, and `init_ranges`), and versioned aliases validate
  `version` with [`match.arg()`](https://rdrr.io/r/base/match.arg.html)
  ([\#350](https://github.com/popov-lab/bmm/issues/350)).
- **Removed the unused `void_mu` field** from all model definitions and
  the template. It was assigned but never read anywhere — response-mean
  suppression is already handled via `fixed_parameters` and the family’s
  `dpars` ([\#350](https://github.com/popov-lab/bmm/issues/350)).

## bmm 1.3.1

CRAN release: 2026-06-05

#### Bug fixes

- Fix `swald_lccdf()` returning incorrect log-survival probability when
  response time equals non-decision time in the **cswald** model.
  Previously returned `-Inf` instead of `0` (log of survival = 1)
  ([\#348](https://github.com/popov-lab/bmm/issues/348)).

#### Other changes

- The **ddm** model supports both `cmdstanr` and `rstan` backends.
  Previously, `cmdstanr` was required.

## bmm 1.3.0

CRAN release: 2026-03-30

#### New models

- Add the **Diffusion Decision Model** (`ddm`) for speeded
  decision-making tasks with trial-level RT and response data. The model
  estimates drift rate, boundary separation, non-decision time, and
  (optionally) relative starting point. Includes distribution functions
  [`dddm()`](https://venpopov.com/bmm/dev/reference/ddm_dist.md) and
  [`rddm()`](https://venpopov.com/bmm/dev/reference/ddm_dist.md)
  ([\#280](https://github.com/popov-lab/bmm/issues/280)).
- Add the **EZ-Diffusion Model** (`ezdm`) for speeded decision-making
  tasks. The model estimates drift rate, boundary separation, and
  non-decision time from aggregated summary statistics (mean RT,
  variance of RT, accuracy) using the closed-form equations derived by
  Wagenmakers et al. (2007). Supports both 3-parameter (symmetric
  starting point) and 4-parameter (asymmetric starting point) versions
  based on Srivastava et al. (2016). Implements Bayesian hierarchical
  estimation following Chavez & Vandekerckhove (2025). Includes
  distribution functions
  [`dezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md) and
  [`rezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md)
  ([\#281](https://github.com/popov-lab/bmm/issues/281)).
- Add the **Censored Shifted Wald Model** (`cswald`) for choice reaction
  time tasks with two response boundaries. The model estimates drift
  rate, boundary separation, and non-decision time from trial-level RT
  and response data. Implements two versions: **simple** (treats errors
  as censored correct responses, appropriate for high-accuracy tasks)
  and **crisk** (competing risks version with separate accumulators for
  each response, suitable for balanced accuracy). Includes distribution
  functions
  [`dcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md),
  [`pcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md),
  [`qcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md),
  and
  [`rcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md).
  Thanks to [@GidonFrischkorn](https://github.com/GidonFrischkorn)

#### New features

- New S3 method **conditional_effects()** for `bmmfit` objects. Provides
  an intuitive interface for visualizing predictor effects on model
  parameters, with automatic routing between distributional and
  non-linear parameters, inverse link transformations to show parameters
  on their natural scale (`scale = "native"`), softmax handling for
  mixture3p weight parameters, and filtering of internal model variables
  ([\#203](https://github.com/popov-lab/bmm/issues/203)).
- New S3 methods for **emmeans** support on `bmmfit` objects. Users can
  now call `emmeans(fit, ~ condition, dpar = "kappa")` for any bmmodel
  ([\#323](https://github.com/popov-lab/bmm/issues/323)).
- New **pp_check()** method for multinomial models (e.g., `m3`). Since
  [`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
  does not support the multinomial family, `bmm` now provides a custom
  method that compares observed and predicted response proportions in
  the `ppc_bars` style from `bayesplot`. The method supports faceting by
  experimental conditions via `group`, configurable credible intervals
  via `probs`, and population-level predictions via `re_formula = NA`.
  For non-multinomial models,
  [`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
  delegates to
  [`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
  and auto-selects the grouped plot variant when `group` is specified
  ([\#324](https://github.com/popov-lab/bmm/issues/324)).
- New function **parameters()** lists all parameters of a `bmmodel` or
  `bmmfit` object with descriptions, link functions, and fixed values
  ([\#329](https://github.com/popov-lab/bmm/issues/329)).
- New function **extract_stan_blocks()** extracts individual program
  blocks (functions, data, parameters, etc.) from compiled Stan code
  ([\#286](https://github.com/popov-lab/bmm/issues/286)).
- New function **extract_parameter_dimensions()** extracts parameter
  names, dimensions, and types from a Stan parameters block.
- New functions **ezdm_summary_stats()** and **adjust_ezdm_accuracy()**
  to compute and pre-process summary statistics from trial-level RT data
  for the EZ-Diffusion Model
  ([\#291](https://github.com/popov-lab/bmm/issues/291)).
- New functions **flag_contaminant_rts()** and
  **validate_fast_guesses()** for trial-level contamination detection in
  RT data. Identifies fast guesses and attention lapses using mixture
  modeling and provides Bayesian validation of fast guess assumptions
  ([\#307](https://github.com/popov-lab/bmm/issues/307)).
- New function **create_initfun()** creates initialization functions for
  models that benefit from or require initial values for MCMC sampling
  ([\#285](https://github.com/popov-lab/bmm/issues/285)).

#### Documentation

- New online
  [article](https://venpopov.com/bmm/dev/articles/bmm_ddm.html) to
  accompany the **ddm** model
- New online
  [article](https://venpopov.com/bmm/dev/articles/bmm_ezdm.html) to
  accompany the **ezdm** model
- New online
  [article](https://venpopov.com/bmm/dev/articles/bmm_cswald.html) to
  accompany the **cswald** model
- New online
  [article](https://venpopov.com/bmm/dev/articles/bmm_rt_contamination.html)
  on pre-processing and contamination detection for reaction time data

#### Other changes

- Improved **rm3()** random generation function for the M3 model
  ([\#279](https://github.com/popov-lab/bmm/issues/279)).
- Removed magrittr dependency; replaced `%>%` with the native pipe `|>`
  ([\#341](https://github.com/popov-lab/bmm/issues/341)).
- Minimum R version is now 4.1.0.

## bmm 1.2.0

CRAN release: 2025-07-24

#### New models

- Add the Memory Measurement Model (Oberauer & Lewandowsky, 2019) and
  its generalization as the Multinomial Measurement Model for
  categorical decision tasks as new model class **m3** with three
  versions: simple span (**ss**), complex span (**cs**), and **custom**.
  For details, see the
  [article](https://venpopov.com/bmm/articles/bmm_m3.html) on the `bmm`
  website ([\#237](https://github.com/popov-lab/bmm/issues/237)). Thanks
  to [@GidonFrischkorn](https://github.com/GidonFrischkorn) and
  [@chenyu-psy](https://github.com/chenyu-psy)

#### New features

- Updates to the `bmf2bf` S3 methods for more flexible translation of
  `bmmformulas` into `brmsformulas`
  ([\#227](https://github.com/popov-lab/bmm/issues/227)).
- New function **apply_links** adds link functions to all non-linear
  formulas in a **bmmformula** object.
- New example data set **oberauer_lewandowsky_2019_e1** for exploring
  the **m3** model.
- The `file_refit` argument of the `bmm` function now accepts character
  strings like `brms`. A warning is given when “on_change” is specified,
  as this is not currently implemented for `bmmodels`
  ([\#228](https://github.com/popov-lab/bmm/issues/228)). (The warning
  was removed and “on_change” implemented in the development version,
  see [\#411](https://github.com/popov-lab/bmm/issues/411).)
- New function **rejection_sampling**

#### Bug fixes

- Fix conflict in setting default priors when model parameters were
  transformed in a non-linear formula
  ([\#232](https://github.com/popov-lab/bmm/issues/232)).
- Allow a NULL formula (`formula(NULL)`) to be added to a bmmformula for
  consistentcy with brms
  ([\#264](https://github.com/popov-lab/bmm/issues/264))
- Improve error messages when attempting to construct bmmformulas
  without a left-hand-side variable

#### Documentation

- Add documentation to the [continuous reproduction
  task](https://venpopov.com/bmm/articles/bmm_vwm_crt.html) article for
  pre-processing half-circular stimulus spaces when using `bmmodels` of
  the `circular` model class
  ([\#229](https://github.com/popov-lab/bmm/issues/229),
  [\#233](https://github.com/popov-lab/bmm/issues/233)).
- New online [article](https://venpopov.com/bmm/articles/bmm_m3.html) to
  accompany the m3 model

#### Other changes

- vectorize [`k2sd()`](https://venpopov.com/bmm/dev/reference/k2sd.md)
  function for improved performance
- various internal refactorings
  ([\#246](https://github.com/popov-lab/bmm/issues/246),
  [\#242](https://github.com/popov-lab/bmm/issues/242))
- dplyr, magrittr and tidyr dependencies are now optional
  ([\#240](https://github.com/popov-lab/bmm/issues/240))
- new contributor - Chenyun Li (chenyu-psy) for his work on the m3 model

## bmm 1.0.0

First version of the package on published on CRAN!

#### New features

- you can now specify to save the **bmmfit** object generated by
  **bmm()** to a file with the **file** argument, similarly to
  **brms::brm()** ([\#190](https://github.com/popov-lab/bmm/issues/190))
- the parameterization of the **imm** was adapted to accurately reflect
  the model as implemented by Oberauer et al. (2017)
- prepare package for CRAN submission

#### Bug fixes

- fix incorrect specification of default priors when only an interaction
  is specified ([\#201](https://github.com/popov-lab/bmm/issues/201))
- the random generation function for the **mixture3p** and **imm**
  returned incorrect samples for some rare parameter combinations, this
  has now been fixed, so that the functions now return correct samples
  for all parameter combinations.

#### Deprecated functions and arguments

- BREAKING CHANGE: the arguments for the distribution functions of the
  **mixture2p** and **mixture3p** model have been change to match the
  snake_case coding scheme. Instead of **pMem** and **pNT** these are
  now **p_mem** and **p_nt**. The old names are deprecated and are no
  longer supported

## bmm 0.5.1

#### Bug fixes

- fix the display of the model call in the summary method for bmm models

## bmm 0.5.0

#### New features

- add a **summary()** method for **bmmfit** objects
  ([\#144](https://github.com/popov-lab/bmm/issues/144))
- add a global option **bmm.summary_backend** to control the backend
  used for the **summary()** method (choices are *“bmm”* and *“brms”*)
- function **restructure()** now allows to apply methods introduced in
  newer **bmm** versions to **bmmfit** objects created by older **bmm**
  versions
- you can now specify any model parameter to be a constant by using an
  equal sign in the **bmmformula**
  ([\#142](https://github.com/popov-lab/bmm/issues/142))
- you can now choose to estimate parameters that are fixed to a constant
  by default for all models
  ([\#145](https://github.com/popov-lab/bmm/issues/145))
- default priors for all models are now specified via the
  **configure_prior()** S3 method
  ([\#145](https://github.com/popov-lab/bmm/issues/145))
- **cmdstanr** will be used as the default backend for **brms** if the
  user has it installed
  ([\#145](https://github.com/popov-lab/bmm/issues/145))
- various updates to the documentation and data sets

#### Documentation

- two new online articles that [introduce the **bmmformula**
  syntax](https://venpopov.com/bmm/articles/bmm_bmmformula.html) and
  explain [how to extract information from
  **bmmodels**](https://venpopov.com/bmm/articles/bmm_extract_info.html)
  such as the generated Stan code and Stan data for each model

#### Bug fixes

- fix a bug preventing the **sort_data** check from being executed
  ([\#72](https://github.com/popov-lab/bmm/issues/72))
- fix bugs with the **summary()** function not displaying implicit
  parameters ([\#152](https://github.com/popov-lab/bmm/issues/152)) and
  not working properly with some hierarchical designs
  ([\#173](https://github.com/popov-lab/bmm/issues/173))
- fix a bug in which the **sort_data** check occurred in cases where it
  shouldn’t ([\#158](https://github.com/popov-lab/bmm/issues/158))

#### Deprecated functions and arguments

- BREAKING CHANGE: remove **get_model_prior(), get_stancode() and
  get_standata()**. Due to [recent
  changes](https://github.com/paul-buerkner/brms/pull/1604) in *brms*
  version 2.21.0, you can now use the *brms* functions
  **default_prior**, **stancode** and **standata** directly with *bmm*
  models.
- the function **fit_model()** is deprecated in favor of **bmm()** and
  will be removed in a future version
  ([\#163](https://github.com/popov-lab/bmm/issues/163))
- the argument **setsize** for the **mixture3p** and **IMM** models is
  now called **set_size** for consistency. The old argument name is
  deprecated and will be removed in a future version
  ([\#163](https://github.com/popov-lab/bmm/issues/163))
- the distributions functions for the imm model are renamed from
  **dIMM**, **pIMM**, **rIMM** and **qIMM** to **dimm**, **pimm**,
  **rimm** and **qimm**
  ([\#163](https://github.com/popov-lab/bmm/issues/163))
- the argument parallel for the **bmm()** function is deprecated and
  will be removed in a future version. Use **cores** instead, as for
  **brms::brm()** ([\#163](https://github.com/popov-lab/bmm/issues/163))
- the models **IMMfull()**, **IMMabc()** and **IMMbsc()** are now called
  via **imm()**, **imm(version = “abc”)** or **imm(version = “bsc”)**.
  The old names are deprecated and will be removed in a future version
  ([\#163](https://github.com/popov-lab/bmm/issues/163))
- the **sdmSimple()** model is now called **sdm()**. The old name is
  deprecated and will be removed in a future version
  ([\#163](https://github.com/popov-lab/bmm/issues/163))

#### Other changes

- **bmm** now requires the at least version 2.21.0 of **brms**.

## bmm 0.4.0

#### New features

- add a check for the **sdmSimple** model if the data is sorted by
  predictors. This leads to much faster sampling. The user can control
  the default behavior with the **sort_data** argument
  ([\#72](https://github.com/popov-lab/bmm/issues/72))
- the **mixture3p** and **IMM** models now require that the intercept
  must be suppressed when set size is used as a predictor
  ([\#96](https://github.com/popov-lab/bmm/issues/96)).
- add postprocessing methods for **sdmSimple** to allow the use of
  **pp_check()**, **conditional_effects** and **bridgesampling** with
  the model ([\#30](https://github.com/popov-lab/bmm/issues/30))
- add informed default priors for all models. You can always use the
  **get_model_prior()** function to see the default priors for a model
- add a new function **set_default_prior** for developers, which allows
  them to more easily set default priors on new models regardless of the
  user-specified formula
- you can now specify variables for models via regular expressions
  rather than character vectors
  ([\#102](https://github.com/popov-lab/bmm/issues/102))
- you can now view and set all **bmm** global options via
  **bmm_options()**. See **?bmm_options** for more information
- add a start-up message upon loading the package

#### Bug fixes

- fix a bug in the **mixture3p** and **IMM** models which caused an
  error when intercept was not suppressed and set size was used as
  predictor
- **update()** now works properly with **bmmfit** objects
  ([\#95](https://github.com/popov-lab/bmm/issues/95))
- fix a bug in the **sort_data** check which caused an error when using
  grouped covariance structure in random effects across different
  parameters

#### Other changes

- **brms** is now loaded automatically when loading **bmm** with
  **library(bmm)**

## bmm 0.3.0

#### New features

- BREAKING CHANGE: The **fit_model** function now requires a
  **bmmformula** to be passed. The syntax of the **bmmformula** or its
  short form **bmf** is equal to specifying a **brmsformula**. However,
  as of this version the **bmmformula** only specifies how parameters of
  a **bmmodel** change across experimental conditions or continuous
  predictors. The response variables that the model is fit to now have
  to be specified when the model is defined using **model = bmmodel()**.
  ([\#79](https://github.com/popov-lab/bmm/issues/79))
- BREAKING CHANGE: The **non_target** and **spaPos** variables for the
  **mixture3p** and **IMM** models were relabeled to **nt_features** and
  **nt_distances** for consistency. This is also to communicate that
  distance is not limited to spatial distance but distances on any
  feature dimensions of the retrieval cues. Currently, still only a
  single generalization gradient for the cue features is possible.
- This release includes reference fits for all implemented models to
  ensure that future changes to the package do not compromise the
  included models and change the results that their implementations
  produce.
- The **check_formula** methods have been adapted to match the new
  **bmmformula** syntax. It now evaluates if formulas have been
  specified using the **bmmformula** function, if formulas for all
  parameters of a **bmmodel** have been specified and warns the user
  that only a fixed intercept will be estimated if no formula for one of
  the parameters was provided. Additionally, **check_formula** throws an
  error should formulas be provided that do not match a parameter of the
  called **bmmodel** unless they are part of a non-linear
  transformation.
- You can now specify formulas for internally fixed parameters such as
  **mu** in all visual working memory models. This allows you to predict
  if there is response bias in the data. If a formula is not provided
  for **mu**, the model will assume that the mean of the response
  distribution is fixed to zero.
- there is now an option **bmm.silent** that allows to suppress messages
- the baseline activation **b** was removed from the **IMM** models, as
  this is internally fixed to zero for scaling and as of now cannot be
  predicted by independent variables because the model would be
  unidentifiable.
- the arguments used to fit the **bmmodel** are now accessible in the
  **bmmfit** object via the `fit$bmm$fit_args` list.
- add class(‘bmmfit’) to the object returned from fit_model() allowing
  for more flexible postprocessing of the underlying **brmsfit** object.
  The object is now of class(‘bmmfit’, ‘brmsfit’)
- changes to column names of datasets **zhang_luck_2008** and
  **oberauer_lin_2017** to make them more consistent

#### Bug Fixes

- an error with the treatment of distances in the **IMMfull** and the
  **IMMbsc** has been corrected. This versions ensures that only
  positive distances can be passed to any of the two models.
- removed a warning regarding the scaling of the distances in the
  **IMMfull** and the **IMMbsc** that was specific only for circular
  distances.

#### Documentation

- All articles have been update to the new **bmmformula** syntax.

## bmm 0.2.2

#### Bug Fixes

- fixed a bug where passing a character vector or negative values to
  set_size argument of visual working memory models caused an error or
  incorrect behavior
  ([\#97](https://github.com/popov-lab/bmm/issues/97))

## bmm 0.2.1

#### Bug Fixes

- Minor change to sdmSimple Stan helper functions to avoid a harmless
  warning message in the Stan output

## bmm 0.2.0

#### New features

- New model available - The Signal Discrimination Model by
  Oberauer (2023) for visual working memory continuous reproduction
  tasks. See ?sdmSimple. The current version does not take into account
  non-target activation
- Add ability to extract information about the default priors in **bmm**
  models with **get_model_prior()**
  ([\#53](https://github.com/popov-lab/bmm/issues/53))
- Add ability to generate stan code and stan data for each model with
  **get_model_stancode()** and **get_model_standata()**
  ([\#81](https://github.com/popov-lab/bmm/issues/81))
- BREAKING CHANGE: Add distribution functions for likelihood
  (e.g. **dIMM()**) and random variate generation **rIMM()**) for all
  models in the package. Remove deprecated **gen_3p_data()** and
  **gen_imm_data()** functions
  ([\#69](https://github.com/popov-lab/bmm/issues/69))
- Two new data sets are available: **zhang_luck_2008** and
  **oberauer_lin_2017**
  ([\#22](https://github.com/popov-lab/bmm/issues/22))

#### Documentation

- Website for the development version of the package is now available at
  <https://venpopov.com/bmm/dev/>
  ([\#18](https://github.com/popov-lab/bmm/issues/18))
- Add articles for each model to the website at
  <https://venpopov.com/bmm/dev/articles/>
- Add a detailed developer’s guide to the website at
  <https://venpopov.com/bmm/dev/dev-notes>
  ([\#21](https://github.com/popov-lab/bmm/issues/21))
- Improve README with more detailed information about the package’s
  goals and its models
  ([\#21](https://github.com/popov-lab/bmm/issues/21))

#### Other changes

- Save **bmm** package version in the **brmsfit** object for
  reproducibility - e.g. `fit$version$bmm`
  ([\#88](https://github.com/popov-lab/bmm/issues/88))

## bmm 0.1.1

#### New features

- BREAKING CHANGE: Improve user interface to fit_model() ensures package
  stability and future development. Model specific arguments are now
  passed to the model functions as named arguments
  (e.g. **mixture3p(non_targets, setsize)**). This allows for a more
  flexible and intuitive way to specify model arguments. Passing model
  specific arguments directly to the **fit_model()** function is now
  deprecated ([\#43](https://github.com/popov-lab/bmm/issues/43)).
- Add information about each model such as domain, task, name, version,
  citation, requirements and parameters
  ([\#42](https://github.com/popov-lab/bmm/issues/42))
- Add ability to generate a template file for adding new models to the
  package with **use_model_template()** (for developers)
  ([\#39](https://github.com/popov-lab/bmm/issues/39))

#### Other changes

- Improve documentation of model functions. You can now get help on each
  model by typing **?model_name** into your console. For example,
  calling the information on the full version of the Interference
  Measurement Model would look like this: **?IMMfull**

## bmm 0.1.0

A major restructuring of the package to support stable and generalizable
development of future models
([\#41](https://github.com/popov-lab/bmm/issues/41)).

#### New Features

- Refactor the **fit_model()** function to be generic and independent of
  the model being fit
  ([\#20](https://github.com/popov-lab/bmm/issues/20))
- Transform models to be S3 objects.
  ([\#41](https://github.com/popov-lab/bmm/issues/41)).
- View currently supported models with new function
  **supported_models()**. Currently supported models are:
  **mixture2p()**, **mixture3p()**, **IMMabc()**, **IMMbsc()**,
  **IMMfull()**
- Add S3 methods for checking the data, formula, model and priors
  ([\#41](https://github.com/popov-lab/bmm/issues/41))
- Add distribution functions for the Signal Discrimination Model. See
  **?SDM** for usage
  ([\#27](https://github.com/popov-lab/bmm/issues/27))
- Add **softmax** and **softmaxinv** functions

#### Bug Fixes

- Change default prior on log(kappa) to Normal(2,1) for the
  **mixture3p()** model
  ([\#15](https://github.com/popov-lab/bmm/issues/15))

#### Other changes

- BREAKING CHANGE: deprecate **model_type** argument in **fit_model()**.
  Models must now be specified with S3 functions passed to argument
  **model** rather than model names as strings passed to argument
  **model_type** ([\#41](https://github.com/popov-lab/bmm/issues/41))
- Add extensive unit testing

## bmm 0.0.1

- Initial release version
