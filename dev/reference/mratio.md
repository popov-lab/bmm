# Metacognitive efficiency (M-ratio) from a fitted meta-d' SDT model

Extracts the posterior of the M-ratio (\\\mathrm{meta\text{-}d'}/d'\\)
and of meta-d' from a
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
model fit with `version = "metad"`. That model estimates the log M-ratio
(`logmratio`) directly, so the M-ratio is `exp(logmratio)` and meta-d'
is `exp(logmratio) * d`. The M-ratio is the standard measure of
metacognitive efficiency (Maniscalco & Lau, 2012; Fleming, 2017): 1 is
ideal metacognition, below 1 is inefficiency, and above 1 is
hyper-efficiency. Type-1 and type-2 sensitivity are reported on one
scale, \\d'\\ or, once `sdratio` is estimated, \\d_a\\ (see
[`sdt_sensitivity()`](https://popov-lab.github.io/bmm/dev/reference/sdt_sensitivity.md)),
so the ratio is unaffected by `sdratio`. The type-1 criterion keeps its
location in the meta-d' space, as in the HMeta-d model equations
(Fleming, 2017, Appendix); maximum-likelihood meta-d' constrains meta-c'
= c' (Maniscalco & Lau, 2014), so the two differ when the criterion is
far from the midpoint.

## Usage

``` r
mratio(fit, conditions = NULL, probs = c(0.025, 0.975), ...)
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) from
  an SDT model.

- conditions:

  Optional data frame of predictor values at which to evaluate the
  M-ratio. Its columns must be population-level predictors of
  `logmratio` or `d` in the formula; grouping variables of random
  effects are refused. If `NULL` (default), unique predictor
  combinations are derived from the data, and a fit without such a
  predictor gives one row.

- probs:

  Numeric vector of length 2. Lower and upper quantiles for the credible
  interval (default `c(0.025, 0.975)`).

- ...:

  Additional arguments passed to
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html),
  such as `draw_ids` to use a subset of the posterior draws. `ndraws` is
  refused: each parameter is drawn by its own call, so random subsets
  would not match across parameters.

## Value

A data frame of class `"bmm_sdt_mratio"` summarising the posterior:
columns `parameter` (`"mratio"` or `"metad"`), `mean`, `median`,
`lower`, `upper` (the credible-interval bounds at `probs`), and any
condition columns. The full per-draw posteriors are kept in the `draws`
attribute (long format: `parameter`, `value`, `.draw`, plus condition
columns) for downstream computation such as plotting or condition
contrasts.

## References

Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach
for estimating metacognitive sensitivity from confidence ratings.
*Consciousness and Cognition*, *21*(1), 422–430.
[doi:10.1016/j.concog.2011.09.021](https://doi.org/10.1016/j.concog.2011.09.021)

Maniscalco, B., & Lau, H. (2014). Signal detection theory analysis of
type 1 and type 2 data: meta-d', response-specific meta-d', and the
unequal variance SDT model. In S. M. Fleming & C. D. Frith (Eds.), *The
cognitive neuroscience of metacognition* (pp. 25–66). Springer.

Fleming, S. M. (2017). HMeta-d: hierarchical Bayesian estimation of
metacognitive efficiency from confidence ratings. *Neuroscience of
Consciousness*, *2017*(1), nix007.
[doi:10.1093/nc/nix007](https://doi.org/10.1093/nc/nix007)

## See also

[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md),
[`latent_sdt()`](https://popov-lab.github.io/bmm/dev/reference/latent_sdt.md)
