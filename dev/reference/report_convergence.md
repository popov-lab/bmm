# Convergence diagnostics of a fitted model

Returns the convergence diagnostics of a model fitted with
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md), per parameter
or as the sentence that
[`report_methods()`](https://venpopov.com/bmm/dev/reference/report_methods.md)
uses.

## Usage

``` r
report_convergence(fit, format = c("table", "text"))
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md). It must
  contain posterior draws, so mock fits and fits made with `chains = 0`
  or `empty = TRUE` are refused.

- format:

  `"table"` for one row per parameter, `"text"` for the sentence.

## Value

An object of class `bmm_report_convergence`: a `data.frame` with columns
`parameter`, `rhat`, `ess_bulk`, `ess_tail` and `in_summary` for
`format = "table"`, a character string for `format = "text"`. Subsetting
the table with `[` returns a plain `data.frame`.

## Details

The table holds the rank-normalized R-hat and the bulk and tail
effective sample sizes (ESS) of every estimated parameter, computed with
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html).
Constants, `lp__`, `lprior` and the draws from the prior (`prior_*`) are
left out. The column `in_summary` marks the parameters that
[`summary()`](https://rdrr.io/r/base/summary.html) shows; brms warns
when one of them has an R-hat above 1.05. The algorithm, the number of
divergent transitions after warmup and the number of transitions that
reached the maximum tree depth are attributes of the table, printed
below it.

The sentence gives the largest R-hat, the smallest bulk and tail ESS and
the number of divergent transitions, without verdicts. For algorithms
other than `"sampling"` the table is empty and the text says that there
are no diagnostics.

## See also

[`report_methods()`](https://venpopov.com/bmm/dev/reference/report_methods.md)

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
fit <- bmm(
  bmf(c ~ 0 + set_size, kappa ~ 1),
  data = oberauer_lin_2017,
  model = sdm(resp_error = "dev_rad")
)
report_convergence(fit)
report_convergence(fit, format = "text")
}
```
