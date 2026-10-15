# Posterior predictive p-values for count models

Computes the goodness-of-fit statistics \\T_1\\ and \\T_2\\ of Klauer
(2010) and their posterior predictive p-values (PPP) for models that fit
aggregated response counts:
[`mpt()`](https://venpopov.com/bmm/reference/mpt.md),
[`m3()`](https://venpopov.com/bmm/reference/m3.md),
[`sdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn.md),
[`sdt_rating()`](https://venpopov.com/bmm/reference/sdt_rating.md),
[`sdt_mafc()`](https://venpopov.com/bmm/reference/sdt_mafc.md),
[`sdt_ranking()`](https://venpopov.com/bmm/reference/sdt_ranking.md) and
[`sdt_cdp()`](https://venpopov.com/bmm/reference/sdt_cdp.md). TreeBUGS'
`PPP()` reports the same two statistics; Details says where their values
differ.

## Usage

``` r
ppp(fit, ...)

# S3 method for class 'bmmfit'
ppp(
  fit,
  statistic = c("T1", "T2"),
  by = NULL,
  ndraws = NULL,
  level = c("group", "participant"),
  participant = NULL,
  ...
)
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/reference/bmm.md).

- ...:

  Unused.

- statistic:

  Character. `"T1"`, `"T2"`, or both (the default). If the participants
  cannot be identified, the default computes `"T1"` only.

- by:

  Character vector naming the columns of the model data that define the
  cells. `NULL` (default) uses the cells described in Details;
  `character(0)` pools all rows into one cell. For
  [`mpt()`](https://venpopov.com/bmm/reference/mpt.md) fits, the tree
  column can be named even if it is not a predictor.

- ndraws:

  Integer. Number of posterior draws to use. `NULL` (default) uses all
  draws.

- level:

  Character. `"group"` (default) for the statistics of the whole data
  set; `"participant"` for one \\T_1\\ and PPP per participant, as
  TreeBUGS reports them, to find the participants the model fits poorly.

- participant:

  Character. The column that identifies participants, needed for \\T_2\\
  and for `level = "participant"`. It must be a variable of the model
  formula, because the model data contain only those; a model without
  random effects has one only if a predictor such as `0 + id` names the
  participants. `NULL` (default) uses the grouping variable of the
  random effects; when the model has more than one (e.g. participants
  and items), name it here.

## Value

A data frame of class `bmm_ppp` with one row per statistic (and per
participant for `level = "participant"`) and the columns `statistic`,
`observed` and `replicated` (the means of the statistic over draws, for
the observed and the simulated counts, as reported by Klauer, 2010, and
TreeBUGS), `ppp` and `ndraws`. The draws of the statistics are stored in
the attribute `"draws"`.

## Details

Both statistics are computed once per posterior draw, for the observed
counts and for counts simulated from that draw with
[`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html),
each against the expected counts of the same draw from
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html).
The simulated data are new data from the same participants. The PPP is
the proportion of draws in which the simulated counts deviate more from
the expectation than the observed counts do. A PPP near 0 flags misfit;
a PPP away from 0 means that the model reproduces the aspect of the data
the statistic measures.

**\\T_1\\: mean structure.** The rows of the data are summed into cells,
and \\T_1\\ is Pearson's \\X^2\\ over all cells and response categories,
\$\$T_1 = \sum\_{c,k} \frac{(n\_{ck} -
\hat{n}\_{ck})^2}{\hat{n}\_{ck}},\$\$ where \\n\_{ck}\\ is the observed
and \\\hat{n}\_{ck}\\ the expected count of category \\k\\ in cell
\\c\\. Klauer (2010, Eq. 17) and TreeBUGS average the counts over
participants instead of summing them. With the same number \\T\\ of
participants in every cell, the two give the same PPP and \\T_1\\ here
is \\T\\ times theirs; otherwise the mean gives the cells with fewer
participants more weight.

**\\T_2\\: covariance structure.** For every participant, the counts are
summed per cell, giving one vector per participant over all cells and
categories. \\T_2\\ compares the covariance matrix \\s\\ of these
vectors across participants with the covariance matrix \\\sigma\\ the
model expects (Klauer, 2010, Eq. 18): \$\$T_2 = \sum\_{i,j}
\frac{(s\_{ij} - \sigma\_{ij})^2}{\sqrt{\sigma\_{ii}\sigma\_{jj}}}.\$\$
\\\sigma\\ is the covariance of the expected counts across participants
plus the average multinomial covariance within a participant, which is
the expected value of \\s\\ given the draw. Both covariances divide by
\\T - 1\\, where \\T\\ is the number of participants. Klauer (2010)
multiplies the multinomial part by \\(T - 1) / T\\, the expected value
for covariances that divide by \\T\\; that version of \\T_2\\ is
\\(T - 1) / T\\ times this one and has the same PPP. TreeBUGS combines
the factor with covariances that divide by \\T - 1\\, so its \\\sigma\\
falls short of the expected \\s\\ by \\1 / T\\ of the multinomial part.
The sum runs over the full matrix, so the covariances between trees or
conditions, which carry the correlations between the person parameters,
are part of the check. When participants differ in the cells they have,
as in a between-participants design, they are grouped by their set of
cells, \\T_2\\ is computed per group, and the results are summed. A
participant who shares the set of cells with nobody else has no
covariance and is left out, with a message.

**Cells.** By default, a cell is a combination of the tree
([`mpt()`](https://venpopov.com/bmm/reference/mpt.md)), the stimulus
(the SDT models) or the set size
([`sdt_ranking()`](https://venpopov.com/bmm/reference/sdt_ranking.md)),
and the categorical population-level predictors of the model formula.
Numeric predictors are not part of the default, because every distinct
value would form a cell of its own. Rows of a cell are summed, so the
trees or stimuli must be part of `by` for the statistic to compare each
tree's own categories.

**Impossible categories.** Categories with an expected count of about 0
in every draw, such as the impossible categories of an
[`mpt()`](https://venpopov.com/bmm/reference/mpt.md) or
[`m3()`](https://venpopov.com/bmm/reference/m3.md) model or the ranks
beyond the set size of an
[`sdt_ranking()`](https://venpopov.com/bmm/reference/sdt_ranking.md)
model, are left out of both statistics. A cell with a single possible
category, such as an [`m3()`](https://venpopov.com/bmm/reference/m3.md)
set size without other candidates, is left out of \\T_2\\: the number of
trials fixes its count.

## References

Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A
latent-trait approach. *Psychometrika, 75*(1), 70–98.
https://doi.org/10.1007/s11336-009-9141-0

## See also

[`pp_check.bmmfit()`](https://venpopov.com/bmm/reference/pp_check.bmmfit.md)
for the graphical check.

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
trees <- list(
  mpt_tree("target", list(old = "D + (1 - D) * g", new = "(1 - D) * (1 - g)")),
  mpt_tree("lure", list(old = "(1 - D) * g", new = "D + (1 - D) * (1 - g)"))
)
dat <- expand.grid(id = factor(1:20), item_type = c("target", "lure"),
                   stringsAsFactors = FALSE)
dat$old <- rbinom(nrow(dat), 50, ifelse(dat$item_type == "target", 0.8, 0.2))
dat$new <- 50 - dat$old

fit <- bmm(
  bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id)), dat, mpt(trees, "item_type"),
  cores = 4, backend = "cmdstanr"
)
ppp(fit)
ppp(fit, level = "participant")
}
```
