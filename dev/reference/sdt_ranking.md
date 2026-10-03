# Ranking Signal Detection Theory Model

Ranking Signal Detection Theory Model

## Usage

``` r
sdt_ranking(response, m, dist = c("gumbel_min", "normal"), links = NULL, ...)
```

## Arguments

- response:

  A character vector of column names with the target rank-count columns,
  ordered from rank 1 (most likely target) to rank `m` (least). With a
  constant `m`, supply exactly `m` columns; with a varying set size,
  supply `max(m)` columns.

- m:

  Either a single integer \>= 2 giving the number of ranked items
  (constant across all rows), or a single string naming a data column
  that gives the number of ranked items per row.

- dist:

  Character. The distribution assumed for the latent evidence, given
  here by its cumulative distribution function:

  - "gumbel_min" (default): smallest extreme value, \\1 -
    \exp(-\exp(x))\\ (complementary log-log). Closed form via
    gamma-function ratios.

  - "normal": Gaussian, \\\Phi(x)\\ (supports unequal variance via
    `sdratio`)

- links:

  A named list of link functions for the parameters. Only the link of
  `d` can be changed, and only for `dist = "gumbel_min"`: with
  `dist = "normal"` the quadrature loses accuracy at the large `d` a log
  link reaches, so `d` keeps the identity link. `sdratio` always keeps
  the identity link.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Perception & Recognition Memory

- **Task:** Ranking Task

- **Name:** Signal Detection Theory (Ranking)

- **Citation:**

  - Meyer-Grant, C. G., Kellen, D., Harding, S. M., & Singmann, H.
    (2026). Extreme-value signal detection theory for recognition
    memory: The parametric road not taken. Psychological Review. Advance
    online publication. https://doi.org/10.1037/rev0000615

- **Requirements:**

  Provide pre-aggregated ranking counts in wide format:

&nbsp;

- Rank-count columns (response): one column per rank position, each
  giving the number of trials in which the target received that rank
  (column 1 = most likely target, column m = least)

- Set size (m): a constant or a column giving the number of ranked items
  per row; rows with fewer ranks leave the surplus columns at 0 No
  stimulus column needed (all trials include exactly one target)

&nbsp;

- **Parameters:**

  - `d`: Sensitivity: d' under equal variance (g' for gumbel_min). When
    sdratio is estimated, d is d_a, the distance between the target and
    lure distributions in units of their root-mean-square SD

- **Fixed parameters:**

- **Default parameter links:**

  - d = identity

- **Default priors:**

  - `d`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

Models the rank ordering of `m` items by perceived strength. Only `d` is
estimated (no criterion). Supports `dist = "gumbel_min"` (closed-form
via lgamma ratios) and `dist = "normal"` (Gauss-Hermite quadrature).

The model uses the native brms multinomial family: each rank position is
a multinomial category whose logit is set to `log p(rank)`, so `softmax`
recovers the rank distribution exactly. This means `log_lik`,
`posterior_predict`, `posterior_epred`, and `pp_check` come from brms as
proper joint multinomial draws.

## Sensitivity is on the same scale as [`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)

`d` is \\d'\\ whenever the target and lure distributions share an SD:
always for `dist = "gumbel_min"`, where it is the \\g'\\ of Meyer-Grant
et al. (2026), and for `dist = "normal"` unless you give `sdratio` a
formula. With `sdratio` estimated, `d` is the balanced index \\d_a\\
that
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
reports, the separation divided by the root-mean-square of the two SDs;
the noise-standardized separation is then
`d * sqrt((1 + exp(sdratio)^2) / 2)`.

Ranking has no criterion, so \\d_a\\ is not read off an ROC here. It is
adopted to keep `d` on the same scale across the SDT family, and it has
a direct meaning for rankings: it fixes the two-item accuracy at
\\\Phi(d/\sqrt{2})\\ whatever `sdratio` is (see below).

Ranking is the one SDT design that identifies the variance ratio from a
single condition. The rank distribution supplies `m - 1` free
probabilities per set size, so with `m >= 3` there is enough information
to separate `d` from `sdratio` without the criterion sweep that
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
needs – the *shape* of the rank distribution, not just its mean, carries
the ratio.

At `m = 2` the model reduces to 2AFC and the two parameters are no
longer separable: the probability of ranking the target first is the
area under the yes/no ROC, which for Gaussian noise is
\\\Phi(d/\sqrt{2})\\ whatever `sdratio` is. Keep `sdratio` fixed for
two-item designs.

For Gaussian ranking (`dist = "normal"`), `sdratio` is fixed to 0 by
default. Add `sdratio ~ 1` to the formula for unequal-variance ranking.

The set size `m` may be a constant or the name of a data column. Supply
a column to fit trials with different set sizes in a single model: the
response has `max(m)` columns, and rows with a smaller set size switch
off the surplus rank categories (the multinomial then renormalizes over
the valid ranks).

## References

Meyer-Grant, C. G., Kellen, D., Harding, S. M., & Singmann, H. (2026).
Extreme-value signal detection theory for recognition memory: The
parametric road not taken. *Psychological Review*. Advance online
publication.
[doi:10.1037/rev0000615](https://doi.org/10.1037/rev0000615)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- data.frame(id = 1:20)
dat <- cbind(dat, rsdt_ranking(20, 200, m = 4, d = 1.5))

model <- sdt_ranking(
  response = c("rank1", "rank2", "rank3", "rank4"),
  m = 4
)

fit <- bmm(
  formula = bmf(d ~ 1),
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# Mixed set sizes in one model: response has max(m) columns, m is a column
model_mixed <- sdt_ranking(
  response = c("rank1", "rank2", "rank3", "rank4", "rank5"),
  m = "set_size"
)
} # }
```
