# Ranking SDT rank log-probability (multinomial logit)

R companion to the Stan `sdt_ranking_logmu` function. It returns
`log(p(rank))` for rank position `cat`, which the
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md)
multinomial formula uses as the category logit (so `softmax` recovers
the ranking probabilities). `brms` evaluates the non-linear formula in R
for `posterior_predict()` and `posterior_epred()`, so this function must
be on the search path; it is exported for that reason and is not called
directly.

## Usage

``` r
sdt_ranking_logmu(cat, max_rank, d, sdratio = 0, dist = 2L)
```

## Arguments

- cat:

  Integer rank position index.

- max_rank:

  Set size (number of ranked items) for the observation; ranks above it
  return the finite `-100` sentinel so the category switches off.

- d:

  Sensitivity (draws-by-observation matrix supplied by brms).

- sdratio:

  Log SD ratio; `0` for the gumbel_min distribution.

- dist:

  Integer noise-distribution id (2 = gumbel_min, 1 = normal).

## Value

`log(p(cat))`, matching the shape of `d`.
