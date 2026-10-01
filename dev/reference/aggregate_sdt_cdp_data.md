# Aggregate long-format Remember/Know data for [`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md)

Reshapes long-format Remember/Know recognition data – one row per trial,
or one row per response category with a count column – into the
aggregated wide format
[`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md)
is fit to: one row per cell and one count column per response category,
named `new1 ... new<n_new>`, `know<k>` / `remember<k>` (and `guess<k>`
when the data contain "guess" judgments), with `k` on the unified
confidence scale `n_new + 1 ... n_new + n_old`.

## Usage

``` r
aggregate_sdt_cdp_data(data, judgment, confidence, count = NULL, response = "")
```

## Arguments

- data:

  A data frame in long format. Cells are defined by all columns other
  than `judgment`, `confidence`, and `count` (e.g. participant, stimulus
  class, and any condition variables), which are carried over to the
  output unchanged.

- judgment:

  The name of the column coding the memory judgment, with values
  `"new"`, `"remember"`, `"know"`, and optionally `"guess"`.

- confidence:

  The name of the column coding confidence on the unified old/new scale
  (integer, 1 = most confident "new" to K = most confident "old"; "new"
  judgments occupy the low levels, "old" judgments the high levels).

- count:

  Optional name of a column of response counts. If `NULL` (default),
  each row is treated as a single trial. Non-integer counts are rounded
  with a warning; `NA` counts and `NA` confidences are refused.

- response:

  Optional common prefix for the generated count columns (default `""`).
  Pass the same value to the `response` argument of
  [`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md).

## Value

A data frame with one row per cell: the cell-defining columns followed
by the response count columns, ready to pass to
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) with an
[`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md)
model. The numbers of "new" and "old" confidence levels are inferred
from the data.

## Examples

``` r
dat <- data.frame(
  stimulus = c(1L, 1L, 1L, 0L, 0L, 0L),
  judgment = c("remember", "know", "new", "new", "know", "new"),
  confidence = c(3L, 2L, 1L, 1L, 3L, 1L),
  count = c(40L, 25L, 15L, 55L, 10L, 20L)
)
aggregate_sdt_cdp_data(dat, "judgment", "confidence", "count")
#>   stimulus new1 know2 know3 remember2 remember3
#> 1        1   15    25     0         0        40
#> 2        0   75     0    10         0         0
```
