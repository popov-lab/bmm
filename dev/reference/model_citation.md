# References for a measurement model

Returns the published source(s) of a model, for the reference list of a
paper that uses it.

## Usage

``` r
model_citation(x, ...)
```

## Arguments

- x:

  A model object (for example `sdm(resp_error = "y")`) or a `bmmfit`
  object returned by
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md).

- ...:

  Currently ignored.

## Value

A character vector with one reference in APA style per element, or
`character(0)` if no reference is available.

## Details

The references come from the model constructor of the installed version
of bmm, so a fit made with an older version gets the current, corrected
references. When the installed version no longer has a constructor for
the model, or its references are empty, the copy stored in the fit or
model object is returned instead.

Some versions of a model cite the papers that introduced them instead,
for example `sdt_rating(version = "dpsdt")`.

A custom [`m3()`](https://popov-lab.github.io/bmm/dev/reference/m3.md)
model is cited with the paper that introduced the M3 framework. The
model structure itself is defined by the user, so its source cannot be
known to bmm. The returned vector then has the attribute `uncited`,
naming the part of the model you have to cite yourself.

bmm itself is cited separately, see `citation("bmm")`.

## See also

[`citation()`](https://rdrr.io/r/utils/citation.html)

## Examples

``` r
model_citation(sdm(resp_error = "y"))
#> [1] "Oberauer, K. (2023). Measurement models for visual working memory—A factorial model comparison. Psychological Review, 130(3), 841-852. https://doi.org/10.1037/rev0000328"
model_citation(
  sdt_rating(response = paste0("r", 1:4), stimulus = "old", version = "dpsdt")
)
#> [1] "Yonelinas, A. P. (1994). Receiver-operating characteristics in recognition memory: Evidence for a dual-process model. Journal of Experimental Psychology: Learning, Memory, and Cognition, 20(6), 1341-1354. https://doi.org/10.1037/0278-7393.20.6.1341"
#> [2] "Yonelinas, A. P. (2024). The role of recollection and familiarity in visual working memory: A mixture of threshold and signal detection processes. Psychological Review, 131(2), 321-348. https://doi.org/10.1037/rev0000432"                            
```
