# Extract information from a brmsfit object

Extract information from a brmsfit object

## Usage

``` r
fit_info(fit, what)
```

## Arguments

- fit:

  A brmsfit object, or a list of brmsfit objects

- what:

  String. What to return:

  - "time" for the sampling time per chain

  - "time_mean" for the mean sampling time

## Value

Depends on `what` and the class of `fit`. For `brmsfit` objects,
information about the single fit is returned. For `brmsfit_list`
objects, a list or data.frame with the information for each fit is
returned.

- "time": A data.frame with the sampling time per chain

- "time_mean": A named numeric vector with the mean sampling time

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
fit <- bmm(
  formula = bmmformula(c ~ 1, kappa ~ 1),
  data = data.frame(y = rsdm(1000)),
  model = sdm(resp_error = "y")
)

fit_info(fit, "time")
}
```
