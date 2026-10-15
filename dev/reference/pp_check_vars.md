# List the posterior predictive checks available for a model fit

Lists the values `resp_var` accepts in
[`pp_check.bmmfit()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
for a fitted model whose likelihood involves several observables.

## Usage

``` r
pp_check_vars(fit)
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md).

## Value

A `data.frame` with one row per available check (columns `resp_var`,
`label`, `default_type`, `slot`, listing the brms standata slots the
check reads, and `default`, flagging the observable that
[`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
plots when `resp_var` is not specified), or `NULL` invisibly for models
without multi-observable support.

## See also

[`pp_check.bmmfit()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bmm(bmmformula(drift ~ condition), data, ddm(rt = "rt", response = "response"))
pp_check_vars(fit)
pp_check(fit, resp_var = "response")
} # }
```
