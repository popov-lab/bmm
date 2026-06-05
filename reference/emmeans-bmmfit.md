# emmeans support for bmmfit objects

These S3 methods allow `emmeans::emmeans()` to work seamlessly with
bmmfit objects by automatically resolving whether a parameter is a
distributional (`dpar`) or non-linear (`nlpar`) parameter in brms. Users
can always use `dpar = "parameter_name"` regardless of internal
classification.

## Usage

``` r
# S3 method for class 'bmmfit'
recover_data(object, ..., dpar = NULL, nlpar = NULL)

# S3 method for class 'bmmfit'
emm_basis(object, trms, xlev, grid, ..., dpar = NULL, nlpar = NULL)
```

## Arguments

- object:

  A bmmfit object (created by
  [`bmm()`](https://venpopov.com/bmm/reference/bmm.md))

- ...:

  Additional arguments passed to the brmsfit methods.

- dpar:

  Character string. Name of the model parameter (e.g., `"kappa"`, `"c"`,
  `"thetat"`). Automatically resolved to the correct brms parameter
  type.

- nlpar:

  Character string. Explicit non-linear parameter name. If provided,
  takes precedence over `dpar`.

- trms, xlev, grid:

  Arguments passed by emmeans internally.

## Value

See `emmeans::recover_data()` and `emmeans::emm_basis()` for return
values.

## Details

bmm models use two types of parameters internally in brms:

- **Distributional parameters (`dpar`)**: Used in models with
  `custom_family(dpars = ...)` (e.g., SDM, EZDM)

- **Non-linear parameters (`nlpar`)**: Used in models with
  [`bmf2bf()`](https://venpopov.com/bmm/reference/bmf2bf.md) + `nlf()`
  (e.g., mixture2p, mixture3p, IMM, M3)

Users should not need to know this distinction. These methods intercept
the `dpar` argument and, if the parameter is actually an nlpar, silently
re-route it to `nlpar`.

## See also

`emmeans::emmeans()`

## Examples

``` r
if (FALSE) { # \dontrun{
# Works for any bmm model — no need to know dpar vs nlpar
em <- emmeans(fit, ~ condition, dpar = "kappa")
pairs(em)
confint(em)
} # }
```
