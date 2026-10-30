# Transform values between the native and the parameter scale of a link

`link_transform()` applies one of the links bmm uses for model
parameters to a numeric vector. With `inverse = FALSE` it maps values
from the native scale of a parameter, the scale of the model equations,
to the parameter scale, the scale on which
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) samples the
intercept and coefficients of the parameter's formula. With
`inverse = TRUE` it maps values from the parameter scale back to the
native scale.

## Usage

``` r
link_transform(values, link, inverse = FALSE)
```

## Arguments

- values:

  A numeric vector, matrix or array. With `inverse = FALSE` the values
  are on the native scale, with `inverse = TRUE` on the parameter scale.
  `NA` stays `NA` at its position.

- link:

  A single string naming the link: `"identity"`, `"log"`, `"softplus"`,
  `"log1p"`, `"logm1"`, `"inverse"`, `"sqrt"`, `"logit"`, `"probit"`,
  `"tan_half"`, `"loglog"` or `"cloglog"`. A `NULL` link is read as the
  identity link, as for a model parameter that declares no link.

- inverse:

  Logical. `FALSE` (the default) transforms from the native to the
  parameter scale; `TRUE` transforms from the parameter scale to the
  native scale.

## Value

The transformed values, with the length and dimensions of `values` (a
matrix stays a matrix).

## Details

The links, with the native values each accepts. `x` is a value on the
native scale, `eta` one on the parameter scale.

|  |  |  |  |
|----|----|----|----|
| `link` | `inverse = FALSE`: `eta` from `x` | `inverse = TRUE`: `x` from `eta` | Native values `x` |
| `"identity"` | `x` | `eta` | all real numbers |
| `"log"` | `log(x)` | `exp(eta)` | `x > 0` |
| `"softplus"` | `log(exp(x) - 1)` | `log(1 + exp(eta))` | `x > 0` |
| `"log1p"` | `log(1 + x)` | `exp(eta) - 1` | `x > -1` |
| `"logm1"` | `log(x - 1)` | `exp(eta) + 1` | `x > 1` |
| `"inverse"` | `1 / x` | `1 / eta` | `x != 0` |
| `"sqrt"` | `sqrt(x)` | `eta^2` | `x >= 0` |
| `"logit"` | `log(x / (1 - x))` | `1 / (1 + exp(-eta))` | `0 < x < 1` |
| `"probit"` | `qnorm(x)` | `pnorm(eta)` | `0 < x < 1` |
| `"tan_half"` | `tan(x / 2)` | `2 * atan(eta)` | `-pi < x < pi` |
| `"loglog"` | `log(-log(x))` | `exp(-exp(eta))` | `0 < x < 1` |
| `"cloglog"` | `log(-log(1 - x))` | `1 - exp(-exp(eta))` | `0 < x < 1` |

Values outside these ranges are not checked. For the links other than
`"tan_half"`, a value at an edge the range excludes, such as `0` for
`"log"` and `"inverse"` or `1` for `"logit"`, gives `-Inf` or `Inf`
without a warning, and a value beyond it gives `NaN` (`NA` for
`"softplus"`) with a warning. `"tan_half"` reads `x` as an angle: a
value outside `-pi` to `pi` gives the `eta` of the same angle inside
that range, and returns that angle on the way back. `"sqrt"` squares on
the way back, so a negative `eta` comes back as its absolute value: `-3`
gives `9`, and `9` gives `3`.

The table holds the links that transform one value at a time. It does
not include `"softmax"`, the link of the mixture weights of models such
as [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md):
a softmax needs all weights of a group at once and a reference position,
so `link_transform()` stops for it, as for any other name that is not in
the table.
[`parameter_info()`](https://venpopov.com/bmm/dev/reference/parameter_info.md)
shows the link of each parameter of a model in its `link` column; a
parameter with `"softmax"` there is not covered by `link_transform()`.
For the weights on the native scale, use
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md).
How the weights of a softmax group combine depends on the model, for
[`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md) on
the set size (at set size 1 the non-target weight is zero), and
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
applies that.

## See also

[`parameter_info()`](https://venpopov.com/bmm/dev/reference/parameter_info.md)
for the link of each parameter of a model, and
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
for the parameters of a fitted model on their native scale.

## Examples

``` r
# a probability on the parameter scale of a logit link, and back
eta <- link_transform(c(0.2, 0.5, 0.8), "logit")
eta
#> [1] -1.386294  0.000000  1.386294
link_transform(eta, "logit", inverse = TRUE)
#> [1] 0.2 0.5 0.8

# the native value of a parameter with a log link, from its parameter scale
link_transform(c(-1, 0, 1), "log", inverse = TRUE)
#> [1] 0.3678794 1.0000000 2.7182818

# the link of each parameter of a model
parameter_info(sdm(resp_error = "y"))$link
#> [1] "tan_half" "log"      "log"     
```
