# Expected response of a fitted bmm model

`posterior_epred()` returns draws of the expected value of the response
for each observation: the mean of what
[`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
simulates for it. What that is depends on the model:

|  |  |
|----|----|
| Model | Expected response |
| [`ddm()`](https://popov-lab.github.io/bmm/dev/reference/ddm.md), [`cswald()`](https://popov-lab.github.io/bmm/dev/reference/cswald.md) | Mean response time, averaged over both responses, under the diffusion process that `posterior_predict()` simulates from |
| [`ezdm()`](https://popov-lab.github.io/bmm/dev/reference/ezdm.md), version `"3par"` | Mean response time of the cell |
| [`ezdm()`](https://popov-lab.github.io/bmm/dev/reference/ezdm.md), version `"4par"` | Mean response time of the cell's upper-boundary responses |
| [`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md) | Number of "old"/"signal" responses, the number of trials times their probability |
| [`sdt_mafc()`](https://popov-lab.github.io/bmm/dev/reference/sdt_mafc.md) | Number of correct responses, the number of trials times their probability |
| [`m3()`](https://popov-lab.github.io/bmm/dev/reference/m3.md), [`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md), [`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md), [`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md) | Expected count of each response category, from the multinomial family of brms |
| [`sdm()`](https://popov-lab.github.io/bmm/dev/reference/sdm.md), [`mixture2p()`](https://popov-lab.github.io/bmm/dev/reference/mixture2p.md), [`mixture3p()`](https://popov-lab.github.io/bmm/dev/reference/mixture3p.md), [`imm()`](https://popov-lab.github.io/bmm/dev/reference/imm.md) | Not defined: the mean of a circular response error is not a useful quantity, so these models stop with an error |

With `dpar` or `nlpar`, `posterior_epred()` returns draws of that model
parameter for every model, as in brms: a `dpar` on its native scale, an
`nlpar` on the scale of its link.
[`native_parameters()`](https://popov-lab.github.io/bmm/dev/reference/native_parameters.md)
returns the model parameters on their native scale over a grid of
predictor values.

## Usage

``` r
# S3 method for class 'bmmfit'
posterior_epred(object, ..., dpar = NULL, nlpar = NULL)

# S3 method for class 'bmmfit'
fitted(object, ..., scale = c("response", "linear"), dpar = NULL, nlpar = NULL)
```

## Arguments

- object:

  A `bmmfit` object.

- ...:

  Further arguments passed to
  [`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
  such as `newdata` or `ndraws`.

- dpar, nlpar:

  Name of a distributional or non-linear parameter whose draws are
  returned instead of the expected response.

- scale:

  As in
  [`brms::fitted.brmsfit()`](https://paulbuerkner.com/brms/reference/fitted.brmsfit.html):
  `"response"` is the expected response, `"linear"` the linear predictor
  of `mu`.

## Value

A draws by observations matrix (an array with a third dimension for the
response categories of the multinomial models).

## See also

[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
[`native_parameters()`](https://popov-lab.github.io/bmm/dev/reference/native_parameters.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bmm(
  bmf(drift ~ 1, bound ~ 1, ndt ~ 1),
  data = rddm(200, drift = 1.5, bound = 1.2, ndt = 0.3),
  model = ddm(rt = "rt", response = "response"),
  backend = "cmdstanr"
)
# expected response time of each observation, one row per draw
epred <- posterior_epred(fit)
} # }
```
