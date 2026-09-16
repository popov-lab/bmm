# Report the priors used in a fitted bmm model

For each parameter of a fitted bmm model, reports the link function, the
prior actually applied on the sampling (link) scale, and where that
prior came from: a bmm default, a brms default, or a user-specified
prior. Parameters without a proper prior are flagged as flat.

## Usage

``` r
report_priors(fit, format = "table")
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)

- format:

  Character. `"table"` (default) prints the report as a table; `"text"`
  prints sentences ready for a methods section.

## Value

A `data.frame` of class `bmm_report_priors` with columns `parameter`,
`link`, `class`, `coef`, `group`, `prior` and `source`. Subsetting the
report with `[` returns a plain `data.frame`, as the report-specific
printing depends on columns and attributes that subsetting drops.

## Details

All priors and constants apply on the sampling scale set by each
parameter's link function: a prior for a log-link parameter describes
the log of that parameter, and a constant fixes the parameter on that
scale (e.g. `constant(0)` with a log link fixes the parameter to 1 on
the native scale). Constants are therefore printed together with their
exact native-scale value whenever the two differ.

The provenance of each prior is determined by re-deriving the default
priors from the model, formula and data stored in the fit. A
user-specified prior that is identical to the bmm default is therefore
reported as a default. Coefficients that inherit their prior from a more
general class are collapsed into the row of the prior they inherit from.

The one exception is a parameter fixed to a constant in the formula: bmm
folds such a constant into the model object at fit time, overwriting the
default it replaced, so the original default is no longer recoverable
from the fit. Constants named in the formula are therefore always
reported as `"user"`, even when the value restates the bmm default.

Re-derivation uses the `default_priors` stored inside the fit, so later
changes to a model's default *values* do not affect the report of an
older fit. It does use the currently installed code that turns those
values into prior rows, so if that construction changes between
versions, the reported provenance of a fit made with an earlier version
may differ from what was actually applied. No check against
`fit$version$bmm` is performed.

Flat priors are flagged because they are improper: Bayes factors via
bridge sampling are undefined when any parameter has an improper prior.

Parameters that exist only because the family machinery requires them
are omitted: a fixed parameter the model does not declare (the `mu` brms
forces on the custom families of the response-time models, or the `mu2`
and `kappa2` of a two-component
[`brms::mixture()`](https://paulbuerkner.com/brms/reference/mixture.html)
family), and the mixture-weight reference component `theta2`, which brms
lists with a default prior although the sampler holds it at zero. Fixed
parameters the model does declare are reported, since the user can
estimate them: `mu` for
[`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md), `mu1` for the
circular mixture models, `zr` for
[`ddm()`](https://venpopov.com/bmm/dev/reference/ddm.md).

## See also

[`default_prior()`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md),
[`parameters()`](https://venpopov.com/bmm/dev/reference/parameters.md)

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
fit <- bmm(
  bmf(c ~ 0 + set_size, kappa ~ 1),
  data = oberauer_lin_2017,
  model = sdm(resp_error = "dev_rad")
)
report_priors(fit)
report_priors(fit, format = "text")
}
```
