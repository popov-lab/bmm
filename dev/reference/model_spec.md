# The specification of a measurement model

`model_spec()` returns everything bmm states about a model in one list
of plain R data: its name and version, its parameters with their links
and fixed values, the default priors, the ranges initial values are
drawn from, the data columns it reads, and the names of its density and
generator functions. Packages that build on bmm read a model here, not
from the fields of the model object.

## Usage

``` r
model_spec(x, ...)

# S3 method for class 'bmmodel'
model_spec(x, ...)

# S3 method for class 'bmmfit'
model_spec(x, ...)
```

## Arguments

- x:

  A `bmmodel` object, such as `sdm(resp_error = "y")`, or a `bmmfit`
  object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md).

- ...:

  Unused.

## Value

A list of class `bmm_model_spec` with these elements, in this order:

- `name`:

  The model's name.

- `package`:

  The package that defines the model: `"bmm"` for the models bmm ships.

- `version`:

  The model's version, such as `"3par"` for
  [`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md), or `NA`
  for a model without versions. This is not the version of bmm.

- `domain`, `task`:

  The research domains and the tasks the model is for, as listed by
  [`bmm_models()`](https://venpopov.com/bmm/dev/reference/bmm_models.md).

- `citation`:

  The references for the model, one per element.

- `parameters`:

  A data frame with one row per parameter and the columns `parameter`,
  `description`, `link` (the link function, see
  [`link_transform()`](https://venpopov.com/bmm/dev/reference/link_transform.md)),
  `fixed` (logical), `value` (the value of a fixed parameter on the link
  scale, else `NA`), `alternative_group` (the parameters of which a
  formula may use only one, joined by `"|"`, such as `"criterion|bias"`,
  else `NA`) and `softmax` (`TRUE` for mixture weights that are
  transformed together through a softmax, whose `link` then does not
  apply on its own). For
  [`m3()`](https://venpopov.com/bmm/dev/reference/m3.md) with
  `version = "custom"` the formula defines the parameters, so the
  specification of the constructor lists `b` only and the specification
  of a fit lists them all. A model with variable precision
  ([`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md),
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md)
  or [`imm()`](https://venpopov.com/bmm/dev/reference/imm.md) with
  `variable_precision = TRUE`) is recognised by `tau` among its
  parameters; `variable_precision` and `vp_nodes` are not items of the
  specification.

- `default_priors`:

  A data frame with the columns `parameter`, `main`, `effects` and `sd`:
  the default prior on a parameter's intercept, on its other regression
  coefficients, and on its random effect standard deviations, as Stan
  strings.

- `init_ranges`:

  A data frame with the columns `parameter`, `lower` and `upper`: the
  range initial values are drawn from, on the native scale (softmax
  weights on the sampling scale). It has no rows for models that build
  their initial values otherwise
  ([`m3()`](https://venpopov.com/bmm/dev/reference/m3.md),
  [`mpt()`](https://venpopov.com/bmm/dev/reference/mpt.md)).

- `resp_vars`, `other_vars`:

  Named lists of the arguments of the model's constructor, as the model
  was created with them: the response variables, and the other variables
  and settings. Values the constructor stores beyond its arguments are
  not listed. Expressions in them, such as the branches of
  [`mpt()`](https://venpopov.com/bmm/dev/reference/mpt.md) trees, are
  given as strings.

- `data_columns`:

  A data frame with the columns `argument` and `role`: what the data
  column named by each argument of the constructor must hold.

- `functions`:

  A list with the names of the model's exported `density` and
  `generator` (random generation) functions, such as `"dsdm"` and
  `"rsdm"`; `NA` where a model has none. They are the model's exported
  d/r pair and may not cover every version:
  [`dmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md)
  and
  [`rmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md)
  do not express the item-wise precision of
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md)
  `version = "slot_averaging"` above `K = set_size`.

- `contract_version`:

  The version of this specification, `"0.9"` until bmm 1.6.0 fixes it as
  `"1.0"`.

For a `bmmfit` the list ends with one more element, `bmm_version`: the
version of bmm the model was fitted with. Its `name`, `domain` and
`task` are those of the installed bmm (the vocabulary of
[`bmm_models()`](https://venpopov.com/bmm/dev/reference/bmm_models.md));
the rest is what the fit was built with.

Every element is a character, logical or numeric vector, a data frame,
or a list of these, so the list converts to JSON as it is. `citation`,
`domain`, `task` and the entries of `resp_vars` and `other_vars` are
vectors, which `auto_unbox = TRUE` writes as a single value when they
have one element and as an array otherwise. An `NA` in a data frame row
drops the key unless `na = "null"` is given. Use
`jsonlite::toJSON(spec, auto_unbox = TRUE, na = "null")`.

## Stability

The specification is the stable surface of a bmm model for other
packages (decision record 0005, "The model contract for companion and
extension packages"): "The contract is exactly what `model_spec()`
returns." and "The lifecycle rules apply to these items as they apply to
exported functions. Reading the model object with `$`, its attributes,
and internal functions stay outside the contract." A change to an
element is announced in the NEWS under "Companion packages" before the
release that makes it. Fields of a `bmmodel` object that `model_spec()`
does not return can change in any release.

## See also

Other model information:
[`bmm_models()`](https://venpopov.com/bmm/dev/reference/bmm_models.md),
[`parameter_info()`](https://venpopov.com/bmm/dev/reference/parameter_info.md),
[`simulate_data()`](https://venpopov.com/bmm/dev/reference/simulate_data.md)

## Examples

``` r
spec <- model_spec(sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials"))
spec
#> Yes/No Signal Detection Theory Model
#> Package: bmm
#> Version: none
#> Parameters:
#>  parameter     link fixed value alternative_group
#>          d identity FALSE    NA              <NA>
#>  criterion identity FALSE    NA    criterion|bias
#>       bias identity FALSE    NA    criterion|bias
#>    sdratio      log  TRUE     0              <NA>
#> Use `$` for the default priors, initial value ranges, data columns and functions.
spec$default_priors
#>   parameter             main        effects             sd
#> 1         d     normal(1, 1) normal(0, 0.5) exponential(1)
#> 2 criterion normal(0.5, 1.5) normal(0, 0.5) exponential(1)
#> 3      bias   normal(0, 1.5) normal(0, 0.5) exponential(2)
#> 4   sdratio   normal(0, 0.3) normal(0, 0.3) exponential(2)
spec$functions
#> $density
#> [1] "dsdt_yn"
#> 
#> $generator
#> [1] "rsdt_yn"
#> 

# mixture weights transformed together through a softmax
model_spec(mixture3p(resp_error = "y", nt_features = "nt", set_size = 4))$parameters
#>   parameter
#> 1        mu
#> 2     kappa
#> 3    thetat
#> 4   thetant
#>                                                                                                                                                                                                                                                               description
#> 1                                                                                                                                                              Location of the von Mises distribution of memory responses (in radians). Fixed internally to 0 by default.
#> 2                                                                                                                                                                                                         Concentration of the von Mises distribution of memory responses
#> 3 Mixture weight for target responses. The weights of the target, the non-targets and guessing are normalised across the components present on a trial, with guessing as the reference, so thetat maps to the probability of a target response only jointly with thetant.
#> 4                                                                                                                              Mixture weight for non-target responses, summed over the non-targets, which share it equally. Like thetat, it is weighed against guessing.
#>       link fixed value alternative_group softmax
#> 1 tan_half  TRUE     0              <NA>   FALSE
#> 2      log FALSE    NA              <NA>   FALSE
#> 3  softmax FALSE    NA              <NA>    TRUE
#> 4  softmax FALSE    NA              <NA>    TRUE
```
