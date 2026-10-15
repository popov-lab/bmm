# Transform parameter draws from the sampling scale to the native scale

Applies a `bmmodel`'s inverse link functions to posterior draws of its
parameters. This is the extension point used by
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md):
the transformation is a property of the *model*, not of an individual
parameter, because some models map several parameters jointly (e.g.
mixture weights through a softmax).

**A model whose parameters transform elementwise, or through one softmax
group that is active on every row of the data, needs no method** —
declaring `links` in its `.model_*()` constructor is sufficient, and
that covers every model `bmm` currently ships. Write a method when
either condition fails. Two cases that are known to fail it: a model
with several independent multinomial branches, because a `"softmax"`
link carries no group identity and all such parameters are softmaxed
together as one group; and a transformation that depends on the design
rather than only on the parameter's own value.

## Usage

``` r
native_transform(model, linpred, data, ...)

# Default S3 method
native_transform(model, linpred, data, ...)

# S3 method for class 'non_targets'
native_transform(model, linpred, data, ...)
```

## Arguments

- model:

  A `bmmodel` object.

- linpred:

  A named list of matrices of linear predictor draws, one per model
  parameter, each with draws in rows and prediction grid cells in
  columns. These are
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html)
  output, so they are on the *link* scale for distributional parameters
  and are the raw value for non-linear parameters. For `bmm`'s
  `nlpar`-based models the link is a fiction maintained by the model's
  own `nlf()` expression — `imm` declares `c = "log"` because
  `configure_model.imm_abc()` writes `exp(c)` into the formula, not
  because `brms` applies a link — so a model that declares a link it
  does not actually apply in its `nlf()` will get silently wrong output
  here.

- data:

  The prediction grid the draws were computed on, with one row per
  column of the `linpred` matrices. These are rows of the fitted model's
  own data, so they carry the `bmm`-internal columns a design-dependent
  transformation keys on (`LureIdx*`, `inv_ss`, `Idx_*`, matrix
  columns). `native_transform.non_targets()` is the worked example.

- ...:

  Currently unused.

## Value

A named list of matrices on the native scale. A method **must** return
one element for every name it was given and preserve each matrix's
dimensions;
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
checks both and errors, naming the model, if either is broken. The order
of the list is free.

## Details

Methods for specific models follow the `bmm` S3 chain, which runs from
the general classes to the specific ones. A method should therefore
transform the parameters it owns, remove them from `linpred`, and call
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) so the remaining
parameters reach the default method:

    native_transform.mymodel <- function(model, linpred, data, ...) {
      own <- .my_joint_transform(linpred[c("p1", "p2")])
      linpred <- linpred[not_in(names(linpred), c("p1", "p2"))]
      c(NextMethod(), own)
    }

Removing them before
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) is what keeps a
method idempotent with respect to the parameters it does not own: a
parameter left in `linpred` is transformed a second time by the next
method in the chain, which does not error and does not warn.

## See also

[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
