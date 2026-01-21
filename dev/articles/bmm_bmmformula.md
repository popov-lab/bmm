# Introduction to the bmmformula syntax

## 1 The `bmmformula` syntax

The `bmmformula` syntax generally follows the same principles as the
formula syntax for [`brms`
models](https://rdrr.io/cran/brms/man/brmsformula.html). A `bmmformula`
follows the general structure:

``` r
parameter ~ peffects + ( geffects | group)
```

However, instead of predicting a `response` - as is typically done in
`brms` - in `bmm` a `bmmformula` predicts a `parameter` from a
`bmmodel`[¹](#fn1). The `peffects` part specifies all effects that are
assumed to be the same across observations. In typical vocabulary such
effects are called ‘population-level’ or ‘overall’ effects. The
`geffects` part specifies effects that are assumed to vary across
grouping variables specified in `group`. Typically, these effects are
called ‘group-level’ or ‘varying’ effects.

As you will typically provide a `bmmformula` for each parameter of a
`bmmodel`, we recommend to set up the formula in a separate object
before passing it to the
[`bmm()`](https://venpopov.github.io/bmm/dev/reference/bmm.md) function.
There are two ways to set up formulas, both work equally easy and well:

1.  a single call to the function `bmmformula` or its short form `bmf`,
    separating formulas for different parameters by commas.
2.  add a `bmmformula` for each parameter by using the `+` operator.

Ηere are two examples for each method that result in the same
`bmmformula` object:

``` r
my_formula <- bmmformula(
  thetat ~ 1 + set_size,
  thetant ~ 1 + set_size,
  kappa ~ 1
)

my_formula <- bmf(thetat ~ 1 + set_size) +
  bmf(thetant ~ 1 + set_size) +
  bmf(kappa ~ 1)
  
my_formula
#> thetat ~ 1 + set_size
#> thetant ~ 1 + set_size
#> kappa ~ 1
```

As noted above, `bmm` generally assumes that you pass a `bmmformula` for
each parameter of a model. If you do not pass a `bmmformula` for one or
more parameters, `bmm` will throw a warning and only estimate a fixed
intercept for the parameters for which no `bmmformula` was provided. If
you are unsure about the parameters of a model, call the help via
`?bmmodel`, replacing ‘bmmodel’ with the name of the model
(e.g. [`?sdm`](https://venpopov.github.io/bmm/dev/reference/sdm.md)).

As `bmmformula` syntax builds upon `brmsformula` syntax, you can use any
functionality you might use in `brms` for specifying a `bmmformula`. One
difference is, that you *do not* have to explicitly specify if a formula
is a non-linear formula, as `bmm` recognizes a formula as non-linear as
soon as one of the left-hand side arguments of a formula is also used as
a right-hand side argument in the whole formula. So the argumetn `nl`
does not exist for `bmmformula`. Similarly, adding grouping statements
around group-level effects to estimate them separately for different
sub-groups in your design can be added. So for example, a more
complicated `bmmformula` could look like this:

``` r
complex_formula <- bmf(
  # a non-linear function on thetat
  thetat ~ end_theta + (start_theta - end_theta) * exp(-R*ss_num),
  # other model parameters
  thetant ~ 1 + (1 | gr(ID, by = age_group)),
  kappa ~ 1 + (1 | gr(ID, by = age_group)),
  # parameters of the non-linear function
  start_theta ~ 1 + (1 | gr(ID, by = age_group)),
  end_theta ~ 1 + (1 | gr(ID, by = age_group)),
  R ~ 1 + (1 | gr(ID, by = age_group))
)

complex_formula
#> thetat ~ end_theta + (start_theta - end_theta) * exp(-R * ss_num)
#> thetant ~ 1 + (1 | gr(ID, by = age_group))
#> kappa ~ 1 + (1 | gr(ID, by = age_group))
#> start_theta ~ 1 + (1 | gr(ID, by = age_group))
#> end_theta ~ 1 + (1 | gr(ID, by = age_group))
#> R ~ 1 + (1 | gr(ID, by = age_group))
```

This formula implements an exponential reduction in `thetat` from a
`start_theta` at `ss_num = 0` towards the lower asymptote `end_theta`
that will be reached with the slope `R`. In addition, the group level
effects have been grouped by `age_group` so that the variation in the
`Intercept` would be estimated separately for different levels of the
`age_group` variable. This formula just serves as an example how
`bmmformula` syntax allows to use the complex formula syntax implemented
in `brms` and passes the information on to `brms` without losing any of
the functionality.

## 2 Seperating response variables from model parameters

As the response is not provided to the `bmmformula`, in `bmm` the
response variables, and also other variables relevant for a `bmmodel`,
are linked to the model when setting up the `bmmodel` object. For
example, when fitting the `mixture3p` model, you have to provide the
names of the variables that reflect the response error (e.g., `y`), the
location of the non-target features (e.g. `nt_col1`, …, `nt_col4`), and
the set-size (e.g., `ss`) in our data:

``` r
my_model <- mixture3p(resp_error = "y", 
                      nt_features = paste0("nt_col",1:4), 
                      set_size = "ss")
```

Internally, `bmm` will then link the response variables to the
`bmmformula` given for the parameters of a `bmmodel`.

## 3 Fixing parameters to constant values

In some cases it can be reasonable to fix parameters to a specific
value. This can either be done via the `bmmformula` or via a [constant
prior](https://github.com/paul-buerkner/brms/issues/783#issuecomment-577909409)
specified for the parameter to be fixed. Here we will only outline
fixing parameters via the `bmmformula`.

To fix a parameter via the `bmmformula` you have to set the parameter to
the value it should be fixed to by using `=` instead of `~`. For
example, the following formula would fix `kappa` to the value `1`.

``` r
fixed_formula <- bmf(kappa = 1)
```

When fixing a parameter to a constant value, you need to keep in mind
that the value you fix it to, is the value on the parameter scale.
Oftentimes, `bmm` uses link functions to convert parameters from a
bounded native scale (e.g. only positive values for precision
parameters, such as `kappa`) to an unbounded parameter scale. In the
case of `kappa`, `bmm` uses a `log` link function. Thus, when fixing the
parameter to `1`, this will effectively fix the parameter to \\e^1\\ on
the native scale.

## 4 Estimating & prediciting parameters internally fixed by `bmm`

`bmm` fixes some parameters internally, as they are usually not an
integral part of a measurement model, or need to be fixed to properly
identify the model. For example in models for continuous reproduction
visual working memory tasks, the distribution of target responses is
usually assumed to be centered around the target. Consequently, the mean
of this distribution `mu` is internally fixed to `0`. Such parameters
are listed in the model object. So, you can see if a model has such
parameters by accessing `my_model$fixed_parameters`.

``` r
my_model <- sdm(resp_error = "error")
my_model$fixed_parameters
#> $mu
#> [1] 0
```

Fixed parameters that can be freely estimated are also listed in the
`parameters` of a model, that you can access using
`my_model$parameters`. Parameters that are only listed in the
`fixed_parameters` but not part of the `parameters` need to be fixed for
identification purposes and thus cannot be estimated freely. For example
in the `mixture` models for visual working memory, the location `mu2`
and precision `kappa2` of the guessing distributions needs to be fixed
for model identification. These parameters are listed in the
`fixed_parameters` but not in the `parameters` and therefore cannot be
estimated. Only, the internally fixed location of the target responses
`mu1` can be estimated freely if needed.

``` r
my_model <- mixture2p(resp_error = "error")
names(my_model$fixed_parameters)
#> [1] "mu1"    "mu2"    "kappa2"

names(my_model$parameters)
#> [1] "mu1"    "kappa"  "thetat"
```

If you want to freely estimate an internally fixed parameters, then you
can freely estimate this parameter by providing a `bmmformula` for it.
In the case of the above mentioned visual working memory models, you
could for example be interested in seeing if subjects were biased (e.g.,
by distractors). If you want to estimate the overall bias in your data
including variations between subjects indicated by `ID`, then you would
additionally specify a `bmmformula` for `mu` in addition to the other
formulas for the core mode parameters:

``` r
my_model <- mixture2p(resp_error = "error")
my_formula <- bmf(
  mu1 ~ 1 + (1 | ID),
  thetat ~ 0 + set_size + (0 + set_size | ID),
  kappa ~ 0 + set_size + (0 + set_size | ID)
)
```

## 5 Accessing the `brmsformula` generated by `bmm`

Let’s assume you have fit a `mixture2p` model to a data set `my_data` as
described in the [the mixture models
article](https://venpopov.github.io/bmm/articles/bmm_mixture_models.html%22).
For this, you had to specify a `bmmformula` and model object:

``` r
user_formula <- bmf(
  thetat ~ 0 + set_size + (0 + set_size | id),
  kappa ~ 0 + set_size + (0 + set_size | id)
)

my_model <- mixture2p(resp_error = "error")

bmmfit <- bmm(
  formula = user_formula,
  data = my_data,
  model = my_model,
  file = "assets/bmmfit_mixture2p_vignette"
)
```

If you were now interested to see how the `bmmformula` that you passed
to `bmm` is converted to a `brmsformula`, you can access the
`brmsformula` via the `bmmfit` object that will be returned by the
[`bmm()`](https://venpopov.github.io/bmm/dev/reference/bmm.md) function:

``` r
bmmfit$formula
#> error ~ 1 
#> mu1 ~ 1
#> kappa ~ 0 + set_size + (0 + set_size || id)
#> thetat ~ 0 + set_size + (0 + set_size || id)
#> kappa2 ~ 1
#> mu2 ~ 1
#> kappa1 ~ kappa
#> theta1 ~ thetat
```

Similarly, you can access the distributional family that was used to
implement the specified model, so `bmmfit$family` will return the family
object that was generated by `bmm` and then passed to `brms` for
fitting.

``` r
bmmfit$family
#> 
#> Mixture
#> 
#> Family: von_mises 
#> Link function: tan_half 
#> 
#> Family: von_mises 
#> Link function: tan_half
```

Finally, you can also access the data that was used by `bmm` to fit the
model via `bmmfit$data`. In many cases this data will be equal to the
data provided by the user. But sometimes `bmm` internally computes
additional index variables for specifying the models adequately. The
data stored in the `bmmfit` object contains these additional variables.

``` r
head(bmmfit$data)
#>        error set_size id
#> 1 -0.0130000        1  1
#> 2 -0.0660000        1  1
#> 3  0.1891853        1  1
#> 4 -0.4500000        1  1
#> 5 -0.0210000        1  1
#> 6 -0.0740000        1  1
```

This way, should you be interested to customize models to fit them in
`brms` without `bmm`, you are able to obtain the most important
information that is essential for specifying the models implemented in
`bmm`.

Fit objects from `bmm` use a custom summary function to format the
output. However, you can call a `brms` summary instead of a `bmm`
summary by setting the `backend` option in the summary function to
`brms`. This way the summary function for `brmsfit` objects will be used
instead of the function for `bmmfit` objects. The `brmsfit` summary also
contains the `brmsformula` that is created by `bmm`.

``` r
summary(bmmfit, backend = "brms")
#> Warning: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: mixture(von_mises, von_mises) 
#>   Links: mu1 = tan_half; mu2 = tan_half; kappa2 = log; kappa1 = log; theta1 = identity 
#> Formula: error ~ 1 
#>          mu1 ~ 1
#>          kappa ~ 0 + set_size + (0 + set_size || id)
#>          thetat ~ 0 + set_size + (0 + set_size || id)
#>          kappa2 ~ 1
#>          mu2 ~ 1
#>          kappa1 ~ kappa
#>          theta1 ~ thetat
#>    Data: dat_preprocessed (Number of observations: 7271) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Multilevel Hyperparameters:
#> ~id (Number of levels: 12) 
#>                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sd(kappa_set_size1)      0.36      0.10     0.22     0.60 1.00     1267     2370
#> sd(kappa_set_size2)      0.21      0.08     0.08     0.39 1.00     1302     1020
#> sd(kappa_set_size4)      0.34      0.11     0.18     0.60 1.00     1924     2711
#> sd(kappa_set_size6)      0.43      0.15     0.21     0.79 1.00     1779     2444
#> sd(thetat_set_size1)     0.57      0.45     0.03     1.73 1.00     1265     1452
#> sd(thetat_set_size2)     0.93      0.29     0.51     1.64 1.00     1747     2568
#> sd(thetat_set_size4)     0.95      0.27     0.55     1.63 1.00     1412     2109
#> sd(thetat_set_size6)     0.67      0.20     0.39     1.15 1.00     1643     2310
#> 
#> Regression Coefficients:
#>                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> mu1_Intercept        0.00      0.00     0.00     0.00   NA       NA       NA
#> mu2_Intercept        0.00      0.00     0.00     0.00   NA       NA       NA
#> kappa2_Intercept  -100.00      0.00  -100.00  -100.00   NA       NA       NA
#> kappa_set_size1      2.89      0.12     2.66     3.12 1.00     1406     1716
#> kappa_set_size2      2.40      0.08     2.24     2.56 1.00     2632     2527
#> kappa_set_size4      2.08      0.12     1.84     2.32 1.00     2426     2690
#> kappa_set_size6      1.94      0.15     1.65     2.25 1.00     2517     2627
#> thetat_set_size1     4.53      0.40     3.89     5.43 1.00     2316     1088
#> thetat_set_size2     2.56      0.31     1.94     3.17 1.00     1658     1999
#> thetat_set_size4     1.08      0.29     0.50     1.67 1.00     1312     1735
#> thetat_set_size6     0.32      0.21    -0.08     0.75 1.00     1502     2120
#> 
#> Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
```

Compare this with the default summary method for `bmmfit` objects:

``` r
summary(bmmfit)
Warning: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
```

``` fansi
  Model: mixture2p(resp_error = "error") 
  Links: mu1 = tan_half; kappa = log; thetat = identity 
Formula: mu1 = 0
         kappa ~ 0 + set_size + (0 + set_size || id)
         thetat ~ 0 + set_size + (0 + set_size || id) 
   Data: dat_preprocessed (Number of observations: 7271)
  Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup draws = 4000

Multilevel Hyperparameters:
~id (Number of levels: 12) 
                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(kappa_set_size1)      0.36      0.10     0.22     0.60 1.00     1267     2370
sd(kappa_set_size2)      0.21      0.08     0.08     0.39 1.00     1302     1020
sd(kappa_set_size4)      0.34      0.11     0.18     0.60 1.00     1924     2711
sd(kappa_set_size6)      0.43      0.15     0.21     0.79 1.00     1779     2444
sd(thetat_set_size1)     0.57      0.45     0.03     1.73 1.00     1265     1452
sd(thetat_set_size2)     0.93      0.29     0.51     1.64 1.00     1747     2568
sd(thetat_set_size4)     0.95      0.27     0.55     1.63 1.00     1412     2109
sd(thetat_set_size6)     0.67      0.20     0.39     1.15 1.00     1643     2310

Regression Coefficients:
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
kappa_set_size1      2.89      0.12     2.66     3.12 1.00     1406     1716
kappa_set_size2      2.40      0.08     2.24     2.56 1.00     2632     2527
kappa_set_size4      2.08      0.12     1.84     2.32 1.00     2426     2690
kappa_set_size6      1.94      0.15     1.65     2.25 1.00     2517     2627
thetat_set_size1     4.53      0.40     3.89     5.43 1.00     2316     1088
thetat_set_size2     2.56      0.31     1.94     3.17 1.00     1658     1999
thetat_set_size4     1.08      0.29     0.50     1.67 1.00     1312     1735
thetat_set_size6     0.32      0.21    -0.08     0.75 1.00     1502     2120

Constant Parameters:
                  Value
mu1_Intercept      0.00

Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

------------------------------------------------------------------------

1.  You can always find information on the parameters of a `bmmodel` by
    calling the help for the model you want to fit: `?bmmodel`
