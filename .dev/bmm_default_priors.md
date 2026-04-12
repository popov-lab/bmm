Default priors for bmm models
================

## Introduction

Here I explore options for how to set up default bmm prior flexibly
based on the formula.

## How it is done in `brms`

First, we can see what happens by default in `brms`. Let’s look at a
simple linear regression.

<details>
<summary>
Click to expand
</summary>

### Intercept only

It puts a default prior on all parameters that have only an intercept

``` r
dat <- oberauer_lin_2017
dat$cond <- factor(rep(1:4, each=nrow(dat)/4))   # fake condition for testing formulas
get_prior(bf(dev_rad ~ 1, sigma ~ 1), dat)
```

    ##                 prior     class coef group resp  dpar nlpar lb ub  source
    ##  student_t(3, 0, 2.5) Intercept                                   default
    ##  student_t(3, 0, 2.5) Intercept                 sigma             default

### Intercept and a predictor

The intercept gets a default prior, but the effects get a flat prior

``` r
get_prior(bf(dev_rad ~ 1, sigma ~ set_size), dat)
```

    ##                 prior     class      coef group resp  dpar nlpar lb ub       source
    ##  student_t(3, 0, 2.5) Intercept                                             default
    ##                (flat)         b                      sigma                  default
    ##                (flat)         b set_size2            sigma             (vectorized)
    ##                (flat)         b set_size3            sigma             (vectorized)
    ##                (flat)         b set_size4            sigma             (vectorized)
    ##                (flat)         b set_size5            sigma             (vectorized)
    ##                (flat)         b set_size6            sigma             (vectorized)
    ##                (flat)         b set_size7            sigma             (vectorized)
    ##                (flat)         b set_size8            sigma             (vectorized)
    ##  student_t(3, 0, 2.5) Intercept                      sigma                  default

### Intercept supressed

Just flat prior on all coefficients

``` r
get_prior(bf(dev_rad ~ 1, sigma ~ 0+set_size), dat)
```

    ##                 prior     class      coef group resp  dpar nlpar lb ub       source
    ##  student_t(3, 0, 2.5) Intercept                                             default
    ##                (flat)         b                      sigma                  default
    ##                (flat)         b set_size1            sigma             (vectorized)
    ##                (flat)         b set_size2            sigma             (vectorized)
    ##                (flat)         b set_size3            sigma             (vectorized)
    ##                (flat)         b set_size4            sigma             (vectorized)
    ##                (flat)         b set_size5            sigma             (vectorized)
    ##                (flat)         b set_size6            sigma             (vectorized)
    ##                (flat)         b set_size7            sigma             (vectorized)
    ##                (flat)         b set_size8            sigma             (vectorized)

### Intercept supressed and random effects

``` r
get_prior(bf(dev_rad ~ 1, sigma ~ 0+set_size + (0+set_size|ID)), dat)
```

    ##                 prior     class      coef group resp  dpar nlpar lb ub       source
    ##                lkj(1)       cor                                             default
    ##                lkj(1)       cor              ID                        (vectorized)
    ##  student_t(3, 0, 2.5) Intercept                                             default
    ##                (flat)         b                      sigma                  default
    ##                (flat)         b set_size1            sigma             (vectorized)
    ##                (flat)         b set_size2            sigma             (vectorized)
    ##                (flat)         b set_size3            sigma             (vectorized)
    ##                (flat)         b set_size4            sigma             (vectorized)
    ##                (flat)         b set_size5            sigma             (vectorized)
    ##                (flat)         b set_size6            sigma             (vectorized)
    ##                (flat)         b set_size7            sigma             (vectorized)
    ##                (flat)         b set_size8            sigma             (vectorized)
    ##  student_t(3, 0, 2.5)        sd                      sigma        0         default
    ##  student_t(3, 0, 2.5)        sd              ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size1    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size2    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size3    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size4    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size5    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size6    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size7    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size8    ID      sigma        0    (vectorized)

### Intercept predictor and random effects

``` r
get_prior(bf(dev_rad ~ 1, sigma ~ set_size + (set_size|ID)), dat)
```

    ##                 prior     class      coef group resp  dpar nlpar lb ub       source
    ##                lkj(1)       cor                                             default
    ##                lkj(1)       cor              ID                        (vectorized)
    ##  student_t(3, 0, 2.5) Intercept                                             default
    ##                (flat)         b                      sigma                  default
    ##                (flat)         b set_size2            sigma             (vectorized)
    ##                (flat)         b set_size3            sigma             (vectorized)
    ##                (flat)         b set_size4            sigma             (vectorized)
    ##                (flat)         b set_size5            sigma             (vectorized)
    ##                (flat)         b set_size6            sigma             (vectorized)
    ##                (flat)         b set_size7            sigma             (vectorized)
    ##                (flat)         b set_size8            sigma             (vectorized)
    ##  student_t(3, 0, 2.5) Intercept                      sigma                  default
    ##  student_t(3, 0, 2.5)        sd                      sigma        0         default
    ##  student_t(3, 0, 2.5)        sd              ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd Intercept    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size2    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size3    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size4    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size5    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size6    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size7    ID      sigma        0    (vectorized)
    ##  student_t(3, 0, 2.5)        sd set_size8    ID      sigma        0    (vectorized)

</details>

## How should we do it in `bmm`

Currently we put a default prior on all parameters, assuming naively
that they are specified with a suppressed intercept (except for the
sdm, for which we have no default priors yet). We don’t want to
follow `brms` approach, because these measurement models need meaningful
default priors to help with sampling.

Options: 1) we do it like `brms` and let the user specify priors for all
parameters 2) we put reasonable priors on intercepts and on all effects
when an intercept is suprresed, but leave flat priors for effects when
an intercept is present 3) we put reasonable priors on all parameters,
regardless of the formula

I don’t like option 1, because it is not user friendly. For my own
models I don’t want to have to specify the priors every time, unless I
want priors that are different from the defaults. Options 2 and 3 are
both ok, but we need to figure out when to put the priors and when not.

## Prior structure for models with non-linear parameters (e.g. mixture3p)

In this section I just print the brms prior structure for a variety of
different formulas. I turn off the bmm priors, and also print the
structure of `terms(formula)`, which can help identify what information
we need to figure out how to set default priors.

<details>
<summary>
Click to expand
</summary>

Disable currently used default priors:

``` r
options(list(bmm.default_priors = FALSE))
```

### Intercept only

All model parameters are `nlpar` so they get class `b` with coef
`Intercept`

``` r
model <- mixture3p('dev_rad', nt_features = paste0('col_nt',1:7), set_size='set_size')
formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##           prior     class      coef group resp   dpar   nlpar   lb   ub       source
    ##  logistic(0, 1)    theta9                                     -Inf  Inf      default
    ##          (flat)         b                               kappa                default
    ##          (flat)         b Intercept                     kappa           (vectorized)
    ##          (flat)         b                             thetant                default
    ##          (flat)         b Intercept                   thetant           (vectorized)
    ##          (flat)         b                              thetat                default
    ##          (flat)         b Intercept                    thetat           (vectorized)
    ##     constant(0) Intercept                         mu1         <NA> <NA>         user
    ##  constant(-100) Intercept                      kappa9         <NA> <NA>         user
    ##     constant(0) Intercept                         mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ 1
    ## attr(,"variables")
    ## list(kappa)
    ## attr(,"factors")
    ## integer(0)
    ## attr(,"term.labels")
    ## character(0)
    ## attr(,"order")
    ## integer(0)
    ## attr(,"intercept")
    ## [1] 1
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

### Intercept and a predictor

For kappa, which we include and intercept and a predictor, we get class
`b` with coef `Intercept` and `session2`

``` r
formula <- bmf(kappa ~ session, thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##           prior     class      coef group resp   dpar   nlpar   lb   ub       source
    ##  logistic(0, 1)    theta9                                     -Inf  Inf      default
    ##          (flat)         b                               kappa                default
    ##          (flat)         b Intercept                     kappa           (vectorized)
    ##          (flat)         b  session2                     kappa           (vectorized)
    ##          (flat)         b                             thetant                default
    ##          (flat)         b Intercept                   thetant           (vectorized)
    ##          (flat)         b                              thetat                default
    ##          (flat)         b Intercept                    thetat           (vectorized)
    ##     constant(0) Intercept                         mu1         <NA> <NA>         user
    ##  constant(-100) Intercept                      kappa9         <NA> <NA>         user
    ##     constant(0) Intercept                         mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ session
    ## attr(,"variables")
    ## list(kappa, session)
    ## attr(,"factors")
    ##         session
    ## kappa         0
    ## session       1
    ## attr(,"term.labels")
    ## [1] "session"
    ## attr(,"order")
    ## [1] 1
    ## attr(,"intercept")
    ## [1] 1
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

### Intercept supressed

For kappa, which we include and intercept and a predictor, we get class
`b` with coef `session1` and `session2`

``` r
formula <- bmf(kappa ~ 0+session, thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##           prior     class      coef group resp   dpar   nlpar   lb   ub       source
    ##  logistic(0, 1)    theta9                                     -Inf  Inf      default
    ##          (flat)         b                               kappa                default
    ##          (flat)         b  session1                     kappa           (vectorized)
    ##          (flat)         b  session2                     kappa           (vectorized)
    ##          (flat)         b                             thetant                default
    ##          (flat)         b Intercept                   thetant           (vectorized)
    ##          (flat)         b                              thetat                default
    ##          (flat)         b Intercept                    thetat           (vectorized)
    ##     constant(0) Intercept                         mu1         <NA> <NA>         user
    ##  constant(-100) Intercept                      kappa9         <NA> <NA>         user
    ##     constant(0) Intercept                         mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ 0 + session
    ## attr(,"variables")
    ## list(kappa, session)
    ## attr(,"factors")
    ##         session
    ## kappa         0
    ## session       1
    ## attr(,"term.labels")
    ## [1] "session"
    ## attr(,"order")
    ## [1] 1
    ## attr(,"intercept")
    ## [1] 0
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

### Intercept supressed and random effects

``` r
formula <- bmf(kappa ~ 0+session + ( 0+session|ID), thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##                 prior     class      coef group resp   dpar   nlpar   lb   ub       source
    ##                lkj(1)       cor                                                    default
    ##                lkj(1)       cor              ID                               (vectorized)
    ##        logistic(0, 1)    theta9                                     -Inf  Inf      default
    ##                (flat)         b                               kappa                default
    ##                (flat)         b  session1                     kappa           (vectorized)
    ##                (flat)         b  session2                     kappa           (vectorized)
    ##  student_t(3, 0, 2.5)        sd                               kappa    0           default
    ##  student_t(3, 0, 2.5)        sd              ID               kappa    0      (vectorized)
    ##  student_t(3, 0, 2.5)        sd  session1    ID               kappa    0      (vectorized)
    ##  student_t(3, 0, 2.5)        sd  session2    ID               kappa    0      (vectorized)
    ##                (flat)         b                             thetant                default
    ##                (flat)         b Intercept                   thetant           (vectorized)
    ##                (flat)         b                              thetat                default
    ##                (flat)         b Intercept                    thetat           (vectorized)
    ##           constant(0) Intercept                         mu1         <NA> <NA>         user
    ##        constant(-100) Intercept                      kappa9         <NA> <NA>         user
    ##           constant(0) Intercept                         mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ 0 + session + (0 + session | ID)
    ## attr(,"variables")
    ## list(kappa, session, 0 + session | ID)
    ## attr(,"factors")
    ##                  session 0 + session | ID
    ## kappa                  0                0
    ## session                1                0
    ## 0 + session | ID       0                1
    ## attr(,"term.labels")
    ## [1] "session"          "0 + session | ID"
    ## attr(,"order")
    ## [1] 1 1
    ## attr(,"intercept")
    ## [1] 0
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

### Two factors plus intercept

``` r
formula <- bmf(kappa ~ session + cond, thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##           prior     class      coef group resp   dpar   nlpar   lb   ub       source
    ##  logistic(0, 1)    theta9                                     -Inf  Inf      default
    ##          (flat)         b                               kappa                default
    ##          (flat)         b     cond2                     kappa           (vectorized)
    ##          (flat)         b     cond3                     kappa           (vectorized)
    ##          (flat)         b     cond4                     kappa           (vectorized)
    ##          (flat)         b Intercept                     kappa           (vectorized)
    ##          (flat)         b  session2                     kappa           (vectorized)
    ##          (flat)         b                             thetant                default
    ##          (flat)         b Intercept                   thetant           (vectorized)
    ##          (flat)         b                              thetat                default
    ##          (flat)         b Intercept                    thetat           (vectorized)
    ##     constant(0) Intercept                         mu1         <NA> <NA>         user
    ##  constant(-100) Intercept                      kappa9         <NA> <NA>         user
    ##     constant(0) Intercept                         mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ session + cond
    ## attr(,"variables")
    ## list(kappa, session, cond)
    ## attr(,"factors")
    ##         session cond
    ## kappa         0    0
    ## session       1    0
    ## cond          0    1
    ## attr(,"term.labels")
    ## [1] "session" "cond"   
    ## attr(,"order")
    ## [1] 1 1
    ## attr(,"intercept")
    ## [1] 1
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

### Two factors no intercept

Ok, so this is tricky. coefs session1 and session2 are main levels while
cond=1, but cond2,cond3 and cond4 are additive effects. So we can’t just
check if an intercept is present in coef or not.

What we can do, is for each predictor, check how many levels appear.

``` r
formula <- bmf(kappa ~ 0 + session + cond, thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##           prior     class      coef group resp   dpar   nlpar   lb   ub       source
    ##  logistic(0, 1)    theta9                                     -Inf  Inf      default
    ##          (flat)         b                               kappa                default
    ##          (flat)         b     cond2                     kappa           (vectorized)
    ##          (flat)         b     cond3                     kappa           (vectorized)
    ##          (flat)         b     cond4                     kappa           (vectorized)
    ##          (flat)         b  session1                     kappa           (vectorized)
    ##          (flat)         b  session2                     kappa           (vectorized)
    ##          (flat)         b                             thetant                default
    ##          (flat)         b Intercept                   thetant           (vectorized)
    ##          (flat)         b                              thetat                default
    ##          (flat)         b Intercept                    thetat           (vectorized)
    ##     constant(0) Intercept                         mu1         <NA> <NA>         user
    ##  constant(-100) Intercept                      kappa9         <NA> <NA>         user
    ##     constant(0) Intercept                         mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ 0 + session + cond
    ## attr(,"variables")
    ## list(kappa, session, cond)
    ## attr(,"factors")
    ##         session cond
    ## kappa         0    0
    ## session       1    0
    ## cond          0    1
    ## attr(,"term.labels")
    ## [1] "session" "cond"   
    ## attr(,"order")
    ## [1] 1 1
    ## attr(,"intercept")
    ## [1] 0
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

### Two factors no intercept, yes interaction

``` r
formula <- bmf(kappa ~ 0 + session * cond, thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##           prior     class           coef group resp   dpar   nlpar   lb   ub       source
    ##  logistic(0, 1)    theta9                                          -Inf  Inf      default
    ##          (flat)         b                                    kappa                default
    ##          (flat)         b          cond2                     kappa           (vectorized)
    ##          (flat)         b          cond3                     kappa           (vectorized)
    ##          (flat)         b          cond4                     kappa           (vectorized)
    ##          (flat)         b       session1                     kappa           (vectorized)
    ##          (flat)         b       session2                     kappa           (vectorized)
    ##          (flat)         b session2:cond2                     kappa           (vectorized)
    ##          (flat)         b session2:cond3                     kappa           (vectorized)
    ##          (flat)         b session2:cond4                     kappa           (vectorized)
    ##          (flat)         b                                  thetant                default
    ##          (flat)         b      Intercept                   thetant           (vectorized)
    ##          (flat)         b                                   thetat                default
    ##          (flat)         b      Intercept                    thetat           (vectorized)
    ##     constant(0) Intercept                              mu1         <NA> <NA>         user
    ##  constant(-100) Intercept                           kappa9         <NA> <NA>         user
    ##     constant(0) Intercept                              mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ 0 + session * cond
    ## attr(,"variables")
    ## list(kappa, session, cond)
    ## attr(,"factors")
    ##         session cond session:cond
    ## kappa         0    0            0
    ## session       1    0            1
    ## cond          0    1            1
    ## attr(,"term.labels")
    ## [1] "session"      "cond"         "session:cond"
    ## attr(,"order")
    ## [1] 1 1 2
    ## attr(,"intercept")
    ## [1] 0
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

### Two factors no intercept, yes interaction

``` r
formula <- bmf(kappa ~ 0 + session:cond + cond + session, thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##           prior     class           coef group resp   dpar   nlpar   lb   ub       source
    ##  logistic(0, 1)    theta9                                          -Inf  Inf      default
    ##          (flat)         b                                    kappa                default
    ##          (flat)         b          cond1                     kappa           (vectorized)
    ##          (flat)         b          cond2                     kappa           (vectorized)
    ##          (flat)         b cond2:session2                     kappa           (vectorized)
    ##          (flat)         b          cond3                     kappa           (vectorized)
    ##          (flat)         b cond3:session2                     kappa           (vectorized)
    ##          (flat)         b          cond4                     kappa           (vectorized)
    ##          (flat)         b cond4:session2                     kappa           (vectorized)
    ##          (flat)         b       session2                     kappa           (vectorized)
    ##          (flat)         b                                  thetant                default
    ##          (flat)         b      Intercept                   thetant           (vectorized)
    ##          (flat)         b                                   thetat                default
    ##          (flat)         b      Intercept                    thetat           (vectorized)
    ##     constant(0) Intercept                              mu1         <NA> <NA>         user
    ##  constant(-100) Intercept                           kappa9         <NA> <NA>         user
    ##     constant(0) Intercept                              mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ 0 + session:cond + cond + session
    ## attr(,"variables")
    ## list(kappa, session, cond)
    ## attr(,"factors")
    ##         cond session session:cond
    ## kappa      0       0            0
    ## session    0       1            1
    ## cond       1       0            1
    ## attr(,"term.labels")
    ## [1] "cond"         "session"      "session:cond"
    ## attr(,"order")
    ## [1] 1 1 2
    ## attr(,"intercept")
    ## [1] 0
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

### Two factors no intercept, no main effects, all levels explicit

``` r
formula <- bmf(kappa ~ 0 + session:cond, thetat ~ 1, thetant ~ 1)
get_model_prior(formula, dat, model)
```

    ##           prior     class           coef group resp   dpar   nlpar   lb   ub       source
    ##  logistic(0, 1)    theta9                                          -Inf  Inf      default
    ##          (flat)         b                                    kappa                default
    ##          (flat)         b session1:cond1                     kappa           (vectorized)
    ##          (flat)         b session1:cond2                     kappa           (vectorized)
    ##          (flat)         b session1:cond3                     kappa           (vectorized)
    ##          (flat)         b session1:cond4                     kappa           (vectorized)
    ##          (flat)         b session2:cond1                     kappa           (vectorized)
    ##          (flat)         b session2:cond2                     kappa           (vectorized)
    ##          (flat)         b session2:cond3                     kappa           (vectorized)
    ##          (flat)         b session2:cond4                     kappa           (vectorized)
    ##          (flat)         b                                  thetant                default
    ##          (flat)         b      Intercept                   thetant           (vectorized)
    ##          (flat)         b                                   thetat                default
    ##          (flat)         b      Intercept                    thetat           (vectorized)
    ##     constant(0) Intercept                              mu1         <NA> <NA>         user
    ##  constant(-100) Intercept                           kappa9         <NA> <NA>         user
    ##     constant(0) Intercept                              mu9         <NA> <NA>         user

``` r
terms(formula$kappa)
```

    ## kappa ~ 0 + session:cond
    ## attr(,"variables")
    ## list(kappa, session, cond)
    ## attr(,"factors")
    ##         session:cond
    ## kappa              0
    ## session            2
    ## cond               2
    ## attr(,"term.labels")
    ## [1] "session:cond"
    ## attr(,"order")
    ## [1] 2
    ## attr(,"intercept")
    ## [1] 0
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>

</details>

## Algorithm

Looking at the structure of the priors and the terms(formula\$kappa),
here is a suggested algorithm for the initial `set_default_prior`
function:

1)  Get dpar names of bmmformula

2)  Loop over dpar names

3)  Get terms from the formula component with stats::terms

4)  Get the value of \`attr(terms, “intercept”)

    - if 1, then:
      - set the default prior on class `b`, coef `Intercept`
      - set wide priors centered on 0 for the remaining coefficients
      - go to step 2 and continue with next parameter
    - if 0, continue to step 5

5)  check number of predictors via rhs_vars

    - if 1, then:
      - set the default prior on class `b`
      - go to step 2 and continue with next parameter
    - if \>1, continue to step 6

6)  check if there is only interaction (e.g. <session:cond>) via
    stats::terms, length(attr(,“order”)) == 1 & attr(,“order”) == 2

    - if true, then:
      - set the default prior on class `b`
      - go to step 2 and continue with next parameter
    - if false, continue to step 7

7)  find the first term in stats::terms attr(,“term.labels”):

    - set the default prior on class `b` with coefs equal to that
      predictor with the level appended

## Mixture3p: Comparison with new function “set_default_prior”

A helper function to compare the priors with and without the default
priors.

``` r
compare_priors <- function(formula,dat,model) {
  p1 <- withr::with_options(list('bmm.default_priors'=FALSE), get_model_prior(formula, dat, model)) |> 
    brms:::prepare_print_prior() |> 
    rename(priorBRMS=prior) |> 
    select(priorBRMS,class,coef,group,resp,dpar,nlpar)
  p2 <- withr::with_options(list('bmm.default_priors'=TRUE), get_model_prior(formula, dat, model)) |> 
    brms:::prepare_print_prior() |> 
    rename(priorBMM=prior) |> 
    select(priorBMM,class,coef,group,resp,dpar,nlpar)
  p <- left_join(p1,p2, by=c('class','coef','group','resp','dpar','nlpar')) |>
    select(priorBRMS,priorBMM,class,coef,group,resp,dpar,nlpar) |> 
    arrange(nlpar, dpar, class, priorBMM, coef)
  p
}
```

``` r
model3p <- mixture3p('dev_rad', nt_features = paste0('col_nt',1:7), set_size='set_size')
```

<details>
<summary>
Click to expand
</summary>

### Intercept only

All model parameters are `nlpar` so they get class `b` with coef
`Intercept`

``` r
formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class      coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                     (unknown)
    ##  constant(-100) constant(-100) Intercept                      kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu9         (unknown)
    ##          (flat)         (flat)         b                               kappa (unknown)
    ##          (flat)    normal(2,1)         b Intercept                     kappa (unknown)
    ##          (flat)         (flat)         b                             thetant (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                              thetat (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                    thetat (unknown)

### Intercept and a predictor

For kappa, which we include and intercept and a predictor, we get class
`b` with coef `Intercept` and `session2`

``` r
formula <- bmf(kappa ~ session, thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class      coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                     (unknown)
    ##  constant(-100) constant(-100) Intercept                      kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu9         (unknown)
    ##          (flat)    normal(0,1)         b                               kappa (unknown)
    ##          (flat)    normal(0,1)         b  session2                     kappa (unknown)
    ##          (flat)    normal(2,1)         b Intercept                     kappa (unknown)
    ##          (flat)         (flat)         b                             thetant (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                              thetat (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                    thetat (unknown)

### Intercept supressed

For kappa, which we include and intercept and a predictor, we get class
`b` with coef `session1` and `session2`

``` r
formula <- bmf(kappa ~ 0+session, thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class      coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                     (unknown)
    ##  constant(-100) constant(-100) Intercept                      kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu9         (unknown)
    ##          (flat)    normal(2,1)         b                               kappa (unknown)
    ##          (flat)    normal(2,1)         b  session1                     kappa (unknown)
    ##          (flat)    normal(2,1)         b  session2                     kappa (unknown)
    ##          (flat)         (flat)         b                             thetant (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                              thetat (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                    thetat (unknown)

### Intercept supressed and random effects

``` r
formula <- bmf(kappa ~ 0+session + ( 0+session|ID), thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##             priorBRMS             priorBMM     class      coef group resp   dpar   nlpar    source
    ##                lkj(1)               lkj(1)       cor                                     (unknown)
    ##                lkj(1)               lkj(1)       cor              ID                     (unknown)
    ##        logistic(0, 1)       logistic(0, 1)    theta9                                     (unknown)
    ##        constant(-100)       constant(-100) Intercept                      kappa9         (unknown)
    ##           constant(0)          constant(0) Intercept                         mu1         (unknown)
    ##           constant(0)          constant(0) Intercept                         mu9         (unknown)
    ##                (flat)          normal(2,1)         b                               kappa (unknown)
    ##                (flat)          normal(2,1)         b  session1                     kappa (unknown)
    ##                (flat)          normal(2,1)         b  session2                     kappa (unknown)
    ##  student_t(3, 0, 2.5) student_t(3, 0, 2.5)        sd                               kappa (unknown)
    ##  student_t(3, 0, 2.5) student_t(3, 0, 2.5)        sd              ID               kappa (unknown)
    ##  student_t(3, 0, 2.5) student_t(3, 0, 2.5)        sd  session1    ID               kappa (unknown)
    ##  student_t(3, 0, 2.5) student_t(3, 0, 2.5)        sd  session2    ID               kappa (unknown)
    ##                (flat)               (flat)         b                             thetant (unknown)
    ##                (flat)       logistic(0, 1)         b Intercept                   thetant (unknown)
    ##                (flat)               (flat)         b                              thetat (unknown)
    ##                (flat)       logistic(0, 1)         b Intercept                    thetat (unknown)

### Two factors plus intercept

``` r
formula <- bmf(kappa ~ session + cond, thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class      coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                     (unknown)
    ##  constant(-100) constant(-100) Intercept                      kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu9         (unknown)
    ##          (flat)    normal(0,1)         b                               kappa (unknown)
    ##          (flat)    normal(0,1)         b     cond2                     kappa (unknown)
    ##          (flat)    normal(0,1)         b     cond3                     kappa (unknown)
    ##          (flat)    normal(0,1)         b     cond4                     kappa (unknown)
    ##          (flat)    normal(0,1)         b  session2                     kappa (unknown)
    ##          (flat)    normal(2,1)         b Intercept                     kappa (unknown)
    ##          (flat)         (flat)         b                             thetant (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                              thetat (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                    thetat (unknown)

### Two factors no intercept

``` r
formula <- bmf(kappa ~ 0 + session + cond, thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class      coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                     (unknown)
    ##  constant(-100) constant(-100) Intercept                      kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu9         (unknown)
    ##          (flat)    normal(0,1)         b                               kappa (unknown)
    ##          (flat)    normal(0,1)         b     cond2                     kappa (unknown)
    ##          (flat)    normal(0,1)         b     cond3                     kappa (unknown)
    ##          (flat)    normal(0,1)         b     cond4                     kappa (unknown)
    ##          (flat)    normal(2,1)         b  session1                     kappa (unknown)
    ##          (flat)    normal(2,1)         b  session2                     kappa (unknown)
    ##          (flat)         (flat)         b                             thetant (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                              thetat (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                    thetat (unknown)

### Two factors no intercept, yes interaction

``` r
formula <- bmf(kappa ~ 0 + session * cond, thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class           coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                          (unknown)
    ##  constant(-100) constant(-100) Intercept                           kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                              mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                              mu9         (unknown)
    ##          (flat)    normal(0,1)         b                                    kappa (unknown)
    ##          (flat)    normal(0,1)         b          cond2                     kappa (unknown)
    ##          (flat)    normal(0,1)         b          cond3                     kappa (unknown)
    ##          (flat)    normal(0,1)         b          cond4                     kappa (unknown)
    ##          (flat)    normal(0,1)         b session2:cond2                     kappa (unknown)
    ##          (flat)    normal(0,1)         b session2:cond3                     kappa (unknown)
    ##          (flat)    normal(0,1)         b session2:cond4                     kappa (unknown)
    ##          (flat)    normal(2,1)         b       session1                     kappa (unknown)
    ##          (flat)    normal(2,1)         b       session2                     kappa (unknown)
    ##          (flat)         (flat)         b                                  thetant (unknown)
    ##          (flat) logistic(0, 1)         b      Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                                   thetat (unknown)
    ##          (flat) logistic(0, 1)         b      Intercept                    thetat (unknown)

### Two factors no intercept, yes interaction

``` r
formula <- bmf(kappa ~ 0 + session:cond + cond + session, thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class           coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                          (unknown)
    ##  constant(-100) constant(-100) Intercept                           kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                              mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                              mu9         (unknown)
    ##          (flat)    normal(0,1)         b                                    kappa (unknown)
    ##          (flat)    normal(0,1)         b cond2:session2                     kappa (unknown)
    ##          (flat)    normal(0,1)         b cond3:session2                     kappa (unknown)
    ##          (flat)    normal(0,1)         b cond4:session2                     kappa (unknown)
    ##          (flat)    normal(0,1)         b       session2                     kappa (unknown)
    ##          (flat)    normal(2,1)         b          cond1                     kappa (unknown)
    ##          (flat)    normal(2,1)         b          cond2                     kappa (unknown)
    ##          (flat)    normal(2,1)         b          cond3                     kappa (unknown)
    ##          (flat)    normal(2,1)         b          cond4                     kappa (unknown)
    ##          (flat)         (flat)         b                                  thetant (unknown)
    ##          (flat) logistic(0, 1)         b      Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                                   thetat (unknown)
    ##          (flat) logistic(0, 1)         b      Intercept                    thetat (unknown)

### Two factors no intercept, no main effects, all levels explicit

``` r
formula <- bmf(kappa ~ 0 + session:cond, thetat ~ 1, thetant ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class           coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                          (unknown)
    ##  constant(-100) constant(-100) Intercept                           kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                              mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                              mu9         (unknown)
    ##          (flat)    normal(2,1)         b                                    kappa (unknown)
    ##          (flat)    normal(2,1)         b session1:cond1                     kappa (unknown)
    ##          (flat)    normal(2,1)         b session1:cond2                     kappa (unknown)
    ##          (flat)    normal(2,1)         b session1:cond3                     kappa (unknown)
    ##          (flat)    normal(2,1)         b session1:cond4                     kappa (unknown)
    ##          (flat)    normal(2,1)         b session2:cond1                     kappa (unknown)
    ##          (flat)    normal(2,1)         b session2:cond2                     kappa (unknown)
    ##          (flat)    normal(2,1)         b session2:cond3                     kappa (unknown)
    ##          (flat)    normal(2,1)         b session2:cond4                     kappa (unknown)
    ##          (flat)         (flat)         b                                  thetant (unknown)
    ##          (flat) logistic(0, 1)         b      Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                                   thetat (unknown)
    ##          (flat) logistic(0, 1)         b      Intercept                    thetat (unknown)

### with non-linear predictors

``` r
formula <- bmf(kappa ~ exp(a)*trial, thetat ~ 1, thetant ~ 1, a ~ 1)
compare_priors(formula, dat, model3p)
```

    ##       priorBRMS       priorBMM     class      coef group resp   dpar   nlpar    source
    ##  logistic(0, 1) logistic(0, 1)    theta9                                     (unknown)
    ##  constant(-100) constant(-100) Intercept                      kappa9         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu1         (unknown)
    ##     constant(0)    constant(0) Intercept                         mu9         (unknown)
    ##          (flat)         (flat)         b                                   a (unknown)
    ##          (flat)         (flat)         b Intercept                         a (unknown)
    ##          (flat)         (flat)         b                             thetant (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                   thetant (unknown)
    ##          (flat)         (flat)         b                              thetat (unknown)
    ##          (flat) logistic(0, 1)         b Intercept                    thetat (unknown)

</details>

## Check the same for the sdm

``` r
options(bmm.sort_data = F)
modelSDM <- sdm('dev_rad')
```

<details>
<summary>
Click to expand
</summary>

### Intercept only

All model parameters are `nlpar` so they get class `b` with coef
`Intercept`

``` r
formula <- bmf(kappa ~ 1, c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##         priorBRMS               priorBMM     class coef group resp  dpar nlpar    source
    ##       constant(0)            constant(0) Intercept                             (unknown)
    ##            (flat)    student_t(5,2,0.75) Intercept                     c       (unknown)
    ##  normal(5.0, 0.8) student_t(5,1.75,0.75) Intercept                 kappa       (unknown)

### Intercept and a predictor

For kappa, which we include and intercept and a predictor, we get class
`b` with coef `Intercept` and `session2`

``` r
formula <- bmf(kappa ~ session, c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##         priorBRMS               priorBMM     class     coef group resp  dpar nlpar    source
    ##       constant(0)            constant(0) Intercept                                 (unknown)
    ##            (flat)    student_t(5,2,0.75) Intercept                         c       (unknown)
    ##  normal(5.0, 0.8) student_t(5,1.75,0.75) Intercept                     kappa       (unknown)
    ##            (flat)            normal(0,1)         b                     kappa       (unknown)
    ##            (flat)            normal(0,1)         b session2            kappa       (unknown)

### Intercept supressed

For kappa, which we include and intercept and a predictor, we get class
`b` with coef `session1` and `session2`

``` r
formula <- bmf(kappa ~ 0+session, c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##    priorBRMS               priorBMM     class     coef group resp  dpar nlpar    source
    ##  constant(0)            constant(0) Intercept                                 (unknown)
    ##       (flat)    student_t(5,2,0.75) Intercept                         c       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b                     kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session1            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session2            kappa       (unknown)

### Intercept supressed and random effects

``` r
formula <- bmf(kappa ~ 0+session + ( 0+session|ID), c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##             priorBRMS               priorBMM     class     coef group resp  dpar nlpar    source
    ##           constant(0)            constant(0) Intercept                                 (unknown)
    ##                lkj(1)                 lkj(1)       cor                                 (unknown)
    ##                lkj(1)                 lkj(1)       cor             ID                  (unknown)
    ##                (flat)    student_t(5,2,0.75) Intercept                         c       (unknown)
    ##                (flat) student_t(5,1.75,0.75)         b                     kappa       (unknown)
    ##                (flat) student_t(5,1.75,0.75)         b session1            kappa       (unknown)
    ##                (flat) student_t(5,1.75,0.75)         b session2            kappa       (unknown)
    ##  student_t(3, 0, 2.5)   student_t(3, 0, 2.5)        sd                     kappa       (unknown)
    ##  student_t(3, 0, 2.5)   student_t(3, 0, 2.5)        sd             ID      kappa       (unknown)
    ##  student_t(3, 0, 2.5)   student_t(3, 0, 2.5)        sd session1    ID      kappa       (unknown)
    ##  student_t(3, 0, 2.5)   student_t(3, 0, 2.5)        sd session2    ID      kappa       (unknown)

### Two factors plus intercept

``` r
formula <- bmf(kappa ~ session + cond, c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##         priorBRMS               priorBMM     class     coef group resp  dpar nlpar    source
    ##       constant(0)            constant(0) Intercept                                 (unknown)
    ##            (flat)    student_t(5,2,0.75) Intercept                         c       (unknown)
    ##  normal(5.0, 0.8) student_t(5,1.75,0.75) Intercept                     kappa       (unknown)
    ##            (flat)            normal(0,1)         b                     kappa       (unknown)
    ##            (flat)            normal(0,1)         b    cond2            kappa       (unknown)
    ##            (flat)            normal(0,1)         b    cond3            kappa       (unknown)
    ##            (flat)            normal(0,1)         b    cond4            kappa       (unknown)
    ##            (flat)            normal(0,1)         b session2            kappa       (unknown)

### Two factors no intercept

``` r
formula <- bmf(kappa ~ 0 + session + cond, c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##    priorBRMS               priorBMM     class     coef group resp  dpar nlpar    source
    ##  constant(0)            constant(0) Intercept                                 (unknown)
    ##       (flat)    student_t(5,2,0.75) Intercept                         c       (unknown)
    ##       (flat)            normal(0,1)         b                     kappa       (unknown)
    ##       (flat)            normal(0,1)         b    cond2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b    cond3            kappa       (unknown)
    ##       (flat)            normal(0,1)         b    cond4            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session1            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session2            kappa       (unknown)

### Two factors no intercept, yes interaction

``` r
formula <- bmf(kappa ~ 0 + session * cond, c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##    priorBRMS               priorBMM     class           coef group resp  dpar nlpar    source
    ##  constant(0)            constant(0) Intercept                                       (unknown)
    ##       (flat)    student_t(5,2,0.75) Intercept                               c       (unknown)
    ##       (flat)            normal(0,1)         b                           kappa       (unknown)
    ##       (flat)            normal(0,1)         b          cond2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b          cond3            kappa       (unknown)
    ##       (flat)            normal(0,1)         b          cond4            kappa       (unknown)
    ##       (flat)            normal(0,1)         b session2:cond2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b session2:cond3            kappa       (unknown)
    ##       (flat)            normal(0,1)         b session2:cond4            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b       session1            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b       session2            kappa       (unknown)

### Supressed intercept, 2 factors, only maineffect for facto1 and interaction

``` r
formula <- bmf(kappa ~ 0 + set_size + set_size:session, c~0 + set_size + set_size:session)
compare_priors(formula, dat, modelSDM)
```

    ##    priorBRMS               priorBMM     class               coef group resp  dpar nlpar    source
    ##  constant(0)            constant(0) Intercept                                           (unknown)
    ##       (flat)            normal(0,1)         b                                   c       (unknown)
    ##       (flat)            normal(0,1)         b set_size1:session2                c       (unknown)
    ##       (flat)            normal(0,1)         b set_size2:session2                c       (unknown)
    ##       (flat)            normal(0,1)         b set_size3:session2                c       (unknown)
    ##       (flat)            normal(0,1)         b set_size4:session2                c       (unknown)
    ##       (flat)            normal(0,1)         b set_size5:session2                c       (unknown)
    ##       (flat)            normal(0,1)         b set_size6:session2                c       (unknown)
    ##       (flat)            normal(0,1)         b set_size7:session2                c       (unknown)
    ##       (flat)            normal(0,1)         b set_size8:session2                c       (unknown)
    ##       (flat)    student_t(5,2,0.75)         b          set_size1                c       (unknown)
    ##       (flat)    student_t(5,2,0.75)         b          set_size2                c       (unknown)
    ##       (flat)    student_t(5,2,0.75)         b          set_size3                c       (unknown)
    ##       (flat)    student_t(5,2,0.75)         b          set_size4                c       (unknown)
    ##       (flat)    student_t(5,2,0.75)         b          set_size5                c       (unknown)
    ##       (flat)    student_t(5,2,0.75)         b          set_size6                c       (unknown)
    ##       (flat)    student_t(5,2,0.75)         b          set_size7                c       (unknown)
    ##       (flat)    student_t(5,2,0.75)         b          set_size8                c       (unknown)
    ##       (flat)            normal(0,1)         b                               kappa       (unknown)
    ##       (flat)            normal(0,1)         b set_size1:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b set_size2:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b set_size3:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b set_size4:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b set_size5:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b set_size6:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b set_size7:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b set_size8:session2            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          set_size1            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          set_size2            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          set_size3            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          set_size4            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          set_size5            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          set_size6            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          set_size7            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          set_size8            kappa       (unknown)

### Supressed intercept, 3-way interaction only

``` r
formula <- bmf(kappa ~ 0 + set_size:session:cond, c~1)
compare_priors(formula, dat, modelSDM)
```

    ##    priorBRMS               priorBMM     class                     coef group resp  dpar nlpar
    ##  constant(0)            constant(0) Intercept                                                
    ##       (flat)    student_t(5,2,0.75) Intercept                                         c      
    ##       (flat) student_t(5,1.75,0.75)         b                                     kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size1:session1:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size1:session1:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size1:session1:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size1:session1:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size1:session2:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size1:session2:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size1:session2:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size1:session2:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size2:session1:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size2:session1:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size2:session1:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size2:session1:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size2:session2:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size2:session2:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size2:session2:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size2:session2:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size3:session1:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size3:session1:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size3:session1:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size3:session1:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size3:session2:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size3:session2:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size3:session2:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size3:session2:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size4:session1:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size4:session1:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size4:session1:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size4:session1:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size4:session2:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size4:session2:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size4:session2:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size4:session2:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size5:session1:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size5:session1:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size5:session1:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size5:session1:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size5:session2:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size5:session2:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size5:session2:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size5:session2:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size6:session1:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size6:session1:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size6:session1:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size6:session1:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size6:session2:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size6:session2:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size6:session2:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size6:session2:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size7:session1:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size7:session1:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size7:session1:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size7:session1:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size7:session2:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size7:session2:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size7:session2:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size7:session2:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size8:session1:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size8:session1:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size8:session1:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size8:session1:cond4            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size8:session2:cond1            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size8:session2:cond2            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size8:session2:cond3            kappa      
    ##       (flat) student_t(5,1.75,0.75)         b set_size8:session2:cond4            kappa      
    ##     source
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)
    ##  (unknown)

### Two factors no intercept, yes interaction

``` r
formula <- bmf(kappa ~ 0 + session:cond + cond + session, c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##    priorBRMS               priorBMM     class           coef group resp  dpar nlpar    source
    ##  constant(0)            constant(0) Intercept                                       (unknown)
    ##       (flat)    student_t(5,2,0.75) Intercept                               c       (unknown)
    ##       (flat)            normal(0,1)         b                           kappa       (unknown)
    ##       (flat)            normal(0,1)         b cond2:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b cond3:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b cond4:session2            kappa       (unknown)
    ##       (flat)            normal(0,1)         b       session2            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          cond1            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          cond2            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          cond3            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b          cond4            kappa       (unknown)

### Two factors no intercept, no main effects, all levels explicit

``` r
formula <- bmf(kappa ~ 0 + session:cond, c ~ 1)
compare_priors(formula, dat, modelSDM)
```

    ##    priorBRMS               priorBMM     class           coef group resp  dpar nlpar    source
    ##  constant(0)            constant(0) Intercept                                       (unknown)
    ##       (flat)    student_t(5,2,0.75) Intercept                               c       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b                           kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session1:cond1            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session1:cond2            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session1:cond3            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session1:cond4            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session2:cond1            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session2:cond2            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session2:cond3            kappa       (unknown)
    ##       (flat) student_t(5,1.75,0.75)         b session2:cond4            kappa       (unknown)

</details>
