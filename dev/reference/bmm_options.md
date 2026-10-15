# View or change global bmm options

View or change global bmm options

## Usage

``` r
bmm_options(
  sort_data,
  parallel,
  default_priors,
  silent,
  color_summary,
  file_refit,
  step_size,
  reset_options = FALSE
)
```

## Arguments

- sort_data:

  logical. If TRUE, the data will be sorted by the predictors. If FALSE,
  the data will not be sorted, but sampling will be slower. If "check"
  (the default),
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) will check if
  the data is sorted, and ask you via a console prompt if it should be
  sorted. **Default: "check"**

- parallel:

  logical. If TRUE, chains will be run in parallel. If FALSE, chains
  will be run sequentially. You can also set these value for each model
  separately via the argument `parallel` in
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md). **Default:
  FALSE**

- default_priors:

  logical. If TRUE (default), the default bmm priors will be used. If
  FALSE, only the basic `brms` priors will be used. **Default: TRUE**

- silent:

  numeric. Verbosity level between 0 and 2. If 1 ( the default), most of
  the informational messages of compiler and sampler are suppressed. If
  2, even more messages are suppressed. The actual sampling progress is
  still printed. **Default: 1**

- color_summary:

  logical. If TRUE, the summary of the model will be printed in color.
  **Default: TRUE**

- file_refit:

  logical or character. Controls when
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) re-uses a fit
  saved via its `file` argument. `TRUE` or "always" always refits,
  `FALSE` or "never" always re-uses the saved fit, and "on_change"
  re-uses it only while the Stan code and data are unchanged. See
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) for details.
  **Default: FALSE**

- step_size:

  numeric or `FALSE`. The step size at which
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) and
  [`update()`](https://rdrr.io/r/stats/update.html) start Stan's
  step-size search, passed as `control = list(step_size = )` (`stepsize`
  for the rstan backend). Stan's own starting value of 1 is far above
  the step sizes the hierarchical models in bmm adapt to, and the
  oversized trial steps of the search print `lkj_corr_cholesky_lpdf` and
  `von_mises_lpdf` exceptions at the start of warmup. The adapted step
  size and the posterior do not depend on the starting value. `FALSE`
  leaves the starting step size to Stan; a `step_size` (or `stepsize`)
  in the `control` list of
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) always wins.
  **Default: 0.01**

- reset_options:

  logical. If TRUE, the options will be reset to their default values
  **Default: FALSE**

## Value

A message with the current bmm options and their values, and invisibly
returns the old options for use with on.exit() and friends.

## Details

The `bmm_options` function is used to view or change the current bmm
options. If no arguments are provided, the function will return the
current options. If arguments are provided, the function will change the
options and return the old options invisibly. If you provide only some
of the arguments, the other options will not be changed. The options are
stored in the global options list and will be used by
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) and other
functions in the `bmm` package. Each of these options can also be set
manually using the built-in
[`options()`](https://rdrr.io/r/base/options.html) function, by setting
the `bmm.sort_data`, `bmm.default_priors`, and `bmm.silent` options.

## Examples

``` r

# view the current options
bmm_options()
#> Current bmm options:
#>   sort_data = check
#>   parallel = FALSE
#>   default_priors = TRUE
#>   silent = 1
#>   file_refit = FALSE
#>   step_size = 0.01
#>   color_summary = TRUE
#> For more information on these options or how to change them, see help(bmm_options).

# change the options to always sort the data and to use parallel sampling
bmm_options(sort_data = TRUE, parallel = TRUE)
#> Current bmm options:
#>   sort_data = TRUE
#>   parallel = TRUE
#>   default_priors = TRUE
#>   silent = 1
#>   file_refit = FALSE
#>   step_size = 0.01
#>   color_summary = TRUE
#> For more information on these options or how to change them, see help(bmm_options).

# restore the default options
bmm_options(reset_options = TRUE)
#> Current bmm options:
#>   sort_data = check
#>   parallel = FALSE
#>   default_priors = TRUE
#>   silent = 1
#>   file_refit = FALSE
#>   step_size = 0.01
#>   color_summary = TRUE
#> For more information on these options or how to change them, see help(bmm_options).

# you can change the options using the options() function as well
options(bmm.sort_data = TRUE, bmm.parallel = TRUE)
bmm_options()
#> Current bmm options:
#>   sort_data = TRUE
#>   parallel = TRUE
#>   default_priors = TRUE
#>   silent = 1
#>   file_refit = FALSE
#>   step_size = 0.01
#>   color_summary = TRUE
#> For more information on these options or how to change them, see help(bmm_options).

# reset the options to their default values
bmm_options(reset_options = TRUE)
#> Current bmm options:
#>   sort_data = check
#>   parallel = FALSE
#>   default_priors = TRUE
#>   silent = 1
#>   file_refit = FALSE
#>   step_size = 0.01
#>   color_summary = TRUE
#> For more information on these options or how to change them, see help(bmm_options).

# bmm_options(sort_data = TRUE, parallel = TRUE) will also return the old options
# so you can use it with on.exit()
old_op <- bmm_options(sort_data = TRUE, parallel = TRUE)
#> Current bmm options:
#>   sort_data = TRUE
#>   parallel = TRUE
#>   default_priors = TRUE
#>   silent = 1
#>   file_refit = FALSE
#>   step_size = 0.01
#>   color_summary = TRUE
#> For more information on these options or how to change them, see help(bmm_options).
on.exit(bmm_options(old_op))
#> Error in !missing(sort_data) && sort_data != "check": 'length = 2' in coercion to 'logical(1)'

bmm_options(reset_options = TRUE)
#> Current bmm options:
#>   sort_data = check
#>   parallel = FALSE
#>   default_priors = TRUE
#>   silent = 1
#>   file_refit = FALSE
#>   step_size = 0.01
#>   color_summary = TRUE
#> For more information on these options or how to change them, see help(bmm_options).
```
