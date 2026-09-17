# Human-readable pre-fit data report for bmm models

Inspect whether your data is coded the way a `bmmodel` expects *before*
committing to Stan compilation and sampling. The function runs the same
validation pipeline that
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) runs internally
([`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md)
and
[`check_formula()`](https://venpopov.com/bmm/dev/reference/check_formula.md)),
but captures all errors, warnings and messages instead of throwing them,
and combines them with a readable summary of:

- the response variables, their observed values and the coding the model
  expects

- the data columns each parameter formula uses, and how they are coded
  (factor levels, character vectors, numeric variables with few unique
  values)

- the number of observations per design cell (crossing the grouping
  variables from random effects terms with the categorical predictors),
  flagging sparse or empty cells

- model-specific diagnostics for common data mistakes (e.g. responses in
  degrees rather than radians for circular models, or misplaced `NA`
  values in the non-target features of set size varying designs). These
  are provided by
  [`data_check_findings()`](https://venpopov.com/bmm/dev/reference/data_check_findings.md)
  methods, which new models can extend

Problems are reported as findings rather than errors, so a single call
shows everything that needs fixing at once.

The report covers the data and the formula only. It stops after
[`check_formula()`](https://venpopov.com/bmm/dev/reference/check_formula.md)
and takes no `prior` argument, so a malformed custom prior or a problem
with the initial values is out of its reach and will only surface when
you call [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md).

## Usage

``` r
bmm_data_check(formula, data, model, min_trials = 10)
```

## Arguments

- formula:

  An object of class `bmmformula`. A symbolic description of the model
  to be fitted.

- data:

  An object of class data.frame, containing data of all variables used
  in the model. The names of the variables must match the variable names
  passed to the `bmmodel` object for required argurments.

- model:

  A description of the model to be fitted. This is a call to a `bmmodel`
  such as
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md)
  function. Every model function has a number of required arguments
  which need to be specified within the function call. Call
  [`supported_models()`](https://venpopov.com/bmm/dev/reference/supported_models.md)
  to see the list of supported models and their required arguments

- min_trials:

  Numeric. Design cells with fewer observations are flagged in the
  report. Defaults to 10. Ignored for models fit to aggregated data
  (e.g. `m3`, `ezdm`), where one row summarizes many trials.

## Value

An object of class `bmm_data_check` with a print method. The object is a
list containing the summarized `response`, `predictors` and `cells`
information (including the full cell count table in `$cells$counts`),
all `findings`, and the captured results of the validation `pipeline`.

## See also

[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md),
[`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md),
[`data_check_findings()`](https://venpopov.com/bmm/dev/reference/data_check_findings.md)

## Examples

``` r
# three-parameter mixture model: radians check, NA padding of the
# non-target features across set sizes, trials per participant x set size
bmm_data_check(
  bmf(kappa ~ 0 + set_size + (0 + set_size | ID), thetat ~ 0 + set_size),
  data = oberauer_lin_2017,
  model = mixture3p(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size"
  )
)
#> Pre-fit data check
#> Model: mixture3p(resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size") 
#> Data:  15200 observations ('oberauer_lin_2017') 
#> 
#> Response variables:
#>   dev_rad  numeric, range [-3.14, 3.12] (expected: radians in [-pi, pi])
#> 
#> Data columns used by each parameter formula:
#>   kappa   set_size, ID
#>   thetat  set_size
#> 
#> Predictor coding:
#>   set_size  [predictor] factor with 8 levels: '1', '2', '3', '4', '5', '6',
#>                         '7', '8'
#>   ID        [grouping]  numeric, range [1, 19]
#> 
#> Observations per design cell:
#>   ID x set_size: 152 non-empty cell(s)
#>   observations per cell: min 100 / median 100 / max 100
#> 
#> Findings: No issues detected.
#> 
#> Hard checks: passed
#>   - No formula for parameter thetant provided. Only a fixed intercept will be
#>     estimated.

# interference measurement model: additionally checks the NA padding
# of the non-target distances
bmm_data_check(
  bmf(c ~ 1, a ~ 1, kappa ~ 1, s ~ 1),
  data = oberauer_lin_2017,
  model = imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "full"
  )
)
#> Pre-fit data check
#> Model: imm(resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7), nt_distances = paste0("dist_nt", 1:7), set_size = "set_size", version = "full") 
#> Data:  15200 observations ('oberauer_lin_2017') 
#> 
#> Response variables:
#>   dev_rad  numeric, range [-3.14, 3.12] (expected: radians in [-pi, pi])
#> 
#> Findings: No issues detected.
#> 
#> Hard checks: passed

# signal discrimination model
bmm_data_check(
  bmf(c ~ 0 + setsize, kappa ~ 0 + setsize),
  data = zhang_luck_2008,
  model = sdm(resp_error = "response_error")
)
#> Pre-fit data check
#> Model: sdm(resp_error = "response_error") 
#> Data:  4000 observations ('zhang_luck_2008') 
#> 
#> Response variables:
#>   response_error  numeric, range [-3.14, 3.14] (expected: radians in [-pi, pi])
#> 
#> Data columns used by each parameter formula:
#>   c      setsize
#>   kappa  setsize
#> 
#> Predictor coding:
#>   setsize  [predictor] factor with 4 levels: '1', '2', '3', '6'
#> 
#> Observations per design cell:
#>   setsize: 4 non-empty cell(s)
#>   observations per cell: min 1000 / median 1000 / max 1000
#> 
#> Findings: No issues detected.
#> 
#> Hard checks: passed

# memory measurement model: response category counts as response variables
bmm_data_check(
  bmf(c ~ 1 + (1 | ID), a ~ 1 + (1 | ID)),
  data = oberauer_lewandowsky_2019_e1,
  model = m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "simple",
    version = "ss"
  )
)
#> Pre-fit data check
#> Model: m3(resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"), choice_rule = "simple", version = "ss") 
#> Data:  120 observations ('oberauer_lewandowsky_2019_e1') 
#> 
#> Response variables:
#>   corr   numeric, range [5, 95] (expected: counts per response category)
#>   other  numeric, range [4, 66] (expected: counts per response category)
#>   npl    numeric, range [0, 70] (expected: counts per response category)
#> 
#> Data columns used by each parameter formula:
#>   c  ID
#>   a  ID
#> 
#> Predictor coding:
#>   ID  [grouping] numeric, range [1, 40]
#> 
#> Observations per design cell:
#>   ID: 40 non-empty cell(s)
#>   observations per cell: min 3 / median 3 / max 3
#> 
#> Findings: No issues detected.
#> 
#> Hard checks: passed

# diffusion decision model: reaction time plausibility and response coding
bmm_data_check(
  bmf(drift ~ condition + (condition | ID)),
  data = data_color_judgement_task,
  model = ddm(rt = "rt", response = "response_correct")
)
#> Pre-fit data check
#> Model: ddm(rt = "rt", response = "response_correct") 
#> Data:  9941 observations ('data_color_judgement_task') 
#> 
#> Response variables:
#>   rt                numeric, range [0.124, 4.99] (expected: seconds)
#>   response_correct  logical, 8601 TRUE / 1340 FALSE (expected: 0/1 or logical;
#>                     1 = upper boundary)
#> 
#> Data columns used by each parameter formula:
#>   drift  condition, ID
#> 
#> Predictor coding:
#>   condition  [predictor] character with 2 levels: 'easy', 'hard'
#>   ID         [grouping]  character with 50 levels: 'AAUE02', 'AAUI51',
#>                          'AEHE22', 'AGOE08', 'AIEO12', 'AMAO80', 'ANAA57',
#>                          'ANAO01', 'ANRR94', 'AOAE64', ...
#> 
#> Observations per design cell:
#>   ID x condition: 100 non-empty cell(s)
#>   observations per cell: min 84 / median 100 / max 100
#> 
#> Findings:
#>   - Character predictor(s) 'condition' will be coerced to factors with
#>     alphabetical level order. Use factor() to control the reference level.
#> 
#> Hard checks: passed with 1 warning(s):
#>   ! The response variable is boolean and will be internally transformed to an
#>     integer variable with values 0 for FALSE and 1 for TRUE.
#>   - No formula for parameter bound provided. Only a fixed intercept will be
#>     estimated.
#>   - No formula for parameter ndt provided. Only a fixed intercept will be
#>     estimated.
```
