# Calculate response error relative to non-target values

Given a vector of responses, and the values of non-targets, this
function computes the error relative to each of the non-targets.

## Usage

``` r
calc_error_relative_to_nontargets(data, response, nt_features)
```

## Arguments

- data:

  A `data.frame` object where each row is a single observation

- response:

  Character. The name of the column in `data` which contains the
  response

- nt_features:

  Character vector. The names of the columns in `data` which contain the
  values of the non-targets

## Value

A `data.frame` with n\*m rows, where n is the number of rows of `data`
and m is the number of non-target variables. It preserves all other
columns of `data`, except for the non-target locations, and adds a
column `y_nt`, which contains the transformed response error relative to
the non-targets

## Examples

``` r
data <- oberauer_lin_2017
data <- calc_error_relative_to_nontargets(data, "dev_rad", paste0("col_nt", 1:7))
hist(data$y_nt, breaks = 100)

```
