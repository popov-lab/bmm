# Softmax function and its inverse

`softmax` returns the value of the softmax function `softmaxinv` returns
the value of the inverse-softmax function

## Usage

``` r
softmax(eta, lambda = 1)

softmaxinv(p, lambda = 1, ref_position = length(p), ref_value = 0)
```

## Arguments

- eta:

  A numeric vector input

- lambda:

  Tuning parameter (a single positive value)

- p:

  A probability vector (i.e., numeric vector of non-negative values that
  sum to one)

- ref_position:

  The reference position that should be used to calculate the inverse
  softmax function. The default is the last position.

- ref_value:

  The value the reference position will be set to. The default is 0.

## Value

Value of the softmax function or its inverse

## Details

The softmax function is a bijective function that maps a real vector
with length `m` to a probability vector with length `m` with all
non-zero probabilities. The present functions define the softmax
function and its inverse, both with a tuning parameter.

The current functions define the softmax as:

\$\$\Large P(\eta_i) = \frac{e^{\lambda \eta_i}}{\sum\_{j=1}^m
e^{\lambda \eta_j}}\$\$

## Examples

``` r
softmax(5:7)
#> [1] 0.09003057 0.24472847 0.66524096
softmaxinv(softmax(5:7), ref_position = 1, ref_value = 5)
#> [1] 5 6 7
```
