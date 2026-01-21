# Softmax function and its inverse

`softmax` returns the value of the softmax function `softmaxinv` returns
the value of the inverse-softmax function

## Usage

``` r
softmax(eta, lambda = 1)

softmaxinv(p, lambda = 1)
```

## Arguments

- eta:

  A numeric vector input

- lambda:

  Tuning parameter (a single positive value)

- p:

  A probability vector (i.e., numeric vector of non-negative values that
  sum to one)

## Value

Value of the softmax function or its inverse

## Details

The softmax function is a bijective function that maps a real vector
with length `m-1` to a probability vector with length `m` with all
non-zero probabilities. The present functions define the softmax
function and its inverse, both with a tuning parameter.

The current functions define the softmax as:

\$\$\Large P(\eta_i) = \frac{e^{\lambda \eta_i}}{1+ \sum\_{j=1}^m
e^{\lambda \eta_j}}\$\$

Code adapted from the
[utilities](https://github.com/ben-oneill/utilities/) package

## Examples

``` r
softmax(5:7)
#> [1] 0.0899759918 0.2445801036 0.6648376511 0.0006062535
softmaxinv(softmax(5:7))
#> [1] 5 6 7
```
