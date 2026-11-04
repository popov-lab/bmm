# Distribution function for the Diffusion Decision Model (`ddm`)

Density and random generation function for the Diffusion Decision Model.

## Usage

``` r
dddm(rt, response, drift, bound, ndt, zr = 0.5, log = TRUE)

rddm(n, drift, bound, ndt, zr = 0.5)
```

## Arguments

- rt:

  Vector of response times for which the density should be returned

- response:

  Vector of responses for which the density should be returned

- drift:

  Drift rates of the ddm

- bound:

  Boundary separation of the ddm

- ndt:

  Non-decision time of the ddm

- zr:

  relative starting point of the ddm

- log:

  Logical; if `TRUE` (default until bmm 1.6.0), log-densities are
  returned. From bmm 1.6.0 on the default is `FALSE`, as in base R and
  the other density functions of bmm. Pass `log` to get the same result
  in every version.

- n:

  Number of random samples to generate
