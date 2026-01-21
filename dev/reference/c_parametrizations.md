# Convert between parametrizations of the c parameter of the SDM distribution

Convert between parametrizations of the c parameter of the SDM
distribution

## Usage

``` r
c_sqrtexp2bessel(c, kappa)

c_bessel2sqrtexp(c, kappa)
```

## Arguments

- c:

  Vector of memory strength values

- kappa:

  Vector of precision values

## Value

A numeric vector of the same length as `c` and `kappa`.

## Details

`c_bessel2sqrtexp` converts the memory strength parameter (c) from the
bessel parametrization to the sqrtexp parametrization,
`c_sqrtexp2bessel` converts from the sqrtexp parametrization to the
bessel parametrization.

See [the online
article](https://venpopov.github.io/bmm/articles/bmm_sdm_simple.html)
for details on the parameterization. The sqrtexp parametrization is the
default in the `bmm` package.

## Examples

``` r
c_bessel <- c_sqrtexp2bessel(c = 4, kappa = 3)
c_sqrtexp <- c_bessel2sqrtexp(c = c_bessel, kappa = 3)
```
