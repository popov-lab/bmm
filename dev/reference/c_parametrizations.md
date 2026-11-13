# Convert between parametrizations of the c parameter of the SDM distribution

Convert between parametrizations of the c parameter of the SDM
distribution

## Usage

``` r
c_sqrtexp2bessel(c, kappa)

c_bessel2sqrtexp(c, kappa)

c_laplace_circle2line(c, kappa)

c_laplace_line2circle(c, kappa)
```

## Arguments

- c:

  Vector of memory strength values

- kappa:

  Vector of precision values

## Value

A numeric vector of the same length as `c` and `kappa`.

## Details

With `scaling = "density"`,
[`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md) divides the
memory strength `c` by the integral of the activation function. Two
pairs of functions convert `c` to and from the forms used elsewhere.

**Von Mises activation.** `bmm` approximates the integral by
\\\sqrt{2\pi/\kappa}\\ (the sqrtexp parametrization, the default in
`bmm`); Oberauer (2023) uses the exact von Mises density (the bessel
parametrization). `c_sqrtexp2bessel` converts from the sqrtexp to the
bessel parametrization, `c_bessel2sqrtexp` back. See [the online
article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html).

**Laplace activation.** `bmm` divides by the exact integral over the
circle, \\2(1 - e^{-\kappa \pi}) / \kappa\\; Oberauer (2023, eq. 2)
multiplies the activation by \\\kappa / 2\\, that is, divides by the
integral over the real line, \\2 / \kappa\\. `c_laplace_circle2line`
converts the `c` of `bmm` to that of eq. 2, \\c\_{line} = c\_{circle} /
(1 - e^{-\kappa \pi})\\, and `c_laplace_line2circle` back. The two
differ by less than 1% once \\\kappa\\ exceeds 1.5. At \\\kappa = 0\\
the real-line form has no finite `c`; `c_laplace_circle2line()` returns
`Inf` there, and `NaN` when `c` is 0 as well.

## References

Oberauer, K. (2023). Measurement models for visual working memory - A
factorial model comparison. Psychological Review, 130(3), 841–852

## Examples

``` r
c_bessel <- c_sqrtexp2bessel(c = 4, kappa = 3)
c_sqrtexp <- c_bessel2sqrtexp(c = c_bessel, kappa = 3)

# c of a Laplace fit with scaling = "density", on the scale of eq. 2
c_line <- c_laplace_circle2line(c = 4, kappa = 1)
c_laplace_line2circle(c = c_line, kappa = 1)
#> [1] 4
```
