# Signal Discrimination Model (SDM)

Signal Discrimination Model (SDM)

## Usage

``` r
sdm(
  resp_error,
  version = "simple",
  activation = c("vonmises", "laplace"),
  scaling = NULL,
  ...
)

sdmSimple(resp_error, version = "simple", ...)
```

## Arguments

- resp_error:

  The name of the variable in the dataset containing the response error.
  The response error should code the response relative to the
  to-be-recalled target in radians. You can transform the response error
  in degrees to radians using the `deg2rad` function.

- version:

  Character. The version of the model to use. Currently only "simple" is
  supported.

- activation:

  Character. The activation kernel of the SDM, `"vonmises"` (the
  default, as in Oberauer, 2023) or `"laplace"`. See the section
  "Activation functions".

- scaling:

  Character. What the memory strength `c` measures: `"density"`,
  `"peak"` or `"range"`. `NULL` (the default) means `"density"` for the
  von Mises and `"peak"` for the Laplace activation. See the sections
  "Scaling of the activation" and "What `c` means".

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

see [the online
article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for a
detailed description of the model and how to use it. \* **Domain:**
Working memory

- **Task:** Continuous reproduction

- **Name:** Signal Discrimination Model (SDM)

- **Citation:**

  - Oberauer, K. (2023). Measurement models for visual working memory—A
    factorial model comparison. Psychological Review, 130(3), 841-852.
    https://doi.org/10.1037/rev0000328

- **Version:** simple

- **Requirements:**

  - The response variable should be in radians and represent the angular
    error relative to the target

- **Parameters:**

  - `mu`: Location parameter of the SDM distribution (in radians; by
    default fixed internally to 0)

  - `c`: Memory strength parameter of the SDM distribution

  - `kappa`: Precision parameter of the SDM distribution

- **Fixed parameters:**

  - `mu` = 0

- **Default parameter links:**

  - mu = tan_half; c = log; kappa = log

- **Default priors:**

  - `mu`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `kappa`:

    - `main`: student_t(5, 1.75, 0.75)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `c`:

    - `main`: student_t(5, 2, 0.75)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

## Activation functions

Each response \\x\\ is activated by \\c^\* s(x - \mu)\\, where the
activation function \\s\\ is 1 at the target and \\c^\*\\ is the
activation of the target; the response density is \\\exp(c^\* s(x -
\mu)) / Z\\, with \\Z\\ its integral over the circle.

- `activation = "vonmises"`: \\s(x) = \exp(\kappa (\cos x - 1))\\, the
  von Mises activation kernel of the signal discrimination model of
  Oberauer (2023).

- `activation = "laplace"`: \\s(x) = \exp(-\kappa \|x\|)\\, with
  \\\|x\|\\ the angular distance and \\\kappa\\ the rate of the Laplace
  activation in Oberauer (2023, eq. 2). This is the SDM with the Laplace
  activation kernel of Oberauer's model space (Table 2): the same Gumbel
  noise and decision rule, with \\\kappa\\ estimated freely. The density
  has a kink at \\\mu\\, so estimating `mu` puts a kink at every
  response, which gradient-based sampling handles poorly;
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) warns when
  the formula estimates `mu`.

## Scaling of the activation

`scaling` decides how `c` translates into the activation of the target
\\c^\*\\; the table in "What `c` means" gives the formulas. `"density"`
divides `c` by the integral of \\s\\ over the circle, so `c` weighs an
activation of unit area. For the von Mises activation that integral is
approximated by \\\sqrt{2\pi/\kappa}\\, the form `bmm` has always used
(the "sqrtexp" parametrization;
[`c_sqrtexp2bessel()`](https://venpopov.com/bmm/dev/reference/c_parametrizations.md)
converts to Oberauer's exact Bessel form), and fits from before
`scaling` existed are this cell. For the Laplace activation the integral
is exact, \\2(1 - e^{-\kappa \pi}) / \kappa\\; Oberauer (2023, eq. 2)
multiplies by \\\kappa / 2\\, which divides by the integral over the
real line, \\2 / \kappa\\. `bmm` keeps the exact integral over the
circle, so its `c` is that of eq. 2 times \\1 - e^{-\kappa \pi}\\ (0.998
at \\\kappa = 2\\, 0.79 at \\\kappa = 0.5\\);
[`c_laplace_circle2line()`](https://venpopov.com/bmm/dev/reference/c_parametrizations.md)
converts it to the `c` of eq. 2 and
[`c_laplace_line2circle()`](https://venpopov.com/bmm/dev/reference/c_parametrizations.md)
back. `"peak"` takes `c` as the activation of the target, and `"range"`
as the difference between the activation of the target and that of the
location opposite it.

## Default priors

`kappa` and `c` mean different things in each cell, so each cell has its
own `main` priors (log scale; the `effects` and `sd` priors are the same
everywhere). They were chosen in a prior predictive check so that the
implied distributions of the circular SD of the response error are
comparable across cells, not exactly matched. The last column is the
share of prior predictive distributions with a circular SD below 0.05
radians; the Laplace cells still produce more of these spikes. The von
Mises peak and range cells keep the priors of the default cell, which
already give them comparable error distributions (Kolmogorov-Smirnov
distance 0.06 to the default cell); the Laplace cells needed their own.
Set your own priors where that matters.

|  |  |  |  |  |
|----|----|----|----|----|
| `activation` | `scaling` | `c` | `kappa` | circular SD \< 0.05 |
| `"vonmises"` | `"density"` | `student_t(5, 2, 0.75)` | `student_t(5, 1.75, 0.75)` | 0.07 |
| `"vonmises"` | `"peak"` | `student_t(5, 2, 0.75)` | `student_t(5, 1.75, 0.75)` | 0.03 |
| `"vonmises"` | `"range"` | `student_t(5, 2, 0.75)` | `student_t(5, 1.75, 0.75)` | 0.03 |
| `"laplace"` | `"density"` | `student_t(5, 3.25, 0.75)` | `student_t(5, -0.5, 0.75)` | 0.16 |
| `"laplace"` | `"peak"` | `student_t(5, 2.5, 0.75)` | `student_t(5, -0.5, 0.75)` | 0.11 |
| `"laplace"` | `"range"` | `student_t(5, 2.25, 0.75)` | `student_t(5, -0.5, 0.75)` | 0.09 |

## What `c` means

|  |  |  |  |
|----|----|----|----|
| `scaling` | \\c^\*\\, von Mises | \\c^\*\\, Laplace | `c` is |
| `"density"` | \\c \sqrt{\kappa / 2\pi}\\ | \\c \kappa / (2 (1 - e^{-\kappa \pi}))\\ | the weight of an activation of unit area |
| `"peak"` | \\c\\ | \\c\\ | the activation of the target |
| `"range"` | \\c / (1 - e^{-2\kappa})\\ | \\c / (1 - e^{-\kappa \pi})\\ | the target-opposite contrast in activation |

For a given data set the three scalings give the same likelihood at
values of `c` related by these formulas; what differs is which quantity
the regression formula for `c` and its prior act on.

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
# simulate data from the model
dat <- data.frame(y = rsdm(n = 1000, c = 4, kappa = 3))

# specify formula
ff <- bmf(
  c ~ 1,
  kappa ~ 1
)

# specify the model
fit <- bmm(
  formula = ff,
  data = dat,
  model = sdm(resp_error = "y"),
  cores = 4,
  backend = "cmdstanr"
)
}
```
