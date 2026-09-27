# Check whether this machine can fit bmm models

Fitting a model with
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) needs a
C++ toolchain (a compiler and `make`) and a Stan backend, the software
that turns the model into a program and samples from it. There are two
backends: the R package `cmdstanr`, which runs CmdStan, a separate
program installed with
[`cmdstanr::install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html);
and the R package `rstan`. `rstan` must load with either backend,
because `brms` uses it to store every fit. `bmm_setup()` checks each of
these, prints whether it passed, and gives one fix for every check that
failed. It does not install or change anything.

Two checks concern the toolchain. "C++ toolchain" lets R compile a small
test file; "CmdStan toolchain" is `cmdstanr`'s own check, which only
looks for `make` and a compiler on the `PATH`. The second can pass while
the first fails, e.g. on a Mac whose command line tools are missing, and
the first is the one that decides.

With `smoke_test = TRUE`, it finally compiles and samples a small
two-parameter mixture model
([`mixture2p()`](https://popov-lab.github.io/bmm/dev/reference/mixture2p.md))
through [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md),
on the backend
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) will
use. This is the only check that shows the whole chain works. On an
Apple Silicon Mac it takes about 10 seconds with `cmdstanr` and 45
seconds with `rstan`, and it is skipped when a check it needs has
failed. Once every check passes, rerun the
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) call
that failed.

## Usage

``` r
bmm_setup(smoke_test = TRUE, backend = getOption("brms.backend", NULL))
```

## Arguments

- smoke_test:

  Logical. Compile and sample a small test model? Defaults to `TRUE`.

- backend:

  The backend to check, `"cmdstanr"` or `"rstan"`. The default checks
  the backend
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) would
  choose: the `brms.backend` option if it is set, otherwise `cmdstanr`
  whenever the `cmdstanr` package is installed, even when CmdStan is
  not, and `rstan` otherwise.

## Value

A data frame of class `bmm_setup` with one row per check and the columns
`check`, `status` (`"pass"`, `"fail"` or `"skip"`), `detail` and `fix`
(`NA` unless the check failed). It prints as a report.

## See also

[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md),
[`bmm_data_check()`](https://popov-lab.github.io/bmm/dev/reference/bmm_data_check.md)

## Examples

``` r
if (FALSE) { # interactive()
# checks without compiling anything
bmm_setup(smoke_test = FALSE)

# compiles and samples a small model as well
bmm_setup()
}
```
