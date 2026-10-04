# Deprecated: use `bmm_models()`

`supported_models()` is deprecated as of bmm 1.4.0 and will be removed
in bmm 1.6.0. It shares its name with `insight::supported_models()`,
which the **parameters** package re-exports, so whichever package is
attached last decides what `supported_models()` returns. Replace both
`supported_models()` and `supported_models(print_call = FALSE)` with
[`bmm_models()`](https://popov-lab.github.io/bmm/dev/reference/bmm_models.md);
wrap it in [`as.character()`](https://rdrr.io/r/base/character.html) if
you need a plain character vector.

## Usage

``` r
supported_models(print_call = TRUE)
```

## Arguments

- print_call:

  Logical. If `TRUE` (default), returns the output of
  [`bmm_models()`](https://popov-lab.github.io/bmm/dev/reference/bmm_models.md),
  which prints the models grouped by task. If `FALSE`, returns the model
  names as a plain character vector.

## Value

The output of
[`bmm_models()`](https://popov-lab.github.io/bmm/dev/reference/bmm_models.md),
or the model names as a plain character vector if `print_call = FALSE`.
Before bmm 1.4.0 the default returned the printed list as one string;
use `capture.output(print(bmm_models()))` for the printed list as text,
one line per element.
