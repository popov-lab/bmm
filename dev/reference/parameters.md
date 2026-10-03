# Deprecated: use `parameter_info()`

`parameters()` is deprecated as of bmm 1.4.0 and will be removed in bmm
1.6.0. It shares its name with `parameters::parameters()`, so whichever
of the two packages is attached last decides what `parameters()`
returns. Use
[`parameter_info()`](https://popov-lab.github.io/bmm/dev/reference/parameter_info.md)
instead; it takes the same arguments and returns the same table.

## Usage

``` r
parameters(x, ...)
```

## Arguments

- x:

  A `bmmodel` object (e.g., `sdm(resp_error = "y")`) or a `bmmfit`
  object (a fitted model returned by
  [`bmm`](https://popov-lab.github.io/bmm/dev/reference/bmm.md)).

- ...:

  Passed on to
  [`parameter_info()`](https://popov-lab.github.io/bmm/dev/reference/parameter_info.md),
  e.g. `formula` for a custom M3 model.

## Value

The output of
[`parameter_info()`](https://popov-lab.github.io/bmm/dev/reference/parameter_info.md).
