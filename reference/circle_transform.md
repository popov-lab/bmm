# Convert degrees to radians or radians to degrees.

The helper functions `deg2rad` and `rad2deg` should add convenience in
transforming data from degrees to radians and from radians to degrees.

## Usage

``` r
deg2rad(deg)

rad2deg(rad)
```

## Arguments

- deg:

  A numeric vector of values in degrees.

- rad:

  A numeric vector of values in radians.

## Value

A numeric vector of the same length as `deg` or `rad`.

## Examples

``` r
degrees <- runif(100, min = 0, max = 360)
radians <- deg2rad(degrees)
degrees_again <- rad2deg(radians)
```
