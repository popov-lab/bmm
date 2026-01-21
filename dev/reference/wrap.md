# Wrap angles that extend beyond (-pi;pi)

On the circular space, angles can be only in the range (-pi;pi or
-180;180). When subtracting angles, this can result in values outside of
this range. For example, when calculating the difference between a value
of 10 degrees minus 340 degrees, this results in a difference of 330
degrees. However, the true difference between these two values is -30
degrees. This function wraps such values, so that they occur in the
circle

## Usage

``` r
wrap(x, radians = TRUE)
```

## Arguments

- x:

  A numeric vector, matrix or data.frame of angles to be wrapped. In
  radians (default) or degrees.

- radians:

  Logical. Is x in radians (default=TRUE) or degrees (FALSE)

## Value

An object of the same type as x

## Examples

``` r
x <- runif(1000, -pi, pi)
y <- runif(1000, -pi, pi)
diff <- x - y
hist(diff)

wrapped_diff <- wrap(x - y)
hist(wrapped_diff)

```
