# Transform kappa of the von Mises distribution to the circular standard deviation

This function transforms the precision parameter kappa of the von Mises
distribution to the circular standard deviation. Adapted from Matlab
code by Paul Bays (https://www.paulbays.com/code.php)

## Usage

``` r
k2sd(K)
```

## Arguments

- K:

  numeric. A vector of kappa values.

## Value

A vector of sd values.

## Examples

``` r
kappas <- runif(1000, 0.01, 100)

# calcualte SD (in radians)
SDs <- k2sd(kappas)

# transform SDs from radians to degrees
SDs_degress <- SDs * 180 / pi

# plot the relationship between kappa and circular SD
plot(kappas, SDs)

plot(kappas, SDs_degress)
```
