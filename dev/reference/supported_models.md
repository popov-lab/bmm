# Measurement models available in `bmm`

Measurement models available in `bmm`

## Usage

``` r
supported_models(print_call = TRUE)
```

## Arguments

- print_call:

  Logical; If TRUE (default), the function will print information about
  how each model function should be called and its required arguments.
  If FALSE, the function will return a character vector with the names
  of the available models

## Value

A character vector of measurement models available in `bmm`

## Examples

``` r
supported_models()
#> The following models are supported:
#> 
#> -  ezdm(mean_rt, var_rt, n_upper, n_trials, links, version) 
#> -  imm(resp_error, nt_features, nt_distances, set_size, regex, version) 
#> -  m3(resp_cats, num_options, choice_rule, version) 
#> -  mixture2p(resp_error) 
#> -  mixture3p(resp_error, nt_features, set_size, regex) 
#> -  sdm(resp_error, version) 
#> 
#> Type  ?modelname  to get information about a specific model, e.g.  ?imm 
```
