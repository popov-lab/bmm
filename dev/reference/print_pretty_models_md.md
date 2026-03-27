# Generate a markdown list of the measurement models available in `bmm`

Used internally to automatically populate information in the README file

## Usage

``` r
print_pretty_models_md()
```

## Value

Markdown code for printing the list of measurement models available in
`bmm`

## Examples

``` r
print_pretty_models_md()
#> **Processing Speed, Decision Making**
#> 
#> * Censored-Shifted Wald Model 
#> * EZ-Diffusion Model 
#> 
#> **Visual working memory**
#> 
#> * Interference measurement model by Oberauer and Lin (2017). 
#> * Two-parameter mixture model by Zhang and Luck (2008). 
#> * Three-parameter mixture model by Bays et al (2009). 
#> * Signal Discrimination Model (SDM) by Oberauer (2023) 
#> 
#> **Working Memory (categorical), Categorical Decision Making**
#> 
#> * The Multinomial / Memory Measurement Model 
#> 
```
