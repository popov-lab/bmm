# Measurement models available in `bmm`

Measurement models available in `bmm`

## Usage

``` r
bmm_models()
```

## Value

A character vector of model names with class `bmm_models`, which prints
as the grouped list. Use it as a character vector in base R or dplyr,
e.g. `"imm" %in% bmm_models()`;
[`as.character()`](https://rdrr.io/r/base/character.html) gives a plain
one where a function refuses the class.

## Details

Printed, the result lists the models grouped by the task they are meant
for, one line per model with its constructor and name. The groups are:
continuous reproduction; categorical recall and n-AFC decisions;
detection, recognition and confidence judgments; choices and response
times. A model that fits several tasks, such as
[`mpt()`](https://venpopov.com/bmm/dev/reference/mpt.md), is listed
under each. Type `?modelname` (for example
[`?imm`](https://venpopov.com/bmm/dev/reference/imm.md)) for the
arguments of a model.

## Examples

``` r
bmm_models()
#> The following models are supported:
#> 
#> Continuous reproduction
#> 
#> - imm(): Interference Measurement Model
#> - mixture2p(): Two-Parameter Mixture Model
#> - mixture3p(): Three-Parameter Mixture Model
#> - sdm(): Signal Discrimination Model (SDM)
#> 
#> Categorical recall and n-AFC decisions
#> 
#> - m3(): Memory Measurement Model
#> - mpt(): Multinomial Processing Tree (MPT) Model
#> 
#> Detection, recognition and confidence judgments
#> 
#> - mpt(): Multinomial Processing Tree (MPT) Model
#> - sdt_cdp(): Continuous Dual-Process Signal Detection Theory Model
#> - sdt_mafc(): m-Alternative Forced Choice Signal Detection Theory Model
#> - sdt_ranking(): Ranking Signal Detection Theory Model
#> - sdt_rating(): Confidence Rating Signal Detection Theory Model
#> - sdt_yn(): Yes/No Signal Detection Theory Model
#> 
#> Choices and response times
#> 
#> - cswald(): Censored-Shifted Wald Model
#> - ddm(): Diffusion Decision Model
#> - ezdm(): EZ-Diffusion Model
#> 
#> Type  ?modelname  to get information about a specific model, e.g.  ?imm 
"imm" %in% bmm_models()
#> [1] TRUE
```
