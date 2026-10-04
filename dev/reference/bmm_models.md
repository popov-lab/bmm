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
for, one line per model with its constructor and full name. The groups
are: continuous reproduction; categorical recall and n-AFC decisions;
detection, recognition and confidence judgments; choices and response
times. Type `?modelname` (for example
[`?imm`](https://popov-lab.github.io/bmm/dev/reference/imm.md)) for the
arguments of a model.

## Examples

``` r
bmm_models()
#> The following models are supported:
#> 
#> Continuous reproduction
#> 
#> - imm(): Interference measurement model by Oberauer and Lin (2017)
#> - mixture2p(): Two-parameter mixture model by Zhang and Luck (2008)
#> - mixture3p(): Three-parameter mixture model by Bays et al (2009)
#> - sdm(): Signal Discrimination Model (SDM) by Oberauer (2023)
#> 
#> Categorical recall and n-AFC decisions
#> 
#> - m3(): The Multinomial / Memory Measurement Model
#> 
#> Detection, recognition and confidence judgments
#> 
#> - sdt_cdp(): Continuous Dual-Process Signal Detection Theory (CDP)
#> - sdt_mafc(): Signal Detection Theory (m-AFC)
#> - sdt_ranking(): Signal Detection Theory (Ranking)
#> - sdt_rating(): Signal Detection Theory (Confidence Rating)
#> - sdt_yn(): Signal Detection Theory (Yes/No)
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
