# Generate a markdown list of the measurement models available in `bmm`

Used internally to populate the README and the "Get started" article.
Models are grouped as in
[`supported_models()`](https://popov-lab.github.io/bmm/dev/reference/supported_models.md),
and every model links to its reference page on the website.

## Usage

``` r
print_pretty_models_md(group = NULL)
```

## Arguments

- group:

  Optional character vector of group labels as printed by
  [`supported_models()`](https://popov-lab.github.io/bmm/dev/reference/supported_models.md).
  Only those groups are listed and the group headers are omitted, so a
  document can add its own text per group.

## Value

Markdown code for printing the list of measurement models available in
`bmm`

## Examples

``` r
print_pretty_models_md()
#> **Continuous reproduction**
#> 
#> - [`imm()`](https://popov-lab.github.io/bmm/reference/imm.html): Interference measurement model by Oberauer and Lin (2017)
#> - [`mixture2p()`](https://popov-lab.github.io/bmm/reference/mixture2p.html): Two-parameter mixture model by Zhang and Luck (2008)
#> - [`mixture3p()`](https://popov-lab.github.io/bmm/reference/mixture3p.html): Three-parameter mixture model by Bays et al (2009)
#> - [`sdm()`](https://popov-lab.github.io/bmm/reference/sdm.html): Signal Discrimination Model (SDM) by Oberauer (2023)
#> 
#> **Categorical recall and n-AFC decisions**
#> 
#> - [`m3()`](https://popov-lab.github.io/bmm/reference/m3.html): The Multinomial / Memory Measurement Model
#> 
#> **Detection, recognition and confidence judgments**
#> 
#> - [`sdt_cdp()`](https://popov-lab.github.io/bmm/reference/sdt_cdp.html): Continuous Dual-Process Signal Detection Theory (CDP)
#> - [`sdt_mafc()`](https://popov-lab.github.io/bmm/reference/sdt_mafc.html): Signal Detection Theory (m-AFC)
#> - [`sdt_ranking()`](https://popov-lab.github.io/bmm/reference/sdt_ranking.html): Signal Detection Theory (Ranking)
#> - [`sdt_rating()`](https://popov-lab.github.io/bmm/reference/sdt_rating.html): Signal Detection Theory (Confidence Rating)
#> - [`sdt_yn()`](https://popov-lab.github.io/bmm/reference/sdt_yn.html): Signal Detection Theory (Yes/No)
#> 
#> **Choices and response times**
#> 
#> - [`cswald()`](https://popov-lab.github.io/bmm/reference/cswald.html): Censored-Shifted Wald Model
#> - [`ddm()`](https://popov-lab.github.io/bmm/reference/ddm.html): Diffusion Decision Model
#> - [`ezdm()`](https://popov-lab.github.io/bmm/reference/ezdm.html): EZ-Diffusion Model
#> 
print_pretty_models_md(group = "Continuous reproduction")
#> - [`imm()`](https://popov-lab.github.io/bmm/reference/imm.html): Interference measurement model by Oberauer and Lin (2017)
#> - [`mixture2p()`](https://popov-lab.github.io/bmm/reference/mixture2p.html): Two-parameter mixture model by Zhang and Luck (2008)
#> - [`mixture3p()`](https://popov-lab.github.io/bmm/reference/mixture3p.html): Three-parameter mixture model by Bays et al (2009)
#> - [`sdm()`](https://popov-lab.github.io/bmm/reference/sdm.html): Signal Discrimination Model (SDM) by Oberauer (2023)
#> 
```
