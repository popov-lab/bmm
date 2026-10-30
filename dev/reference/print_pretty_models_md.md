# Generate a markdown list of the measurement models available in `bmm`

Used internally to populate the README and the "Get started" article.
Models are grouped as in
[`bmm_models()`](https://venpopov.com/bmm/dev/reference/bmm_models.md),
and every model links to its reference page on the website.

## Usage

``` r
print_pretty_models_md(group = NULL)
```

## Arguments

- group:

  Optional character vector of group labels as printed by
  [`bmm_models()`](https://venpopov.com/bmm/dev/reference/bmm_models.md).
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
#> - [`imm()`](https://venpopov.com/bmm/reference/imm.html): Interference Measurement Model
#> - [`mixture2p()`](https://venpopov.com/bmm/reference/mixture2p.html): Two-Parameter Mixture Model
#> - [`mixture3p()`](https://venpopov.com/bmm/reference/mixture3p.html): Three-Parameter Mixture Model
#> - [`sdm()`](https://venpopov.com/bmm/reference/sdm.html): Signal Discrimination Model (SDM)
#> 
#> **Categorical recall and n-AFC decisions**
#> 
#> - [`m3()`](https://venpopov.com/bmm/reference/m3.html): Memory Measurement Model
#> - [`mpt()`](https://venpopov.com/bmm/reference/mpt.html): Multinomial Processing Tree (MPT) Model
#> 
#> **Detection, recognition and confidence judgments**
#> 
#> - [`mpt()`](https://venpopov.com/bmm/reference/mpt.html): Multinomial Processing Tree (MPT) Model
#> - [`sdt_cdp()`](https://venpopov.com/bmm/reference/sdt_cdp.html): Continuous Dual-Process Signal Detection Theory Model
#> - [`sdt_mafc()`](https://venpopov.com/bmm/reference/sdt_mafc.html): m-Alternative Forced Choice Signal Detection Theory Model
#> - [`sdt_ranking()`](https://venpopov.com/bmm/reference/sdt_ranking.html): Ranking Signal Detection Theory Model
#> - [`sdt_rating()`](https://venpopov.com/bmm/reference/sdt_rating.html): Confidence Rating Signal Detection Theory Model
#> - [`sdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn.html): Yes/No Signal Detection Theory Model
#> 
#> **Choices and response times**
#> 
#> - [`cswald()`](https://venpopov.com/bmm/reference/cswald.html): Censored-Shifted Wald Model
#> - [`ddm()`](https://venpopov.com/bmm/reference/ddm.html): Diffusion Decision Model
#> - [`ezdm()`](https://venpopov.com/bmm/reference/ezdm.html): EZ-Diffusion Model
#> 
print_pretty_models_md(group = "Continuous reproduction")
#> - [`imm()`](https://venpopov.com/bmm/reference/imm.html): Interference Measurement Model
#> - [`mixture2p()`](https://venpopov.com/bmm/reference/mixture2p.html): Two-Parameter Mixture Model
#> - [`mixture3p()`](https://venpopov.com/bmm/reference/mixture3p.html): Three-Parameter Mixture Model
#> - [`sdm()`](https://venpopov.com/bmm/reference/sdm.html): Signal Discrimination Model (SDM)
#> 
```
