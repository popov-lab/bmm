# Recognition ROC data from Broeder & Schuetz (2009, Experiment 3)

Binary old/new recognition data from 40 subjects, aggregated to response
counts. Each subject was tested under five base-rate conditions: the
proportion of old items shifts the decision criterion from conservative
(`br1`) to liberal (`br5`) while leaving sensitivity unchanged, tracing
a five-point binary ROC per subject. That is the cleanest design for
estimating the unequal-variance ratio (`sdratio`) of
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md), which
needs more than one operating point: a single condition with no other
varying predictor yields only one hit/false-alarm pair and cannot
separate a wider signal distribution from a larger d'. Counts were
digitised from the frequencies reported in the original article.

## Usage

``` r
broeder_schuetz_2009_e3
```

## Format

### `broeder_schuetz_2009_e3`

A data frame with 400 rows (40 subjects x 5 conditions x 2 stimulus
types) and 5 columns:

- id:

  Integer uniquely identifying each subject

- condition:

  Factor with five base-rate conditions, ordered from the most
  conservative (`br1`) to the most liberal (`br5`) induced criterion

- stimulus:

  Integer stimulus type: 0 = new/lure, 1 = old/target

- n_old:

  Integer count of "old" responses in that cell: hits for old items
  (`stimulus == 1`) and false alarms for new items (`stimulus == 0`)

- n_trials:

  Integer number of items presented in that cell

## Source

Broeder, A., & Schuetz, J. (2009). Recognition ROCs are curvilinear—or
are they? On premature arguments against the two-high-threshold model of
recognition. *Journal of Experimental Psychology: Learning, Memory, and
Cognition*, 35(3), 587–606.
[doi:10.1037/a0015279](https://doi.org/10.1037/a0015279)

## Examples

``` r
if (FALSE) { # \dontrun{
# Unequal-variance yes/no SDT: the criterion varies across base-rate
# conditions, while sensitivity (d) and the signal/noise SD ratio (sdratio)
# are held constant across conditions.
model <- sdt_yn(
  response = "n_old", stimulus = "stimulus", n_trials = "n_trials"
)
fit <- bmm(
  formula = bmf(
    d ~ 1 + (1 | id),
    criterion ~ 0 + condition + (1 | id),
    sdratio ~ 1
  ),
  data = broeder_schuetz_2009_e3,
  model = model,
  backend = "cmdstanr"
)
} # }
```
