# Draft a methods paragraph for a fitted model

Drafts the text that describes a model fitted with
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) in the methods
section of a paper: the model, the software and sampler settings, the
priors, and the convergence diagnostics. The draft is meant to be
edited, not pasted unread.

## Usage

``` r
report_methods(fit, sections = c("model", "priors", "convergence"))
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md). It must
  contain posterior draws, so mock fits and fits made with `chains = 0`
  or `empty = TRUE` are refused.

- sections:

  The parts of the paragraph to draft, in the order they should appear:
  `"model"` (the model, the versions of bmm, brms, Stan, the backend and
  R, and the sampler settings), `"priors"` (the sentences of
  [`prior_info()`](https://venpopov.com/bmm/dev/reference/prior_info.md)
  with `format = "text"`) and `"convergence"` (the sentence of
  [`report_convergence()`](https://venpopov.com/bmm/dev/reference/report_convergence.md)).

## Value

An object of class `bmm_report_methods`, a list with

- `paragraph`: the drafted text, one string;

- `formula`: the bmm formula with its constants, one line per parameter
  (empty without the `"model"` section);

- `references`: the references, numbered by position;

- `notes`: what you have to check or add, not meant for the paper;

- `sections`: the sections drafted, in order.

[`print()`](https://rdrr.io/r/base/print.html) shows all of it;
[`format()`](https://rdrr.io/r/base/format.html) returns the paragraph,
the formula and the numbered references as one string to paste into a
document. The references keep the Markdown emphasis marks that
[`citation()`](https://rdrr.io/r/utils/citation.html) writes
(`_journal_`, `*volume*`): Quarto renders them as italics; in Word,
remove them.

## Details

Every value in the text is read from the fit, and every reference comes
from the model's constructor (see
[`model_citation()`](https://venpopov.com/bmm/dev/reference/model_citation.md))
or from [`citation()`](https://rdrr.io/r/utils/citation.html) of an
installed package. Nothing else enters the text:

- The versions are those the fit stores. The R version is stored by fits
  made with bmm 1.5.0 or later; for older fits it is left out and a note
  asks you to add it.

- The references are those of the packages installed now. A note says
  when an installed version differs from the version the fit was made
  with.

- The convergence sentence gives numbers, not verdicts: the largest
  R-hat and the smallest bulk and tail effective sample sizes over all
  estimated parameters (constants, `lp__` and `lprior` excluded), and
  the number of divergent transitions after warmup. Where the summary of
  brms would warn (an R-hat above 1.05 among the parameters it shows, or
  a divergent transition), a note says so.

- For algorithms other than `"sampling"`, the convergence sentence is
  replaced by a note, as there are no Markov chains to diagnose.

In the text, the numbers in square brackets refer to the reference list.
A model without a reference gets `[CITATION NEEDED]` in that list.

## See also

[`report_convergence()`](https://venpopov.com/bmm/dev/reference/report_convergence.md),
[`prior_info()`](https://venpopov.com/bmm/dev/reference/prior_info.md),
[`model_citation()`](https://venpopov.com/bmm/dev/reference/model_citation.md)

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
fit <- bmm(
  bmf(c ~ 0 + set_size, kappa ~ 1),
  data = oberauer_lin_2017,
  model = sdm(resp_error = "dev_rad")
)
report_methods(fit)
report_methods(fit, sections = c("model", "convergence"))
cat(format(report_methods(fit)))
}
```
