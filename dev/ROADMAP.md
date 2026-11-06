# bmm roadmap

Updated 2026-11-06. The reasoning is in
`.dev/decisions/0008-roadmap-to-2027-03.md`; the live list of work is in
the GitHub milestones
[1.5.0](https://github.com/popov-lab/bmm/milestone/8),
[1.6.0](https://github.com/popov-lab/bmm/milestone/9) and Future.

## What ships when

**bmm 1.4.2** is on CRAN (6 November 2026). It fixes wrong
[`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md) estimates in
some designs and a handful of smaller bugs; install it with
`install.packages("bmm")`.

**bmm 1.5.0**, CRAN submission 1 February 2027, on CRAN by mid-February.
This is the release the TeaP 2027 pre-conference workshop (Mannheim, 21
March 2027) teaches from.

- A documented description of every model, `model_spec()`, and
  [`simulate_data()`](https://venpopov.com/bmm/dev/reference/simulate_data.md)
  to simulate data from any model, so that packages built on bmm
  (bmmtools, the model explorer) no longer read bmm’s internals.
- The circular mixture models
  ([`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md),
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md),
  [`imm()`](https://venpopov.com/bmm/dev/reference/imm.md)) rebuilt as
  custom families, with variable-precision and capacity versions. Your
  1.4 code and fits keep working; two parameter names change (`mu1`
  becomes `mu`, `pnt` becomes `ptarget`) and the old names warn.
- [`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md) with a
  Laplace activation, peak and range scalings, and discrete response
  options.
- `psychometric_function()`, a new model for detection and
  discrimination data.
- One row per trial for the signal detection and MPT models, so item
  random effects are possible.
- The target confusability competition model `tcc()` and the EZ circular
  diffusion model `ezcdm()`, if their reviews finish before 18 December
  2026; otherwise they move to 1.6.0.
- A written policy on default priors: they change only in a minor
  release and are listed in NEWS.

**bmm 1.6.0**, CRAN submission 1 June 2027, before the MathPsych
meeting.

- Other packages can define bmm models: `register_bmm_model()` and a
  versioned model contract.
- Change detection versions of the visual working memory models.
- The racing models: linear ballistic accumulator, log-normal race,
  racing diffusion.
- A cumulative response model for Likert items.
- Fits record the default-prior set they used, and
  `bmm_options(prior_version = )` reproduces an older set.
- Multivariate models, if their interface is settled by mid-January
  2027; otherwise 1.7.0.
- The forms deprecated in bmm 0.4 and 1.0 are removed
  ([`fit_model()`](https://venpopov.com/bmm/dev/reference/bmm.md),
  [`IMMfull()`](https://venpopov.com/bmm/dev/reference/imm.md),
  [`sdmSimple()`](https://venpopov.com/bmm/dev/reference/sdm.md), …), as
  are
  [`supported_models()`](https://venpopov.com/bmm/dev/reference/supported_models.md)
  and
  [`parameters()`](https://venpopov.com/bmm/dev/reference/parameters.md).
  The `log` default of
  [`dm3()`](https://venpopov.com/bmm/dev/reference/m3dist.md),
  [`dmpt()`](https://venpopov.com/bmm/dev/reference/mptdist.md),
  [`dddm()`](https://venpopov.com/bmm/dev/reference/ddm_dist.md),
  [`dcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md)
  and [`dezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md)
  becomes `FALSE`.

## Not planned

Models without a tractable likelihood (diffusion models for conflict
tasks), models with recursion over trial sequences (reinforcement
learning, trial-by-trial category learning), free-recall architectures,
and simulation-based likelihoods are not planned for bmm itself. A
reinforcement-learning family is a candidate for a separate package
built on the extension API.

## How to contribute a model

Read `.github/CONTRIBUTING.md` and the developer notes
(<https://venpopov.com/bmm/dev-notes/>). In short: open an issue with
the “New model” template, build the model with
[`use_model_template()`](https://venpopov.com/bmm/dev/reference/use_model_template.md),
add a recovery plan under `tests/recovery/` and attach its bmmtools
validation card, and open the pull request against `develop` (or
`dev-1.6.0` for a model planned for 1.6.0). A new model merges with one
approving review and a green card.
