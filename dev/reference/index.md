# Package index

## About bmm

These pages provide a summary of the functionality available in bmm

- [`bmm-package`](https://venpopov.com/bmm/dev/reference/bmm-package.md)
  : Easy and Accesible Bayesian Measurement Models Using 'brms'

- [`supported_models()`](https://venpopov.com/bmm/dev/reference/supported_models.md)
  :

  Measurement models available in `bmm`

- [`bmm_options()`](https://venpopov.com/bmm/dev/reference/bmm_options.md)
  : View or change global bmm options

## Fitting models

Main functions for model fitting

- [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)
  [`fit_model()`](https://venpopov.com/bmm/dev/reference/bmm.md) : Fit
  Bayesian Measurement Models

- [`bmmformula()`](https://venpopov.com/bmm/dev/reference/bmmformula.md)
  [`bmf()`](https://venpopov.com/bmm/dev/reference/bmmformula.md) :

  Create formula for predicting parameters of a `bmmodel`

- [`update(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/update.bmmfit.md)
  : Update a bmm model

- [`summary(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/summary.bmmfit.md)
  :

  Create a summary of a fitted model represented by a `bmmfit` object

- [`recover_data(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/emmeans-bmmfit.md)
  [`emm_basis(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/emmeans-bmmfit.md)
  : emmeans support for bmmfit objects

- [`pp_check(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
  : Posterior predictive check for bmmfit objects

## Extract model information

Utility functions for extracting default priors, generate stan code,
stan data, etc.

- [`default_prior(`*`<bmmformula>`*`)`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
  : Get Default priors for Measurement Models specified in BMM

- [`extract_parameter_dimensions()`](https://venpopov.com/bmm/dev/reference/extract_parameter_dimensions.md)
  : Extract dimension from parameters in STAN parameter block

- [`extract_stan_blocks()`](https://venpopov.com/bmm/dev/reference/extract_stan_blocks.md)
  : Extract code from different STAN program blocks

- [`fit_info()`](https://venpopov.com/bmm/dev/reference/fit_info.md) :
  Extract information from a brmsfit object

- [`stancode(`*`<bmmformula>`*`)`](https://venpopov.com/bmm/dev/reference/stancode.bmmformula.md)
  : Generate Stan code for bmm models

- [`standata(`*`<bmmformula>`*`)`](https://venpopov.com/bmm/dev/reference/standata.bmmformula.md)
  :

  Stan data for `bmm` models

- [`parameters()`](https://venpopov.com/bmm/dev/reference/parameters.md)
  : Get parameter information for a bmm model

## Specifying models

Functions for specifying which model to fit

- [`cswald()`](https://venpopov.com/bmm/dev/reference/cswald.md) :
  Censored-Shifted Wald Model
- [`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md) :
  EZ-Diffusion Model
- [`imm()`](https://venpopov.com/bmm/dev/reference/imm.md)
  [`IMMfull()`](https://venpopov.com/bmm/dev/reference/imm.md)
  [`IMMbsc()`](https://venpopov.com/bmm/dev/reference/imm.md)
  [`IMMabc()`](https://venpopov.com/bmm/dev/reference/imm.md) :
  Interference measurement model by Oberauer and Lin (2017).
- [`m3()`](https://venpopov.com/bmm/dev/reference/m3.md) : The
  Multinomial / Memory Measurement Model
- [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md) :
  Two-parameter mixture model by Zhang and Luck (2008).
- [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md) :
  Three-parameter mixture model by Bays et al (2009).
- [`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md)
  [`sdmSimple()`](https://venpopov.com/bmm/dev/reference/sdm.md) :
  Signal Discrimination Model (SDM) by Oberauer (2023)

## Distributions

Functions for special distributions

- [`dimm()`](https://venpopov.com/bmm/dev/reference/IMMdist.md)
  [`pimm()`](https://venpopov.com/bmm/dev/reference/IMMdist.md)
  [`qimm()`](https://venpopov.com/bmm/dev/reference/IMMdist.md)
  [`rimm()`](https://venpopov.com/bmm/dev/reference/IMMdist.md) :
  Distribution functions for the Interference Measurement Model (IMM)

- [`dsdm()`](https://venpopov.com/bmm/dev/reference/SDMdist.md)
  [`psdm()`](https://venpopov.com/bmm/dev/reference/SDMdist.md)
  [`qsdm()`](https://venpopov.com/bmm/dev/reference/SDMdist.md)
  [`rsdm()`](https://venpopov.com/bmm/dev/reference/SDMdist.md) :
  Distribution functions for the Signal Discrimination Model (SDM)

- [`dcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md)
  [`rcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md)
  [`pcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md)
  [`qcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md) :

  Distribution functions for the censored shifted Wald model (`cswald`)

- [`dezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md)
  [`rezdm()`](https://venpopov.com/bmm/dev/reference/ezdm_dist.md) :
  Distribution functions for the EZ-Diffusion Model (ezdm)

- [`dm3()`](https://venpopov.com/bmm/dev/reference/m3dist.md)
  [`rm3()`](https://venpopov.com/bmm/dev/reference/m3dist.md) :
  Distribution functions for the Memory Measurement Model (M3)

- [`dmixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p_dist.md)
  [`pmixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p_dist.md)
  [`qmixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p_dist.md)
  [`rmixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p_dist.md)
  : Distribution functions for the two-parameter mixture model
  (mixture2p)

- [`dmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md)
  [`pmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md)
  [`qmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md)
  [`rmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md)
  : Distribution functions for the three-parameter mixture model
  (mixture3p)

- [`rejection_sampling()`](https://venpopov.com/bmm/dev/reference/rejection_sampling.md)
  : Rejection Sampling

## Data, model and parameter transformations

Utility functions for transforming data and parameters

- [`adjust_ezdm_accuracy()`](https://venpopov.com/bmm/dev/reference/adjust_ezdm_accuracy.md)
  : Adjust Accuracy Counts for Contamination

- [`c_sqrtexp2bessel()`](https://venpopov.com/bmm/dev/reference/c_parametrizations.md)
  [`c_bessel2sqrtexp()`](https://venpopov.com/bmm/dev/reference/c_parametrizations.md)
  : Convert between parametrizations of the c parameter of the SDM
  distribution

- [`calc_error_relative_to_nontargets()`](https://venpopov.com/bmm/dev/reference/calc_error_relative_to_nontargets.md)
  : Calculate response error relative to non-target values

- [`deg2rad()`](https://venpopov.com/bmm/dev/reference/circle_transform.md)
  [`rad2deg()`](https://venpopov.com/bmm/dev/reference/circle_transform.md)
  : Convert degrees to radians or radians to degrees.

- [`construct_m3_act_funs()`](https://venpopov.com/bmm/dev/reference/construct_m3_act_funs.md)
  : Get Activation Functions for different M3 versions

- [`ezdm_summary_stats()`](https://venpopov.com/bmm/dev/reference/ezdm_summary_stats.md)
  : Compute Robust Summary Statistics for EZ-Diffusion Model

- [`flag_contaminant_rts()`](https://venpopov.com/bmm/dev/reference/flag_contaminant_rts.md)
  : Flag contaminant reaction times using mixture modeling

- [`k2sd()`](https://venpopov.com/bmm/dev/reference/k2sd.md) : Transform
  kappa of the von Mises distribution to the circular standard deviation

- [`restructure(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/restructure.bmmfit.md)
  :

  Restructure Old `bmmfit` Objects

- [`softmax()`](https://venpopov.com/bmm/dev/reference/softmax.md)
  [`softmaxinv()`](https://venpopov.com/bmm/dev/reference/softmax.md) :
  Softmax function and its inverse

- [`validate_fast_guesses()`](https://venpopov.com/bmm/dev/reference/validate_fast_guesses.md)
  : Test if fast contaminants show random guessing behavior

- [`wrap()`](https://venpopov.com/bmm/dev/reference/wrap.md) : Wrap
  angles that extend beyond (-pi;pi)

## Datasets

Available datasets for fitting the models

- [`oberauer_lewandowsky_2019_e1`](https://venpopov.com/bmm/dev/reference/oberauer_lewandowsky_2019_e1.md)
  : Data from Experiment 1 reported by Oberauer & Lewandowsky (2019)
- [`oberauer_lin_2017`](https://venpopov.com/bmm/dev/reference/oberauer_lin_2017.md)
  : Data from Experiment 1 reported by Oberauer & Lin (2017)
- [`zhang_luck_2008`](https://venpopov.com/bmm/dev/reference/zhang_luck_2008.md)
  : Data from Experiment 2 reported by Zhang & Luck (2008)

## Developers’ corner

Functions to assist in developing new models

- [`apply_links()`](https://venpopov.com/bmm/dev/reference/apply_links.md)
  : Apply link functions for parameters in a formula or bmmformula

- [`bmf2bf()`](https://venpopov.com/bmm/dev/reference/bmf2bf.md) :

  Convert `bmmformula` objects to `brmsformula` objects

- [`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md)
  : Generic S3 method for checking data based on model type

- [`check_formula()`](https://venpopov.com/bmm/dev/reference/check_formula.md)
  : Generic S3 method for checking if the formula is valid for the
  specified model

- [`check_model()`](https://venpopov.com/bmm/dev/reference/check_model.md)
  : Generic S3 method for checking if the model is supported and model
  preprocessing

- [`configure_model()`](https://venpopov.com/bmm/dev/reference/configure_model.md)
  : Generic S3 method for configuring the model to be fit by brms

- [`configure_prior()`](https://venpopov.com/bmm/dev/reference/configure_prior.md)
  : Generic S3 method for configuring the default prior for a bmmodel

- [`create_initfun()`](https://venpopov.com/bmm/dev/reference/create_initfun.md)
  : Generic S3 method for creating an initial values function

- [`postprocess_brm()`](https://venpopov.com/bmm/dev/reference/postprocess_brm.md)
  : Generic S3 method for postprocessing the fitted brm model

- [`revert_postprocess_brm()`](https://venpopov.com/bmm/dev/reference/revert_postprocess_brm.md)
  : Generic S3 method for reverting any postprocessing of the fitted brm
  model

- [`use_model_template()`](https://venpopov.com/bmm/dev/reference/use_model_template.md)
  : Create a file with a template for adding a new model (for
  developers)
