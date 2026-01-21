# Package index

## About bmm

These pages provide a summary of the functionality available in bmm

- [`bmm-package`](https://venpopov.github.io/bmm/dev/reference/bmm-package.md)
  : Easy and Accesible Bayesian Measurement Models Using 'brms'

- [`supported_models()`](https://venpopov.github.io/bmm/dev/reference/supported_models.md)
  :

  Measurement models available in `bmm`

- [`bmm_options()`](https://venpopov.github.io/bmm/dev/reference/bmm_options.md)
  : View or change global bmm options

## Fitting models

Main functions for model fitting

- [`bmm()`](https://venpopov.github.io/bmm/dev/reference/bmm.md)
  [`fit_model()`](https://venpopov.github.io/bmm/dev/reference/bmm.md) :
  Fit Bayesian Measurement Models

- [`bmmformula()`](https://venpopov.github.io/bmm/dev/reference/bmmformula.md)
  [`bmf()`](https://venpopov.github.io/bmm/dev/reference/bmmformula.md)
  :

  Create formula for predicting parameters of a `bmmodel`

- [`update(`*`<bmmfit>`*`)`](https://venpopov.github.io/bmm/dev/reference/update.bmmfit.md)
  : Update a bmm model

- [`summary(`*`<bmmfit>`*`)`](https://venpopov.github.io/bmm/dev/reference/summary.bmmfit.md)
  :

  Create a summary of a fitted model represented by a `bmmfit` object

## Extract model information

Utility functions for extracting default priors, generate stan code,
stan data, etc.

- [`default_prior(`*`<bmmformula>`*`)`](https://venpopov.github.io/bmm/dev/reference/default_prior.bmmformula.md)
  : Get Default priors for Measurement Models specified in BMM

- [`extract_parameter_dimensions()`](https://venpopov.github.io/bmm/dev/reference/extract_parameter_dimensions.md)
  : Extract dimension from parameters in STAN parameter block

- [`extract_stan_blocks()`](https://venpopov.github.io/bmm/dev/reference/extract_stan_blocks.md)
  : Extract code from different STAN program blocks

- [`fit_info()`](https://venpopov.github.io/bmm/dev/reference/fit_info.md)
  : Extract information from a brmsfit object

- [`stancode(`*`<bmmformula>`*`)`](https://venpopov.github.io/bmm/dev/reference/stancode.bmmformula.md)
  : Generate Stan code for bmm models

- [`standata(`*`<bmmformula>`*`)`](https://venpopov.github.io/bmm/dev/reference/standata.bmmformula.md)
  :

  Stan data for `bmm` models

## Specifying models

Functions for specifying which model to fit

- [`imm()`](https://venpopov.github.io/bmm/dev/reference/imm.md)
  [`IMMfull()`](https://venpopov.github.io/bmm/dev/reference/imm.md)
  [`IMMbsc()`](https://venpopov.github.io/bmm/dev/reference/imm.md)
  [`IMMabc()`](https://venpopov.github.io/bmm/dev/reference/imm.md) :
  Interference measurement model by Oberauer and Lin (2017).
- [`m3()`](https://venpopov.github.io/bmm/dev/reference/m3.md) : The
  Multinomial / Memory Measurement Model
- [`mixture2p()`](https://venpopov.github.io/bmm/dev/reference/mixture2p.md)
  : Two-parameter mixture model by Zhang and Luck (2008).
- [`mixture3p()`](https://venpopov.github.io/bmm/dev/reference/mixture3p.md)
  : Three-parameter mixture model by Bays et al (2009).
- [`sdm()`](https://venpopov.github.io/bmm/dev/reference/sdm.md)
  [`sdmSimple()`](https://venpopov.github.io/bmm/dev/reference/sdm.md) :
  Signal Discrimination Model (SDM) by Oberauer (2023)

## Distributions

Functions for special distributions

- [`dimm()`](https://venpopov.github.io/bmm/dev/reference/IMMdist.md)
  [`pimm()`](https://venpopov.github.io/bmm/dev/reference/IMMdist.md)
  [`qimm()`](https://venpopov.github.io/bmm/dev/reference/IMMdist.md)
  [`rimm()`](https://venpopov.github.io/bmm/dev/reference/IMMdist.md) :
  Distribution functions for the Interference Measurement Model (IMM)
- [`dsdm()`](https://venpopov.github.io/bmm/dev/reference/SDMdist.md)
  [`psdm()`](https://venpopov.github.io/bmm/dev/reference/SDMdist.md)
  [`qsdm()`](https://venpopov.github.io/bmm/dev/reference/SDMdist.md)
  [`rsdm()`](https://venpopov.github.io/bmm/dev/reference/SDMdist.md) :
  Distribution functions for the Signal Discrimination Model (SDM)
- [`dm3()`](https://venpopov.github.io/bmm/dev/reference/m3dist.md)
  [`rm3()`](https://venpopov.github.io/bmm/dev/reference/m3dist.md) :
  Distribution functions for the Memory Measurement Model (M3)
- [`dmixture2p()`](https://venpopov.github.io/bmm/dev/reference/mixture2p_dist.md)
  [`pmixture2p()`](https://venpopov.github.io/bmm/dev/reference/mixture2p_dist.md)
  [`qmixture2p()`](https://venpopov.github.io/bmm/dev/reference/mixture2p_dist.md)
  [`rmixture2p()`](https://venpopov.github.io/bmm/dev/reference/mixture2p_dist.md)
  : Distribution functions for the two-parameter mixture model
  (mixture2p)
- [`dmixture3p()`](https://venpopov.github.io/bmm/dev/reference/mixture3p_dist.md)
  [`pmixture3p()`](https://venpopov.github.io/bmm/dev/reference/mixture3p_dist.md)
  [`qmixture3p()`](https://venpopov.github.io/bmm/dev/reference/mixture3p_dist.md)
  [`rmixture3p()`](https://venpopov.github.io/bmm/dev/reference/mixture3p_dist.md)
  : Distribution functions for the three-parameter mixture model
  (mixture3p)
- [`rejection_sampling()`](https://venpopov.github.io/bmm/dev/reference/rejection_sampling.md)
  : Rejection Sampling

## Data, model and parameter transformations

Utility functions for transforming data and parameters

- [`c_sqrtexp2bessel()`](https://venpopov.github.io/bmm/dev/reference/c_parametrizations.md)
  [`c_bessel2sqrtexp()`](https://venpopov.github.io/bmm/dev/reference/c_parametrizations.md)
  : Convert between parametrizations of the c parameter of the SDM
  distribution

- [`calc_error_relative_to_nontargets()`](https://venpopov.github.io/bmm/dev/reference/calc_error_relative_to_nontargets.md)
  : Calculate response error relative to non-target values

- [`deg2rad()`](https://venpopov.github.io/bmm/dev/reference/circle_transform.md)
  [`rad2deg()`](https://venpopov.github.io/bmm/dev/reference/circle_transform.md)
  : Convert degrees to radians or radians to degrees.

- [`construct_m3_act_funs()`](https://venpopov.github.io/bmm/dev/reference/construct_m3_act_funs.md)
  : Get Activation Functions for different M3 versions

- [`k2sd()`](https://venpopov.github.io/bmm/dev/reference/k2sd.md) :
  Transform kappa of the von Mises distribution to the circular standard
  deviation

- [`restructure(`*`<bmmfit>`*`)`](https://venpopov.github.io/bmm/dev/reference/restructure.bmmfit.md)
  :

  Restructure Old `bmmfit` Objects

- [`softmax()`](https://venpopov.github.io/bmm/dev/reference/softmax.md)
  [`softmaxinv()`](https://venpopov.github.io/bmm/dev/reference/softmax.md)
  : Softmax function and its inverse

- [`wrap()`](https://venpopov.github.io/bmm/dev/reference/wrap.md) :
  Wrap angles that extend beyond (-pi;pi)

## Datasets

Available datasets for fitting the models

- [`oberauer_lewandowsky_2019_e1`](https://venpopov.github.io/bmm/dev/reference/oberauer_lewandowsky_2019_e1.md)
  : Data from Experiment 1 reported by Oberauer & Lewandowsky (2019)
- [`oberauer_lin_2017`](https://venpopov.github.io/bmm/dev/reference/oberauer_lin_2017.md)
  : Data from Experiment 1 reported by Oberauer & Lin (2017)
- [`zhang_luck_2008`](https://venpopov.github.io/bmm/dev/reference/zhang_luck_2008.md)
  : Data from Experiment 2 reported by Zhang & Luck (2008)

## Developers’ corner

Functions to assist in developing new models

- [`apply_links()`](https://venpopov.github.io/bmm/dev/reference/apply_links.md)
  : Apply link functions for parameters in a formula or bmmformula

- [`bmf2bf()`](https://venpopov.github.io/bmm/dev/reference/bmf2bf.md) :

  Convert `bmmformula` objects to `brmsformula` objects

- [`check_data()`](https://venpopov.github.io/bmm/dev/reference/check_data.md)
  : Generic S3 method for checking data based on model type

- [`check_formula()`](https://venpopov.github.io/bmm/dev/reference/check_formula.md)
  : Generic S3 method for checking if the formula is valid for the
  specified model

- [`check_model()`](https://venpopov.github.io/bmm/dev/reference/check_model.md)
  : Generic S3 method for checking if the model is supported and model
  preprocessing

- [`configure_model()`](https://venpopov.github.io/bmm/dev/reference/configure_model.md)
  : Generic S3 method for configuring the model to be fit by brms

- [`configure_prior()`](https://venpopov.github.io/bmm/dev/reference/configure_prior.md)
  : Generic S3 method for configuring the default prior for a bmmodel

- [`postprocess_brm()`](https://venpopov.github.io/bmm/dev/reference/postprocess_brm.md)
  : Generic S3 method for postprocessing the fitted brm model

- [`revert_postprocess_brm()`](https://venpopov.github.io/bmm/dev/reference/revert_postprocess_brm.md)
  : Generic S3 method for reverting any postprocessing of the fitted brm
  model

- [`use_model_template()`](https://venpopov.github.io/bmm/dev/reference/use_model_template.md)
  : Create a file with a template for adding a new model (for
  developers)
