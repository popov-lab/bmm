# Package index

## About bmm

These pages provide a summary of the functionality available in bmm

- [`bmm-package`](https://venpopov.com/bmm/dev/reference/bmm-package.md)
  : Easy and Accesible Bayesian Measurement Models Using 'brms'

- [`bmm_models()`](https://venpopov.com/bmm/dev/reference/bmm_models.md)
  :

  Measurement models available in `bmm`

- [`bmm_options()`](https://venpopov.com/bmm/dev/reference/bmm_options.md)
  : View or change global bmm options

- [`bmm_setup()`](https://venpopov.com/bmm/dev/reference/bmm_setup.md) :
  Check whether this machine can fit bmm models

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

- [`conditional_effects(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/conditional_effects.bmmfit.md)
  : Conditional Effects for BMM Models

- [`recover_data(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/emmeans-bmmfit.md)
  [`emm_basis(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/emmeans-bmmfit.md)
  : emmeans support for bmmfit objects

- [`pp_check(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
  : Posterior predictive check for bmmfit objects

- [`ppp()`](https://venpopov.com/bmm/dev/reference/ppp.md) : Posterior
  predictive p-values for count models

- [`posterior_epred(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/posterior_epred.bmmfit.md)
  [`fitted(`*`<bmmfit>`*`)`](https://venpopov.com/bmm/dev/reference/posterior_epred.bmmfit.md)
  : Expected response of a fitted bmm model

## Model evaluation & post-processing

ROC/AUC, quantile-transformed ROC, and latent-distribution analysis for
signal detection models

- [`roc_sdt()`](https://venpopov.com/bmm/dev/reference/roc_sdt.md) :
  Model-implied ROC curve from a fitted SDT model
- [`roc_observed()`](https://venpopov.com/bmm/dev/reference/roc_observed.md)
  : Observed ROC points from a fitted SDT model
- [`auc_sdt()`](https://venpopov.com/bmm/dev/reference/auc_sdt.md) :
  Area under the ROC curve from a fitted SDT model
- [`latent_sdt()`](https://venpopov.com/bmm/dev/reference/latent_sdt.md)
  : Model-implied latent decision-variable distributions
- [`sdt_thresholds()`](https://venpopov.com/bmm/dev/reference/sdt_thresholds.md)
  : Latent decision thresholds from a fitted rating SDT model
- [`mratio()`](https://venpopov.com/bmm/dev/reference/mratio.md) :
  Metacognitive efficiency (M-ratio) from a fitted meta-d' SDT model
- [`sdt_sensitivity()`](https://venpopov.com/bmm/dev/reference/sdt_sensitivity.md)
  : Sensitivity on the noise, signal, or root-mean-square scale
- [`plot(`*`<bmm_sdt_roc>`*`)`](https://venpopov.com/bmm/dev/reference/plot.bmm_sdt_roc.md)
  : Plot a model-implied SDT ROC curve
- [`plot(`*`<bmm_sdt_auc>`*`)`](https://venpopov.com/bmm/dev/reference/plot.bmm_sdt_auc.md)
  : Plot the posterior AUC distribution from a fitted SDT model
- [`plot(`*`<bmm_sdt_latent>`*`)`](https://venpopov.com/bmm/dev/reference/plot.bmm_sdt_latent.md)
  : Plot the latent decision-variable distributions of an SDT model

## Extract model information

Utility functions for extracting default priors, generate stan code,
stan data, etc.

- [`bmm_data_check()`](https://venpopov.com/bmm/dev/reference/bmm_data_check.md)
  : Human-readable pre-fit data report for bmm models

- [`default_prior(`*`<bmmformula>`*`)`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
  : Get Default priors for Measurement Models specified in BMM

- [`extract_parameter_dimensions()`](https://venpopov.com/bmm/dev/reference/extract_parameter_dimensions.md)
  : Extract dimension from parameters in STAN parameter block

- [`extract_stan_blocks()`](https://venpopov.com/bmm/dev/reference/extract_stan_blocks.md)
  : Extract code from different STAN program blocks

- [`fit_info()`](https://venpopov.com/bmm/dev/reference/fit_info.md) :
  Extract information from a brmsfit object

- [`model_citation()`](https://venpopov.com/bmm/dev/reference/model_citation.md)
  : References for a measurement model

- [`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
  : Posterior draws of model parameters on the native scale

- [`parameter_info()`](https://venpopov.com/bmm/dev/reference/parameter_info.md)
  : Get parameter information for a bmm model

- [`pp_check_vars()`](https://venpopov.com/bmm/dev/reference/pp_check_vars.md)
  : List the posterior predictive checks available for a model fit

- [`prior_info()`](https://venpopov.com/bmm/dev/reference/prior_info.md)
  : Priors of a fitted bmm model and where they came from

- [`stancode(`*`<bmmformula>`*`)`](https://venpopov.com/bmm/dev/reference/stancode.bmmformula.md)
  : Generate Stan code for bmm models

- [`standata(`*`<bmmformula>`*`)`](https://venpopov.com/bmm/dev/reference/standata.bmmformula.md)
  :

  Stan data for `bmm` models

## Specifying models

Functions for specifying which model to fit

- [`cswald()`](https://venpopov.com/bmm/dev/reference/cswald.md) :
  Censored-Shifted Wald Model
- [`ddm()`](https://venpopov.com/bmm/dev/reference/ddm.md) : Diffusion
  Decision Model
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
- [`mpt()`](https://venpopov.com/bmm/dev/reference/mpt.md) : Multinomial
  Processing Tree (MPT) models
- [`mpt_from_eqn()`](https://venpopov.com/bmm/dev/reference/mpt_from_eqn.md)
  : Create an MPT model from an EQN model file
- [`mpt_from_string()`](https://venpopov.com/bmm/dev/reference/mpt_from_string.md)
  : Create an MPT model from an MPTinR-style model string
- [`mpt_tree()`](https://venpopov.com/bmm/dev/reference/mpt_tree.md) :
  Create a tree for a Multinomial Processing Tree (MPT) model
- [`plot(`*`<mpt_tree>`*`)`](https://venpopov.com/bmm/dev/reference/plot.mpt_tree.md)
  [`plot(`*`<mpt>`*`)`](https://venpopov.com/bmm/dev/reference/plot.mpt_tree.md)
  : Plot the structure of MPT trees
- [`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md)
  [`sdmSimple()`](https://venpopov.com/bmm/dev/reference/sdm.md) :
  Signal Discrimination Model (SDM) by Oberauer (2023)
- [`sdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp.md) :
  Continuous Dual-Process Signal Detection Theory Model
- [`sdt_mafc()`](https://venpopov.com/bmm/dev/reference/sdt_mafc.md) :
  m-Alternative Forced Choice Signal Detection Theory Model
- [`sdt_ranking()`](https://venpopov.com/bmm/dev/reference/sdt_ranking.md)
  : Ranking Signal Detection Theory Model
- [`sdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating.md)
  : Confidence Rating Signal Detection Theory Model
- [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) :
  Yes/No Signal Detection Theory Model

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

- [`sdt_d()`](https://venpopov.com/bmm/dev/reference/SDTdist.md)
  [`sdt_criterion()`](https://venpopov.com/bmm/dev/reference/SDTdist.md)
  [`sdt_bias()`](https://venpopov.com/bmm/dev/reference/SDTdist.md) :
  Utility functions for Signal Detection Theory

- [`dcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md)
  [`rcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md)
  [`pcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md)
  [`qcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md) :

  Distribution functions for the censored shifted Wald model (`cswald`)

- [`dddm()`](https://venpopov.com/bmm/dev/reference/ddm_dist.md)
  [`rddm()`](https://venpopov.com/bmm/dev/reference/ddm_dist.md) :

  Distribution function for the Diffusion Decision Model (`ddm`)

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

- [`dmpt()`](https://venpopov.com/bmm/dev/reference/mptdist.md)
  [`rmpt()`](https://venpopov.com/bmm/dev/reference/mptdist.md) :
  Distribution functions for Multinomial Processing Tree (MPT) models

- [`rejection_sampling()`](https://venpopov.com/bmm/dev/reference/rejection_sampling.md)
  : Rejection Sampling

- [`dsdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp_dist.md)
  [`rsdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp_dist.md)
  : Distribution functions for Continuous Dual-Process SDT (CDP)

- [`dsdt_dpsdt()`](https://venpopov.com/bmm/dev/reference/sdt_dpsdt_dist.md)
  [`rsdt_dpsdt()`](https://venpopov.com/bmm/dev/reference/sdt_dpsdt_dist.md)
  : Distribution functions for dual-process SDT (DPSDT)

- [`dsdt_mafc()`](https://venpopov.com/bmm/dev/reference/sdt_mafc_dist.md)
  [`rsdt_mafc()`](https://venpopov.com/bmm/dev/reference/sdt_mafc_dist.md)
  : Distribution functions for m-AFC SDT

- [`dsdt_metad()`](https://venpopov.com/bmm/dev/reference/sdt_metad_dist.md)
  [`rsdt_metad()`](https://venpopov.com/bmm/dev/reference/sdt_metad_dist.md)
  : Distribution functions for meta-d' SDT

- [`dsdt_ranking()`](https://venpopov.com/bmm/dev/reference/sdt_ranking_dist.md)
  [`rsdt_ranking()`](https://venpopov.com/bmm/dev/reference/sdt_ranking_dist.md)
  : Distribution functions for Ranking SDT

- [`dsdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating_dist.md)
  [`rsdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating_dist.md)
  : Distribution functions for Confidence Rating SDT

- [`dsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md)
  [`rsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md) :
  Distribution functions for Yes/No SDT

## Data, model and parameter transformations

Utility functions for transforming data and parameters

- [`adjust_ezdm_accuracy()`](https://venpopov.com/bmm/dev/reference/adjust_ezdm_accuracy.md)
  : Adjust Accuracy Counts for Contamination (deprecated)

- [`aggregate_sdt_cdp_data()`](https://venpopov.com/bmm/dev/reference/aggregate_sdt_cdp_data.md)
  :

  Aggregate long-format Remember/Know data for
  [`sdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp.md)

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

- [`mpt_data_reshape()`](https://venpopov.com/bmm/dev/reference/mpt_data_reshape.md)
  :

  Reshape wide MPT frequency tables to the long format of
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)

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

- [`broeder_schuetz_2009_e3`](https://venpopov.com/bmm/dev/reference/broeder_schuetz_2009_e3.md)
  : Recognition ROC data from Broeder & Schuetz (2009, Experiment 3)
- [`data_color_judgement_task`](https://venpopov.com/bmm/dev/reference/data_color_judgement_task.md)
  : Example Data from a Color Judgement Task
- [`meyer_grant_jakob_2025`](https://venpopov.com/bmm/dev/reference/meyer_grant_jakob_2025.md)
  : Ranking recognition data from Meyer-Grant & Jakob (2025)
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

- [`data_check_finding()`](https://venpopov.com/bmm/dev/reference/data_check_finding.md)
  : Construct a single finding for a pre-fit data report

- [`data_check_findings()`](https://venpopov.com/bmm/dev/reference/data_check_findings.md)
  : Generic S3 method for model-specific pre-fit data findings

- [`native_transform()`](https://venpopov.com/bmm/dev/reference/native_transform.md)
  : Transform parameter draws from the sampling scale to the native
  scale

- [`postprocess_brm()`](https://venpopov.com/bmm/dev/reference/postprocess_brm.md)
  : Generic S3 method for postprocessing the fitted brm model

- [`pp_observables()`](https://venpopov.com/bmm/dev/reference/pp_observables.md)
  [`pp_simulate()`](https://venpopov.com/bmm/dev/reference/pp_observables.md)
  : Declare the observables of a bmm model for posterior predictive
  checks

- [`print_model_details()`](https://venpopov.com/bmm/dev/reference/print_model_details.md)
  :

  Generic S3 method for model-specific lines in
  [`print()`](https://rdrr.io/r/base/print.html) of a model

- [`revert_postprocess_brm()`](https://venpopov.com/bmm/dev/reference/revert_postprocess_brm.md)
  : Generic S3 method for reverting any postprocessing of the fitted brm
  model

- [`use_model_template()`](https://venpopov.com/bmm/dev/reference/use_model_template.md)
  : Create a file with a template for adding a new model (for
  developers)
