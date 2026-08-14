# bmm (development version)

### New features
* Add `softplus` as an opt-in link function for positively-bounded parameters, as an alternative to the default `log` link. `softplus(x) = log(1 + exp(x))` keeps parameters positive while growing linearly for large values, avoiding the numerical blow-up of `exp()` and giving predictor effects an additive (rather than multiplicative) interpretation on the natural scale. Enable it per parameter via the model's `links` list, e.g. `m3(...)$links <- list(c = "softplus", a = "softplus")` or `ddm(rt, response, links = list(bound = "softplus"))` (#363).

### Bug fixes
* Fix initial values being set in two places, where the `init` returned by `configure_model()` was silently overwritten by `create_initfun()`. This caused the **m3** model's intended `init = 0` (needed for stable sampling with the `simple` choice rule and an `identity` link) to be lost, and left dead `init` code in the **sdm** model. `create_initfun()` is now the single source of truth for initial values, with model-specific behaviour expressed through S3 methods (#375).
* Fix `.pwald()` returning `NaN`/`-Inf` in the upper tail of the shifted-Wald survival function, which propagated to `dcswald()` (and therefore `log_lik`/`posterior_predict`) for the **cswald** model at extreme reaction times. The R-side survival now uses the stable `log_diff_exp` form already used by the Stan likelihood (`swald_lccdf`) (#376).
* Fix `print()` for model summaries selecting regression-coefficient rows by an unanchored substring match, so a parameter such as `a` could pull in rows of another parameter like `kappa` (e.g. in the **imm** model). Rows are now matched on the exact parameter prefix. This also fixes a crash when only a single coefficient row is shown (#379, #369).
* `create_initfun()` now matches Stan parameters to model parameters with a word-boundary regex (`(^|_)param(_|$)`) instead of a substring match, preventing collisions in models with short parameter names (e.g. `s`, `c`, `a`) that are substrings of longer ones (`sim`, `correct`, `activation`); the longest (most specific) match is selected when several apply (#354, #355).
* `create_initfun()` now resolves initialization terms from `nlpars` when a model parameter is not a distributional parameter, so models built as non-linear brms formulas (e.g. native-multinomial models whose parameters live in `bterms$nlpars`) no longer error with `no applicable method for 'has_intercept' applied to an object of class "NULL"` (#362).

### Other changes
* Added an internal consistency check in `configure_prior()`: if a model's `fixed_parameters` includes a parameter that its `configure_model()` never wires into the formula (neither a dpar nor an nlpar), bmm now fails with a clear model-definition error instead of letting a malformed `b_Intercept ~ constant()` prior reach `brm()`. This is a safety net for model development; it cannot be reached through the normal `bmm()` interface, where an unrecognized parameter is already caught earlier by `check_formula()` (#377).
* `print.bmmodel()` now also displays the required response variables (one per line, annotated with the expected coding — e.g. radians in [-pi, pi] for the circular models, seconds and 0/1 responses for the trial-wise RT models) and the default parameter links, answering "what does my data frame need to look like?" directly at the console (#392).
* Recalibrated the default priors for the **m3** activation parameters (`a`, `c`) so that the `simple` and `softmax` choice rules imply a comparable, broad prior-predictive range of average performance, and so that the `softmax` defaults place equal prior means on general (`a`) and context (`c`) activation — centering the implied `c - a` prior at zero for fair comparisons under cell-means coding. The mis-scaled `normal(0, 2)` effect prior on `c` is replaced by the shared `normal(0, 0.5)`. Because the `simple` rule requires `c > a` to predict accurate recall, its defaults remain asymmetric; direct comparisons of context and general activation should use the `softmax` choice rule (#364).

### Developer-facing changes
* **`use_model_template()` now scaffolds the current model-specification patterns.** It generates a flat `.{model}_defaults` block for unversioned models (like `ddm`) or, with the new `versions` argument, a `.{model}_version_table` block for versioned models (like `cswald`). The generated constructor spells out every field of the model object inline (referencing the defaults/version table for `parameters`, `links`, `fixed_parameters`, `default_priors`, and `init_ranges`), and versioned aliases validate `version` with `match.arg()` (#350).
* **Removed the unused `void_mu` field** from all model definitions and the template. It was assigned but never read anywhere — response-mean suppression is already handled via `fixed_parameters` and the family's `dpars` (#350).

# bmm 1.3.1

### Bug fixes
* Fix `swald_lccdf()` returning incorrect log-survival probability when response time equals non-decision time in the **cswald** model. Previously returned `-Inf` instead of `0` (log of survival = 1) (#348).

### Other changes
* The **ddm** model supports both `cmdstanr` and `rstan` backends. Previously, `cmdstanr` was required.

# bmm 1.3.0

### New models
* Add the **Diffusion Decision Model** (`ddm`) for speeded decision-making tasks with trial-level RT and response data. The model estimates drift rate, boundary separation, non-decision time, and (optionally) relative starting point. Includes distribution functions `dddm()` and `rddm()` (#280).
* Add the **EZ-Diffusion Model** (`ezdm`) for speeded decision-making tasks. The model estimates drift rate, boundary separation, and non-decision time from aggregated summary statistics (mean RT, variance of RT, accuracy) using the closed-form equations derived by Wagenmakers et al. (2007). Supports both 3-parameter (symmetric starting point) and 4-parameter (asymmetric starting point) versions based on Srivastava et al. (2016). Implements Bayesian hierarchical estimation following Chavez & Vandekerckhove (2025). Includes distribution functions `dezdm()` and `rezdm()` (#281).
* Add the **Censored Shifted Wald Model** (`cswald`) for choice reaction time tasks with two response boundaries. The model estimates drift rate, boundary separation, and non-decision time from trial-level RT and response data. Implements two versions: **simple** (treats errors as censored correct responses, appropriate for high-accuracy tasks) and **crisk** (competing risks version with separate accumulators for each response, suitable for balanced accuracy). Includes distribution functions `dcswald()`, `pcswald()`, `qcswald()`, and `rcswald()`. Thanks to @GidonFrischkorn

### New features
* New S3 method **conditional_effects()** for `bmmfit` objects. Provides an intuitive interface for visualizing predictor effects on model parameters, with automatic routing between distributional and non-linear parameters, inverse link transformations to show parameters on their natural scale (`scale = "native"`), softmax handling for mixture3p weight parameters, and filtering of internal model variables (#203).
* New S3 methods for **emmeans** support on `bmmfit` objects. Users can now call `emmeans(fit, ~ condition, dpar = "kappa")` for any bmmodel (#323).
* New **pp_check()** method for multinomial models (e.g., `m3`). Since `brms::pp_check()` does not support the multinomial family, `bmm` now provides a custom method that compares observed and predicted response proportions in the `ppc_bars` style from `bayesplot`. The method supports faceting by experimental conditions via `group`, configurable credible intervals via `probs`, and population-level predictions via `re_formula = NA`. For non-multinomial models, `pp_check()` delegates to `brms::pp_check()` and auto-selects the grouped plot variant when `group` is specified (#324).
* New function **parameters()** lists all parameters of a `bmmodel` or `bmmfit` object with descriptions, link functions, and fixed values (#329).
* New function **extract_stan_blocks()** extracts individual program blocks (functions, data, parameters, etc.) from compiled Stan code (#286).
* New function **extract_parameter_dimensions()** extracts parameter names, dimensions, and types from a Stan parameters block.
* New functions **ezdm_summary_stats()** and **adjust_ezdm_accuracy()** to compute and pre-process summary statistics from trial-level RT data for the EZ-Diffusion Model (#291).
* New functions **flag_contaminant_rts()** and **validate_fast_guesses()** for trial-level contamination detection in RT data. Identifies fast guesses and attention lapses using mixture modeling and provides Bayesian validation of fast guess assumptions (#307).
* New function **create_initfun()** creates initialization functions for models that benefit from or require initial values for MCMC sampling (#285).

### Documentation
* New online [article](https://venpopov.com/bmm/dev/articles/bmm_ddm.html) to accompany the **ddm** model
* New online [article](https://venpopov.com/bmm/dev/articles/bmm_ezdm.html) to accompany the **ezdm** model
* New online [article](https://venpopov.com/bmm/dev/articles/bmm_cswald.html) to accompany the **cswald** model
* New online [article](https://venpopov.com/bmm/dev/articles/bmm_rt_contamination.html) on pre-processing and contamination detection for reaction time data

### Other changes
* Improved **rm3()** random generation function for the M3 model (#279).
* Removed magrittr dependency; replaced `%>%` with the native pipe `|>` (#341).
* Minimum R version is now 4.1.0.

# bmm 1.2.0

### New models
* Add the Memory Measurement Model (Oberauer & Lewandowsky, 2019) and its generalization as the Multinomial Measurement Model for categorical decision tasks as new model class **m3** with three versions: simple span (**ss**), complex span (**cs**), and **custom**. For details, see the [article](https://venpopov.com/bmm/articles/bmm_m3.html) on the `bmm` website (#237). Thanks to @GidonFrischkorn and @chenyu-psy

### New features
* Updates to the `bmf2bf` S3 methods for more flexible translation of `bmmformulas` into `brmsformulas` (#227).
* New function **apply_links** adds link functions to all non-linear formulas in a **bmmformula** object.
* New example data set **oberauer_lewandowsky_2019_e1** for exploring the **m3** model.
* The `file_refit` argument of the `bmm` function now accepts character strings like `brms`. A warning is given when "on_change" is specified, as this is not currently implemented for `bmmodels` (#228).
* New function **rejection_sampling**

### Bug fixes
* Fix conflict in setting default priors when model parameters were transformed in a non-linear formula (#232).
* Allow a NULL formula (`formula(NULL)`) to be added to a bmmformula for consistentcy with brms (#264)
* Improve error messages when attempting to construct bmmformulas without a left-hand-side variable

### Documentation
* Add documentation to the [continuous reproduction task](https://venpopov.com/bmm/articles/bmm_vwm_crt.html) article for pre-processing half-circular stimulus spaces when using `bmmodels` of the `circular` model class (#229, #233).
* New online [article](https://venpopov.com/bmm/articles/bmm_m3.html) to accompany the m3 model

### Other changes
* vectorize `k2sd()` function for improved performance
* various internal refactorings (#246, #242)
* dplyr, magrittr and tidyr dependencies are now optional (#240)
* new contributor - Chenyun Li (chenyu-psy) for his work on the m3 model

# bmm 1.0.0

First version of the package on published on CRAN!

### New features
* you can now specify to save the **bmmfit** object generated by **bmm()** to a file with the **file** argument, similarly to **brms::brm()** (#190)
* the parameterization of the **imm** was adapted to accurately reflect the model as implemented by Oberauer et al. (2017)
* prepare package for CRAN submission

### Bug fixes
* fix incorrect specification of default priors when only an interaction is specified (#201)
* the random generation function for the **mixture3p** and **imm** returned incorrect samples for some rare parameter combinations, this has now been fixed, so that the functions now return correct samples for all parameter combinations.

### Deprecated functions and arguments
* BREAKING CHANGE: the arguments for the distribution functions of the  **mixture2p** and **mixture3p** model have been change to match the snake_case coding scheme. Instead of **pMem** and **pNT** these are now **p_mem** and **p_nt**. The old names are deprecated and are no longer supported

# bmm 0.5.1

### Bug fixes
* fix the display of the model call in the summary method for bmm models

# bmm 0.5.0

### New features
* add a **summary()** method for **bmmfit** objects (#144) 
* add a global option **bmm.summary_backend** to control the backend used for the **summary()** method (choices are *"bmm"* and *"brms"*)
* function **restructure()** now allows to apply methods introduced in newer **bmm** versions to **bmmfit** objects created by older **bmm** versions
* you can now specify any model parameter to be a constant by using an equal sign in the **bmmformula** (#142)
* you can now choose to estimate parameters that are fixed to a constant by default for all models (#145)
* default priors for all models are now specified via the **configure_prior()** S3 method (#145)
* **cmdstanr** will be used as the default backend for **brms** if the user has it installed (#145)
* various updates to the documentation and data sets

### Documentation
* two new online articles that [introduce the **bmmformula** syntax](https://venpopov.com/bmm/articles/bmm_bmmformula.html) and explain [how to extract information from **bmmodels**](https://venpopov.com/bmm/articles/bmm_extract_info.html) such as the generated Stan code and Stan data for each model

### Bug fixes
* fix a bug preventing the **sort_data** check from being executed (#72)
* fix bugs with the **summary()** function not displaying implicit parameters (#152) and not working properly with some hierarchical designs (#173)
* fix a bug in which the **sort_data** check occurred in cases where it shouldn't (#158)

### Deprecated functions and arguments
* BREAKING CHANGE: remove **get_model_prior(), get_stancode() and get_standata()**. Due to [recent changes](https://github.com/paul-buerkner/brms/pull/1604) in *brms* version 2.21.0, you can now use the *brms* functions **default_prior**, **stancode** and **standata** directly with *bmm* models.
* the function **fit_model()** is deprecated in favor of **bmm()** and will be removed in a future version (#163)
* the argument **setsize** for the **mixture3p** and **IMM** models is now called **set_size** for consistency. The old argument name is deprecated and will be removed in a future version (#163)
* the distributions functions for the imm model are renamed from **dIMM**, **pIMM**, **rIMM** and **qIMM** to **dimm**, **pimm**, **rimm** and **qimm** (#163)
* the argument parallel for the **bmm()** function is deprecated and will be removed in a future version. Use **cores** instead, as for **brms::brm()** (#163)
* the models **IMMfull()**, **IMMabc()** and **IMMbsc()** are now called via **imm()**, **imm(version = "abc")** or **imm(version = "bsc")**. The old names are deprecated and will be removed in a future version (#163)
* the **sdmSimple()** model is now called **sdm()**. The old name is deprecated and will be removed in a future version (#163)

### Other changes
* **bmm** now requires the at least version 2.21.0 of **brms**.

# bmm 0.4.0

### New features

* add a check for the **sdmSimple** model if the data is sorted by predictors. This leads to much faster sampling. The user can control the default behavior with the **sort_data** argument (#72)
* the **mixture3p** and **IMM** models now require that the intercept must be suppressed when set size is used as a predictor (#96).
* add postprocessing methods for **sdmSimple** to allow the use of **pp_check()**, **conditional_effects** and **bridgesampling** with the model (#30)
* add informed default priors for all models. You can always use the **get_model_prior()** function to see the default priors for a model
* add a new function **set_default_prior** for developers, which allows them to more easily set default priors on new models regardless of the user-specified formula
* you can now specify variables for models via regular expressions rather than character vectors (#102)
* you can now view and set all **bmm** global options via **bmm_options()**. See **?bmm_options** for more information
* add a start-up message upon loading the package 

### Bug fixes
* fix a bug in the **mixture3p** and **IMM** models which caused an error when intercept was not suppressed and set size was used as predictor
* **update()** now works properly with **bmmfit** objects (#95)
* fix a bug in the **sort_data** check which caused an error when using grouped covariance structure in random effects across different parameters

### Other changes
* **brms** is now loaded automatically when loading **bmm** with **library(bmm)**

# bmm 0.3.0

### New features

* BREAKING CHANGE: The **fit_model** function now requires a **bmmformula** to be passed. The syntax of the **bmmformula** or its short form **bmf** is equal to specifying a **brmsformula**. However, as of this version the **bmmformula** only specifies how parameters of a **bmmodel** change across experimental conditions or continuous predictors. The response variables that the model is fit to now have to be specified when the model is defined using **model = bmmodel()**. (#79)
* BREAKING CHANGE: The **non_target** and **spaPos** variables for the **mixture3p** and **IMM** models were relabeled to **nt_features** and **nt_distances** for consistency. This is also to communicate that distance is not limited to spatial distance but distances on any feature dimensions of the retrieval cues. Currently, still only a single generalization gradient for the cue features is possible. 
* This release includes reference fits for all implemented models to ensure that future changes to the package do not compromise the included models and change the results that their implementations produce.
* The **check_formula** methods have been adapted to match the new **bmmformula** syntax. It now evaluates if formulas have been specified using the **bmmformula** function, if formulas for all parameters of a **bmmodel** have been specified and warns the user that only a fixed intercept will be estimated if no formula for one of the parameters was provided. Additionally, **check_formula** throws an error should formulas be provided that do not match a parameter of the called **bmmodel** unless they are part of a non-linear transformation.
* You can now specify formulas for internally fixed parameters such as **mu** in all visual working memory models. This allows you to predict if there is response bias in the data. If a formula is not provided for **mu**, the model will assume that the mean of the response distribution is fixed to zero.
* there is now an option **bmm.silent** that allows to suppress messages
* the baseline activation **b** was removed from the **IMM** models, as this is internally fixed 
to zero for scaling and as of now cannot be predicted by independent variables because the model would be unidentifiable.
* the arguments used to fit the **bmmodel** are now accessible in the **bmmfit** object via the `fit$bmm$fit_args` list.
* add class('bmmfit') to the object returned from fit_model() allowing for more flexible postprocessing of the underlying **brmsfit** object. The object is now of class('bmmfit', 'brmsfit')
* changes to column names of datasets **zhang_luck_2008** and **oberauer_lin_2017** to make them more consistent

### Bug Fixes
* an error with the treatment of distances in the **IMMfull** and the **IMMbsc** has been corrected. This versions ensures that only positive distances can be passed to any of the two models.
* removed a warning regarding the scaling of the distances in the **IMMfull** and the **IMMbsc** that was specific only for circular distances.

### Documentation
* All articles have been update to the new **bmmformula** syntax.

# bmm 0.2.2

### Bug Fixes
* fixed a bug where passing a character vector or negative values to set_size argument of visual working memory models caused an error or incorrect behavior (#97)

# bmm 0.2.1

### Bug Fixes
* Minor change to sdmSimple Stan helper functions to avoid a harmless warning message in the Stan output


# bmm 0.2.0

### New features

* New model available - The Signal Discrimination Model by Oberauer (2023) for visual working memory continuous reproduction tasks. See ?sdmSimple. The current version does not take into account non-target activation
* Add ability to extract information about the default priors in **bmm** models with **get_model_prior()** (#53)
* Add ability to generate stan code and stan data for each model with **get_model_stancode()** and **get_model_standata()** (#81)
* BREAKING CHANGE: Add distribution functions for likelihood (e.g. **dIMM()**) and random variate generation **rIMM()**) for all models in the package. Remove deprecated **gen_3p_data()** and **gen_imm_data()** functions (#69)
* Two new data sets are available: **zhang_luck_2008** and **oberauer_lin_2017** (#22)

### Documentation

* Website for the development version of the package is now available at https://venpopov.com/bmm/dev/ (#18)
* Add articles for each model to the website at https://venpopov.com/bmm/dev/articles/
* Add a detailed developer's guide to the website at https://venpopov.com/bmm/dev/dev-notes (#21)
* Improve README with more detailed information about the package's goals and its models (#21)

### Other changes

* Save **bmm** package version in the **brmsfit** object for reproducibility - e.g. `fit$version$bmm` (#88)


# bmm 0.1.1

### New features

* BREAKING CHANGE: Improve user interface to fit_model() ensures package stability and future development. Model specific arguments are now passed to the model functions as named arguments (e.g. **mixture3p(non_targets, setsize)**). This allows for a more flexible and intuitive way to specify model arguments. Passing model specific arguments directly to the **fit_model()** function is now deprecated (#43).
* Add information about each model such as domain, task, name, version, citation, requirements and parameters (#42)
* Add ability to generate a template file for adding new models to the package with **use_model_template()** (for developers) (#39)

### Other changes

* Improve documentation of model functions. You can now get help on each model by typing **?model_name** into your console. For example, calling the information on the full version of the Interference Measurement Model would look like this: **?IMMfull**


# bmm 0.1.0

A major restructuring of the package to support stable and generalizable development of future models (#41). 

### New Features

* Refactor the **fit_model()** function to be generic and independent of the model being fit (#20)
* Transform models to be S3 objects. (#41). 
* View currently supported models with new function **supported_models()**. Currently supported models are: **mixture2p()**, **mixture3p()**, **IMMabc()**, **IMMbsc()**, **IMMfull()** 
* Add S3 methods for checking the data, formula, model and priors (#41)
* Add distribution functions for the Signal Discrimination Model. See **?SDM** for usage (#27)
* Add **softmax** and **softmaxinv** functions 

### Bug Fixes

* Change default prior on log(kappa) to Normal(2,1) for the **mixture3p()** model (#15)

### Other changes

* BREAKING CHANGE: deprecate **model_type** argument in **fit_model()**. Models must now be specified with S3 functions passed to argument **model** rather than model names as strings passed to argument **model_type** (#41)
* Add extensive unit testing


# bmm 0.0.1

* Initial release version
