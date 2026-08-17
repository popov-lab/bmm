############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.imm_parameters <- list(
  mu = glue(
    "Location of the von Mises distribution of memory responses (in radians). \\
    Fixed internally to 0 by default."
  ),
  kappa = "Concentration of the von Mises distribution of memory responses",
  a = "General activation of memory items",
  c = "Context activation",
  s = "Spatial similarity gradient",
  b = glue(
    "Background activation, the reference against which the activations of \\
    the target and the non-targets are weighed. Fixed internally to 1; \\
    freeing it leaves the remaining activations unidentified, since only \\
    their ratio to the background enters the likelihood."
  )
)

.imm_priors <- list(
  mu = list(main = "student_t(1, 0, 1)"),
  kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
  a = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
  c = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
  s = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
  b = list(main = "normal(0, 1)")
)

.imm_links <- list(
  mu = "tan_half", kappa = "log", a = "log", c = "log", s = "log", b = "log"
)

.imm_init_ranges <- list(
  mu = c(-0.1, 0.1), kappa = c(3, 8), a = c(0.2, 1), c = c(1, 4),
  s = c(0.5, 2), b = c(0.8, 1.2)
)

.imm_version_spec <- function(version) {
  kept <- switch(version,
    full = c("mu", "kappa", "a", "c", "s", "b"),
    bsc = c("mu", "kappa", "c", "s", "b"),
    abc = c("mu", "kappa", "a", "c", "b")
  )
  list(
    weight_parameters = switch(version,
      full = c("c", "a", "s", "b"),
      bsc = c("c", "s", "b"),
      abc = c("c", "a", "b")
    ),
    needs_distances = version %in% c("full", "bsc"),
    parameters = .imm_parameters[kept],
    links = .imm_links[kept],
    priors = .imm_priors[kept],
    init_ranges = .imm_init_ranges[kept]
  )
}

.model_imm <- function(resp_error = NULL, nt_features = NULL, nt_distances = NULL,
                       set_size = NULL, regex = FALSE, version = "full", links = NULL,
                       variable_precision = FALSE, vp_nodes = 41L,
                       call = NULL, ...) {
  spec <- .imm_version_spec(version)
  if (variable_precision) {
    spec <- .circmix_add_variable_precision(spec)
  }

  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(nt_features, nt_distances, set_size),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Interference measurement model by Oberauer and Lin (2017).",
      version = version,
      citation = glue(
        "Oberauer, K., & Lin, H.Y. (2017). An interference model \\
          of visual working memory. Psychological Review, 124(1), 21-59"
      ),
      requirements = glue(
        "- The response vairable should be in radians and \\
          represent the angular error relative to the target
          - The non-target features should be in radians and be \\
          centered relative to the target"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = list(mu = 0, b = 0),
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      deprecated_parameters = list(mu1 = "mu"),
      variable_precision = variable_precision,
      vp_nodes = as.integer(vp_nodes)
    ),
    # attributes
    regex = regex,
    regex_vars = if (version == "abc") "nt_features" else c("nt_features", "nt_distances"),
    class = c("bmmodel", "circular", "non_targets", "imm", paste0("imm_", version)),
    call = call
  )

  out$links[names(links)] <- links
  out
}

# user facing alias

#' @title `r .model_imm()$name`
#' @description Three versions of the `r .model_imm()$name` - the full, bsc, and abc.
#' `IMMfull()`, `IMMbsc()`, and `IMMabc()` are deprecated and will be removed in the future.
#' Please use `imm(version = 'full')`, `imm(version = 'bsc')`, or `imm(version = 'abc')` instead.
#'
#' @name imm
#' @details `r model_info(.model_imm(), components =c('domain', 'task', 'name', 'citation'))`
#' #### Version: `full`
#' `r model_info(.model_imm(version = "full"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `bsc`
#' `r model_info(.model_imm(version = "bsc"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `abc`
#' `r model_info(.model_imm(version = "abc"), components =c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#'
#' Additionally, all imm models have an internal parameter that is fixed to 0 to
#' allow the model to be identifiable. This parameter is not estimated and is not
#' included in the model formula. The parameter is:
#'
#'   - b = "Background activation (internally fixed to 0)"
#'
#' @param resp_error The name of the variable in the provided dataset containing
#'   the response error. The response Error should code the response relative to
#'   the to-be-recalled target in radians. You can transform the response error
#'   in degrees to radian using the `deg2rad` function.
#' @param nt_features A character vector with the names of the non-target
#'   variables. The non_target variables should be in radians and be centered
#'   relative to the target. Alternatively, if regex=TRUE, a regular
#'   expression can be used to match the non-target feature columns in the
#'   dataset.
#' @param nt_distances A vector of names of the columns containing the distances
#'   of non-target items to the target item. Alternatively, if regex=TRUE, a regular
#'   expression can be used to match the non-target distances columns in the
#'   dataset. Only necessary for the `bsc` and `full` versions.
#' @param set_size Name of the column containing the set size variable (if
#'   set_size varies) or a numeric value for the set_size, if the set_size is
#'   fixed.
#' @param regex Logical. If TRUE, the `nt_features` and `nt_distances` arguments
#'   are interpreted as a regular expression to match the non-target feature
#'   columns in the dataset.
#' @param version Character. The version of the IMM model to use. Can be one of
#'  `full`, `bsc`, or `abc`. The default is `full`.
#' @param variable_precision Logical; if `TRUE`, the precision of memory varies
#'   from trial to trial (van den Berg et al., 2012). This adds the parameter
#'   `tau`; see [mixture2p()] for the parameterisation. The capacity models of
#'   Zhang and Luck (2008) do not apply to the IMM, which has no discrete slots
#'   -- its bottleneck is the ratio of context to general activation -- so
#'   variable precision is the only additional version offered here.
#' @param vp_nodes Number of quadrature nodes used when
#'   `variable_precision = TRUE`; must be an odd number of at least 41. See
#'   [mixture2p()].
#' @param links A named list of link functions for the model parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # load data
#' data <- oberauer_lin_2017
#'
#' # define formula
#' ff <- bmmformula(
#'   kappa ~ 0 + set_size,
#'   c ~ 0 + set_size,
#'   a ~ 0 + set_size,
#'   s ~ 0 + set_size
#' )
#'
#' # specify the full IMM model with explicit column names for non-target features and distances
#' # by default this fits the full version of the model
#' model1 <- imm(
#'   resp_error = "dev_rad",
#'   nt_features = paste0("col_nt", 1:7),
#'   nt_distances = paste0("dist_nt", 1:7),
#'   set_size = "set_size"
#' )
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = data,
#'   model = model1,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # alternatively specify the IMM model with a regular expression to match non-target features
#' # this is equivalent to the previous call, but more concise
#' model2 <- imm(
#'   resp_error = "dev_rad",
#'   nt_features = "col_nt",
#'   nt_distances = "dist_nt",
#'   set_size = "set_size",
#'   regex = TRUE
#' )
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = data,
#'   model = model2,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # you can also specify the `bsc` or `abc` versions of the model to fit a reduced version
#' model3 <- imm(
#'   resp_error = "dev_rad",
#'   nt_features = "col_nt",
#'   set_size = "set_size",
#'   regex = TRUE,
#'   version = "abc"
#' )
#' fit <- bmm(
#'   formula = ff,
#'   data = data,
#'   model = model3,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#' @export
imm <- function(resp_error, nt_features, nt_distances, set_size, regex = FALSE,
                version = "full", variable_precision = FALSE, vp_nodes = 41L,
                links = NULL, ...) {
  call <- match.call()
  dots <- list(...)
  if ("setsize" %in% names(dots)) {
    set_size <- dots$setsize
    warning2("The argument 'setsize' is deprecated. Please use 'set_size' instead.")
  }
  if (version == "abc") nt_distances <- NULL
  stop_missing_args()
  version <- match.arg(version, c("full", "bsc", "abc"))
  .circmix_check_variable_precision(variable_precision, vp_nodes)

  .model_imm(
    resp_error = resp_error, nt_features = nt_features,
    nt_distances = nt_distances, set_size = set_size, regex = regex,
    version = version, variable_precision = variable_precision,
    vp_nodes = vp_nodes, links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_DATA S3 methods                                                  ####
############################################################################# !
# A check_data.* function should be defined for each class of the model.
# If a model shares methods with other models, the shared methods should be
# defined in data-helpers.R. Put here only the methods that are specific to
# the model. See ?check_data for details

#' @export
check_data.imm_bsc <- function(model, data, formula) {
  data <- .check_data_imm_dist(model, data, formula)
  NextMethod("check_data")
}

#' @export
check_data.imm_full <- function(model, data, formula) {
  data <- .check_data_imm_dist(model, data, formula)
  NextMethod("check_data")
}

.check_data_imm_dist <- function(model, data, formula) {
  nt_distances <- model$other_vars$nt_distances
  max_set_size <- attr(data, "max_set_size")

  stopif(
    !isTRUE(all.equal(length(nt_distances), max_set_size - 1)),
    "The number of columns for non-target distances in the argument \\
    'nt_distances' should equal max(set_size)-1})"
  )

  # replace NA values with 999 so they have 0 effect through the distance formula
  data[, nt_distances][is.na(data[, nt_distances])] <- 999

  stopif(
    any(data[, nt_distances] < 0),
    "All non-target distances to the target need to be postive."
  )

  data
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
bmf2bf.imm <- function(model, formula) {
  spec <- .imm_version_spec(model$version)
  covariates <- model$other_vars$nt_features
  if (spec$needs_distances) {
    covariates <- c(covariates, model$other_vars$nt_distances)
  }
  brms::bf(.circmix_aterm(
    model$resp_vars$resp_error,
    vint = "ss_numeric", vreal = covariates
  ))
}

#' @export
configure_model.imm <- function(model, data, formula) {
  spec <- .imm_version_spec(model$version)
  n_nt <- attr(data, "max_set_size") - 1
  groups <- if (spec$needs_distances) list(nt = n_nt, dist = n_nt) else list(nt = n_nt)

  formula <- bmf2bf(model, formula)
  formula$family <- .circmix_custom_family(
    model,
    family = paste0("imm_", model$version),
    weight_parameters = spec$weight_parameters,
    vint = TRUE, n_vreal = sum(unlist(groups)),
    log_lik = .imm_log_lik(model$version),
    posterior_predict = .imm_posterior_predict(model$version)
  )

  nlist(
    formula, data,
    stanvars = .circmix_model_stanvars(model, formula$family, "imm_funs.stan",
      vint = "ss", vreal = groups
    )
  )
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

.imm_log_lik <- function(version) {
  switch(version,
    full = log_lik_imm_full, bsc = log_lik_imm_bsc, abc = log_lik_imm_abc
  )
}

.imm_posterior_predict <- function(version) {
  switch(version,
    full = posterior_predict_imm_full, bsc = posterior_predict_imm_bsc,
    abc = posterior_predict_imm_abc
  )
}

# the non-target features come first in the vreal block, the distances after
.imm_covariates <- function(prep, i, with_distances) {
  covariates <- .circmix_prep_nt(prep, i)
  if (!with_distances) {
    return(list(nt = covariates, dist = NULL))
  }
  half <- length(covariates) / 2
  list(nt = covariates[seq_len(half)], dist = covariates[half + seq_len(half)])
}

# the density takes the response and the node count on top of what the
# sampler needs
.imm_args <- function(i, prep, version, density = TRUE) {
  covariates <- .imm_covariates(prep, i, version %in% c("full", "bsc"))
  c(
    if (density) list(x = prep$data$Y[i]),
    list(
      mu = brms::get_dpar(prep, "mu", i = i),
      kappa = brms::get_dpar(prep, "kappa", i = i),
      c = brms::get_dpar(prep, "c", i = i),
      a = if (version != "bsc") brms::get_dpar(prep, "a", i = i),
      s = if (version != "abc") brms::get_dpar(prep, "s", i = i),
      b = brms::get_dpar(prep, "b", i = i),
      set_size = prep$data$vint1[i], nt = covariates$nt, dist = covariates$dist,
      tau = .circmix_prep_tau(prep, i), version = version
    ),
    if (density) list(nodes = .circmix_prep_nodes(prep))
  )
}

log_lik_imm_full <- function(i, prep) {
  brms::do_call(.dimm_version, .imm_args(i, prep, "full"))
}

log_lik_imm_bsc <- function(i, prep) {
  brms::do_call(.dimm_version, .imm_args(i, prep, "bsc"))
}

log_lik_imm_abc <- function(i, prep) {
  brms::do_call(.dimm_version, .imm_args(i, prep, "abc"))
}

posterior_predict_imm_full <- function(i, prep, ...) {
  brms::do_call(.rimm_version, .imm_args(i, prep, "full", density = FALSE))
}

posterior_predict_imm_bsc <- function(i, prep, ...) {
  brms::do_call(.rimm_version, .imm_args(i, prep, "bsc", density = FALSE))
}

posterior_predict_imm_abc <- function(i, prep, ...) {
  brms::do_call(.rimm_version, .imm_args(i, prep, "abc", density = FALSE))
}

# ---- deprecated calls for specific versions ----

#' @rdname imm
#' @keywords deprecated
#' @export
IMMfull <- function(resp_error, nt_features, nt_distances, set_size, regex = FALSE, ...) {
  call <- match.call()
  dots <- list(...)
  warning2("The function `IMMfull()` is deprecated. Please use `imm(version = 'full')` instead.")
  if ("setsize" %in% names(dots)) {
    set_size <- dots$setsize
    warning2("The argument 'setsize' is deprecated. Please use 'set_size' instead.")
  }
  stop_missing_args()
  .model_imm(
    resp_error = resp_error, nt_features = nt_features,
    nt_distances = nt_distances, set_size = set_size, regex = regex,
    version = "full", call = call, ...
  )
}

#' @rdname imm
#' @keywords deprecated
#' @export
IMMbsc <- function(resp_error, nt_features, nt_distances, set_size, regex = FALSE, ...) {
  call <- match.call()
  dots <- list(...)
  warning2("The function `IMMbsc()` is deprecated. Please use `imm(version = 'bsc')` instead.")
  if ("setsize" %in% names(dots)) {
    set_size <- dots$setsize
    warning2("The argument 'setsize' is deprecated. Please use 'set_size' instead.")
  }
  stop_missing_args()
  .model_imm(
    resp_error = resp_error, nt_features = nt_features,
    nt_distances = nt_distances, set_size = set_size, regex = regex,
    version = "bsc", call = call, ...
  )
}

#' @rdname imm
#' @keywords deprecated
#' @export
IMMabc <- function(resp_error, nt_features, set_size, regex = FALSE, ...) {
  call <- match.call()
  dots <- list(...)
  warning2("The function `IMMabc()` is deprecated. Please use `imm(version = 'abc')` instead.")
  if ("setsize" %in% names(dots)) {
    set_size <- dots$setsize
    warning2("The argument 'setsize' is deprecated. Please use 'set_size' instead.")
  }
  stop_missing_args()
  .model_imm(
    resp_error = resp_error, nt_features = nt_features, set_size = set_size,
    regex = regex, version = "abc", call = call, ...
  )
}
