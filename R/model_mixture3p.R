############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.mixture3p_location <- glue(
  "Location of the von Mises distribution of memory responses (in radians). \\
  Fixed internally to 0 by default."
)

.mixture3p_version_table <- list(
  simple = list(
    weight_parameters = c("thetat", "thetant"),
    swap_parameter = "thetant",
    parameters = list(
      mu = .mixture3p_location,
      kappa = "Concentration of the von Mises distribution of memory responses",
      thetat = glue(
        "Mixture weight for target responses. The weights of the target, the \\
        non-targets and guessing are normalised across the components present \\
        on a trial, with guessing as the reference, so thetat maps to the \\
        probability of a target response only jointly with thetant."
      ),
      thetant = "Mixture weight for non-target responses"
    ),
    links = list(
      mu = "tan_half", kappa = "log", thetat = "identity", thetant = "identity"
    ),
    priors = list(
      mu = list(main = "student_t(1, 0, 1)"),
      kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
      thetat = list(main = "logistic(0, 1)"),
      thetant = list(main = "logistic(0, 1)")
    ),
    init_ranges = list(
      mu = c(-0.1, 0.1), kappa = c(3, 8), thetat = c(1, 2), thetant = c(-1, 0)
    )
  ),
  slot = list(
    weight_parameters = c("K", "pnt"),
    swap_parameter = "pnt",
    parameters = list(
      mu = .mixture3p_location,
      kappa = "Concentration of the von Mises distribution of memory responses",
      K = glue(
        "Capacity: the number of items that can be held in memory. An item is \\
        either held at full precision or not held at all, so the probability \\
        that a response comes from memory is min(1, K / set_size)."
      ),
      pnt = glue(
        "Probability that a response reports a non-target rather than the \\
        target, given that it comes from memory."
      )
    ),
    links = list(mu = "tan_half", kappa = "log", K = "log", pnt = "logit"),
    priors = list(
      mu = list(main = "student_t(1, 0, 1)"),
      kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
      K = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)"),
      pnt = list(main = "logistic(-1, 1)", effects = "normal(0, 0.5)")
    ),
    init_ranges = list(
      mu = c(-0.1, 0.1), kappa = c(3, 8), K = c(2, 4), pnt = c(0.05, 0.2)
    )
  ),
  slot_averaging = list(
    weight_parameters = c("K", "pnt"),
    swap_parameter = "pnt",
    parameters = list(
      mu = .mixture3p_location,
      kappa = glue(
        "Concentration of the von Mises distribution for an item held in a \\
        single slot. An item holding several slots is more precise, because \\
        averaging independent samples adds their Fisher information."
      ),
      K = glue(
        "Capacity: the number of memory slots distributed over the items in \\
        the array. An item receives floor(K / set_size) or one more slot, and \\
        an item that receives none is guessed."
      ),
      pnt = glue(
        "Probability that a response reports a non-target rather than the \\
        target, given that it comes from memory."
      )
    ),
    links = list(mu = "tan_half", kappa = "log", K = "log", pnt = "logit"),
    priors = list(
      mu = list(main = "student_t(1, 0, 1)"),
      kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
      K = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)"),
      pnt = list(main = "logistic(-1, 1)", effects = "normal(0, 0.5)")
    ),
    init_ranges = list(
      mu = c(-0.1, 0.1), kappa = c(3, 8), K = c(2, 4), pnt = c(0.05, 0.2)
    )
  )
)

.model_mixture3p <- function(resp_error = NULL, nt_features = NULL, set_size = NULL,
                             regex = FALSE, links = NULL, version = "simple",
                             variable_precision = FALSE, vp_nodes = 41L,
                             call = NULL, ...) {
  spec <- .mixture3p_version_table[[version]]
  if (variable_precision) {
    spec <- .circmix_add_variable_precision(spec)
  }

  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(nt_features, set_size),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Three-parameter mixture model by Bays et al (2009).",
      version = version,
      citation = glue(
        "Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). \\
        The precision of visual working memory is set by allocation \\
        of a shared resource. Journal of Vision, 9(10), 1-11"
      ),
      requirements = glue(
        "- The response vairable should be in radians and \\
        represent the angular error relative to the target
        - The non-target features should be in radians and be \\
        centered relative to the target"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = list(mu = 0),
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      deprecated_parameters = list(mu1 = "mu"),
      variable_precision = variable_precision,
      vp_nodes = as.integer(vp_nodes)
    ),
    # attributes
    regex = regex,
    regex_vars = c("nt_features"),
    class = c(
      "bmmodel", "circular", "non_targets", "mixture3p",
      paste0("mixture3p_", version)
    ),
    call = call
  )
  out$links[names(links)] <- links
  out
}


# user facing alias
#' @title `r .model_mixture3p()$name`
#' @details `r model_info(.model_mixture3p())`
#' @param resp_error The name of the variable in the dataset containing
#'   the response error. The response error should code the response relative to
#'   the to-be-recalled target in radians. You can transform the response error
#'   in degrees to radians using the `deg2rad` function.
#' @param nt_features A character vector with the names of the non-target
#'   feature values. The non_target feature values should be in radians and
#'   centered relative to the target. Alternatively, if regex=TRUE, a regular
#'   expression can be used to match the non-target feature columns in the
#'   dataset.
#' @param set_size Name of the column containing the set size variable (if
#'   set_size varies) or a numeric value for the set_size, if the set_size is
#'   fixed.
#' @param regex Logical. If TRUE, the `nt_features` argument is interpreted as
#'  a regular expression to match the non-target feature columns in the dataset.
#' @param version A character string specifying which version of the model to
#'   use. Options are:
#'   \itemize{
#'     \item `"simple"` (default): the standard three-parameter mixture model.
#'       The weights of the target (`thetat`), the non-targets (`thetant`) and
#'       guessing are normalised across the components a trial actually has,
#'       with guessing as the reference.
#'     \item `"slot"`: the probability that a response comes from memory follows
#'       the fixed-resolution slot model of Zhang and Luck (2008),
#'       `min(1, K / set_size)`, and `pnt` is the probability that a remembered
#'       response reports a non-target.
#'     \item `"slot_averaging"`: as `"slot"`, but `K` slots are distributed over
#'       the items and precision grows with the number of slots an item holds.
#'   }
#'   Storage and swapping are separate in the capacity versions: the capacity
#'   rule says how often a response comes from memory and `pnt` says which
#'   stored item is reported. In `"simple"` the two are jointly constrained by
#'   the shared normalisation, so neither weight is separately interpretable.
#' @param variable_precision Logical; if `TRUE`, the precision of memory varies
#'   from trial to trial (van den Berg et al., 2012). See [mixture2p()] for
#'   details of the parameterisation.
#' @param vp_nodes Number of quadrature nodes used when
#'   `variable_precision = TRUE`; must be an odd number of at least 41. See
#'   [mixture2p()].
#' @param links A named list of link functions for the model parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @references Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The
#'   precision of visual working memory is set by allocation of a shared
#'   resource. Journal of Vision, 9(10), 1-11.
#'
#'   Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution representations
#'   in visual working memory. Nature, 453(7192), 233-235.
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # generate artificial data from the Bays et al (2009) 3-parameter mixture model
#' dat <- data.frame(
#'   y = rmixture3p(n = 2000, mu = c(0, 1, -1.5, 2)),
#'   nt1_loc = 1,
#'   nt2_loc = -1.5,
#'   nt3_loc = 2
#' )
#'
#' # define formula
#' ff <- bmmformula(
#'   kappa ~ 1,
#'   thetat ~ 1,
#'   thetant ~ 1
#' )
#'
#' # specify the 3-parameter model with explicit column names for non-target features
#' model1 <- mixture3p(resp_error = "y", nt_features = paste0("nt", 1:3, "_loc"), set_size = 4)
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = model1,
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
#'
#' # alternatively specify the 3-parameter model with a regular expression to match non-target features
#' # this is equivalent to the previous call, but more concise
#' model2 <- mixture3p(resp_error = "y", nt_features = "nt.*_loc", set_size = 4, regex = TRUE)
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = model2,
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
mixture3p <- function(resp_error, nt_features, set_size, regex = FALSE,
                      version = c("simple", "slot", "slot_averaging"),
                      variable_precision = FALSE, vp_nodes = 41L, links = NULL,
                      ...) {
  call <- match.call()
  dots <- list(...)
  if ("setsize" %in% names(dots)) {
    set_size <- dots$setsize
    warning2("The argument 'setsize' is deprecated. Please use 'set_size' instead.")
  }
  stop_missing_args()
  version <- match.arg(version)
  .circmix_check_variable_precision(variable_precision, vp_nodes)

  .model_mixture3p(
    resp_error = resp_error, nt_features = nt_features, set_size = set_size,
    regex = regex, version = version, variable_precision = variable_precision,
    vp_nodes = vp_nodes, links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_FORMULA METHODS                                                  ####
############################################################################# !

#' @export
check_formula.mixture3p <- function(model, data, formula) {
  swap <- .mixture3p_version_table[[model$version]]$swap_parameter
  set_size_var <- model$other_vars$set_size

  # A trial with one item has no non-target to swap to, so that level of the
  # swap parameter is absent from the likelihood and samples its prior. The
  # previous mixture implementation hid this by pinning the level to a constant,
  # which it had to do anyway to switch off the sentinel weight.
  predicted_by_set_size <- set_size_var %in% rhs_vars(formula[[swap]])
  warnif(
    predicted_by_set_size && isTRUE(any(data[["ss_numeric"]] == 1)),
    "Your data contain trials with a set size of 1, where no swap can occur, \\
    so the set-size-1 level of '{swap}' is not identified and will sample its \\
    prior. Drop those trials from the formula for '{swap}', or read that level \\
    as prior-only."
  )

  NextMethod("check_formula")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

#' @export
bmf2bf.mixture3p <- function(model, formula) {
  brms::bf(.circmix_aterm(
    model$resp_vars$resp_error,
    vint = "ss_numeric", vreal = model$other_vars$nt_features
  ))
}

#' @export
configure_model.mixture3p <- function(model, data, formula) {
  spec <- .mixture3p_version_table[[model$version]]
  n_nt <- attr(data, "max_set_size") - 1

  formula <- bmf2bf(model, formula)
  formula$family <- .circmix_custom_family(
    model,
    family = paste0("mixture3p_", model$version),
    weight_parameters = spec$weight_parameters,
    vint = TRUE, n_vreal = n_nt,
    log_lik = .mixture3p_log_lik(model$version),
    posterior_predict = .mixture3p_posterior_predict(model$version)
  )

  nlist(
    formula, data,
    stanvars = .circmix_model_stanvars(model, formula$family, "mixture3p_funs.stan",
      vint = "ss", vreal = list(nt = n_nt)
    )
  )
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

.mixture3p_log_lik <- function(version) {
  switch(version,
    simple = log_lik_mixture3p_simple,
    slot = log_lik_mixture3p_slot,
    slot_averaging = log_lik_mixture3p_slot_averaging
  )
}

.mixture3p_posterior_predict <- function(version) {
  switch(version,
    simple = posterior_predict_mixture3p_simple,
    slot = posterior_predict_mixture3p_slot,
    slot_averaging = posterior_predict_mixture3p_slot_averaging
  )
}

log_lik_mixture3p_simple <- function(i, prep) {
  .dmixture3p_simple(
    prep$data$Y[i],
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    thetat = brms::get_dpar(prep, "thetat", i = i),
    thetant = brms::get_dpar(prep, "thetant", i = i),
    set_size = prep$data$vint1[i], nt = .circmix_prep_nt(prep, i),
    tau = .circmix_prep_tau(prep, i), nodes = .circmix_prep_nodes(prep)
  )
}

log_lik_mixture3p_slot <- function(i, prep) {
  .dmixture3p_slot(
    prep$data$Y[i],
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    K = brms::get_dpar(prep, "K", i = i),
    p_nt = brms::get_dpar(prep, "pnt", i = i),
    set_size = prep$data$vint1[i], nt = .circmix_prep_nt(prep, i),
    tau = .circmix_prep_tau(prep, i), nodes = .circmix_prep_nodes(prep)
  )
}

log_lik_mixture3p_slot_averaging <- function(i, prep) {
  .dmixture3p_slot_averaging(
    prep$data$Y[i],
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    K = brms::get_dpar(prep, "K", i = i),
    p_nt = brms::get_dpar(prep, "pnt", i = i),
    set_size = prep$data$vint1[i], nt = .circmix_prep_nt(prep, i),
    tau = .circmix_prep_tau(prep, i), nodes = .circmix_prep_nodes(prep)
  )
}

posterior_predict_mixture3p_simple <- function(i, prep, ...) {
  weights <- .mixture3p_softmax_weights(
    brms::get_dpar(prep, "thetat", i = i), brms::get_dpar(prep, "thetant", i = i),
    prep$data$vint1[i]
  )
  .rmixture3p(
    brms::get_dpar(prep, "mu", i = i), .circmix_prep_nt(prep, i),
    prep$data$vint1[i], weights,
    brms::get_dpar(prep, "kappa", i = i), .circmix_prep_tau(prep, i)
  )
}

posterior_predict_mixture3p_slot <- function(i, prep, ...) {
  set_size <- prep$data$vint1[i]
  weights <- .mixture3p_nested_weights(
    pmin(1, brms::get_dpar(prep, "K", i = i) / set_size),
    brms::get_dpar(prep, "pnt", i = i), set_size
  )
  .rmixture3p(
    brms::get_dpar(prep, "mu", i = i), .circmix_prep_nt(prep, i), set_size,
    weights, brms::get_dpar(prep, "kappa", i = i), .circmix_prep_tau(prep, i)
  )
}

posterior_predict_mixture3p_slot_averaging <- function(i, prep, ...) {
  .rmixture3p_slot_averaging(
    mu = brms::get_dpar(prep, "mu", i = i),
    nt = .circmix_prep_nt(prep, i), set_size = prep$data$vint1[i],
    kappa = brms::get_dpar(prep, "kappa", i = i),
    K = brms::get_dpar(prep, "K", i = i),
    p_nt = brms::get_dpar(prep, "pnt", i = i),
    tau = .circmix_prep_tau(prep, i)
  )
}
