############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.mixture3p_location <- glue(
  "Location of the von Mises distribution of memory responses (in radians). \\
  Fixed internally to 0 by default."
)

.mixture3p_citation <- paste(
  "Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The precision of",
  "visual working memory is set by allocation of a shared resource. Journal",
  "of Vision, 9(10), Article 7. https://doi.org/10.1167/9.10.7"
)

.mixture3p_capacity_citation <- paste(
  "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution representations",
  "in visual working memory. Nature, 453(7192), 233-235.",
  "https://doi.org/10.1038/nature06860"
)

.mixture3p_location_prior <- list(
  main = "normal(0, 0.5)", effects = "normal(0, 0.25)", sd = "exponential(4)"
)

.mixture3p_kappa_prior <- list(
  main = "normal(2, 1)", effects = "normal(0, 1)", sd = "exponential(1)"
)

.mixture3p_capacity_priors <- list(
  mu = .mixture3p_location_prior,
  kappa = .mixture3p_kappa_prior,
  K = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"),
  ptarget = list(main = "logistic(1, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)")
)

# central 50% of the main default prior on the native scale
.mixture3p_capacity_inits <- list(
  mu = c(-0.65, 0.65), kappa = c(3.8, 15), K = c(1.9, 3.8), ptarget = c(0.48, 0.89)
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
      thetant = glue(
        "Mixture weight for non-target responses, summed over the non-targets, \\
        which share it equally. Like thetat, it is weighed against guessing."
      )
    ),
    links = list(
      mu = "tan_half", kappa = "log", thetat = "softmax", thetant = "softmax"
    ),
    priors = list(
      mu = .mixture3p_location_prior,
      kappa = .mixture3p_kappa_prior,
      # two cells with logistic(0, 1) priors differ by SD pi * sqrt(2 / 3) = 2.57;
      # narrower effects priors shrink set-size effects from a reference near
      # ceiling, so 1 + set_size and 0 + set_size disagree (#466)
      thetat = list(main = "logistic(0, 1)", effects = "normal(0, 2.5)", sd = "exponential(1)"),
      thetant = list(main = "logistic(0, 1)", effects = "normal(0, 2.5)", sd = "exponential(1)")
    ),
    # the softmax weights have no native value of their own, so their ranges
    # are on the sampling scale
    init_ranges = list(
      mu = c(-0.65, 0.65), kappa = c(3.8, 15),
      thetat = c(-1.1, 1.1), thetant = c(-1.1, 1.1)
    )
  ),
  slot = list(
    weight_parameters = c("K", "ptarget"),
    swap_parameter = "ptarget",
    parameters = list(
      mu = .mixture3p_location,
      kappa = "Concentration of the von Mises distribution of memory responses",
      K = glue(
        "Capacity: the number of items that can be held in memory. An item is \\
        either held at full precision or not held at all, so the probability \\
        that a response comes from memory is min(1, K / set_size)."
      ),
      ptarget = glue(
        "Probability that a response reports the target rather than a \\
        non-target, given that it comes from memory."
      )
    ),
    links = list(mu = "tan_half", kappa = "log", K = "log", ptarget = "logit"),
    priors = .mixture3p_capacity_priors,
    init_ranges = .mixture3p_capacity_inits
  ),
  slot_averaging = list(
    weight_parameters = c("K", "ptarget"),
    swap_parameter = "ptarget",
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
      ptarget = glue(
        "Probability that a response reports the target rather than a \\
        non-target, given that it comes from memory."
      )
    ),
    links = list(mu = "tan_half", kappa = "log", K = "log", ptarget = "logit"),
    priors = .mixture3p_capacity_priors,
    init_ranges = .mixture3p_capacity_inits
  )
)

.model_mixture3p <- function(resp_error = NULL, nt_features = NULL, set_size = NULL,
                             regex = FALSE, links = NULL, version = "simple",
                             variable_precision = FALSE, vp_nodes = 81L,
                             call = NULL, ...) {
  spec <- .mixture3p_version_table[[version]]
  if (variable_precision) {
    spec <- .circmix_add_variable_precision(spec, per_slot = version == "slot_averaging")
  }

  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(nt_features, set_size),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Three-parameter mixture model by Bays et al (2009).",
      citation = c(
        .mixture3p_citation,
        if (version != "simple") .mixture3p_capacity_citation,
        if (variable_precision) .circmix_vp_citation
      ),
      version = version,
      # paste0() rather than glue(), which would trim the indent that keeps
      # the requirements in one list
      requirements = paste0(
        "- The response variable should be in radians and represent the ",
        "angular error relative to the target\n",
        "  - The non-target features should be in radians and be centered ",
        "relative to the target\n",
        "  - In the 'slot' and 'slot_averaging' versions, the set size should ",
        "vary for the capacity rule to be testable"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = list(mu = 0),
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      deprecated_parameters = list(mu1 = "mu"),
      dpar_since = "1.5.0",
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
  out <- set_links(out, links)
  out
}

# the softmax is computed inside the likelihood, so a softmax weight cannot take
# another link; every other link reaches the custom family
#' @exportS3Method
settable_links.mixture3p <- function(model) {
  names(model$links)[model$links != "softmax"]
}


# user facing alias
#' @title `r .model_mixture3p()$name`
#' @details `r model_docs(.model_mixture3p(version = "slot", variable_precision = TRUE), components = c('domain', 'task', 'name', 'citation', 'requirements'))`
#' #### Version: `simple`
#' `r model_docs(.model_mixture3p(version = "simple"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `slot`
#' `r model_docs(.model_mixture3p(version = "slot"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `slot_averaging`
#' `r model_docs(.model_mixture3p(version = "slot_averaging"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Variable precision
#' With `variable_precision = TRUE` every version gains `tau`, and `kappa`
#' becomes the precision of the average trial. Shown for `"simple"`:
#' `r model_docs(.model_mixture3p(variable_precision = TRUE), components = c('parameters', 'links', 'prior'))`
#' @section Priors from bmm 1.4:
#' `r .circmix_prior_migration_doc()`
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
#'   fixed. In the `"slot"` and `"slot_averaging"` versions a single set size
#'   cannot test the capacity rule: `K` is then estimated from the responses at
#'   that set size, and only below it.
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
#'       `min(1, K / set_size)`, and `ptarget` is the probability that a
#'       remembered response reports the target rather than a non-target.
#'     \item `"slot_averaging"`: as `"slot"`, but `K` slots are distributed over
#'       the items and precision grows with the number of slots an item holds.
#'   }
#'   Storage and swapping are separate in the capacity versions: the capacity
#'   rule says how often a response comes from memory and `ptarget` says which
#'   stored item is reported. In `"simple"` both weights are log odds against
#'   guessing, `thetat` for the target and `thetant` for all non-targets
#'   together, so neither is a probability on its own; [native_parameters()]
#'   returns the probabilities.
#' @param variable_precision Logical; if `TRUE`, the precision of memory varies
#'   from trial to trial (van den Berg et al., 2012). This adds the parameter
#'   `tau`, and `kappa` then describes the average trial: it is the
#'   concentration whose Fisher information `J(kappa)` is the mean of `J`. See
#'   [mixture2p()] for details of the parameterisation. Under
#'   `version = "slot_averaging"`, the precision of an item held in `j` slots
#'   is the sum of `j` independent slot draws,
#'   `gamma(shape = j * J(kappa) / tau, scale = tau)`, so its variance grows
#'   with `j`.
#' @param vp_nodes Number of quadrature nodes used when
#'   `variable_precision = TRUE`; must be an odd number of at least 41, and
#'   defaults to 81. See [mixture2p()] for the accuracy this gives.
#' @param links A named list of link functions for the model parameters. The
#'   links of `thetat` and `thetant` cannot be changed, because the softmax over
#'   them is part of the likelihood.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @references Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The
#'   precision of visual working memory is set by allocation of a shared
#'   resource. Journal of Vision, 9(10), Article 7.
#'
#'   Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution representations
#'   in visual working memory. Nature, 453(7192), 233-235.
#'
#'   van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. (2012).
#'   Variability in encoding precision accounts for visual short-term memory
#'   limitations. PNAS, 109(22), 8780-8785.
#' @export
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
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
                      variable_precision = FALSE, vp_nodes = 81L, links = NULL,
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
# CHECK_DATA METHODS                                                     ####
############################################################################# !

# check_data.non_targets() runs first and adds ss_numeric. The wording is the
# one check_data.mixture2p() uses for the same limit.
#' @export
check_data.mixture3p <- function(model, data, formula) {
  set_sizes <- unique(data$ss_numeric)
  warnif(
    model$version != "simple" && length(set_sizes) == 1,
    "The set size does not vary, so the data cannot test the capacity rule of \\
    the '{model$version}' version. K is estimated from the responses at set \\
    size {set_sizes}, and reflects only its prior where it reaches {set_sizes}."
  )
  NextMethod("check_data")
}

############################################################################# !
# CHECK_FORMULA METHODS                                                  ####
############################################################################# !

#' @export
check_formula.mixture3p <- function(model, data, formula) {
  .mixture3p_warn_set_size_one(model, data, formula)
  NextMethod("check_formula")
}

# A trial with one item has no non-target to swap to, so the swap parameter is
# absent from its likelihood. The previous mixture implementation hid this by
# pinning the set-size-1 level to a constant, which it had to do anyway to
# switch off the sentinel weight.
.mixture3p_warn_set_size_one <- function(model, data, formula) {
  swap <- .mixture3p_version_table[[model$version]]$swap_parameter
  swap_formula <- formula[[swap]]
  if (is.null(swap_formula) || is_constant(swap_formula) || is_nl(swap_formula) ||
    !isTRUE(any(data[["ss_numeric"]] == 1))) {
    return(invisible())
  }
  .warn_set_size_one(
    .set_size_one_undetermined(swap, formula, data), swap,
    reason = glue("say nothing about '{swap}'")
  )
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

  # the softmax over the weights is part of the likelihood, so brms passes the
  # weights through on the scale they are sampled on
  family_model <- model
  family_model$links[family_model$links == "softmax"] <- "identity"

  formula <- bmf2bf(model, formula)
  formula$family <- .circmix_custom_family(
    family_model,
    family = paste0("mixture3p_", model$version),
    weight_parameters = spec$weight_parameters,
    vint = TRUE, n_vreal = n_nt,
    log_lik = .mixture3p_log_lik(model$version),
    posterior_predict = .mixture3p_posterior_predict(model$version),
    posterior_epred = .mixture3p_posterior_epred(model$version)
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

.mixture3p_posterior_epred <- function(version) {
  switch(version,
    simple = posterior_epred_mixture3p_simple,
    slot = posterior_epred_mixture3p_slot,
    slot_averaging = posterior_epred_mixture3p_slot_averaging
  )
}

# add_posterior_epred() looks these up by family name for fits saved without one
posterior_epred_mixture3p_simple <- posterior_epred_undefined("mixture3p")
posterior_epred_mixture3p_slot <- posterior_epred_undefined("mixture3p")
posterior_epred_mixture3p_slot_averaging <- posterior_epred_undefined("mixture3p")

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
    p_target = brms::get_dpar(prep, "ptarget", i = i),
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
    p_target = brms::get_dpar(prep, "ptarget", i = i),
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
    brms::get_dpar(prep, "ptarget", i = i), set_size
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
    p_target = brms::get_dpar(prep, "ptarget", i = i),
    tau = .circmix_prep_tau(prep, i)
  )
}
