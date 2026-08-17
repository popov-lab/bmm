############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.mixture2p_location <- glue(
  "Location of the von Mises distribution of memory responses (in radians). \\
  Fixed internally to 0 by default."
)

.mixture2p_concentration <- "Concentration of the von Mises distribution of memory responses"

.mixture2p_version_table <- list(
  simple = list(
    weight_parameter = "thetat",
    needs_set_size = FALSE,
    parameters = list(
      mu = .mixture2p_location,
      kappa = .mixture2p_concentration,
      thetat = glue(
        "Mixture weight for target responses, i.e. the probability that a \\
        response comes from memory rather than from guessing."
      )
    ),
    links = list(mu = "tan_half", kappa = "log", thetat = "logit"),
    fixed_parameters = list(mu = 0),
    priors = list(
      mu = list(main = "student_t(1, 0, 1)"),
      kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
      thetat = list(main = "logistic(0, 1)")
    ),
    init_ranges = list(mu = c(-0.1, 0.1), kappa = c(3, 8), thetat = c(0.6, 0.9))
  ),
  slot = list(
    weight_parameter = "K",
    needs_set_size = TRUE,
    parameters = list(
      mu = .mixture2p_location,
      kappa = .mixture2p_concentration,
      K = glue(
        "Capacity: the number of items that can be held in memory. An item is \\
        either held at full precision or not held at all, so the probability \\
        of a memory response is min(1, K / set_size)."
      )
    ),
    links = list(mu = "tan_half", kappa = "log", K = "log"),
    fixed_parameters = list(mu = 0),
    priors = list(
      mu = list(main = "student_t(1, 0, 1)"),
      kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
      K = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)")
    ),
    init_ranges = list(mu = c(-0.1, 0.1), kappa = c(3, 8), K = c(2, 4))
  ),
  slot_averaging = list(
    weight_parameter = "K",
    needs_set_size = TRUE,
    parameters = list(
      mu = .mixture2p_location,
      kappa = glue(
        "Concentration of the von Mises distribution for an item held in a \\
        single slot. An item holding several slots is more precise, because \\
        averaging independent samples adds their Fisher information."
      ),
      K = glue(
        "Capacity: the number of memory slots distributed over the items in \\
        the array. An item receives floor(K / set_size) or one more slot, and \\
        an item that receives none is guessed."
      )
    ),
    links = list(mu = "tan_half", kappa = "log", K = "log"),
    fixed_parameters = list(mu = 0),
    priors = list(
      mu = list(main = "student_t(1, 0, 1)"),
      kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
      K = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)")
    ),
    init_ranges = list(mu = c(-0.1, 0.1), kappa = c(3, 8), K = c(2, 4))
  )
)

.model_mixture2p <- function(resp_error = NULL, set_size = NULL, links = NULL,
                             version = "simple", variable_precision = FALSE,
                             vp_nodes = 41L, call = NULL, ...) {
  spec <- .mixture2p_version_table[[version]]
  if (variable_precision) {
    spec <- .circmix_add_variable_precision(spec)
  }

  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(set_size),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Two-parameter mixture model by Zhang and Luck (2008).",
      version = version,
      citation = glue(
        "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution \\
        representations in visual working memory. Nature, 453(7192), 233-235"
      ),
      requirements = glue(
        "- The response variable should be in radians and represent the \\
        angular error relative to the target", "\n",
        "- The 'slot' and 'slot_averaging' versions additionally require a \\
        set_size variable, which has to vary for the capacity K to be identified"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = spec$fixed_parameters,
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      deprecated_parameters = list(mu1 = "mu"),
      variable_precision = variable_precision,
      vp_nodes = as.integer(vp_nodes)
    ),
    class = c("bmmodel", "circular", "mixture2p", paste0("mixture2p_", version)),
    call = call
  )
  out$links[names(links)] <- links
  out
}

# user facing alias

#' @title `r .model_mixture2p()$name`
#' @details `r model_info(.model_mixture2p())`
#' @param resp_error The name of the variable in the provided dataset containing
#'   the response error. The response Error should code the response relative to
#'   the to-be-recalled target in radians. You can transform the response error
#'   in degrees to radian using the `deg2rad` function.
#' @param set_size Name of the variable in the dataset containing the set size,
#'   or a single number if the set size is constant. Required by the `"slot"`
#'   and `"slot_averaging"` versions and ignored by `"simple"`. The capacity `K`
#'   is only identified if the set size varies.
#' @param version A character string specifying which version of the model to
#'   use. Options are:
#'   \itemize{
#'     \item `"simple"` (default): the standard two-parameter mixture model, in
#'       which the probability of a memory response (`thetat`) and the precision
#'       of memory (`kappa`) are free parameters.
#'     \item `"slot"`: the fixed-resolution slot model of Zhang and Luck (2008).
#'       An item is either held at full precision or not held at all, so the
#'       probability of a memory response is `min(1, K / set_size)` and `kappa`
#'       does not vary with set size. Note that `K` is not identified above
#'       `max(set_size)`, where every item is remembered whatever `K` is, and
#'       that the cap introduces a kink in the gradient at `K = set_size`.
#'     \item `"slot_averaging"`: the slots-plus-averaging model of Zhang and
#'       Luck (2008). `K` slots are distributed over the items, an item
#'       receives `floor(K / set_size)` or one more, and precision grows with
#'       the number of slots because averaging independent samples adds their
#'       Fisher information.
#'   }
#' @param variable_precision Logical; if `TRUE`, the precision of memory varies
#'   from trial to trial rather than being constant (van den Berg et al., 2012;
#'   Fougnie et al., 2012). This adds the parameter `tau`, and the likelihood
#'   integrates the Fisher information `J` of a memory representation over
#'   `gamma(mean = J(kappa), scale = tau)`. `kappa` keeps its meaning as the
#'   mean precision, so `tau -> 0` recovers the constant-precision model.
#'   The integral is evaluated by quadrature, which costs about ten times as
#'   much per gradient at the default `vp_nodes`; expect a larger factor in wall
#'   time, since the posterior geometry also needs more leapfrog steps. `tau` is
#'   the parameter that needs the most data, and is often only weakly identified
#'   in small designs.
#' @param vp_nodes Number of quadrature nodes used when
#'   `variable_precision = TRUE`; must be an odd number of at least 41. The
#'   default holds the log likelihood to about 1e-5 as long as the implied gamma
#'   shape `J(kappa) / tau` stays above roughly 1. Below that the model rejects
#'   rather than returning a biased likelihood, and asks for more nodes; 81
#'   nodes reach a shape of about 0.5, and 161 about 0.25.
#' @param links A named list of link functions for the model parameters, e.g.
#'   `list(kappa = "softplus")`. Defaults are `"log"` for `kappa`, `K` and
#'   `tau`, and `"logit"` for `thetat`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @references Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
#'   representations in visual working memory. Nature, 453(7192), 233-235.
#'
#'   van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. (2012).
#'   Variability in encoding precision accounts for visual short-term memory
#'   limitations. PNAS, 109(22), 8780-8785.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # generate artificial data
#' dat <- data.frame(y = rmixture2p(n = 2000))
#'
#' # define formula
#' ff <- bmmformula(kappa ~ 1, thetat ~ 1)
#'
#' model <- mixture2p(resp_error = "y")
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
#' @export
mixture2p <- function(resp_error, set_size = NULL,
                      version = c("simple", "slot", "slot_averaging"),
                      variable_precision = FALSE, vp_nodes = 41L, links = NULL,
                      ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)

  stopif(
    .mixture2p_version_table[[version]]$needs_set_size && is.null(set_size),
    "The '{version}' version of mixture2p predicts the probability of a memory \\
    response from the set size, so the set_size argument is required."
  )
  .circmix_check_variable_precision(variable_precision, vp_nodes)

  .model_mixture2p(
    resp_error = resp_error, set_size = set_size, version = version,
    variable_precision = variable_precision, vp_nodes = vp_nodes,
    links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_DATA METHODS                                                     ####
############################################################################# !

#' @export
check_data.mixture2p <- function(model, data, formula) {
  if (!.mixture2p_version_table[[model$version]]$needs_set_size) {
    return(NextMethod("check_data"))
  }

  ss <- check_var_set_size(model$other_vars$set_size, data)
  warnif(
    length(unique(ss$ss_numeric)) == 1,
    "The set size does not vary, so the capacity K of the '{model$version}' \\
    version is not identified by the data and will only reflect its prior."
  )
  data$ss_numeric <- ss$ss_numeric
  attr(data, "max_set_size") <- ss$max_set_size

  NextMethod("check_data")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

#' @export
bmf2bf.mixture2p <- function(model, formula) {
  if (!.mixture2p_version_table[[model$version]]$needs_set_size) {
    return(NULL)
  }
  brms::bf(.circmix_aterm(model$resp_vars$resp_error, vint = "ss_numeric"))
}

#' @export
configure_model.mixture2p <- function(model, data, formula) {
  spec <- .mixture2p_version_table[[model$version]]

  formula <- bmf2bf(model, formula)
  formula$family <- .circmix_custom_family(
    model,
    family = paste0("mixture2p_", model$version),
    weight_parameters = spec$weight_parameter,
    vint = spec$needs_set_size,
    log_lik = .mixture2p_log_lik(model$version),
    posterior_predict = .mixture2p_posterior_predict(model$version)
  )

  nlist(
    formula, data,
    stanvars = .circmix_model_stanvars(
      model, formula$family, "mixture2p_funs.stan",
      vint = if (spec$needs_set_size) "ss" else NULL
    )
  )
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

.mixture2p_log_lik <- function(version) {
  switch(version,
    simple = log_lik_mixture2p_simple,
    slot = log_lik_mixture2p_slot,
    slot_averaging = log_lik_mixture2p_slot_averaging
  )
}

.mixture2p_posterior_predict <- function(version) {
  switch(version,
    simple = posterior_predict_mixture2p_simple,
    slot = posterior_predict_mixture2p_slot,
    slot_averaging = posterior_predict_mixture2p_slot_averaging
  )
}

log_lik_mixture2p_simple <- function(i, prep) {
  .dmixture2p_simple(
    prep$data$Y[i],
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    p_mem = brms::get_dpar(prep, "thetat", i = i),
    tau = .circmix_prep_tau(prep, i),
    nodes = .circmix_prep_nodes(prep)
  )
}

log_lik_mixture2p_slot <- function(i, prep) {
  .dmixture2p_slot(
    prep$data$Y[i],
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    K = brms::get_dpar(prep, "K", i = i),
    set_size = prep$data$vint1[i],
    tau = .circmix_prep_tau(prep, i),
    nodes = .circmix_prep_nodes(prep)
  )
}

log_lik_mixture2p_slot_averaging <- function(i, prep) {
  .dmixture2p_slot_averaging(
    prep$data$Y[i],
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    K = brms::get_dpar(prep, "K", i = i),
    set_size = prep$data$vint1[i],
    tau = .circmix_prep_tau(prep, i),
    nodes = .circmix_prep_nodes(prep)
  )
}

posterior_predict_mixture2p_simple <- function(i, prep, ...) {
  .rmixture2p_simple(
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    p_mem = brms::get_dpar(prep, "thetat", i = i),
    tau = .circmix_prep_tau(prep, i)
  )
}

posterior_predict_mixture2p_slot <- function(i, prep, ...) {
  K <- brms::get_dpar(prep, "K", i = i)
  .rmixture2p_simple(
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    p_mem = pmin(1, K / prep$data$vint1[i]),
    tau = .circmix_prep_tau(prep, i)
  )
}

posterior_predict_mixture2p_slot_averaging <- function(i, prep, ...) {
  .rmixture2p_slot_averaging(
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    K = brms::get_dpar(prep, "K", i = i),
    set_size = prep$data$vint1[i],
    tau = .circmix_prep_tau(prep, i)
  )
}
