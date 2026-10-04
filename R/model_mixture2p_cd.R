############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The retrieval side is mixture2p's version table, so the two models share
# their parameter names and meanings; a joint fit of both tasks relies on that.
.model_mixture2p_cd <- function(response = NULL, probe = NULL, target = NULL,
                                set_size = NULL, links = NULL, version = "simple",
                                variable_precision = FALSE, vp_nodes = 41L,
                                knowledge = "limited", call = NULL, ...) {
  spec <- .mixture2p_version_table[[version]]
  if (variable_precision) {
    spec <- .circmix_add_variable_precision(spec)
  }
  spec <- .cd_add_criterion(spec)

  out <- structure(
    list(
      resp_vars = nlist(response, probe, target),
      other_vars = nlist(set_size),
      domain = "Visual working memory (change detection)",
      task = "Change detection",
      name = "Two-parameter mixture model for single-probe change detection.",
      citation = c(
        glue(
          "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution \\
          representations in visual working memory. Nature, 453(7192), 233-235. \\
          https://doi.org/10.1038/nature06860"
        ),
        glue(
          "Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual \\
          working memory: Applications to the change detection task. Cognitive \\
          Psychology, 133, 101463. https://doi.org/10.1016/j.cogpsych.2022.101463"
        ),
        if (variable_precision) .circmix_vp_citation
      ),
      version = version,
      requirements = glue(
        "- The response variable should be coded 0 ('same') and 1 ('change')", "\n",
        "- The probe and the target should be in radians: the probed feature, \\
        and the feature shown at the probed location", "\n",
        "- The 'slot' and 'slot_averaging' versions additionally require a \\
        set_size variable, which has to vary for the capacity K to be identified"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = spec$fixed_parameters,
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      variable_precision = variable_precision,
      vp_nodes = as.integer(vp_nodes),
      knowledge = knowledge
    ),
    class = c(
      "bmmodel", "change_detection", "mixture2p_cd",
      paste0("mixture2p_cd_", version)
    ),
    call = call
  )
  out <- set_links(out, links)
  out
}

# user facing alias

#' @title `r .model_mixture2p_cd()$name`
#' @name mixture2p_cd
#' @details `r model_docs(.model_mixture2p_cd())`
#'
#' The observer retrieves a feature from the same memory as in [mixture2p()]
#' and compares it with the probe (Lin & Oberauer, 2022, Eq. 8). "Change" is
#' the answer when the log-likelihood ratio of a change against no change
#' exceeds `criterion`, which puts a "same" region around the probe; the
#' probability of a "same" response is the retrieval mass inside it. The
#' memory parameters have the same names and meanings as in [mixture2p()],
#' plus `criterion`. Whether a person gives the same values in both tasks is
#' an empirical question rather than something the model guarantees: the
#' precision estimated from continuous reproduction includes response noise
#' that change detection does not have, for example.
#'
#' @param response The name of the variable in the dataset containing the
#'   change-detection response, coded 0 ("same") and 1 ("change"). A logical
#'   variable is accepted and coded `TRUE` = "change".
#' @param probe The name of the variable containing the probed feature, in
#'   radians.
#' @param target The name of the variable containing the feature shown at the
#'   probed location, in radians. On "same" trials it equals the probe.
#' @param set_size Name of the variable in the dataset containing the set size,
#'   or a single number if the set size is constant. Required by the `"slot"`
#'   and `"slot_averaging"` versions and ignored by `"simple"`. The capacity `K`
#'   is only identified if the set size varies.
#' @param version Which storage model to use: `"simple"` (default), `"slot"`
#'   or `"slot_averaging"`. They are the versions of [mixture2p()].
#' @param variable_precision Logical; if `TRUE`, the precision of memory varies
#'   from trial to trial, as in [mixture2p()], and the model gains `tau`. The
#'   likelihood then integrates over precision for every trial, which makes it
#'   much slower than with constant precision. Whether `tau` can be
#'   estimated from change-detection responses alone has not been checked;
#'   expect it to be weakly identified.
#' @param vp_nodes Number of quadrature nodes used when
#'   `variable_precision = TRUE`; must be an odd number of at least 41. See
#'   [mixture2p()].
#' @param knowledge What the observer knows about the precision of the current
#'   trial when precision varies across trials (`variable_precision = TRUE`,
#'   or slot averaging with a capacity that is not a multiple of the set
#'   size). With `"limited"` (default), the observer knows only the
#'   distribution of precision and uses one decision boundary for all trials;
#'   with `"rich"`, the observer knows the precision of the trial and sets the
#'   boundary accordingly. Lin and Oberauer (2022) found the limited observer
#'   fit better for both mechanisms. With a single precision state both give
#'   the same model.
#' @param links A named list of link functions for the model parameters, e.g.
#'   `list(kappa = "softplus")`. Defaults are `"log"` for `kappa`, `K` and
#'   `tau`, `"logit"` for `thetat` and `"identity"` for `criterion`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @references Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
#'   representations in visual working memory. Nature, 453(7192), 233-235.
#'
#'   Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual
#'   working memory: Applications to the change detection task. Cognitive
#'   Psychology, 133, 101463. \doi{10.1016/j.cogpsych.2022.101463}
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # half of the trials probe the target itself, the others a random colour
#' dat <- data.frame(target = 0, probe = c(rep(0, 1000), runif(1000, -pi, pi)))
#' dat$change <- rmixture2p_cd(2000, dat$probe - dat$target, kappa = 6, p_mem = 0.7)
#'
#' fit <- bmm(
#'   formula = bmmformula(kappa ~ 1, thetat ~ 1),
#'   data = dat,
#'   model = mixture2p_cd(response = "change", probe = "probe", target = "target"),
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
#' @export
mixture2p_cd <- function(response, probe, target, set_size = NULL,
                         version = c("simple", "slot", "slot_averaging"),
                         variable_precision = FALSE, vp_nodes = 41L,
                         knowledge = c("limited", "rich"), links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  knowledge <- match.arg(knowledge)

  stopif(
    .mixture2p_version_table[[version]]$needs_set_size && is.null(set_size),
    "The '{version}' version of mixture2p_cd predicts the probability of a \\
    memory response from the set size, so the set_size argument is required."
  )
  .circmix_check_variable_precision(variable_precision, vp_nodes)

  .model_mixture2p_cd(
    response = response, probe = probe, target = target, set_size = set_size,
    version = version, variable_precision = variable_precision,
    vp_nodes = vp_nodes, knowledge = knowledge, links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_DATA METHODS                                                     ####
############################################################################# !

#' @export
check_data.mixture2p_cd <- function(model, data, formula) {
  data <- .mixture2p_check_set_size(model, data)
  NextMethod("check_data")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

#' @export
bmf2bf.mixture2p_cd <- function(model, formula) {
  needs_set_size <- .mixture2p_version_table[[model$version]]$needs_set_size
  brms::bf(.circmix_aterm(
    model$resp_vars$response,
    vint = if (needs_set_size) "ss_numeric",
    vreal = "probe_centered"
  ))
}

#' @export
configure_model.mixture2p_cd <- function(model, data, formula) {
  spec <- .mixture2p_version_table[[model$version]]
  family <- paste0("mixture2p_cd_", model$version)

  formula <- bmf2bf(model, formula)
  formula$family <- .cd_custom_family(
    model, family,
    weight_parameters = spec$weight_parameter,
    vint = spec$needs_set_size,
    log_lik = get(paste0("log_lik_", family), mode = "function"),
    posterior_predict = get(paste0("posterior_predict_", family), mode = "function"),
    posterior_epred = get(paste0("posterior_epred_", family), mode = "function")
  )

  nlist(
    formula, data,
    stanvars = .cd_model_stanvars(
      model, formula$family, c("mixture2p_funs.stan", "mixture2p_cd_funs.stan"),
      vint = if (spec$needs_set_size) "ss" else NULL
    )
  )
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

# P("same") of observation i for every posterior draw. The same function serves
# all three versions, because the version only decides which of p_mem and K
# the family estimates.
.mixture2p_cd_prep_psame <- function(prep, i, version) {
  capacity <- version != "simple"
  .mixture2p_cd_psame(
    probe = .cd_prep_probe(prep, i),
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    p_mem = if (!capacity) brms::get_dpar(prep, "thetat", i = i) else 1,
    criterion = brms::get_dpar(prep, "criterion", i = i),
    tau = .circmix_prep_tau(prep, i),
    K = if (capacity) brms::get_dpar(prep, "K", i = i) else 1,
    set_size = if (capacity) prep$data$vint1[i] else 1,
    version = version,
    rich = .cd_prep_rich(prep),
    nodes = .circmix_prep_nodes(prep)
  )
}

.mixture2p_cd_log_lik <- function(version) {
  function(i, prep) {
    .cd_bernoulli_ld(prep$data$Y[i], .mixture2p_cd_prep_psame(prep, i, version))
  }
}

.mixture2p_cd_posterior_predict <- function(version) {
  function(i, prep, ...) {
    p_same <- .mixture2p_cd_prep_psame(prep, i, version)
    stats::rbinom(length(p_same), 1, 1 - p_same)
  }
}

# The expected response is the probability of a "change" response. Every dpar
# arrives as a draws x observations matrix, so the cells are evaluated as one
# long vector and the data columns are lined up with them.
.mixture2p_cd_posterior_epred <- function(version) {
  function(prep) {
    capacity <- version != "simple"
    dpars <- prep$dpars
    p_same <- .mixture2p_cd_psame(
      probe = .epred_data(prep$data$vreal1, prep),
      mu = as.vector(dpars$mu),
      kappa = as.vector(dpars$kappa),
      p_mem = if (!capacity) as.vector(dpars$thetat) else 1,
      criterion = as.vector(dpars$criterion),
      tau = if (is.null(dpars$tau)) 0 else as.vector(dpars$tau),
      K = if (capacity) as.vector(dpars$K) else 1,
      set_size = if (capacity) .epred_data(prep$data$vint1, prep) else 1,
      version = version,
      rich = .cd_prep_rich(prep),
      nodes = .circmix_prep_nodes(prep)
    )
    .epred_matrix(1 - p_same, prep)
  }
}

log_lik_mixture2p_cd_simple <- .mixture2p_cd_log_lik("simple")
log_lik_mixture2p_cd_slot <- .mixture2p_cd_log_lik("slot")
log_lik_mixture2p_cd_slot_averaging <- .mixture2p_cd_log_lik("slot_averaging")

posterior_predict_mixture2p_cd_simple <- .mixture2p_cd_posterior_predict("simple")
posterior_predict_mixture2p_cd_slot <- .mixture2p_cd_posterior_predict("slot")
posterior_predict_mixture2p_cd_slot_averaging <-
  .mixture2p_cd_posterior_predict("slot_averaging")

posterior_epred_mixture2p_cd_simple <- .mixture2p_cd_posterior_epred("simple")
posterior_epred_mixture2p_cd_slot <- .mixture2p_cd_posterior_epred("slot")
posterior_epred_mixture2p_cd_slot_averaging <-
  .mixture2p_cd_posterior_epred("slot_averaging")
