############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.mixture3p_cd_citation <- paste(
  "Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual working",
  "memory: Applications to the change detection task. Cognitive Psychology,",
  "133, Article 101463. https://doi.org/10.1016/j.cogpsych.2022.101463"
)

# The retrieval mixture is mixture3p's, version by version, so the parameters
# (and their links, priors and inits) come from its version table; only
# criterion is added. Joint fits of both tasks rely on the shared names.
.model_mixture3p_cd <- function(response = NULL, probe = NULL, target = NULL,
                                nt_features = NULL, set_size = NULL,
                                regex = FALSE, links = NULL, version = "simple",
                                variable_precision = FALSE, vp_nodes = 41L,
                                knowledge = "limited", call = NULL, ...) {
  spec <- .mixture3p_version_table[[version]]
  if (variable_precision) {
    spec <- .circmix_add_variable_precision(spec)
  }
  spec <- .cd_add_criterion(spec)

  out <- structure(
    list(
      resp_vars = nlist(response, probe, target),
      other_vars = nlist(nt_features, set_size),
      domain = "Visual working memory (change detection)",
      task = "Change detection",
      name = "Three-parameter mixture model for change detection.",
      citation = c(
        .mixture3p_citation,
        if (version != "simple") .mixture3p_capacity_citation,
        .mixture3p_cd_citation,
        if (variable_precision) .circmix_vp_citation
      ),
      version = version,
      requirements = glue(
        "- The response should be coded 0 for 'same' and 1 for 'change' \\
        (logical values are converted)
        - The probe and the target (the feature shown at the probed location) \\
        should be in radians
        - The non-target features should be in radians and be centered \\
        relative to the target"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = c(list(mu = 0), spec$fixed_parameters),
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      variable_precision = variable_precision,
      vp_nodes = as.integer(vp_nodes),
      knowledge = knowledge
    ),
    regex = regex,
    regex_vars = c("nt_features"),
    class = c(
      "bmmodel", "change_detection", "non_targets", "mixture3p_cd",
      paste0("mixture3p_cd_", version)
    ),
    call = call
  )
  out <- set_links(out, links)
  out
}

# as in mixture3p(), the softmax over thetat and thetant is part of the
# likelihood
#' @exportS3Method
settable_links.mixture3p_cd <- function(model) {
  names(model$links)[model$links != "softmax"]
}

# user facing alias
#' @title `r .model_mixture3p_cd()$name`
#' @name mixture3p_cd
#' @description The three-parameter mixture model of [mixture3p()] applied to
#'   single-probe change detection: on each trial one item of the memory array
#'   is probed, and the observer reports whether the probe differs from the item
#'   shown at that location. The observer retrieves a feature from the same
#'   retrieval mixture as in continuous reproduction (target, non-targets,
#'   guessing) and says "change" when the log-likelihood ratio of a change
#'   against no change exceeds `criterion` (Lin & Oberauer, 2022). Non-targets
#'   enter only through retrieval, so a probe that matches a non-target is
#'   harder to reject than a new one.
#' @details `r model_docs(.model_mixture3p_cd())`
#'
#'   The memory parameters have the same names and meaning as in [mixture3p()],
#'   which is what makes the two tasks comparable. Whether the same person
#'   yields the same values in both tasks is an empirical question, not
#'   something the model guarantees: `kappa` in continuous reproduction also
#'   absorbs response noise, which change detection does not have.
#'
#'   As in [mixture3p()], a trial with a single item has no non-target, so the
#'   set-size-1 level of `thetant` (or `pnt`) is not identified and samples its
#'   prior; the model warns when a formula asks for it.
#' @inheritParams mixture3p
#' @param response The name of the variable in the data containing the change
#'   detection response, coded 0 for "same" and 1 for "change". A logical
#'   variable is converted with a warning.
#' @param probe The name of the variable containing the probe feature in
#'   radians.
#' @param target The name of the variable containing the feature shown at the
#'   probed location in radians. On "same" trials it equals the probe.
#' @param nt_features A character vector with the names of the non-target
#'   feature values. The non-target feature values should be in radians and
#'   centered relative to the target. Alternatively, if `regex = TRUE`, a
#'   regular expression can be used to match the non-target feature columns in
#'   the dataset.
#' @param version A character string specifying which version of the retrieval
#'   mixture to use: `"simple"` (default), `"slot"` or `"slot_averaging"`. They
#'   are the versions of [mixture3p()] and have the same parameters.
#' @param knowledge What the observer knows about the precision of the current
#'   trial when it places the decision boundary. `"limited"` (default): only the
#'   distribution of precision across trials, which Lin and Oberauer (2022)
#'   found to describe their data better. `"rich"`: the precision of the current
#'   trial. The two coincide unless precision varies across trials, that is
#'   unless `variable_precision = TRUE` or `version = "slot_averaging"`.
#' @param links A named list of link functions for the model parameters. The
#'   links of `thetat` and `thetant` cannot be changed, because the softmax over
#'   them is part of the likelihood.
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @references Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The
#'   precision of visual working memory is set by allocation of a shared
#'   resource. Journal of Vision, 9(10), 7.
#'
#'   Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual working
#'   memory: Applications to the change detection task. Cognitive Psychology,
#'   133, 101463.
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simulate change detection with three items: a third of the probes are the
#' # target, a third a non-target and a third a new colour
#' n <- 1500
#' dat <- data.frame(
#'   target = 0, nt1 = runif(n, -pi, pi), nt2 = runif(n, -pi, pi),
#'   type = rep(c("same", "intrusion", "new"), length.out = n)
#' )
#' dat$probe <- ifelse(dat$type == "same", 0,
#'   ifelse(dat$type == "intrusion", dat$nt1, runif(n, -pi, pi))
#' )
#' dat$response <- vapply(seq_len(n), function(i) {
#'   rmixture3p_cd(1, dat$probe[i],
#'     mu = c(0, dat$nt1[i], dat$nt2[i]), kappa = 8, p_mem = 0.7, p_nt = 0.1
#'   )
#' }, numeric(1))
#'
#' model <- mixture3p_cd("response", "probe", "target",
#'   nt_features = c("nt1", "nt2"), set_size = 3
#' )
#' fit <- bmm(
#'   bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1), dat, model,
#'   cores = 4, backend = "cmdstanr"
#' )
mixture3p_cd <- function(response, probe, target, nt_features, set_size,
                         regex = FALSE,
                         version = c("simple", "slot", "slot_averaging"),
                         variable_precision = FALSE, vp_nodes = 41L,
                         knowledge = c("limited", "rich"), links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  knowledge <- match.arg(knowledge)
  .circmix_check_variable_precision(variable_precision, vp_nodes)

  .model_mixture3p_cd(
    response = response, probe = probe, target = target,
    nt_features = nt_features, set_size = set_size, regex = regex,
    version = version, variable_precision = variable_precision,
    vp_nodes = vp_nodes, knowledge = knowledge, links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_FORMULA METHODS                                                  ####
############################################################################# !

#' @export
check_formula.mixture3p_cd <- function(model, data, formula) {
  .mixture3p_warn_set_size_one(model, data, formula)
  NextMethod("check_formula")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

# probe_centered comes first, so the wrapper receives it as vreal1
#' @export
bmf2bf.mixture3p_cd <- function(model, formula) {
  brms::bf(.circmix_aterm(
    model$resp_vars$response,
    vint = "ss_numeric",
    vreal = c("probe_centered", model$other_vars$nt_features)
  ))
}

#' @export
configure_model.mixture3p_cd <- function(model, data, formula) {
  spec <- .mixture3p_version_table[[model$version]]
  n_nt <- attr(data, "max_set_size") - 1
  family <- paste0("mixture3p_cd_", model$version)

  family_model <- model
  family_model$links[family_model$links == "softmax"] <- "identity"

  formula <- bmf2bf(model, formula)
  formula$family <- .cd_custom_family(
    family_model, family, spec$weight_parameters,
    vint = TRUE, n_vreal = 1 + n_nt,
    log_lik = get(paste0("log_lik_", family), mode = "function"),
    posterior_predict = get(paste0("posterior_predict_", family), mode = "function"),
    posterior_epred = get(paste0("posterior_epred_", family), mode = "function")
  )

  nlist(
    formula, data,
    stanvars = .cd_model_stanvars(
      model, formula$family, c("mixture3p_funs.stan", "mixture3p_cd_funs.stan"),
      vint = "ss", vreal = list(probe = 1, nt = n_nt)
    )
  )
}

############################################################################# !
# P("SAME") AND POSTPROCESS METHODS                                      ####
############################################################################# !
# R twins of the mixture3p_cd_*_core functions in mixture3p_cd_funs.stan: the
# weights of mixture3p's version, the component locations relative to the
# probe, and the shared change-detection evaluator. One observation (one set
# size and one row of non-targets) per call; the other arguments are vectors
# over posterior draws or trials.

.mixture3p_cd_vp_psame <- function(args, weights, set_size, nt, rich, nodes) {
  .circmix_cd_vp_psame(
    args$probe - .circmix_locations(args$mu, nt, set_size),
    weights$logw, weights$logw_guess, args$kappa, args$tau, args$criterion,
    exp(weights$logw[, 1]), rich, nodes
  )
}

.mixture3p_cd_psame_simple <- function(probe, mu, kappa, thetat, thetant,
                                       criterion, set_size, nt, tau, rich,
                                       nodes) {
  args <- .circmix_recycle(
    probe = probe, mu = mu, kappa = kappa, thetat = thetat, thetant = thetant,
    criterion = criterion, tau = tau
  )
  .mixture3p_cd_vp_psame(
    args, .mixture3p_softmax_weights(args$thetat, args$thetant, set_size),
    set_size, nt, rich, nodes
  )
}

.mixture3p_cd_psame_slot <- function(probe, mu, kappa, K, p_nt, criterion,
                                     set_size, nt, tau, rich, nodes) {
  args <- .circmix_recycle(
    probe = probe, mu = mu, kappa = kappa, K = K, p_nt = p_nt,
    criterion = criterion, tau = tau
  )
  .mixture3p_cd_vp_psame(
    args,
    .mixture3p_nested_weights(pmin(1, args$K / set_size), args$p_nt, set_size),
    set_size, nt, rich, nodes
  )
}

.mixture3p_cd_psame_slot_averaging <- function(probe, mu, kappa, K, p_nt,
                                               criterion, set_size, nt, tau,
                                               rich, nodes) {
  args <- .circmix_recycle(
    probe = probe, mu = mu, kappa = kappa, K = K, p_nt = p_nt,
    criterion = criterion, tau = tau
  )
  .circmix_cd_slot_averaging_psame(
    args$probe - .circmix_locations(args$mu, nt, set_size),
    .mixture3p_swap_weights(args$p_nt, set_size),
    .circmix_slot_averaging_branches(args$K, set_size, args$kappa),
    args$tau, args$criterion, rich, nodes
  )
}

# P("same") for observation i, one value per posterior draw
.mixture3p_cd_prep_psame <- function(version, i, prep) {
  dpar <- function(name) brms::get_dpar(prep, name, i = i)
  common <- list(
    probe = .cd_prep_probe(prep, i), mu = dpar("mu"), kappa = dpar("kappa"),
    criterion = dpar("criterion"), set_size = prep$data$vint1[i],
    nt = .circmix_prep_nt(prep, i, skip = 1L), tau = .circmix_prep_tau(prep, i),
    rich = .cd_prep_rich(prep), nodes = .circmix_prep_nodes(prep)
  )
  switch(version,
    simple = do.call(.mixture3p_cd_psame_simple, c(common, list(
      thetat = dpar("thetat"), thetant = dpar("thetant")
    ))),
    slot = do.call(.mixture3p_cd_psame_slot, c(common, list(
      K = dpar("K"), p_nt = dpar("pnt")
    ))),
    slot_averaging = do.call(.mixture3p_cd_psame_slot_averaging, c(common, list(
      K = dpar("K"), p_nt = dpar("pnt")
    )))
  )
}

.mixture3p_cd_log_lik <- function(version) {
  function(i, prep) {
    # .cd_bernoulli_ld() decides with ifelse(), whose length is that of y
    .cd_bernoulli_ld(rep(prep$data$Y[i], prep$ndraws), .mixture3p_cd_prep_psame(version, i, prep))
  }
}

.mixture3p_cd_posterior_predict <- function(version) {
  function(i, prep, ...) {
    stats::rbinom(prep$ndraws, 1, 1 - .mixture3p_cd_prep_psame(version, i, prep))
  }
}

# the expected response is the probability of a "change" response
.mixture3p_cd_posterior_epred <- function(version) {
  function(prep) {
    .epred_matrix(vapply(seq_len(prep$nobs), function(i) {
      1 - .mixture3p_cd_prep_psame(version, i, prep)
    }, numeric(prep$ndraws)), prep)
  }
}

# named per family, because brms and add_posterior_epred() look them up by name
log_lik_mixture3p_cd_simple <- .mixture3p_cd_log_lik("simple")
log_lik_mixture3p_cd_slot <- .mixture3p_cd_log_lik("slot")
log_lik_mixture3p_cd_slot_averaging <- .mixture3p_cd_log_lik("slot_averaging")
posterior_predict_mixture3p_cd_simple <- .mixture3p_cd_posterior_predict("simple")
posterior_predict_mixture3p_cd_slot <- .mixture3p_cd_posterior_predict("slot")
posterior_predict_mixture3p_cd_slot_averaging <- .mixture3p_cd_posterior_predict("slot_averaging")
posterior_epred_mixture3p_cd_simple <- .mixture3p_cd_posterior_epred("simple")
posterior_epred_mixture3p_cd_slot <- .mixture3p_cd_posterior_epred("slot")
posterior_epred_mixture3p_cd_slot_averaging <- .mixture3p_cd_posterior_epred("slot_averaging")
