############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The same retrieval mixture as imm(), judged against a probe instead of
# reported: the parameters are imm's, plus the decision criterion
.model_imm_cd <- function(response = NULL, probe = NULL, target = NULL,
                          nt_features = NULL, nt_distances = NULL,
                          set_size = NULL, regex = FALSE, version = "full",
                          links = NULL, variable_precision = FALSE,
                          vp_nodes = 41L, knowledge = "limited", call = NULL,
                          ...) {
  spec <- .imm_version_spec(version)
  if (variable_precision) {
    spec <- .circmix_add_variable_precision(spec)
  }
  spec$fixed_parameters <- list(mu = 0, b = 0)
  spec <- .cd_add_criterion(spec)

  out <- structure(
    list(
      resp_vars = nlist(response, probe, target),
      other_vars = nlist(nt_features, nt_distances, set_size),
      domain = "Visual working memory (change detection)",
      task = "Change detection",
      name = "Interference measurement model for change detection by Lin and Oberauer (2022).",
      citation = c(
        glue(
          "Oberauer, K., & Lin, H.-Y. (2017). An interference model \\
          of visual working memory. Psychological Review, 124(1), 21-59. \\
          https://doi.org/10.1037/rev0000044"
        ),
        glue(
          "Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual \\
          working memory: Applications to the change detection task. \\
          Cognitive Psychology, 133, 101463. \\
          https://doi.org/10.1016/j.cogpsych.2022.101463"
        ),
        if (variable_precision) .circmix_vp_citation
      ),
      version = version,
      requirements = glue(
        "- The response should be coded 0 for 'same' and 1 for 'change' \\
          (or be logical, TRUE = 'change')
          - The probe and the target (the feature shown at the probed location) \\
          should be in radians
          - The non-target features should be in radians and be centered \\
          relative to the target"
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
    regex = regex,
    regex_vars = if (version == "abc") "nt_features" else c("nt_features", "nt_distances"),
    class = c(
      "bmmodel", "change_detection", "non_targets", "imm_cd",
      paste0("imm_cd_", version)
    ),
    call = call
  )

  out <- set_links(out, links)
  out
}

# user facing alias

#' @title `r .model_imm_cd()$name`
#' @description The `r .model_imm_cd()$name` The observer retrieves a feature
#'   from the retrieval distribution of [imm()] and compares it with the probe:
#'   "change" is answered when the retrieved feature is unlikely under the
#'   hypothesis that the probe is the stored item (Lin & Oberauer, 2022, Eq. 8).
#'   The memory parameters are those of [imm()], so the versions `full`, `bsc`
#'   and `abc` and `variable_precision` are the same; the decision adds the
#'   criterion.
#'
#' @name imm_cd
#' @details `r model_docs(.model_imm_cd(), components = c('domain', 'task', 'name', 'citation'))`
#' #### Version: `full`
#' `r model_docs(.model_imm_cd(version = "full"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `bsc`
#' `r model_docs(.model_imm_cd(version = "bsc"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `abc`
#' `r model_docs(.model_imm_cd(version = "abc"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#'
#' #### Comparing with continuous reproduction
#'
#' `imm_cd()` and [imm()] share their parameter names, but that does not make
#' their estimates comparable at face value. In continuous reproduction `kappa`
#' also absorbs noise in reporting on the response wheel, which change detection
#' does not have, and participants may encode differently when they expect one
#' test or the other. Treat equal parameters across the two tasks as a
#' hypothesis to test, not as a given.
#'
#' @param response The name of the variable in the dataset containing the
#'   binary response, 0 for "same" and 1 for "change". A logical variable is
#'   accepted (`TRUE` = "change").
#' @param probe The name of the variable containing the feature of the probe,
#'   in radians.
#' @param target The name of the variable containing the feature the probed item
#'   had in the memory array, in radians. The probe is taken relative to it.
#' @param nt_features A character vector with the names of the non-target
#'   variables. The non-target variables should be in radians and be centered
#'   relative to the target. Alternatively, if `regex = TRUE`, a regular
#'   expression can be used to match the non-target feature columns in the
#'   dataset.
#' @param nt_distances A vector of names of the columns containing the distances
#'   of non-target items to the target item. Alternatively, if `regex = TRUE`, a
#'   regular expression can be used to match the non-target distances columns
#'   in the dataset. Only necessary for the `bsc` and `full` versions.
#' @param set_size Name of the column containing the set size variable (if
#'   set_size varies) or a numeric value for the set_size, if the set_size is
#'   fixed.
#' @param regex Logical. If TRUE, the `nt_features` and `nt_distances` arguments
#'   are interpreted as a regular expression to match the non-target feature
#'   columns in the dataset.
#' @param version Character. The version of the IMM to use. Can be one of
#'   `full`, `bsc`, or `abc`. The default is `full`. See [imm()].
#' @param variable_precision Logical; if `TRUE`, the precision of memory varies
#'   from trial to trial (van den Berg et al., 2012). This adds the parameter
#'   `tau`; see [mixture2p()] for the parameterisation.
#' @param vp_nodes Number of quadrature nodes used when
#'   `variable_precision = TRUE`; must be an odd number of at least 41. See
#'   [mixture2p()].
#' @param knowledge What the observer knows about the precision of the current
#'   trial when deciding. `"limited"` (the default) knows only how precision is
#'   distributed and uses one decision boundary for all trials; `"rich"` knows
#'   the precision of the trial and adapts the boundary to it. Lin and Oberauer
#'   (2022, Table 2) found the limited observer fitted better under variable
#'   precision. Without `variable_precision` the two are the same model.
#' @param links A named list of link functions for the model parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simulate a single-probe change-detection task with 4 items: a third of the
#' # probes repeat the target, a third show a new color, and a third show the
#' # color of a non-target
#' n <- 1200
#' nt <- matrix(runif(n * 3, -pi, pi), ncol = 3, dimnames = list(NULL, paste0("nt", 1:3)))
#' dist <- matrix(runif(n * 3, 0.3, 2), ncol = 3, dimnames = list(NULL, paste0("dist", 1:3)))
#' type <- sample(c("same", "new", "intrusion"), n, replace = TRUE)
#' probe <- ifelse(type == "same", 0, ifelse(type == "new", runif(n, -pi, pi), nt[, 1]))
#' dat <- data.frame(target = 0, probe = probe, nt, dist, set_size = 4)
#' dat$change <- vapply(seq_len(n), function(i) {
#'   rimm_cd(1, dat$probe[i],
#'     mu = c(0, nt[i, ]), dist = c(0, dist[i, ]),
#'     c = 4, a = 0.5, b = 1, s = 2, kappa = 8
#'   )
#' }, numeric(1))
#'
#' model <- imm_cd(
#'   response = "change", probe = "probe", target = "target",
#'   nt_features = paste0("nt", 1:3), nt_distances = paste0("dist", 1:3),
#'   set_size = "set_size"
#' )
#' fit <- bmm(
#'   formula = bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#' @export
imm_cd <- function(response, probe, target, nt_features, nt_distances, set_size,
                   regex = FALSE, version = c("full", "bsc", "abc"),
                   variable_precision = FALSE, vp_nodes = 41L,
                   knowledge = c("limited", "rich"), links = NULL, ...) {
  call <- match.call()
  version <- match.arg(version)
  knowledge <- match.arg(knowledge)
  if (version == "abc") nt_distances <- NULL
  stop_missing_args()
  .circmix_check_variable_precision(variable_precision, vp_nodes)

  .model_imm_cd(
    response = response, probe = probe, target = target,
    nt_features = nt_features, nt_distances = nt_distances,
    set_size = set_size, regex = regex, version = version, links = links,
    variable_precision = variable_precision, vp_nodes = vp_nodes,
    knowledge = knowledge, call = call, ...
  )
}

############################################################################# !
# CHECK_DATA S3 methods                                                  ####
############################################################################# !

#' @export
check_data.imm_cd_bsc <- function(model, data, formula) {
  data <- .check_data_imm_dist(model, data, formula)
  NextMethod("check_data")
}

#' @export
check_data.imm_cd_full <- function(model, data, formula) {
  data <- .check_data_imm_dist(model, data, formula)
  NextMethod("check_data")
}

############################################################################# !
# CHECK_FORMULA METHODS                                                  ####
############################################################################# !

#' @export
check_formula.imm_cd <- function(model, data, formula) {
  .imm_check_set_size_intercept(model, formula)
  NextMethod("check_formula")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

# probe_centered comes first so that it is vreal1, the probe the wrapper passes
# as a bare real; the non-target features and distances follow as for imm()
#' @export
bmf2bf.imm_cd <- function(model, formula) {
  spec <- .imm_version_spec(model$version)
  covariates <- c("probe_centered", model$other_vars$nt_features)
  if (spec$needs_distances) {
    covariates <- c(covariates, model$other_vars$nt_distances)
  }
  brms::bf(.circmix_aterm(
    model$resp_vars$response,
    vint = "ss_numeric", vreal = covariates
  ))
}

#' @export
configure_model.imm_cd <- function(model, data, formula) {
  spec <- .imm_version_spec(model$version)
  n_nt <- attr(data, "max_set_size") - 1
  groups <- c(
    list(probe = 1, nt = n_nt),
    if (spec$needs_distances) list(dist = n_nt)
  )
  family <- paste0("imm_cd_", model$version)

  formula <- bmf2bf(model, formula)
  formula$family <- .cd_custom_family(
    model, family, spec$weight_parameters,
    vint = TRUE, n_vreal = sum(unlist(groups)),
    log_lik = get(paste0("log_lik_", family), mode = "function"),
    posterior_predict = get(paste0("posterior_predict_", family), mode = "function"),
    posterior_epred = get(paste0("posterior_epred_", family), mode = "function")
  )

  nlist(
    formula, data,
    stanvars = .cd_model_stanvars(model, formula$family,
      c("imm_funs.stan", "imm_cd_funs.stan"),
      vint = "ss", vreal = groups
    )
  )
}

############################################################################# !
# CONFIGURE_PRIOR METHODS                                                ####
############################################################################# !

# set-size-1 trials inform a and s no more than they do in imm()
#' @export
configure_prior.imm_cd <- function(model, data, formula, user_prior, ...) {
  .imm_set_size1_prior(model, data, formula)
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

.imm_cd_args <- function(i, prep, version) {
  covariates <- .imm_covariates(prep, i, version != "abc", skip = 1L)
  list(
    probe = .cd_prep_probe(prep, i),
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    c = brms::get_dpar(prep, "c", i = i),
    a = if (version != "bsc") brms::get_dpar(prep, "a", i = i),
    s = if (version != "abc") brms::get_dpar(prep, "s", i = i),
    b = brms::get_dpar(prep, "b", i = i),
    criterion = brms::get_dpar(prep, "criterion", i = i),
    set_size = prep$data$vint1[i], nt = covariates$nt, dist = covariates$dist,
    tau = .circmix_prep_tau(prep, i), nodes = .circmix_prep_nodes(prep),
    rich = .cd_prep_rich(prep), version = version
  )
}

.imm_cd_log_lik <- function(i, prep, version) {
  .cd_bernoulli_ld(
    prep$data$Y[i], brms::do_call(.imm_cd_psame, .imm_cd_args(i, prep, version))
  )
}

.imm_cd_predict <- function(i, prep, version) {
  p_same <- brms::do_call(.imm_cd_psame, .imm_cd_args(i, prep, version))
  stats::rbinom(length(p_same), 1, 1 - p_same)
}

# the expected response is P("change")
.imm_cd_epred <- function(prep, version) {
  p_change <- vapply(seq_len(prep$nobs), function(i) {
    1 - brms::do_call(.imm_cd_psame, .imm_cd_args(i, prep, version))
  }, numeric(prep$ndraws))
  .epred_matrix(p_change, prep)
}

log_lik_imm_cd_full <- function(i, prep) .imm_cd_log_lik(i, prep, "full")
log_lik_imm_cd_bsc <- function(i, prep) .imm_cd_log_lik(i, prep, "bsc")
log_lik_imm_cd_abc <- function(i, prep) .imm_cd_log_lik(i, prep, "abc")

posterior_predict_imm_cd_full <- function(i, prep, ...) .imm_cd_predict(i, prep, "full")
posterior_predict_imm_cd_bsc <- function(i, prep, ...) .imm_cd_predict(i, prep, "bsc")
posterior_predict_imm_cd_abc <- function(i, prep, ...) .imm_cd_predict(i, prep, "abc")

posterior_epred_imm_cd_full <- function(prep) .imm_cd_epred(prep, "full")
posterior_epred_imm_cd_bsc <- function(prep) .imm_cd_epred(prep, "bsc")
posterior_epred_imm_cd_abc <- function(prep) .imm_cd_epred(prep, "abc")
