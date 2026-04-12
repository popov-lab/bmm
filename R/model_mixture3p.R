############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_mixture3p <- function(resp_error = NULL, nt_features = NULL,
                             set_size = NULL, response = NULL,
                             probe = NULL, target = NULL,
                             regex = FALSE, task = "de", links = NULL,
                             call = NULL, ...) {
  domain_class <- if (task == "cd") "change_detection" else "circular"

  out <- structure(
    list(
      resp_vars = if (task == "cd") nlist(response, probe, target) else nlist(resp_error),
      other_vars = nlist(nt_features, set_size),
      domain = "Visual working memory",
      task = if (task == "cd") "Change detection" else "Continuous reproduction",
      name = "Three-parameter mixture model by Bays et al (2009).",
      version = "NA",
      citation = glue(
        "Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). \\
        The precision of visual working memory is set by allocation \\
        of a shared resource. Journal of Vision, 9(10), 1-11"
      ),
      parameters = list(
        kappa = "Concentration parameter of the von Mises distribution",
        thetat = "Mixture weight for target responses",
        thetant = "Mixture weight for non-target responses"
      ),
      links = list(
        mu1 = "tan_half",
        kappa = "log",
        thetat = "softmax",
        thetant = "softmax"
      ),
      fixed_parameters = list(mu1 = 0, mu2 = 0, kappa2 = -100),
      default_priors = list(
        kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
        thetat = list(main = "logistic(0, 1)"),
        thetant = list(main = "logistic(0, 1)")
      ),
      fixed_parameters = list(),
      void_mu = task == "cd"
    ),
    regex = regex,
    regex_vars = c("nt_features"),
    class = c("bmmodel", domain_class, "non_targets", "mixture3p",
              paste0("mixture3p_", task)),
    call = call
  )

  if (task == "cd") {
    out$citation <- glue(
      "{out$citation}; \\
      Lin, H.Y., & Oberauer, K. (2022). An interference model for visual \\
      working memory: Applications to the change detection task. \\
      Cognitive Psychology, 133, 101463."
    )
    out$requirements <- glue(
      "- response: Binary (0='same', 1='change')
      - probe: Probe color in radians
      - target: Target color in radians
      - Non-target features should be in radians and centered relative to the target"
    )
    out$parameters <- c(
      list(mu = glue(
        "Location parameter (bias) of the retrieval distribution \\
        (in radians). Fixed to 0 by default."
      )),
      out$parameters
    )
    out$links <- c(list(mu = "tan_half"), out$links)
    out$fixed_parameters$mu <- 0
    out$default_priors <- c(
      list(mu = list(main = "student_t(1, 0, 1)")),
      out$default_priors
    )
    out$parameters$beta <- glue(
      "Decision criterion (log prior odds). \\
      Fixed to 0 by default for unbiased decision."
    )
    out$links$beta <- "identity"
    out$fixed_parameters$beta <- 0
    out$default_priors$beta <- list(main = "normal(0, 0.5)")
  } else {
    out$requirements <- glue(
      "- The response vairable should be in radians and \\
      represent the angular error relative to the target
      - The non-target features should be in radians and be \\
      centered relative to the target"
    )
    out$parameters <- c(
      list(mu1 = glue(
        "Location parameter of the von Mises distribution for memory responses \\
        (in radians). Fixed internally to 0 by default."
      )),
      out$parameters
    )
    out$links <- c(list(mu1 = "tan_half"), out$links)
    out$fixed_parameters <- list(mu1 = 0, mu2 = 0, kappa2 = -100)
    out$default_priors <- c(
      list(mu1 = list(main = "student_t(1, 0, 1)")),
      out$default_priors
    )
  }

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
#' @param task Character. The experimental task: `"de"` for delayed estimation
#'   (continuous reproduction) or `"cd"` for change detection. Default is `"de"`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
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
#' @param response The name of the variable in the dataset containing the binary
#'   response (0 = "same", 1 = "change"). Required when `task = "cd"`.
#' @param probe The name of the variable containing the probe color in radians.
#'   Required when `task = "cd"`.
#' @param target The name of the variable containing the target color in
#'   radians. Required when `task = "cd"`.
mixture3p <- function(resp_error = NULL, nt_features = NULL, set_size = NULL,
                      response = NULL, probe = NULL, target = NULL,
                      regex = FALSE, task = "de", ...) {
  call <- match.call()
  task <- match.arg(task, c("de", "cd"))
  dots <- list(...)
  if ("setsize" %in% names(dots)) {
    set_size <- dots$setsize
    warning2("The argument 'setsize' is deprecated. Please use 'set_size' instead.")
  }
  if (task == "de") {
    stopif(is.null(resp_error), "Argument 'resp_error' is required for task = 'de'.")
    stopif(is.null(nt_features), "Argument 'nt_features' is required.")
    stopif(is.null(set_size), "Argument 'set_size' is required.")
  } else {
    warning2("`mixture3p(task = \"cd\")` is deprecated. Please use `mixture3p_cd()` instead.")
    stopif(is.null(response), "Argument 'response' is required for task = 'cd'.")
    stopif(is.null(probe), "Argument 'probe' is required for task = 'cd'.")
    stopif(is.null(target), "Argument 'target' is required for task = 'cd'.")
    stopif(is.null(nt_features), "Argument 'nt_features' is required.")
    stopif(is.null(set_size), "Argument 'set_size' is required.")
  }
  .model_mixture3p(
    resp_error = resp_error, nt_features = nt_features,
    set_size = set_size, response = response, probe = probe,
    target = target, regex = regex, task = task, call = call, ...
  )
}

#' @rdname mixture3p
#' @export
mixture3p_cd <- function(response, probe, target, nt_features, set_size,
                         regex = FALSE, ...) {
  call <- match.call()
  stop_missing_args()
  .model_mixture3p(
    nt_features = nt_features,
    set_size = set_size,
    response = response,
    probe = probe,
    target = target,
    regex = regex,
    task = "cd",
    call = call,
    ...
  )
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.mixture3p_de <- function(model, data, formula) {
  # retrieve arguments from the data check
  max_set_size <- attr(data, "max_set_size")
  lure_idx <- attr(data, "lure_idx_vars")
  nt_features <- model$other_vars$nt_features
  set_size_var <- model$other_vars$set_size

  # construct initial brms formula
  formula <- bmf2bf(model, formula) +
    brms::lf(kappa2 ~ 1) +
    brms::lf(mu2 ~ 1) +
    brms::nlf(theta1 ~ thetat) +
    brms::nlf(kappa1 ~ kappa)

  # additional internal terms for the mixture model formula
  kappa_nts <- paste0("kappa", 3:(max_set_size + 1))
  theta_nts <- paste0("theta", 3:(max_set_size + 1))
  mu_nts <- paste0("mu", 3:(max_set_size + 1))

  for (i in 1:(max_set_size - 1)) {
    formula <- formula +
      glue_nlf("{kappa_nts[i]} ~ kappa") +
      glue_nlf(
        "{theta_nts[i]} ~ {lure_idx[i]} * (thetant + log(inv_ss))",
        " + (1 - {lure_idx[i]}) * (-100)"
      ) +
      glue_nlf("{mu_nts[i]} ~ {nt_features[i]}")
  }

  # define mixture family
  formula$family <- brms::mixture(
    brms::von_mises("tan_half"), brms::von_mises("identity"),
    nmix = c(1, max_set_size),
    order = "none"
  )

  nlist(formula, data)
}

#' @export
configure_model.mixture3p_cd <- function(model, data, formula) {
  mixture3p_cd <- brms::custom_family(
    name = "mixture3p_cd",
    dpars = c("mu", "kappa", "thetat", "thetant", "beta"),
    links = c("tan_half", "log", "identity", "identity", "identity"),
    lb = c(NA, 0, NA, NA, NA),
    ub = c(NA, NA, NA, NA, NA),
    type = "int",
    loop = FALSE,
    log_lik = log_lik_mixture3p_cd,
    posterior_predict = posterior_predict_mixture3p_cd
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- .generate_mixture3p_cd_stan()
  stan_tdata <- read_lines2(paste0(sc_path, "/mixture3p_cd_tdata.stan"))
  stan_likelihood <- read_lines2(paste0(sc_path, "/mixture3p_cd_likelihood.stan"))
  stanvars <- brms::stanvar(x = data$probe_centered, name = "probe_cd") +
    brms::stanvar(x = attr(data, "cd_nt_features_matrix"), name = "cd_nt_features") +
    brms::stanvar(x = attr(data, "cd_lure_idx_matrix"), name = "cd_lure_idx") +
    brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(scode = stan_likelihood, block = "likelihood", position = "end")

  formula <- bmf2bf(model, formula)
  formula$family <- mixture3p_cd

  nlist(formula, data, stanvars)
}

#' @export
bmf2bf.mixture3p_cd <- function(model, formula = bmmformula()) {
  resp_name <- model$resp_vars$response
  mu_rhs <- .extract_mu_rhs(formula)
  brms::bf(glue("{resp_name} ~ {mu_rhs}"))
}

.generate_mixture3p_cd_stan <- function() {
  "
  real mixture3p_cd_lpmf(array[] int y, vector mu, vector kappa,
                         vector thetat, vector thetant, vector beta) {
    return 0;
  }

  real mixture3p_cd_log_prob(int y, real probe, real mu, real kappa,
                             real thetat, real thetant, real beta,
                             row_vector nt_features, row_vector lure_idx,
                             vector grid, real dx) {
    int n_quad = size(grid);
    int n_nt = cols(nt_features);
    int n_active = 0;
    real log_uniform = -log(2 * pi());
    real log_vm_norm = log(2 * pi()) + log_modified_bessel_first_kind(0, kappa);
    real log_thetat;
    real log_thetant;
    real log_pguess;
    vector[n_quad] log_p_ret;
    vector[n_quad] log_p_same;
    vector[n_quad] log_integrand;
    real log_p_change;

    for (k in 1:n_nt) {
      if (lure_idx[k] > 0.5) {
        n_active += 1;
      }
    }

    {
      real log_Z = log_sum_exp(log_sum_exp(thetat, thetant), 0);
      log_thetat = thetat - log_Z;
      log_thetant = thetant - log_Z;
      log_pguess = -log_Z;
    }

    for (j in 1:n_quad) {
      real x = grid[j];
      real vm_ret = kappa * cos(x - mu) - log_vm_norm;
      real vm_same = kappa * cos(x - probe - mu) - log_vm_norm;

      log_p_ret[j] = log_sum_exp(log_thetat + vm_ret, log_pguess + log_uniform);
      log_p_same[j] = log_sum_exp(log_thetat + vm_same, log_pguess + log_uniform);

      if (n_active > 0) {
        for (k in 1:n_nt) {
          if (lure_idx[k] > 0.5) {
            real log_nt = log_thetant - log(n_active) +
              kappa * cos(x - nt_features[k]) - log_vm_norm;
            real log_nt_same = log_thetant - log(n_active) +
              kappa * cos(x - nt_features[k] - probe) - log_vm_norm;
            log_p_ret[j] = log_sum_exp(log_p_ret[j], log_nt);
            log_p_same[j] = log_sum_exp(log_p_same[j], log_nt_same);
          }
        }
      }
    }

    for (j in 1:n_quad) {
      log_integrand[j] = log_inv_logit(5 * (log_p_ret[j] - log_p_same[j] - beta)) +
        log_p_ret[j];
    }

    log_p_change = log_sum_exp(log_integrand) + log(dx);
    log_p_change = fmin(fmax(log_p_change, log(1e-10)), log1m(1e-10));
    if (y == 1) return log_p_change;
    return log1m_exp(log_p_change);
  }
  "
}

log_lik_mixture3p_cd <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  thetat <- brms::get_dpar(prep, "thetat", i = i)
  thetant <- brms::get_dpar(prep, "thetant", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- .extract_cd_probe(i, prep)
  y <- prep$data$Y[i]
  nt_data <- .extract_cd_nt_data(i, prep)

  dmixture3p_cd(y, probe, nt_features = nt_data$nt_features,
                lure_idx = nt_data$lure_idx, kappa = kappa, thetat = thetat,
                thetant = thetant, beta = beta, mu = mu, log = TRUE)
}

posterior_predict_mixture3p_cd <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  thetat <- brms::get_dpar(prep, "thetat", i = i)
  thetant <- brms::get_dpar(prep, "thetant", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- .extract_cd_probe(i, prep)
  nt_data <- .extract_cd_nt_data(i, prep)

  rmixture3p_cd(length(kappa), probe, nt_features = nt_data$nt_features,
                lure_idx = nt_data$lure_idx, kappa = kappa, thetat = thetat,
                thetant = thetant, beta = beta, mu = mu)
}

.configure_prior_mixture3p <- function(model, data, formula) {
  prior <- brms::empty_prior()
  set_size_var <- model$other_vars$set_size
  set_size_is_factor_with_level1 <- any(data$ss_numeric == 1) && !is.numeric(data[[set_size_var]])
  if (!set_size_is_factor_with_level1) return(prior)
  prior +
    constrain_set_size1_fixef(formula, "thetant", set_size_var, "constant(-100)") +
    constrain_set_size1_ranef(formula, "thetant", set_size_var, "constant(1e-8)")
}

#' @export
configure_prior.mixture3p_cd <- function(model, data, formula, user_prior, ...) {
  .configure_prior_mixture3p(model, data, formula)
}

#' @export
configure_prior.mixture3p_de <- function(model, data, formula, user_prior, ...) {
  .configure_prior_mixture3p(model, data, formula)
}
