############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_imm <- function(resp_error = NULL, nt_features = NULL, nt_distances = NULL,
                       set_size = NULL, response = NULL, probe = NULL,
                       target = NULL, regex = FALSE, version = "full",
                       task = "de", links = NULL, call = NULL, ...) {
  domain_class <- if (task == "cd") "change_detection" else "circular"

  out <- structure(
    list(
      resp_vars = if (task == "cd") nlist(response, probe, target) else nlist(resp_error),
      other_vars = nlist(nt_features, nt_distances, set_size),
      domain = "Visual working memory",
      task = if (task == "cd") "Change detection" else "Continuous reproduction",
      name = "Interference measurement model by Oberauer and Lin (2017).",
      version = version,
      citation = glue(
        "Oberauer, K., & Lin, H.Y. (2017). An interference model \\
          of visual working memory. Psychological Review, 124(1), 21-59"
      ),
      parameters = list(
        kappa = "Concentration parameter of the von Mises distribution",
        a = "General activation of memory items",
        c = "Context activation",
        s = "Spatial similarity gradient"
      ),
      links = list(kappa = "log", a = "log", c = "log", s = "log"),
      default_priors = list(
        kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
        a = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
        c = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
        s = list(main = "normal(0, 1)", effects = "normal(0, 1)")
      ),
      fixed_parameters = list(),
      void_mu = task == "cd"
    ),
    regex = regex,
    regex_vars = c("nt_features", "nt_distances"),
    class = c("bmmodel", domain_class, "non_targets", "imm",
              paste0("imm_", version), paste0("imm_", version, "_", task)),
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
        - Non-target features should be in radians and centered relative to the target
        - Non-target distances should be positive"
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
        "Location parameter of the von Mises distribution for memory \\
          responses (in radians). Fixed internally to 0 by default."
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

  if (version == "abc") {
    out$parameters$s <- NULL
    out$links$s <- NULL
    out$default_priors$s <- NULL
    attributes(out)$regex_vars <- c("nt_features")
    if (task == "de") out$other_vars$nt_distances <- NULL
  } else if (version == "bsc") {
    out$parameters$a <- NULL
    out$links$a <- NULL
    out$default_priors$a <- NULL
  }

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
#' @param response The name of the variable in the dataset containing the binary
#'   response (0 = "same", 1 = "change"). Required when `task = "cd"`.
#' @param probe The name of the variable containing the probe color in radians.
#'   Required when `task = "cd"`.
#' @param target The name of the variable containing the target color in
#'   radians. Required when `task = "cd"`.
#' @param task Character. The experimental task: `"de"` for delayed estimation
#'   (continuous reproduction) or `"cd"` for change detection. Default is `"de"`.
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
imm <- function(resp_error = NULL, nt_features = NULL, nt_distances = NULL,
                set_size = NULL, response = NULL, probe = NULL, target = NULL,
                regex = FALSE, version = "full", task = "de", ...) {
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
    if (version != "abc") {
      stopif(is.null(nt_distances), "Argument 'nt_distances' is required for version = '{version}'.")
    } else {
      nt_distances <- NULL
    }
  } else {
    warning2("`imm(task = \"cd\")` is deprecated. Please use `imm_cd()` instead.")
    stopif(is.null(response), "Argument 'response' is required for task = 'cd'.")
    stopif(is.null(probe), "Argument 'probe' is required for task = 'cd'.")
    stopif(is.null(target), "Argument 'target' is required for task = 'cd'.")
    stopif(is.null(nt_features), "Argument 'nt_features' is required.")
    stopif(is.null(set_size), "Argument 'set_size' is required.")
    if (version != "abc") {
      stopif(is.null(nt_distances), "Argument 'nt_distances' is required for version = '{version}'.")
    } else {
      nt_distances <- NULL
    }
  }

  .model_imm(
    resp_error = resp_error, nt_features = nt_features,
    nt_distances = nt_distances, set_size = set_size,
    response = response, probe = probe, target = target,
    regex = regex, version = version, task = task, call = call, ...
  )
}

#' @rdname imm
#' @export
imm_cd <- function(response, probe, target, nt_features, nt_distances, set_size,
                   regex = FALSE, version = "full", ...) {
  call <- match.call()
  stop_missing_args()
  stopif(version != "full", "Only version = 'full' is currently supported by `imm_cd()`.")
  .model_imm(
    nt_features = nt_features,
    nt_distances = nt_distances,
    set_size = set_size,
    response = response,
    probe = probe,
    target = target,
    regex = regex,
    version = version,
    task = "cd",
    call = call,
    ...
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

  attr(data, "cd_nt_distances_matrix") <- data.matrix(data[, nt_distances, drop = FALSE])

  data
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

.base_imm_formula <- function(model, formula) {
  bmf2bf(model, formula) +
    brms::lf(kappa2 ~ 1) +
    brms::lf(mu2 ~ 1) +
    brms::nlf(kappa1 ~ kappa)
}

.add_imm_mixture_terms <- function(formula, max_set_size, lure_idx, nt_features, theta_exprs) {
  kappa_nts <- paste0("kappa", 3:(max_set_size + 1))
  theta_nts <- paste0("theta", 3:(max_set_size + 1))
  mu_nts <- paste0("mu", 3:(max_set_size + 1))

  for (i in 1:(max_set_size - 1)) {
    formula <- formula +
      glue_nlf("{kappa_nts[i]} ~ kappa") +
      glue_nlf("{theta_nts[i]} ~ {theta_exprs[i]}") +
      glue_nlf("{mu_nts[i]} ~ {nt_features[i]}")
  }

  formula
}

.imm_mixture_family <- function(max_set_size) {
  brms::mixture(
    brms::von_mises("tan_half"), brms::von_mises("identity"),
    nmix = c(1, max_set_size),
    order = "none"
  )
}

#' @export
configure_model.imm_abc_de <- function(model, data, formula) {
  # retrieve arguments from the data check
  max_set_size <- attr(data, "max_set_size")
  lure_idx <- attr(data, "lure_idx_vars")
  nt_features <- model$other_vars$nt_features

  formula <- .base_imm_formula(model, formula) +
    brms::nlf(theta1 ~ log(exp(c) + exp(a)))

  theta_exprs <- glue("{lure_idx} * a + (1 - {lure_idx}) * (-100)")
  formula <- .add_imm_mixture_terms(formula, max_set_size, lure_idx, nt_features, theta_exprs)
  formula$family <- .imm_mixture_family(max_set_size)

  nlist(formula, data)
}

#' @export
configure_model.imm_bsc_de <- function(model, data, formula) {
  # retrieve arguments from the data check
  max_set_size <- attr(data, "max_set_size")
  lure_idx <- attr(data, "lure_idx_vars")
  nt_features <- model$other_vars$nt_features
  nt_distances <- model$other_vars$nt_distances

  formula <- .base_imm_formula(model, formula) +
    brms::nlf(theta1 ~ c) +
    brms::nlf(expS ~ exp(s))

  theta_exprs <- glue("{lure_idx} * (-expS*{nt_distances} + c) + (1 - {lure_idx}) * (-100)")
  formula <- .add_imm_mixture_terms(formula, max_set_size, lure_idx, nt_features, theta_exprs)
  formula$family <- .imm_mixture_family(max_set_size)

  nlist(formula, data)
}

#' @export
configure_model.imm_full_de <- function(model, data, formula) {
  # retrieve arguments from the data check
  max_set_size <- attr(data, "max_set_size")
  lure_idx <- attr(data, "lure_idx_vars")
  nt_features <- model$other_vars$nt_features
  nt_distances <- model$other_vars$nt_distances

  formula <- .base_imm_formula(model, formula) +
    brms::nlf(theta1 ~ log(exp(c) + exp(a))) +
    brms::nlf(expS ~ exp(s))

  theta_exprs <- glue("{lure_idx} * log(exp(c-expS*{nt_distances}) + exp(a)) + (1 - {lure_idx}) * (-100)")
  formula <- .add_imm_mixture_terms(formula, max_set_size, lure_idx, nt_features, theta_exprs)
  formula$family <- .imm_mixture_family(max_set_size)

  nlist(formula, data)
}

#' @export
configure_model.imm_full_cd <- function(model, data, formula) {
  imm_full_cd <- brms::custom_family(
    name = "imm_full_cd",
    dpars = c("mu", "kappa", "c", "a", "s", "beta"),
    links = c("tan_half", "log", "log", "log", "log", "identity"),
    lb = c(NA, 0, NA, NA, NA, NA),
    ub = c(NA, NA, NA, NA, NA, NA),
    type = "int",
    loop = FALSE,
    log_lik = log_lik_imm_full_cd,
    posterior_predict = posterior_predict_imm_full_cd
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- .generate_imm_full_cd_stan()
  stan_tdata <- read_lines2(paste0(sc_path, "/imm_full_cd_tdata.stan"))
  stan_likelihood <- read_lines2(paste0(sc_path, "/imm_full_cd_likelihood.stan"))
  stanvars <- brms::stanvar(x = data$probe_centered, name = "probe_cd") +
    brms::stanvar(x = attr(data, "cd_nt_features_matrix"), name = "cd_nt_features") +
    brms::stanvar(x = attr(data, "cd_nt_distances_matrix"), name = "cd_nt_distances") +
    brms::stanvar(x = attr(data, "cd_lure_idx_matrix"), name = "cd_lure_idx") +
    brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(scode = stan_likelihood, block = "likelihood", position = "end")

  formula <- bmf2bf(model, formula)
  formula$family <- imm_full_cd

  nlist(formula, data, stanvars)
}

#' @export
bmf2bf.imm_full_cd <- function(model, formula = bmmformula()) {
  resp_name <- model$resp_vars$response
  mu_rhs <- .extract_mu_rhs(formula)
  brms::bf(glue("{resp_name} ~ {mu_rhs}"))
}

.generate_imm_full_cd_stan <- function() {
  "
  real imm_full_cd_lpmf(array[] int y, vector mu, vector kappa, vector c_par,
                        vector a, vector s, vector beta) {
    return 0;
  }

  real imm_full_cd_log_prob(int y, real probe, real mu, real kappa, real c_par,
                            real a, real s, real beta, row_vector nt_features,
                            row_vector nt_distances, row_vector lure_idx,
                            vector grid, real dx) {
    int n_quad = size(grid);
    int n_nt = cols(nt_features);
    int n_active = 0;
    real log_uniform = -log(2 * pi());
    real log_vm_norm = log(2 * pi()) + log_modified_bessel_first_kind(0, kappa);
    real log_w_bg = 0;
    real log_w_target = log_sum_exp(c_par, a);
    vector[n_quad] log_p_ret;
    vector[n_quad] log_p_same;
    vector[n_nt] log_w_nt;
    vector[n_quad] log_integrand;
    real log_total_weight;
    real log_p_change;

    for (k in 1:n_nt) {
      if (lure_idx[k] > 0.5) {
        n_active += 1;
        log_w_nt[k] = log_sum_exp(c_par - exp(s) * nt_distances[k], a);
      } else {
        log_w_nt[k] = negative_infinity();
      }
    }

    log_total_weight = log_sum_exp(log_w_bg, log_w_target);
    if (n_active > 0) {
      for (k in 1:n_nt) {
        if (lure_idx[k] > 0.5) {
          log_total_weight = log_sum_exp(log_total_weight, log_w_nt[k]);
        }
      }
    }

    for (j in 1:n_quad) {
      real x = grid[j];
      real vm_ret = kappa * cos(x - mu) - log_vm_norm;
      real vm_same = kappa * cos(x - probe - mu) - log_vm_norm;

      log_p_ret[j] = log_sum_exp(log_w_target + vm_ret, log_w_bg + log_uniform);
      log_p_same[j] = log_sum_exp(log_w_target + vm_same, log_w_bg + log_uniform);

      if (n_active > 0) {
        for (k in 1:n_nt) {
          if (lure_idx[k] > 0.5) {
            real log_nt = log_w_nt[k] + kappa * cos(x - nt_features[k]) - log_vm_norm;
            real log_nt_same = log_w_nt[k] + kappa * cos(x - nt_features[k] - probe) - log_vm_norm;
            log_p_ret[j] = log_sum_exp(log_p_ret[j], log_nt);
            log_p_same[j] = log_sum_exp(log_p_same[j], log_nt_same);
          }
        }
      }

      log_p_ret[j] -= log_total_weight;
      log_p_same[j] -= log_total_weight;
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

log_lik_imm_full_cd <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  c_par <- brms::get_dpar(prep, "c", i = i)
  a <- brms::get_dpar(prep, "a", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- .extract_cd_probe(i, prep)
  y <- prep$data$Y[i]
  nt_data <- .extract_cd_nt_data(i, prep, has_distances = TRUE)

  dimm_cd(y, probe, nt_features = nt_data$nt_features,
          nt_distances = nt_data$nt_distances, lure_idx = nt_data$lure_idx,
          kappa = kappa, c = c_par, a = a, s = s, beta = beta, mu = mu,
          log = TRUE)
}

posterior_predict_imm_full_cd <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  c_par <- brms::get_dpar(prep, "c", i = i)
  a <- brms::get_dpar(prep, "a", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- .extract_cd_probe(i, prep)
  nt_data <- .extract_cd_nt_data(i, prep, has_distances = TRUE)

  rimm_cd(length(kappa), probe, nt_features = nt_data$nt_features,
          nt_distances = nt_data$nt_distances, lure_idx = nt_data$lure_idx,
          kappa = kappa, c = c_par, a = a, s = s, beta = beta, mu = mu)
}

#' @export
configure_prior.imm_full_cd <- function(model, data, formula, user_prior, ...) {
  .configure_prior_imm(model, data, formula, nlpars = c("a", "s"))
}

############################################################################# !
# CONFIGURE_PRIOR METHODS                                                ####
############################################################################# !

#' @export
configure_prior.imm_abc_de <- function(model, data, formula, user_prior, ...) {
  .configure_prior_imm(model, data, formula, nlpars = "a")
}

#' @export
configure_prior.imm_bsc_de <- function(model, data, formula, user_prior, ...) {
  .configure_prior_imm(model, data, formula, nlpars = "s")
}

#' @export
configure_prior.imm_full_de <- function(model, data, formula, user_prior, ...) {
  .configure_prior_imm(model, data, formula, nlpars = c("a", "s"))
}

.configure_prior_imm <- function(model, data, formula, nlpars, ...) {
  prior <- brms::empty_prior()
  set_size_var <- model$other_vars$set_size

  # Models with non-target errors need constant priors on set-size 1 factors
  set_size_is_factor_with_level1 <- any(data$ss_numeric == 1) && !is.numeric(data[[set_size_var]])
  if (!set_size_is_factor_with_level1) {
    return(prior)
  }

  prior +
    constrain_set_size1_fixef(formula, nlpars, set_size_var, "constant(0)") +
    constrain_set_size1_ranef(formula, nlpars, set_size_var, "constant(1e-8)")
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
