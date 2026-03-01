############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_imm <- function(resp_error = NULL, nt_features = NULL, nt_distances = NULL,
                       set_size = NULL, response = NULL, probe = NULL,
                       target = NULL, regex = FALSE, version = "full",
                       task = "de", links = NULL, call = NULL, ...) {
  domain_class <- if (task == "cd") "change_detection" else "circular"

  if (task == "cd") {
    out <- structure(
      list(
        resp_vars = nlist(response, probe, target),
        other_vars = nlist(nt_features, nt_distances, set_size),
        domain = "Visual working memory",
        task = "Change detection",
        name = "Interference measurement model by Oberauer and Lin (2017).",
        version = version,
        citation = glue(
          "Oberauer, K., & Lin, H.Y. (2017). An interference model \\
            of visual working memory. Psychological Review, 124(1), 21-59; \\
            Lin, H.Y., & Oberauer, K. (2022). An interference model for visual \\
            working memory: Applications to the change detection task. \\
            Cognitive Psychology, 133, 101463."
        ),
        requirements = glue(
          "- response: Binary (0='same', 1='change')
            - probe: Probe color in radians
            - target: Target color in radians
            - Non-target features should be in radians and centered relative to the target
            - Non-target distances should be positive"
        ),
        parameters = list(
          kappa = "Concentration parameter of the von Mises distribution",
          a = "General activation of memory items",
          c = "Context activation",
          s = "Spatial similarity gradient",
          beta = glue(
            "Decision criterion (log prior odds). \\
            Fixed to 0 by default for unbiased decision."
          )
        ),
        links = list(
          kappa = "log",
          a = "log",
          c = "log",
          s = "log",
          beta = "identity"
        ),
        fixed_parameters = list(beta = 0),
        default_priors = list(
          kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
          a = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
          c = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
          s = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
          beta = list(main = "normal(0, 0.5)")
        ),
        void_mu = TRUE
      ),
      regex = regex,
      regex_vars = c("nt_features", "nt_distances"),
      class = c("bmmodel", domain_class, "non_targets", "imm",
                paste0("imm_", version), paste0("imm_", version, "_", task)),
      call = call
    )
  } else {
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
        parameters = list(
          mu1 = glue(
            "Location parameter of the von Mises distribution for memory \\
              responses (in radians). Fixed internally to 0 by default."
          ),
          kappa = "Concentration parameter of the von Mises distribution",
          a = "General activation of memory items",
          c = "Context activation",
          s = "Spatial similarity gradient"
        ),
        links = list(
          mu1 = "tan_half",
          kappa = "log",
          a = "log",
          c = "log",
          s = "log"
        ),
        fixed_parameters = list(mu1 = 0, mu2 = 0, kappa2 = -100),
        default_priors = list(
          mu1 = list(main = "student_t(1, 0, 1)"),
          kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
          a = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
          c = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
          s = list(main = "normal(0, 1)", effects = "normal(0, 1)")
        ),
        void_mu = FALSE
      ),
      regex = regex,
      regex_vars = c("nt_features", "nt_distances"),
      class = c("bmmodel", domain_class, "non_targets", "imm",
                paste0("imm_", version), paste0("imm_", version, "_", task)),
      call = call
    )
  }

  # version-specific parameter removal (applies to both DE and CD)
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
  max_set_size <- attr(data, "max_set_size")
  lure_idx <- attr(data, "lure_idx_vars")
  nt_features <- model$other_vars$nt_features
  nt_distances <- model$other_vars$nt_distances
  n_nt <- max_set_size - 1

  imm_full_cd <- brms::custom_family(
    name = "imm_full_cd",
    dpars = c("mu", "kappa", "c", "a", "s", "beta"),
    links = c("identity", "log", "log", "log", "log", "identity"),
    lb = c(NA, 0, NA, NA, NA, NA),
    ub = c(NA, NA, NA, NA, NA, NA),
    type = "int",
    loop = TRUE,
    vars = c("vreal1[n]",
             paste0("vreal", 2:(n_nt + 1), "[n]"),
             paste0("vreal", (n_nt + 2):(2 * n_nt + 1), "[n]"),
             paste0("vint", 1:n_nt, "[n]")),
    log_lik = log_lik_imm_full_cd,
    posterior_predict = posterior_predict_imm_full_cd
  )

  stan_funs <- .generate_imm_full_cd_stan(n_nt)
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  formula <- bmf2bf(model, formula)
  formula$family <- imm_full_cd

  nlist(formula, data, stanvars)
}

#' @export
bmf2bf.imm_full_cd <- function(model, formula = bmmformula()) {
  resp_name <- model$resp_vars$response
  nt_features <- model$other_vars$nt_features
  nt_distances <- model$other_vars$nt_distances
  lure_idx <- attr(model, "lure_idx_computed") %||%
    paste0("LureIdx", seq_along(nt_features))

  vreal_args <- paste(c("probe_centered", nt_features, nt_distances), collapse = ", ")
  vint_args <- paste(lure_idx, collapse = ", ")

  brms_formula <- brms::bf(
    glue("{resp_name} | vreal({vreal_args}) + vint({vint_args}) ~ 1")
  )
  components <- lapply(formula, function(x) {
    if (is_nl(x)) brms::nlf(x) else brms::lf(x)
  })
  Reduce(`+`, components, init = brms_formula)
}

.generate_imm_full_cd_stan <- function(n_nt) {
  nt_args <- paste0("real nt", seq_len(n_nt), collapse = ", ")
  dist_args <- paste0("real dist", seq_len(n_nt), collapse = ", ")
  lure_args <- paste0("int lure", seq_len(n_nt), collapse = ", ")

  nt_loop <- paste(vapply(seq_len(n_nt), function(i) {
    glue("
      if (lure{i} == 1) {{
        real w_nt{i} = exp(c_par - exp(s) * dist{i}) + exp(a);
        total_weight += w_nt{i};
        real log_nt{i} = log(w_nt{i}) + von_mises_lpdf(x | nt{i}, kappa);
        log_p_retrieve = log_sum_exp(log_p_retrieve, log_nt{i});
        real log_nt{i}_same = log(w_nt{i}) + von_mises_lpdf(x | nt{i} + probe, kappa);
        log_p_x_given_same = log_sum_exp(log_p_x_given_same, log_nt{i}_same);
      }}")
  }, character(1)), collapse = "\n")

  glue("
  #include 'fun_tan_half.stan'

  real imm_full_cd_lpmf(int y, real mu, real kappa, real c_par, real a, real s, real beta, real probe, {nt_args}, {dist_args}, {lure_args}) {{
    int n_quad = 101;
    real dx = 2 * pi() / (n_quad - 1);
    real p_change = 0;
    real log_uniform = -log(2 * pi());

    // background activation (b = 0, so exp(b) = 1)
    real w_bg = 1.0;
    // target activation
    real w_target = exp(c_par) + exp(a);

    for (i in 1:n_quad) {{
      real x = -pi() + (i - 1) * dx;
      real total_weight = w_target + w_bg;

      // Start with target retrieval density (unnormalized)
      real log_p_retrieve = log(w_target) + von_mises_lpdf(x | 0, kappa);

      // Same hypothesis: retrieval density centered at probe
      real log_p_x_given_same = log(w_target) + von_mises_lpdf(x | probe, kappa);

      // Add non-target components
      {nt_loop}

      // Add background uniform component
      log_p_retrieve = log_sum_exp(log_p_retrieve, log(w_bg) + log_uniform);
      log_p_x_given_same = log_sum_exp(log_p_x_given_same, log(w_bg) + log_uniform);

      // Normalize by total weight
      log_p_retrieve -= log(total_weight);
      log_p_x_given_same -= log(total_weight);

      // LLR: compare change vs same hypothesis
      real log_p_same = log_p_x_given_same + log_uniform;
      real log_p_change_hyp = log_p_retrieve + log_uniform;
      real llr = log_p_change_hyp - log_p_same;

      if (llr > beta) {{
        p_change += exp(log_p_retrieve) * dx;
      }}
    }}

    p_change = fmin(fmax(p_change, 1e-10), 1 - 1e-10);
    if (y == 1) return log(p_change);
    return log1m(p_change);
  }}
  ")
}

log_lik_imm_full_cd <- function(i, prep) {
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  c_par <- brms::get_dpar(prep, "c", i = i)
  a <- brms::get_dpar(prep, "a", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- prep$data$vreal1[i]
  y <- prep$data$Y[i]

  n_nt <- (sum(grepl("^vreal[0-9]+$", names(prep$data))) - 1) / 2
  nt_features <- vapply(seq_len(n_nt), function(j) prep$data[[paste0("vreal", j + 1)]][i], numeric(1))
  nt_distances <- vapply(seq_len(n_nt), function(j) prep$data[[paste0("vreal", n_nt + j + 1)]][i], numeric(1))
  lure_idx <- vapply(seq_len(n_nt), function(j) prep$data[[paste0("vint", j)]][i], numeric(1))

  dimm_cd(y, probe, nt_features = nt_features, nt_distances = nt_distances,
          lure_idx = lure_idx, kappa = kappa, c = c_par, a = a, s = s,
          beta = beta, log = TRUE)
}

posterior_predict_imm_full_cd <- function(i, prep, ...) {
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  c_par <- brms::get_dpar(prep, "c", i = i)
  a <- brms::get_dpar(prep, "a", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- prep$data$vreal1[i]

  n_nt <- (sum(grepl("^vreal[0-9]+$", names(prep$data))) - 1) / 2
  nt_features <- vapply(seq_len(n_nt), function(j) prep$data[[paste0("vreal", j + 1)]][i], numeric(1))
  nt_distances <- vapply(seq_len(n_nt), function(j) prep$data[[paste0("vreal", n_nt + j + 1)]][i], numeric(1))
  lure_idx <- vapply(seq_len(n_nt), function(j) prep$data[[paste0("vint", j)]][i], numeric(1))

  rimm_cd(length(kappa), probe, nt_features = nt_features,
          nt_distances = nt_distances, lure_idx = lure_idx,
          kappa = kappa, c = c_par, a = a, s = s, beta = beta)
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
