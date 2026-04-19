############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdm <- function(resp_error = NULL, nt_features = NULL, nt_distances = NULL,
                       set_size = NULL, regex = FALSE, links = NULL,
                       version = "simple", call = NULL, ...) {
  version <- match.arg(version, c("simple", "full", "bsc", "abc"))
  simple_version <- version == "simple"
  uses_distances <- version %in% c("full", "bsc")
  regex_vars <- if (simple_version) {
    character(0)
  } else if (uses_distances) {
    c("nt_features", "nt_distances")
  } else {
    c("nt_features")
  }
  classes <- if (simple_version) {
    c("bmmodel", "circular", "sdm", paste0("sdm_", version))
  } else {
    c("bmmodel", "circular", "non_targets", "sdm", paste0("sdm_", version))
  }
  requirements <- if (simple_version) {
    glue(
      "- The response variable should be in radians and represent the angular \\
      error relative to the target"
    )
  } else if (uses_distances) {
    glue(
      "- The response variable should be in radians and represent the angular \\
      error relative to the target
      - The non-target features should be in radians and be centered relative \\
      to the target
      - The non-target distances should be positive and match the set size"
    )
  } else {
    glue(
      "- The response variable should be in radians and represent the angular \\
      error relative to the target
      - The non-target features should be in radians and be centered relative \\
      to the target"
    )
  }

  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = if (simple_version) nlist() else nlist(nt_features, nt_distances, set_size),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Signal Discrimination Model (SDM) by Oberauer (2023)",
      citation = glue(
        "Oberauer, K. (2023). Measurement models for visual working memory - \\
        A factorial model comparison. Psychological Review, 130(3), 841-852"
      ),
      version = version,
      requirements = requirements,
      parameters = list(
        mu = glue("Location parameter of the SDM distribution (in radians; \\
                  by default fixed internally to 0)"),
        c = "Cue-dependent activation of the target representation",
        kappa = "Precision parameter of the SDM distribution",
        a = "Cue-independent activation of encoded items",
        s = "Spatial similarity gradient"
      ),
      links = list(
        mu = "tan_half",
        c = "log",
        kappa = "log",
        a = "log",
        s = "log"
      ),
      fixed_parameters = list(mu = 0),
      default_priors = list(
        mu = list(main = "student_t(1, 0, 1)"),
        kappa = list(main = "student_t(5, 1.75, 0.75)", effects = "normal(0, 1)"),
        c = list(main = "student_t(5, 2, 0.75)", effects = "normal(0, 1)"),
        a = list(main = "normal(0, 1)", effects = "normal(0, 1)"),
        s = list(main = "normal(0, 1)", effects = "normal(0, 1)")
      ),
      init_ranges = list(
        mu = c(-0.5, 0.5),
        kappa = c(2.5, 3.5),
        c = c(4, 6),
        a = c(0.5, 1.5),
        s = c(0.5, 1.5)
      ),
      void_mu = FALSE
    ),
    regex = regex && !simple_version,
    regex_vars = regex_vars,
    class = classes,
    call = call
  )

  if (simple_version) {
    out$parameters$a <- NULL
    out$parameters$s <- NULL
    out$default_priors$a <- NULL
    out$default_priors$s <- NULL
    out$links$a <- NULL
    out$links$s <- NULL
    out$init_ranges$a <- NULL
    out$init_ranges$s <- NULL
  } else if (version == "abc") {
    out$parameters$s <- NULL
    out$default_priors$s <- NULL
    out$links$s <- NULL
    out$init_ranges$s <- NULL
  } else if (version == "bsc") {
    out$parameters$a <- NULL
    out$default_priors$a <- NULL
    out$links$a <- NULL
    out$init_ranges$a <- NULL
  }

  out$links[names(links)] <- links
  out
}

# user facing alias

#' @title `r .model_sdm()$name`
#' @name sdm
#' @details see [the online article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for a detailed description of the model
#'   and how to use it. `r model_info(.model_sdm())`
#' @param resp_error The name of the variable in the dataset containing the
#'   response error. The response error should code the response relative to the
#'   to-be-recalled target in radians. You can transform the response error in
#'   degrees to radians using the `deg2rad` function.
#' @param nt_features A character vector with the names of the non-target
#'   variables. The non-target variables should be in radians and centered
#'   relative to the target. Only required for the `abc`, `bsc`, and `full`
#'   versions.
#' @param nt_distances A character vector with the names of the columns
#'   containing the distances of non-target items to the cued target. Only
#'   required for the `bsc` and `full` versions.
#' @param set_size Name of the column containing the set size variable (if
#'   set_size varies) or a numeric value for the set_size, if the set_size is
#'   fixed. Only required for the `abc`, `bsc`, and `full` versions.
#' @param regex Logical. If TRUE, the `nt_features` and `nt_distances`
#'   arguments are interpreted as regular expressions to match columns in the
#'   dataset.
#' @param version Character. The version of the model to use. One of
#'   `"simple"`, `"abc"`, `"bsc"`, or `"full"`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @export
#' @keywords bmmodel
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simulate data from the simple SDM
#' dat <- data.frame(y = rsdm(n = 1000, c = 4, kappa = 3))
#'
#' # specify formula
#' ff <- bmf(
#'   c ~ 1,
#'   kappa ~ 1
#' )
#'
#' # specify the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = sdm(resp_error = "y"),
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
sdm <- function(resp_error, nt_features = NULL, nt_distances = NULL, set_size = NULL,
                regex = FALSE, version = "simple", ...) {
  call <- match.call()
  dots <- list(...)
  version <- match.arg(version, c("simple", "full", "bsc", "abc"))

  if ("setsize" %in% names(dots)) {
    set_size <- dots$setsize
    warning2("The argument 'setsize' is deprecated. Please use 'set_size' instead.")
  }

  missing_args <- c()
  if (missing(resp_error)) missing_args <- c(missing_args, "resp_error")
  if (version != "simple" && missing(nt_features)) missing_args <- c(missing_args, "nt_features")
  if (version %in% c("full", "bsc") && missing(nt_distances)) {
    missing_args <- c(missing_args, "nt_distances")
  }
  if (version != "simple" && missing(set_size)) missing_args <- c(missing_args, "set_size")
  stopif(
    length(missing_args) > 0,
    "The following required arguments are missing in sdm(): \\
    {paste(missing_args, collapse = ', ')}"
  )

  if (version == "simple") {
    nt_features <- NULL
    nt_distances <- NULL
    set_size <- NULL
  } else if (version == "abc") {
    nt_distances <- NULL
  }

  .model_sdm(
    resp_error = resp_error,
    nt_features = nt_features,
    nt_distances = nt_distances,
    set_size = set_size,
    regex = regex,
    version = version,
    call = call,
    ...
  )
}

#' @rdname sdm
#' @keywords deprecated
#' @export
sdmSimple <- function(resp_error, version = "simple", ...) {
  warning2("The function `sdmSimple()` is deprecated. Please use `sdm()` instead.")
  call <- match.call()
  dots <- list(...)
  if ("setsize" %in% names(dots)) {
    dots$set_size <- dots$setsize
    dots$setsize <- NULL
  }
  .model_sdm(resp_error = resp_error, version = version, call = call, ...)
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdm_simple <- function(model, data, formula) {
  # data sorted by predictors is necessary for speedy computation of normalizing constant
  data <- order_data_query(model, data, formula)
  NextMethod("check_data")
}

#' @export
check_data.sdm_bsc <- function(model, data, formula) {
  data <- .check_data_sdm_dist(model, data, formula)
  NextMethod("check_data")
}

#' @export
check_data.sdm_full <- function(model, data, formula) {
  data <- .check_data_sdm_dist(model, data, formula)
  NextMethod("check_data")
}

.check_data_sdm_dist <- function(model, data, formula) {
  nt_distances <- model$other_vars$nt_distances
  max_set_size <- attr(data, "max_set_size")

  stopif(
    !isTRUE(all.equal(length(nt_distances), max_set_size - 1)),
    "The number of columns for non-target distances in the argument \\
    'nt_distances' should equal max(set_size)-1"
  )

  data[, nt_distances][is.na(data[, nt_distances])] <- 999

  stopif(
    any(data[, nt_distances] < 0),
    "All non-target distances to the target need to be postive."
  )

  data
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

.sdm_extra_vars <- function(model, data) {
  c(
    model$other_vars$nt_features,
    model$other_vars$nt_distances,
    attr(data, "lure_idx_vars")
  )
}

.sdm_response_formula <- function(model, formula, data, include_distances = FALSE) {
  response <- model$resp_vars[[1]]
  extra_vars <- c(
    model$other_vars$nt_features,
    if (include_distances) model$other_vars$nt_distances,
    attr(data, "lure_idx_vars")
  )
  response_formula <- if (length(extra_vars) == 0) {
    brms::bf(glue("{response} ~ 1"))
  } else {
    brms::bf(glue("{response} | vreal({paste(extra_vars, collapse = ', ')}) ~ 1"))
  }

  components <- lapply(formula, function(x) if (is_nl(x)) brms::nlf(x) else brms::lf(x))
  Reduce(`+`, components, init = response_formula)
}

.sdm_lpdf_signature <- function(n_vreal) {
  if (n_vreal == 0) {
    return("")
  }

  paste0(", array[] real vreal", seq_len(n_vreal), collapse = "")
}

.sdm_lpdf_body <- function(version, n_lures) {
  target_weight <- switch(version,
    abc = "exp(c[n]) + exp(a[n])",
    bsc = "exp(c[n])",
    full = "exp(c[n]) + exp(a[n])"
  )

  lines <- c(
    paste0("      array[", n_lures + 1, "] real item_angles;"),
    paste0("      array[", n_lures + 1, "] real item_weights;"),
    "      int J = 1;",
    "      real activation_y;",
    "      item_angles[1] = 0;",
    paste0("      item_weights[1] = ", target_weight, ";")
  )

  lure_lines <- lapply(seq_len(n_lures), function(k) {
    if (version == "abc") {
      lure_idx_pos <- n_lures + k
      c(
        paste0("      if (vreal", lure_idx_pos, "[n] > 0.5) {"),
        "        J += 1;",
        paste0("        item_angles[J] = vreal", k, "[n];"),
        "        item_weights[J] = exp(a[n]);",
        "      }"
      )
    } else if (version == "bsc") {
      dist_pos <- n_lures + k
      lure_idx_pos <- 2 * n_lures + k
      c(
        paste0("      if (vreal", lure_idx_pos, "[n] > 0.5) {"),
        "        J += 1;",
        paste0("        item_angles[J] = vreal", k, "[n];"),
        paste0("        item_weights[J] = exp(c[n] - s[n] * vreal", dist_pos, "[n]);"),
        "      }"
      )
    } else {
      dist_pos <- n_lures + k
      lure_idx_pos <- 2 * n_lures + k
      c(
        paste0("      if (vreal", lure_idx_pos, "[n] > 0.5) {"),
        "        J += 1;",
        paste0("        item_angles[J] = vreal", k, "[n];"),
        paste0("        item_weights[J] = exp(a[n]) + exp(c[n] - s[n] * vreal", dist_pos, "[n]);"),
        "      }"
      )
    }
  })

  c(
    lines,
    unlist(lure_lines, use.names = FALSE),
    "      activation_y = sdm_spectral_activation_obs(y[n], mu[n], kappa[n], item_angles, item_weights, J);",
    "      out += activation_y;"
  )
}

.sdm_lpdf_wrapper <- function(version, dpars, n_lures, include_distances = FALSE) {
  n_vreal <- n_lures * if (include_distances) 3 else 2
  dpar_signature <- paste0(", vector ", dpars, collapse = "")
  wrapper_name <- paste0("sdm_", version, "_lpdf")
  c(
    paste0("  real ", wrapper_name, "(vector y", dpar_signature, .sdm_lpdf_signature(n_vreal), ") {"),
    "    int N = num_elements(y);",
    "    real out = 0;",
    "    for (n in 1:N) {",
    .sdm_lpdf_body(version, n_lures),
    "    }",
    "    return out;",
    "  }"
  )
}

.sdm_general_stanvars <- function(model, data, version, include_distances = FALSE) {
  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_shared_funs <- read_lines2(paste0(sc_path, "/sdm_spectral_funs.stan"))
  stan_tdata <- read_lines2(paste0(sc_path, "/sdm_spectral_tdata.stan"))
  stan_likelihood <- read_lines2(paste0(sc_path, "/sdm_", version, "_likelihood.stan"))
  wrapper <- .sdm_lpdf_wrapper(
    version = version,
    dpars = switch(version,
      abc = c("mu", "c", "a", "kappa"),
      bsc = c("mu", "c", "s", "kappa"),
      full = c("mu", "c", "a", "s", "kappa")
    ),
    n_lures = length(model$other_vars$nt_features),
    include_distances = include_distances
  )

  stanvars <- brms::stanvar(
    scode = c(stan_shared_funs, "", wrapper),
    block = "functions"
  ) +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(scode = stan_likelihood, block = "likelihood", position = "end") +
    brms::stanvar(x = as.matrix(data[model$other_vars$nt_features]), name = "SDM_NT_FEATURES") +
    brms::stanvar(x = as.matrix(data[attr(data, "lure_idx_vars")]), name = "SDM_LURE_IDX")

  if (include_distances) {
    stanvars <- stanvars +
      brms::stanvar(
        x = as.matrix(data[model$other_vars$nt_distances]),
        name = "SDM_NT_DISTANCES"
      )
  }

  stanvars
}

.sdm_custom_family <- function(name, dpars, links, vars = NULL, log_lik_fun, posterior_predict_fun) {
  brms::custom_family(
    name = name,
    dpars = dpars,
    links = links,
    lb = rep(NA_real_, length(dpars)),
    ub = rep(NA_real_, length(dpars)),
    vars = vars,
    type = "real",
    loop = FALSE,
    log_lik = log_lik_fun,
    posterior_predict = posterior_predict_fun
  )
}

#' @export
configure_model.sdm_simple <- function(model, data, formula) {
  sdm_simple <- .sdm_custom_family(
    name = "sdm_simple",
    dpars = c("mu", "c", "kappa"),
    links = c("tan_half", "identity", "log"),
    log_lik_fun = log_lik_sdm_simple,
    posterior_predict_fun = posterior_predict_sdm_simple
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_shared_funs <- read_lines2(paste0(sc_path, "/sdm_spectral_funs.stan"))
  stan_funs <- read_lines2(paste0(sc_path, "/sdm_simple_funs.stan"))
  stan_tdata <- read_lines2(paste0(sc_path, "/sdm_spectral_tdata.stan"))
  stan_likelihood <- read_lines2(paste0(sc_path, "/sdm_simple_likelihood.stan"))
  stanvars <- brms::stanvar(scode = c(stan_shared_funs, "", stan_funs), block = "functions") +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(scode = stan_likelihood, block = "likelihood", position = "end")

  formula <- bmf2bf(model, formula)
  formula$family <- sdm_simple

  init <- 1

  nlist(formula, data, stanvars, init)
}

.configure_model_sdm_general <- function(model, data, formula, version, dpars, links,
                                         include_distances, log_lik_fun, posterior_predict_fun) {
  extra_vars <- c(
    model$other_vars$nt_features,
    if (include_distances) model$other_vars$nt_distances,
    attr(data, "lure_idx_vars")
  )

  formula <- .sdm_response_formula(model, formula, data, include_distances = include_distances)
  formula$family <- .sdm_custom_family(
    name = paste0("sdm_", version),
    dpars = dpars,
    links = links,
    vars = paste0("vreal", seq_along(extra_vars)),
    log_lik_fun = log_lik_fun,
    posterior_predict_fun = posterior_predict_fun
  )

  stanvars <- .sdm_general_stanvars(model, data, version, include_distances = include_distances)

  nlist(formula, data, stanvars, init = 1)
}

#' @export
configure_model.sdm_abc <- function(model, data, formula) {
  .configure_model_sdm_general(
    model = model,
    data = data,
    formula = formula,
    version = "abc",
    dpars = c("mu", "c", "a", "kappa"),
    links = c("tan_half", "identity", "identity", "log"),
    include_distances = FALSE,
    log_lik_fun = log_lik_sdm_abc,
    posterior_predict_fun = posterior_predict_sdm_abc
  )
}

#' @export
configure_model.sdm_bsc <- function(model, data, formula) {
  .configure_model_sdm_general(
    model = model,
    data = data,
    formula = formula,
    version = "bsc",
    dpars = c("mu", "c", "s", "kappa"),
    links = c("tan_half", "identity", "log", "log"),
    include_distances = TRUE,
    log_lik_fun = log_lik_sdm_bsc,
    posterior_predict_fun = posterior_predict_sdm_bsc
  )
}

#' @export
configure_model.sdm_full <- function(model, data, formula) {
  .configure_model_sdm_general(
    model = model,
    data = data,
    formula = formula,
    version = "full",
    dpars = c("mu", "c", "a", "s", "kappa"),
    links = c("tan_half", "identity", "identity", "log", "log"),
    include_distances = TRUE,
    log_lik_fun = log_lik_sdm_full,
    posterior_predict_fun = posterior_predict_sdm_full
  )
}

############################################################################# !
# CONFIGURE_PRIOR METHODS                                                ####
############################################################################# !

#' @export
configure_prior.sdm_abc <- function(model, data, formula, user_prior, ...) {
  .configure_prior_sdm(model, data, formula, dpars = "a")
}

#' @export
configure_prior.sdm_bsc <- function(model, data, formula, user_prior, ...) {
  .configure_prior_sdm(model, data, formula, dpars = "s")
}

#' @export
configure_prior.sdm_full <- function(model, data, formula, user_prior, ...) {
  .configure_prior_sdm(model, data, formula, dpars = c("a", "s"))
}

.sdm_constrain_set_size1_fixef <- function(formula, dpars, set_size_var, prior_value) {
  dpar_pforms <- formula$pforms[dpars]
  has_setsize <- vapply(dpar_pforms, function(x) set_size_var %in% rhs_vars(x), logical(1))

  if (!any(has_setsize)) {
    return(NULL)
  }

  brms::prior_(
    prior_value,
    class = "b",
    coef = paste0(set_size_var, 1),
    dpar = dpars[has_setsize]
  )
}

.sdm_constrain_set_size1_ranef <- function(formula, dpars, set_size_var, prior_value) {
  re_terms_list <- lapply(dpars, function(dpar) {
    bterms <- brms::brmsterms(formula$pforms[[dpar]])
    re <- bterms$dpars$mu$re
    if (NROW(re) > 0) {
      data.frame(dpar = dpar, group = re$group, form = I(re$form), stringsAsFactors = FALSE)
    }
  })

  re_terms_df <- do.call(rbind, re_terms_list)
  if (is.null(re_terms_df)) {
    return(NULL)
  }

  has_setsize <- vapply(re_terms_df$form, function(x) set_size_var %in% rhs_vars(x), logical(1))
  if (!any(has_setsize)) {
    return(NULL)
  }

  brms::prior_(
    prior_value,
    class = "sd",
    coef = paste0(set_size_var, 1),
    group = re_terms_df$group[has_setsize],
    dpar = re_terms_df$dpar[has_setsize]
  )
}

.configure_prior_sdm <- function(model, data, formula, dpars) {
  prior <- brms::empty_prior()
  set_size_var <- model$other_vars$set_size
  set_size_is_factor_with_level1 <- any(data$ss_numeric == 1) && !is.numeric(data[[set_size_var]])

  if (!set_size_is_factor_with_level1) {
    return(prior)
  }

  prior +
    .sdm_constrain_set_size1_fixef(formula, dpars, set_size_var, "constant(0)") +
    .sdm_constrain_set_size1_ranef(formula, dpars, set_size_var, "constant(1e-8)")
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

.sdm_manual_link_pars <- function(model) {
  intersect(names(model$parameters), c("c", "a"))
}

#' @export
postprocess_brm.sdm <- function(model, fit, ...) {
  for (par in .sdm_manual_link_pars(model)) {
    fit$family[[paste0("link_", par)]] <- "log"
    fit$formula$family[[paste0("link_", par)]] <- "log"
  }
  fit
}

#' @export
revert_postprocess_brm.sdm <- function(model, fit, ...) {
  for (par in .sdm_manual_link_pars(model)) {
    fit$family[[paste0("link_", par)]] <- "identity"
    fit$formula$family[[paste0("link_", par)]] <- "identity"
  }
  fit
}

.sdm_extract_prep_data <- function(i, prep, version) {
  vreal_names <- grep("^vreal", names(prep$data), value = TRUE)
  divisor <- if (version == "abc") 2 else 3
  stopif(
    length(vreal_names) %% divisor != 0,
    "Unexpected number of auxiliary variables found for SDM version '{version}'"
  )

  n_lures <- length(vreal_names) / divisor
  nt_features <- vapply(seq_len(n_lures), function(j) prep$data[[paste0("vreal", j)]][i], numeric(1))
  out <- list(nt_features = nt_features)

  if (version %in% c("full", "bsc")) {
    offset <- n_lures
    out$nt_distances <- vapply(
      seq_len(n_lures),
      function(j) prep$data[[paste0("vreal", offset + j)]][i],
      numeric(1)
    )
    offset <- offset + n_lures
  } else {
    offset <- n_lures
  }

  out$lure_idx <- vapply(
    seq_len(n_lures),
    function(j) prep$data[[paste0("vreal", offset + j)]][i],
    numeric(1)
  )
  out
}

.sdm_dpar_args <- function(i, prep, version) {
  args <- list(
    mu = brms::get_dpar(prep, "mu", i = i),
    c = brms::get_dpar(prep, "c", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    version = version
  )

  if (version %in% c("abc", "full")) {
    args$a <- brms::get_dpar(prep, "a", i = i)
  }
  if (version %in% c("bsc", "full")) {
    args$s <- brms::get_dpar(prep, "s", i = i)
  }
  if (version != "simple") {
    args <- c(args, .sdm_extract_prep_data(i, prep, version))
  }
  args
}

.log_lik_sdm <- function(i, prep, version) {
  args <- .sdm_dpar_args(i, prep, version)
  args$x <- prep$data$Y[i]
  args$log <- TRUE
  do.call(dsdm, args)
}

.posterior_predict_sdm <- function(i, prep, version) {
  args <- .sdm_dpar_args(i, prep, version)
  args$n <- length(args$mu)
  do.call(rsdm, args)
}

log_lik_sdm_simple <- function(i, prep) {
  .log_lik_sdm(i, prep, version = "simple")
}

posterior_predict_sdm_simple <- function(i, prep, ...) {
  .posterior_predict_sdm(i, prep, version = "simple")
}

log_lik_sdm_abc <- function(i, prep) {
  .log_lik_sdm(i, prep, version = "abc")
}

posterior_predict_sdm_abc <- function(i, prep, ...) {
  .posterior_predict_sdm(i, prep, version = "abc")
}

log_lik_sdm_bsc <- function(i, prep) {
  .log_lik_sdm(i, prep, version = "bsc")
}

posterior_predict_sdm_bsc <- function(i, prep, ...) {
  .posterior_predict_sdm(i, prep, version = "bsc")
}

log_lik_sdm_full <- function(i, prep) {
  .log_lik_sdm(i, prep, version = "full")
}

posterior_predict_sdm_full <- function(i, prep, ...) {
  .posterior_predict_sdm(i, prep, version = "full")
}
