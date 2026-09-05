############################################################################# !
# MULTIVARIATE COMPOSITION HELPERS                                       ####
############################################################################# !

# Restrict a component's likelihood to its own rows of the stacked data by
# adding subset() to the response side of the formula. Editing the formula
# call directly keeps brms's attributes (e.g. m3's non-linear flag) and the
# environment intact, which a deparse/re-parse round trip would lose
add_subset <- function(bformula, subset_var) {
  f <- bformula$formula
  subset_call <- call("subset", as.name(subset_var))
  lhs <- f[[2]]
  f[[2]] <- if (is.call(lhs) && identical(lhs[[1]], as.name("|"))) {
    lhs[[3]] <- call("+", lhs[[3]], subset_call)
    lhs
  } else {
    call("|", lhs, subset_call)
  }
  bformula$formula <- f
  bformula
}

# mirror of how brms turns response variable names into multivariate resp
# suffixes (brms:::make_stan_names); priors and post-processing must refer to
# responses by the sanitized name
sanitize_resp_name <- function(x) {
  gsub("\\.|_", "", make.names(x, unique = TRUE))
}

# Stack the checked per-component datasets into one data.frame for brms.
# Subset indicator columns are set on every frame before stacking so they are
# never NA; missing columns are filled with typed NAs so that factor levels
# and matrix response columns (e.g. m3's Y) survive the stacking
stack_component_data <- function(data_list, subset_vars) {
  for (i in seq_along(data_list)) {
    for (j in seq_along(subset_vars)) {
      data_list[[i]][[subset_vars[j]]] <- i == j
    }
  }
  templates <- Reduce(
    function(acc, data) c(acc, as.list(data)[setdiff(names(data), names(acc))]),
    data_list,
    init = list()
  )
  filled <- lapply(data_list, function(data) {
    for (col in setdiff(names(templates), names(data))) {
      data[[col]] <- na_column_like(templates[[col]], nrow(data))
    }
    data[names(templates)]
  })
  stacked <- do.call(rbind, filled)
  rownames(stacked) <- NULL
  stacked
}

# NA column of length n with the same type, class, levels and (for matrix
# columns) dimensions as the template column
na_column_like <- function(template, n) {
  if (is.matrix(template)) {
    template[rep(NA_integer_, n), , drop = FALSE]
  } else {
    template[rep(NA_integer_, n)]
  }
}

# Combine the stanvars of all components. The same Stan code injected by two
# components (e.g. the same model on two tasks) must appear only once, since
# duplicate function definitions fail to compile
merge_stanvars <- function(stanvars_list) {
  stanvars_list <- stanvars_list[!vapply(stanvars_list, is.null, logical(1))]
  if (length(stanvars_list) == 0) {
    return(NULL)
  }
  merged <- Reduce(`+`, stanvars_list)
  keys <- vapply(
    merged,
    function(v) paste(v$block %||% "", v$position %||% "", v$scode %||% "", sep = "\r"),
    character(1)
  )
  structure(unclass(merged)[!duplicated(keys)], class = "stanvars")
}

############################################################################# !
# MULTIVARIATE MODEL CONFIGURATION                                       ####
############################################################################# !

# Runs the standard bmm pipeline for every component, restricts each
# component's likelihood to its own rows, and composes the results into the
# arguments of a single multivariate brms model. The `model` slot returned for
# postprocess_brm() is the list of components, classed so that the
# multivariate fit gets its own postprocessing
#' @export
configure_fit.mvbmmformula <- function(formula, data = NULL, model = NULL, prior = NULL,
                                       init = TRUE) {
  spec <- formula
  stopif(
    length(spec) < 2,
    "A multivariate bmm model requires at least two components. Combine \\
    components created by bmm_component() with `+`."
  )
  components <- lapply(spec, mv_configure_component)
  validate_component_resps(components)
  check_shared_random_effects(spec)
  warn_factor_level_mismatch(spec)

  formula <- Reduce(`+`, lapply(components, function(x) x$cfg$formula)) +
    brms::set_rescor(FALSE)
  data <- stack_component_data(
    lapply(components, function(x) x$data),
    vapply(components, function(x) x$subset_var, character(1))
  )
  data_names <- vapply(components, function(x) x$data_name %||% "data", character(1))
  attr(data, "data_name") <- paste(unique(data_names), collapse = " + ")
  stanvars <- merge_stanvars(lapply(components, function(x) x$cfg$stanvars))
  component_priors <- Reduce(
    combine_prior,
    lapply(components, function(x) x$prior),
    init = brms::empty_prior()
  )
  # the global prior is passed through as-is so users can set joint priors
  # such as prior(lkj(2), class = "cor"); component-targeted rows require an
  # explicit resp
  prior <- combine_prior(component_priors, prior)

  config_args <- nlist(formula, data, stanvars)
  if (init) {
    config_args$init <- mv_create_initfun(components, formula, data, prior, stanvars)
  }
  components <- lapply(components, function(x) {
    x[c("model", "user_formula", "resp", "resp_name", "subset_var", "data_name")]
  })
  nlist(
    config_args, prior,
    model = structure(components, class = "bmm_components"),
    user_formula = spec
  )
}

# The standard univariate pipeline applied to one component, plus the
# multivariate plumbing: subset() injection and resp-tagged priors
mv_configure_component <- function(component) {
  user_formula <- component$formula
  model <- check_model(component$model, component$data, user_formula)
  data <- check_data(model, component$data, user_formula)
  formula <- check_formula(model, data, user_formula)
  cfg <- configure_model(model, data, formula)

  blocks <- vapply(
    unclass(cfg$stanvars %||% list()),
    function(v) v$block %||% "",
    character(1)
  )
  bad_blocks <- setdiff(unique(blocks), "functions")
  stopif(
    length(bad_blocks) > 0,
    "The {class(model)[length(class(model))]} model cannot be used in a \\
    multivariate bmm model: it injects Stan code into the \\
    {collapse_comma(bad_blocks)} block(s), which would conflict with the \\
    other components."
  )

  resp <- all.vars(cfg$formula$formula[[2]])[1]
  # brms refers to responses by their sanitized name in multivariate models,
  # so priors must be tagged with the sanitized name to reach the component
  resp_name <- sanitize_resp_name(resp)
  prior <- configure_prior(model, data, cfg$formula, user_prior = component$prior)
  prior$resp <- rep(resp_name, nrow(prior))

  subset_var <- paste0(".subset_", resp)
  cfg$formula <- add_subset(cfg$formula, subset_var)

  nlist(
    model, user_formula, data, cfg, prior, resp, resp_name, subset_var,
    data_name = component$data_name
  )
}

validate_component_resps <- function(components) {
  resps <- vapply(components, function(x) x$resp, character(1))
  resp_names <- vapply(components, function(x) x$resp_name, character(1))
  duplicates <- unique(resps[duplicated(resp_names)])
  stopif(
    length(duplicates) > 0,
    "Response variables must be unique across the components of a \\
    multivariate bmm model, but the following are duplicated: \\
    {collapse_comma(duplicates)}. Note that fitting the same response \\
    variable in two components (e.g. two m3 components, which both use the \\
    internal response 'Y') is not yet supported, and that variable names are \\
    compared after removing '.' and '_' characters, as done by brms."
  )

  all_cols <- unique(unlist(lapply(components, function(x) names(x$data))))
  clashes <- intersect(vapply(components, function(x) x$subset_var, character(1)), all_cols)
  stopif(
    length(clashes) > 0,
    "The data column(s) {collapse_comma(clashes)} conflict with the subset \\
    indicators bmm creates internally for multivariate models. Please rename \\
    these columns."
  )
}

############################################################################# !
# INITIAL VALUES                                                          ####
############################################################################# !

# Multivariate counterpart of the univariate create_initfun() method. Stan
# parameters belonging to components with init_ranges are initialized from
# those ranges and the shared random-effects and correlation parameters
# generically; the remaining parameters are left to the backend, which accepts
# a partial init list
mv_create_initfun <- function(components, formula, data, prior, stanvars) {
  needs_init <- vapply(
    components,
    function(x) !is.null(x$model$init_ranges),
    logical(1)
  )
  if (!any(needs_init)) {
    return(1)
  }

  # the prior decides which parameters exist: parameters held constant by a
  # fixed_parameters prior are absent from the parameters block
  standata_list <- standata(formula, data, prior = prior, stanvars = stanvars)
  stan_code <- stancode(formula, data, prior = prior, stanvars = stanvars)
  stanpars_list <- extract_parameter_dimensions(extract_stan_blocks(stan_code)$parameters)
  par_matches <- mv_match_stan_parameters(names(stanpars_list), components[needs_init])

  function() {
    force(stanpars_list)
    force(standata_list)
    force(par_matches)
    force(formula)
    force(data)

    bterms <- brms::brmsterms(formula)
    inits <- list()
    for (spar in names(stanpars_list)) {
      spec <- stanpars_list[[spar]]
      type <- spec$type
      dim <- resolve_stan_dim(spec$dims, standata_list)
      match <- par_matches[[spar]]

      inits[[spar]] <- if (!is.null(match)) {
        resp_terms <- bterms$terms[[match$resp_name]]
        par_bterms <- resp_terms$dpars[[match$parameter]] %||%
          resp_terms$nlpars[[match$parameter]]
        switch(type,
          real = init_real_param(spar, match$range, match$link),
          vector = init_vector_param(spar, dim, match$range, match$link, par_bterms, data),
          NULL
        )
      } else {
        switch(type,
          vector = if (grepl("^sd_", spar)) array(runif(prod(dim), 0.05, 0.1), dim = dim),
          matrix = if (grepl("^z_", spar)) matrix(runif(prod(dim), -0.5, 0.5), nrow = dim[1]),
          cholesky_factor_corr = ,
          cholesky_factor_cov = ,
          cov_matrix = ,
          corr_matrix = diag(nrow = dim),
          NULL
        )
      }
    }
    inits[!vapply(inits, is.null, logical(1))]
  }
}

# brms names multivariate Stan parameters b_<resp>_<par> for non-linear
# parameters, and b_<par>_<resp> / Intercept_<par>_<resp> for distributional
# parameters; exact matching against these candidates avoids the substring
# collisions a grepl()-based lookup would produce across components
mv_match_stan_parameters <- function(spar_names, components) {
  matches <- list()
  for (component in components) {
    model <- component$model
    resp_name <- component$resp_name
    for (parameter in names(model$parameters)) {
      range <- model$init_ranges[[parameter]]
      if (is.null(range)) {
        next
      }
      candidates <- c(
        paste0("b_", resp_name, "_", parameter),
        paste0("b_", parameter, "_", resp_name),
        paste0("Intercept_", parameter, "_", resp_name)
      )
      for (spar in intersect(candidates, spar_names)) {
        matches[[spar]] <- nlist(
          resp_name, parameter, range,
          link = model$links[[parameter]]
        )
      }
    }
  }
  matches
}

# Stan declarations carry their dimensions either as a literal (scalars) or as
# the name of a variable in the data block
resolve_stan_dim <- function(dims, standata_list) {
  vapply(dims, function(d) {
    literal <- suppressWarnings(as.numeric(d))
    if (is.na(literal)) as.numeric(standata_list[[d]])[1] else literal
  }, numeric(1), USE.NAMES = FALSE)
}

############################################################################# !
# POSTPROCESSING AND EXTRACTORS                                          ####
############################################################################# !

# fit$bmm deliberately has no `model` or `user_formula` field: methods that
# assume a single model must provide an explicit mvbmmfit method instead of
# silently operating on the wrong structure
#' @export
postprocess_brm.bmm_components <- function(model, fit, ...) {
  dots <- list(...)
  class(fit) <- c("mvbmmfit", "bmmfit", "brmsfit")
  fit <- add_bmm_fit_info(
    fit,
    nlist(components = unclass(model), configure_opts = dots$configure_opts),
    dots$fit_args
  )
  reset_env(fit)
}

#' @rdname stancode.bmmformula
#' @export
stancode.mvbmmformula <- function(object, prior = NULL, ...) {
  withr::local_options(bmm.sort_data = FALSE)
  dots <- list(...)
  local_brms_threads(dots)
  cfg <- configure_fit(object, prior = prior, init = FALSE)
  add_bmm_version_to_stancode(call_brms_extractor(brms::stancode, cfg, dots))
}

#' @rdname standata.bmmformula
#' @export
standata.mvbmmformula <- function(object, ...) {
  withr::local_options(bmm.sort_data = FALSE)
  dots <- list(...)
  local_brms_threads(dots)
  cfg <- configure_fit(object, init = FALSE)
  call_brms_extractor(brms::standata, cfg, dots, prior = NULL)
}

#' @rdname default_prior.bmmformula
#' @export
default_prior.mvbmmformula <- function(object, ...) {
  withr::local_options(bmm.sort_data = FALSE)
  cfg <- configure_fit(object, init = FALSE)
  combine_prior(call_brms_extractor(brms::default_prior, cfg, list(...)), cfg$prior)
}
