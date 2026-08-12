############################################################################# !
# MULTIVARIATE COMPOSITION HELPERS                                       ####
############################################################################# !

# Restrict a component's likelihood to its own rows of the stacked data.
# subset() must be attached with `|` when the response has no addition terms
# yet, and with `+` after existing addition terms; the wrong separator makes
# brms evaluate base::subset() with a cryptic error
add_subset <- function(bformula, subset_var) {
  f <- bformula$formula
  lhs <- deparse1(f[[2]])
  rhs <- deparse1(f[[3]])
  sep <- if (grepl("|", lhs, fixed = TRUE)) " + " else " | "
  new <- stats::as.formula(
    paste0(lhs, sep, "subset(", subset_var, ") ~ ", rhs),
    env = environment(f)
  )
  # as.formula() drops brms's attributes on the formula object (e.g. the
  # non-linear flag on m3's main formula), so they must be restored
  new_attributes <- setdiff(names(attributes(f)), c("class", ".Environment"))
  for (a in new_attributes) {
    attr(new, a) <- attr(f, a)
  }
  bformula$formula <- new
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
# arguments of a single multivariate brms model. Shared engine behind bmm(),
# stancode(), standata() and default_prior() for mvbmmformula objects.
mvbmm_config <- function(spec, prior = NULL) {
  stopif(
    length(spec) < 2,
    "A multivariate bmm model requires at least two components. Combine \\
    components created by bmm_component() with `+`."
  )
  components <- lapply(spec, mv_configure_component)
  validate_component_resps(components)
  warn_missing_shared_re(spec)
  check_grouping_overlap(spec)
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
  init <- mv_create_initfun(components, formula, data, prior)

  nlist(
    config_args = nlist(formula, data, stanvars, init),
    prior,
    components = lapply(components, function(x) {
      x[c("model", "user_formula", "resp", "resp_name", "subset_var", "data_name")]
    })
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

# The point of a multivariate model is the joint correlation matrix; without
# a shared |ID| random-effects label the components are estimated jointly but
# independently, which is almost never what the user wants
warn_missing_shared_re <- function(spec) {
  ids <- lapply(spec, function(component) formula_re_ids(component$formula))
  id_counts <- table(unlist(lapply(ids, unique)))
  warnif(
    !any(id_counts >= 2),
    "No random-effects ID is shared between the components. To estimate the \\
    correlations between subject-level parameters across components, use the \\
    |ID| syntax with the same ID in at least two components, e.g. \\
    (1 | p | id). The model will run, but no cross-component correlations \\
    will be estimated."
  )
}

# the |ID| labels between two bars in random-effects terms, e.g. "p" in
# (1 | p | id)
formula_re_ids <- function(formula) {
  forms <- Filter(is_formula, unclass(formula))
  txt <- vapply(forms, deparse1, character(1))
  ids <- unlist(regmatches(txt, gregexpr("\\|\\s*[^|()~+ ]+\\s*\\|", txt)))
  unique(gsub("[| ]", "", ids))
}

# the (|ID| label, grouping variable) pairs of the random-effects terms in a
# bmmformula, e.g. ("p", "id") for (1 | p | id)
formula_re_pairs <- function(formula) {
  forms <- Filter(is_formula, unclass(formula))
  txt <- vapply(forms, deparse1, character(1))
  matches <- unlist(regmatches(
    txt,
    gregexpr("\\|\\s*[^|()~+ ]+\\s*\\|\\s*[^|()~+ ]+", txt)
  ))
  if (length(matches) == 0) {
    return(NULL)
  }
  parts <- strsplit(gsub(" ", "", matches), "|", fixed = TRUE)
  unique(do.call(rbind, lapply(parts, function(p) {
    data.frame(id = p[2], group = p[3])
  })))
}

# Cross-component correlations are identified only by grouping values (e.g.
# persons) observed in several components. Zero overlap leaves them
# unidentified — the posterior just reproduces the LKJ prior — and arises
# silently when IDs are coded differently across the component datasets
check_grouping_overlap <- function(spec) {
  pair_list <- lapply(spec, function(component) formula_re_pairs(component$formula))
  all_pairs <- unique(do.call(rbind, pair_list))
  for (i in seq_len(NROW(all_pairs))) {
    id_label <- all_pairs$id[i]
    group <- all_pairs$group[i]
    has_pair <- vapply(seq_along(spec), function(j) {
      !is.null(pair_list[[j]]) &&
        any(pair_list[[j]]$id == id_label & pair_list[[j]]$group == group) &&
        group %in% names(spec[[j]]$data)
    }, logical(1))
    if (sum(has_pair) < 2) {
      next
    }
    value_sets <- lapply(spec[has_pair], function(component) {
      unique(as.character(component$data[[group]]))
    })
    shared <- Reduce(intersect, value_sets)
    total <- unique(unlist(value_sets))
    warnif(
      length(shared) == 0,
      "The components tied by the random-effects ID '{id_label}' have no \\
      values of the grouping variable '{group}' in common. The correlations \\
      between their parameters cannot be estimated from these data, and \\
      their posterior will just reproduce the prior. Check that the values \\
      of '{group}' are coded identically across the component datasets."
    )
    if (length(shared) > 0 && length(shared) < length(total)) {
      message2(
        "{length(shared)} of {length(total)} unique values of '{group}' \\
        appear in all components tied by the random-effects ID '{id_label}'. \\
        The correlations between the components' parameters are informed by \\
        these {length(shared)} groups."
      )
    }
  }
}

# Factors used as fixed-effect predictors must have identical levels across
# components: the stacked data unions the levels, and levels absent from a
# component's rows produce all-zero design columns and rank deficiency
warn_factor_level_mismatch <- function(spec) {
  fe_vars <- lapply(spec, function(component) formula_fe_vars(component$formula))
  for (var in unique(unlist(fe_vars))) {
    levels_list <- lapply(seq_along(spec), function(i) {
      column <- spec[[i]]$data[[var]]
      if (var %in% fe_vars[[i]] && is.factor(column)) levels(column)
    })
    levels_list <- Filter(Negate(is.null), levels_list)
    mismatch <- length(levels_list) > 1 &&
      !all(vapply(levels_list, identical, logical(1), y = levels_list[[1]]))
    warnif(
      mismatch,
      "The factor '{var}' is used as a predictor in several components but \\
      has different levels across their datasets. Levels that do not occur \\
      in a component's data will produce inestimable coefficients for that \\
      component."
    )
  }
}

# fixed-effect variables of a bmmformula, excluding random-effects terms
formula_fe_vars <- function(formula) {
  forms <- Filter(is_formula, unclass(formula))
  unique(unlist(lapply(forms, function(f) {
    trms <- try(stats::terms(f), silent = TRUE)
    if (is_try_error(trms)) {
      return(character(0))
    }
    labels <- attr(trms, "term.labels")
    labels <- labels[!grepl("|", labels, fixed = TRUE)]
    unlist(lapply(labels, function(label) all.vars(str2lang(label))))
  })))
}

# Multivariate counterpart of create_initfun.bmmodel. Stan parameters
# belonging to components with init_ranges are initialized from those ranges,
# the shared random-effects and correlation parameters generically, and every
# remaining parameter as brm(init = 1) would initialize it — the value a
# component without init_ranges receives in a univariate bmm() fit. Leaving
# them out instead would make the backend fill the gaps silently and warn
mv_create_initfun <- function(components, formula, data, prior) {
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
  standata_list <- standata(formula, data, prior = prior)
  stan_code <- stancode(formula, data, prior = prior)
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
        ) %||% mv_default_init(spec, dim, standata_list)
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

# Stan draws initial values uniformly on the unconstrained scale and maps them
# through the parameter's declared bounds; with a radius of 1 this is exactly
# what brm(init = 1) does, the init a bmmodel without init_ranges gets in a
# univariate fit. Returns NULL for the types whose transform is not
# reimplemented here (simplex, ordered, ...) and for bounds that cannot be
# resolved to a number, leaving those parameters to the sampler
mv_default_init <- function(spec, dim, standata_list) {
  if (!spec$type %in% c("real", "vector", "row_vector", "matrix")) {
    return(NULL)
  }
  lower <- resolve_stan_bound(spec$bounds$lower, standata_list, -Inf)
  upper <- resolve_stan_bound(spec$bounds$upper, standata_list, Inf)
  if (is.na(lower) || is.na(upper) || anyNA(dim)) {
    return(NULL)
  }

  raw <- runif(prod(dim), min = -1, max = 1)
  value <- if (is.infinite(lower) && is.infinite(upper)) {
    raw
  } else if (is.infinite(upper)) {
    lower + exp(raw)
  } else if (is.infinite(lower)) {
    upper - exp(raw)
  } else {
    lower + (upper - lower) * plogis(raw)
  }
  if (length(dim) > 1) array(value, dim = dim) else value
}

# Stan declarations carry their dimensions either as a literal (scalars) or as
# the name of a variable in the data block
resolve_stan_dim <- function(dims, standata_list) {
  vapply(dims, function(d) {
    literal <- suppressWarnings(as.numeric(d))
    if (is.na(literal)) as.numeric(standata_list[[d]])[1] else literal
  }, numeric(1), USE.NAMES = FALSE)
}

# Bounds are parsed as strings and may be numbers or the names of data
# variables; NA marks a bound that is neither
resolve_stan_bound <- function(value, standata_list, default) {
  if (is.null(value)) {
    return(default)
  }
  number <- suppressWarnings(as.numeric(value))
  if (!is.na(number)) {
    return(number)
  }
  from_data <- standata_list[[value]]
  if (is.numeric(from_data) && length(from_data) == 1) from_data else NA_real_
}

############################################################################# !
# FITTING AND POSTPROCESSING                                             ####
############################################################################# !

# multivariate counterpart of the body of bmm(), mirroring its steps
bmm_mv <- function(spec, prior = NULL,
                   sort_data = getOption("bmm.sort_data", "check"),
                   silent = getOption("bmm.silent", 1),
                   backend = getOption("brms.backend", NULL),
                   file = NULL, file_compress = TRUE,
                   file_refit = getOption("bmm.file_refit", FALSE), ...) {
  deprecated_args(...)
  dots <- list(...)

  x <- try_read_bmmfit(file, file_refit)
  if (!is.null(x)) {
    return(x)
  }

  configure_opts <- nlist(
    sort_data, silent, backend,
    parallel = dots$parallel,
    cores = dots$cores
  )
  opts <- configure_options(configure_opts)
  dots$parallel <- NULL

  cfg <- mvbmm_config(spec, prior = prior)

  fit_args <- combine_args(nlist(
    config_args = cfg$config_args, opts, dots,
    prior = cfg$prior
  ))
  fit <- brms::do_call(brms::brm, fit_args)

  fit <- postprocess_brm_mv(cfg$components, fit,
    fit_args = fit_args,
    configure_opts = configure_opts
  )

  try_save_bmmfit(fit, file, compress = file_compress)
}

# multivariate counterpart of postprocess_brm.bmmodel. Not an S3 method,
# since there is no single model object to dispatch on. fit$bmm deliberately
# has no `model` or `user_formula` field: methods that assume a single model
# must provide an explicit mvbmmfit method instead of silently operating on
# the wrong structure
postprocess_brm_mv <- function(components, fit, fit_args, configure_opts) {
  class(fit) <- c("mvbmmfit", "bmmfit", "brmsfit")
  fit$version$bmm <- utils::packageVersion("bmm")
  fit$bmm <- nlist(components, configure_opts)
  attr(fit$data, "data_name") <- attr(fit_args$data, "data_name")
  fit$model <- add_bmm_version_to_stancode(fit$model)
  reset_env(fit)
}

############################################################################# !
# EXTRACTOR METHODS                                                      ####
############################################################################# !

#' @rdname stancode.bmmformula
#' @export
stancode.mvbmmformula <- function(object, prior = NULL, ...) {
  withr::local_options(bmm.sort_data = FALSE)
  cfg <- mvbmm_config(object, prior = prior)
  dots <- list(...)
  fit_args <- combine_args(nlist(config_args = cfg$config_args, dots, prior = cfg$prior))
  fit_args$object <- fit_args$formula
  fit_args$formula <- NULL
  code <- brms::do_call(brms::stancode, fit_args)
  add_bmm_version_to_stancode(code)
}

#' @rdname standata.bmmformula
#' @export
standata.mvbmmformula <- function(object, ...) {
  withr::local_options(bmm.sort_data = FALSE)
  cfg <- mvbmm_config(object)
  dots <- list(...)
  fit_args <- combine_args(nlist(config_args = cfg$config_args, dots, prior = cfg$prior))
  fit_args$object <- fit_args$formula
  fit_args$formula <- NULL
  brms::do_call(brms::standata, fit_args)
}

#' @rdname default_prior.bmmformula
#' @export
default_prior.mvbmmformula <- function(object, ...) {
  withr::local_options(bmm.sort_data = FALSE)
  cfg <- mvbmm_config(object)
  dots <- list(...)
  prior_args <- combine_args(nlist(config_args = cfg$config_args, dots, prior = cfg$prior))
  prior_args$object <- prior_args$formula
  prior_args$formula <- NULL
  brms_priors <- brms::do_call(brms::default_prior, prior_args)
  combine_prior(brms_priors, prior_args$prior)
}
