############################################################################# !
# LATENT VARIABLE STRUCTURE FOR MULTIVARIATE MODELS                      ####
############################################################################# !

#' @title Latent variable structure for the person-level parameters of a
#'   multivariate bmm model
#' @description `bmm_latent()` describes a structural equation model for the
#'   subject-level deviations (random effects) of the parameters of a
#'   multivariate bmm model, written in `lavaan` syntax. It is added with `+`
#'   to components created by [bmm_component()]. Instead of an unstructured
#'   correlation matrix, the subject-level parameters are then decomposed into
#'   common factors and unique parts, and the factors may be correlated,
#'   regressed on each other, or regressed on person-level covariates.
#'
#' @details The variables of the `lavaan` model are of three kinds:
#'
#'   * **Model parameters**: every random-effects coefficient of the grouping
#'   variable `group` in the components, named `<response>_<parameter>` for
#'   random intercepts and `<response>_<parameter>_<coefficient>` for random
#'   slopes (e.g. `error_kappa`, `rt_drift`, `error_kappa_set_size`). Use
#'   [latent_indicators()] to list the names available for a specification.
#'   Their scale is the link scale of the parameter (e.g. log kappa).
#'   * **Observed person-level variables**: columns of `data` (one row per
#'   level of `group`). A column on the right of `=~` is an observed indicator
#'   with a normal measurement model; a column on the right of `~` is a
#'   covariate measured without error. Both are standardized internally.
#'   * **Latent variables**: everything on the left of `=~`.
#'
#'   Supported operators: `=~` (loadings), `~` (regressions of latent variables
#'   on latent variables or covariates, and of model parameters on
#'   covariates), `~~` (covariances among exogenous latent variables, residual
#'   correlations between two model parameters, and latent variances), fixed
#'   values (`0.5*x`) and equality constraints through labels (`a*x1 + a*x2`).
#'
#'   With `std.lv = TRUE` (the default) the latent (residual) variances are
#'   fixed to 1 and all loadings are free; the sign of each factor is aligned
#'   across draws after sampling, so that the loading largest in absolute
#'   value is positive. With `std.lv = FALSE` the first loading of each factor
#'   is fixed to 1 and the latent variances are estimated.
#'
#'   Random-effects terms that contain at least one indicator are moved into a
#'   block of their own and estimated without a free correlation matrix: all
#'   covariance among them runs through the latent structure, and the
#'   standard deviations brms reports for them are the unique (residual)
#'   standard deviations. Random effects of the same grouping variable that are
#'   not named in the model keep their brms correlations among themselves.
#'
#'   This is an experimental development feature.
#'
#' @param model A character string with a model in `lavaan` syntax.
#' @param group The name of the grouping variable whose random effects the
#'   latent variables describe, e.g. `"id"`.
#' @param data An optional data.frame with one row per level of `group`,
#'   containing observed person-level indicators and covariates.
#' @param std.lv Logical. Identify latent variables by fixing their variances
#'   to 1 (`TRUE`, the default) or by fixing their first loading to 1
#'   (`FALSE`).
#'
#' @return An object of class `bmm_latent`, to be added to a multivariate
#'   specification with `+`.
#'
#' @seealso [bmm_component()], [latent_indicators()], [latent_summary()]
#'
#' @examples
#' \dontrun{
#' spec <- bmm_component(
#'   bmf(thetat ~ 1 + (1 | p | id), kappa ~ 1 + (1 | p | id)),
#'   model = mixture2p(resp_error = "error"), data = dat_vwm
#' ) +
#'   bmm_component(bmf(rt ~ 1 + (1 | p | id)), family = brms::lognormal(), data = dat_rt)
#'
#' latent_indicators(spec)
#'
#' spec <- spec + bmm_latent("
#'   WM =~ error_kappa + error_thetat + span
#'   WM ~ age
#'   rt_mu ~ WM
#' ", group = "id", data = persons)
#' }
#' @export
bmm_latent <- function(model, group, data = NULL, std.lv = TRUE) {
  stop_missing_args()
  stopif(
    !requireNamespace("lavaan", quietly = TRUE),
    "bmm_latent() needs the 'lavaan' package to parse the model syntax. \\
    Please install it with install.packages('lavaan')."
  )
  stopif(
    !is.character(model) || length(model) != 1,
    "The 'model' argument must be a single character string in lavaan syntax."
  )
  stopif(
    !is.character(group) || length(group) != 1,
    "The 'group' argument must be the name of the grouping variable, e.g. \"id\"."
  )
  data_name <- substitute_name(data)
  if (!is.null(data)) {
    stopif(!is.data.frame(data), "The 'data' argument must be a data.frame.")
    stopif(
      not_in(group, names(data)),
      "The grouping variable '{group}' is not a column of 'data'."
    )
    stopif(
      anyDuplicated(data[[group]]) > 0,
      "'data' must have one row per level of '{group}'."
    )
  }
  partable <- lavaan::lavaanify(
    model,
    std.lv = std.lv, auto.fix.first = !std.lv, auto.var = TRUE,
    auto.cov.lv.x = TRUE, fixed.x = TRUE, meanstructure = FALSE
  )
  structure(nlist(model, group, data, std.lv, partable, data_name), class = "bmm_latent")
}

is_bmm_latent <- function(x) {
  inherits(x, "bmm_latent")
}

#' @export
print.bmm_latent <- function(x, ...) {
  cat(style("green")(paste0("Latent structure for the random effects of '", x$group, "':\n")))
  cat(paste0("  ", trimws(strsplit(trimws(x$model), "\n")[[1]]), collapse = "\n"), "\n")
  if (!is.null(x$data)) {
    cat(style("purple1")("  Person-level data:"), x$data_name %||% "", "\n")
  }
  invisible(x)
}

#' @title Names of the model parameters available to a latent structure
#' @description Lists the subject-level parameters (random effects) of a
#'   multivariate bmm specification under the names by which they can be used
#'   as variables in the `lavaan` syntax of [bmm_latent()].
#' @param spec A multivariate specification created by combining
#'   [bmm_component()] objects with `+`, with or without a [bmm_latent()]
#'   structure.
#' @param group Optional name of a grouping variable. Defaults to the group of
#'   the attached [bmm_latent()] structure, or lists all grouping variables if
#'   there is none.
#' @return A data.frame with one row per variable: its `name` for the lavaan
#'   syntax, the `component` (response variable), the model `parameter`, the
#'   random-effects `coefficient`, the `group`, and the `link` that sets its
#'   scale. Observed columns of the `data` of an attached [bmm_latent()] are
#'   listed with `component = NA`.
#' @export
latent_indicators <- function(spec, group = NULL) {
  stopif(
    !is_mvbmmformula(spec),
    "'spec' must be a multivariate specification built from bmm_component() \\
    objects combined with `+`."
  )
  latent <- attr(spec, "latent")
  group <- group %||% latent$group
  components <- suppressMessages(lapply(spec, mv_configure_component))
  table <- latent_indicator_table(components, group)
  observed <- setdiff(names(latent$data), latent$group)
  if (length(observed) == 0) {
    return(table)
  }
  rbind(table, data.frame(
    name = observed, component = NA_character_, parameter = NA_character_,
    coefficient = NA_character_, group = latent$group, link = "identity (standardized)"
  ))
}

# One row per random-effects coefficient of the components, named as the
# variable the lavaan syntax refers to. Parameter names never contain an
# underscore (brms rejects them), so the first two parts of a name are always
# the sanitized response and the parameter
latent_indicator_table <- function(components, group = NULL) {
  rows <- lapply(components, function(x) {
    re <- brms:::frame_re(brms::brmsterms(x$cfg$formula), x$data)
    if (!is.null(group)) {
      re <- re[re$group == group, , drop = FALSE]
    }
    if (nrow(re) == 0) {
      return(NULL)
    }
    parameter <- ifelse(nzchar(re$nlpar), re$nlpar, ifelse(nzchar(re$dpar), re$dpar, "mu"))
    coef_suffix <- ifelse(re$coef == "Intercept", "", paste0("_", gsub("[^A-Za-z0-9_.]", ".", re$coef)))
    data.frame(
      name = paste0(x$resp_name, "_", parameter, coef_suffix),
      component = x$resp,
      parameter = parameter,
      coefficient = re$coef,
      group = re$group,
      link = vapply(parameter, function(p) x$model$links[[p]] %||% "identity", character(1), USE.NAMES = FALSE)
    )
  })
  table <- do.call(rbind, rows)
  if (is.null(table)) {
    return(data.frame(
      name = character(), component = character(), parameter = character(),
      coefficient = character(), group = character(), link = character()
    ))
  }
  rownames(table) <- NULL
  table
}

############################################################################# !
# STRUCTURE: FROM THE LAVAAN PARAMETER TABLE TO STAN EXPRESSIONS         ####
############################################################################# !

LATENT_RE_ID <- "bmmlv"

# Classifies every variable and row of the parameter table and assigns each
# free parameter its Stan expression. Rows that share a label share a
# parameter; this is how lavaan expresses equality constraints
latent_structure <- function(latent, components) {
  group <- latent$group
  indicators <- latent_indicator_table(components, group)
  pt <- latent$partable
  user_eq <- pt$op == "==" & pt$user == 1
  stopif(
    any(user_eq) || any(!pt$op %in% c("=~", "~", "~~", "==")),
    "bmm_latent() supports the operators =~, ~ and ~~, and equality \\
    constraints through shared labels. Other operators and explicit \\
    constraints are not supported yet."
  )
  pt <- pt[pt$op != "==", , drop = FALSE]

  latents <- unique(pt$lhs[pt$op == "=~"])
  observed <- setdiff(names(latent$data), group)
  clash <- intersect(observed, indicators$name)
  stopif(
    length(clash) > 0,
    "The person-level data column(s) {collapse_comma(clash)} have the same \\
    name as model parameters. Please rename them."
  )
  unknown <- setdiff(unique(c(pt$lhs, pt$rhs)), c(latents, indicators$name, observed))
  stopif(
    length(unknown) > 0,
    "The latent structure refers to {collapse_comma(unknown)}, which are \\
    neither latent variables, nor random effects of '{group}' in the \\
    components, nor columns of the person-level data. Use \\
    latent_indicators() to list the available names."
  )
  role <- function(v) {
    ifelse(v %in% latents, "latent", ifelse(v %in% indicators$name, "re", "observed"))
  }
  pt$lhs_role <- role(pt$lhs)
  pt$rhs_role <- role(pt$rhs)
  pt$kind <- latent_row_kind(pt)

  obs_ind <- unique(pt$rhs[pt$kind == "lambda" & pt$rhs_role == "observed"])
  covariates <- unique(pt$rhs[pt$kind == "beta" & pt$rhs_role == "observed"])
  both <- intersect(obs_ind, covariates)
  stopif(
    length(both) > 0,
    "The person-level variable(s) {collapse_comma(both)} are used both as \\
    an indicator (=~) and as a covariate (~)."
  )

  regress <- pt[pt$kind == "beta" & pt$lhs_role == "latent", , drop = FALSE]
  endo <- unique(regress$lhs)
  exo <- setdiff(latents, endo)
  latents <- c(exo, latent_topological_order(endo, regress[regress$rhs_role == "latent", ]))

  cov_rows <- pt[pt$kind == "cov", , drop = FALSE]
  stopif(
    any(cov_rows$lhs %in% endo | cov_rows$rhs %in% endo),
    "Covariances involving endogenous latent variables (those on the left of \\
    ~) are not supported yet."
  )
  exo_cor <- length(exo) > 1 && any(cov_rows$free > 0)
  stopif(
    length(exo) > 1 && !all(cov_rows$free > 0) && !all(cov_rows$free == 0 & cov_rows$ustart == 0),
    "The covariances among exogenous latent variables must be either all \\
    free or all fixed to 0."
  )
  re_names <- unique(c(pt$lhs[pt$lhs_role == "re"], pt$rhs[pt$rhs_role == "re"]))
  rescor <- pt[pt$kind == "rescor", , drop = FALSE]
  stopif(
    anyDuplicated(c(rescor$lhs, rescor$rhs)) > 0,
    "Each model parameter may have a residual correlation (~~) with at most \\
    one other parameter."
  )

  pt <- latent_assign_parameters(pt, latents)
  nlist(
    group, std.lv = latent$std.lv, model = latent$model, partable = pt,
    latents, exo, endo, exo_cor, obs_ind, covariates,
    indicators = indicators[indicators$name %in% re_names, , drop = FALSE],
    person_data = latent$data
  )
}

latent_row_kind <- function(pt) {
  kind <- rep(NA_character_, nrow(pt))
  self <- pt$lhs == pt$rhs
  kind[pt$op == "=~" & pt$rhs_role %in% c("re", "observed")] <- "lambda"
  kind[pt$op == "~" & pt$lhs_role == "latent" & pt$rhs_role %in% c("latent", "observed")] <- "beta"
  kind[pt$op == "~" & pt$lhs_role == "re" & pt$rhs_role == "observed"] <- "beta"
  kind[pt$op == "~~" & self & pt$lhs_role == "latent"] <- "sd"
  kind[pt$op == "~~" & self & pt$lhs_role != "latent" & pt$user == 0] <- "auto"
  kind[pt$op == "~~" & !self & pt$lhs_role == "latent" & pt$rhs_role == "latent"] <- "cov"
  kind[pt$op == "~~" & !self & pt$lhs_role == "re" & pt$rhs_role == "re"] <- "rescor"
  bad <- is.na(kind)
  stopif(
    any(bad),
    "The following parts of the latent structure are not supported yet: \\
    {collapse_comma(paste(pt$lhs[bad], pt$op[bad], pt$rhs[bad]))}. Supported \\
    are loadings of model parameters and person-level variables on latent \\
    variables (=~), regressions of latent variables or model parameters on \\
    covariates and of latent variables on latent variables (~), and \\
    covariances among exogenous latent variables and between two model \\
    parameters (~~)."
  )
  stopif(
    any(kind == "rescor" & (pt$free == 0 | nzchar(pt$label))),
    "Residual correlations between model parameters cannot be fixed or \\
    labelled yet."
  )
  kind
}

latent_topological_order <- function(endo, edges) {
  ordered <- character()
  remaining <- endo
  while (length(remaining) > 0) {
    ready <- remaining[vapply(remaining, function(v) {
      all(setdiff(edges$rhs[edges$lhs == v], ordered) %in% setdiff(edges$rhs, endo))
    }, logical(1))]
    stopif(
      length(ready) == 0,
      "The regressions among latent variables contain a cycle, which is not \\
      supported."
    )
    ordered <- c(ordered, ready)
    remaining <- setdiff(remaining, ready)
  }
  ordered
}

# Stan expression of every row: a literal for fixed values, an element of a
# typed parameter vector for free ones. Loadings are left unconstrained: a
# positivity constraint on one loading traps chains at zero when that indicator
# is weak, so the sign of a factor is aligned per draw after sampling instead
# (see latent_signs())
latent_assign_parameters <- function(pt, latents) {
  pt$key <- ifelse(nzchar(pt$label), pt$label, paste0(".row", pt$id))
  kinds_per_key <- tapply(pt$kind, pt$key, function(k) length(unique(k)))
  stopif(
    any(kinds_per_key > 1),
    "A label is shared by parameters of different kinds (e.g. a loading and \\
    a regression), which is not supported."
  )
  pt$stan <- NA_character_
  fixed <- pt$free == 0
  pt$stan[fixed & pt$kind %in% c("lambda", "beta")] <- as.character(pt$ustart[fixed & pt$kind %in% c("lambda", "beta")])
  pt$stan[fixed & pt$kind == "sd"] <- as.character(sqrt(pt$ustart[fixed & pt$kind == "sd"]))
  vec <- function(select, name) {
    keys <- unique(pt$key[select])
    pt$stan[select] <<- paste0(name, "[", match(pt$key[select], keys), "]")
  }
  vec(!fixed & pt$kind == "lambda", "lv_lambda")
  vec(!fixed & pt$kind == "beta", "lv_beta")
  vec(!fixed & pt$kind == "sd", "lv_sd")
  vec(pt$kind == "rescor", "lv_rescor")
  pt
}

############################################################################# !
# FORMULA REWRITING                                                     ####
############################################################################# !

# Moves every random-effects term of `group` that contains an indicator into
# one uncorrelated block shared by all components, gr(group, cor = FALSE,
# id = "bmmlv"); the latent structure then supplies their covariance
rewrite_latent_terms <- function(components, structure) {
  ind <- structure$indicators
  lapply(components, function(x) {
    pars <- unique(ind$parameter[ind$component == x$resp])
    if (length(pars) == 0) {
      return(x)
    }
    bf <- x$cfg$formula
    for (par in pars) {
      if (par == "mu" && !par %in% names(bf$pforms)) {
        bf$formula[[3]] <- rewrite_re_call(bf$formula[[3]], structure$group)
      } else {
        bf$pforms[[par]][[3]] <- rewrite_re_call(bf$pforms[[par]][[3]], structure$group)
      }
    }
    x$cfg$formula <- bf
    # a correlation matrix the rewrite dissolved has no parameter left for its
    # prior; mv_cor_prior() rebuilds the ones that remain from the joint formula
    x$prior <- x$prior[!x$prior$class %in% c("cor", "L"), ]
    x
  })
}

rewrite_re_call <- function(expr, group) {
  if (!is.call(expr)) {
    return(expr)
  }
  if (identical(expr[[1]], as.name("(")) && is.call(expr[[2]]) &&
    as.character(expr[[2]][[1]]) %in% c("|", "||")) {
    bar <- expr[[2]]
    term_group <- bar[[3]]
    group_name <- if (is.call(term_group) && identical(term_group[[1]], as.name("gr"))) {
      stopif(
        length(term_group) > 2,
        "Random-effects terms written with gr() and further arguments cannot \\
        be part of a latent structure yet."
      )
      as.character(term_group[[2]])
    } else if (is.name(term_group)) {
      as.character(term_group)
    }
    if (identical(group_name, group)) {
      lhs <- bar[[2]]
      if (is.call(lhs) && identical(lhs[[1]], as.name("|"))) {
        lhs <- lhs[[2]]
      }
      new_group <- call("gr", as.name(group), cor = FALSE, id = LATENT_RE_ID)
      return(call("(", call("|", lhs, new_group)))
    }
    return(expr)
  }
  expr[] <- lapply(as.list(expr), rewrite_re_call, group = group)
  expr
}

############################################################################# !
# STAN CODE                                                              ####
############################################################################# !

# Locates the indicators in the brms random-effects block of the composed
# model (its number, the Stan names of their vectors and the order of the
# group levels), and aligns the person-level data with these levels
latent_locate_block <- function(structure, formula, data) {
  re <- brms:::frame_re(brms::brmsterms(formula), data)
  ind <- structure$indicators
  parameter <- ifelse(nzchar(re$nlpar), re$nlpar, ifelse(nzchar(re$dpar), re$dpar, "mu"))
  rows <- match(
    paste(ind$component, ind$parameter, ind$coefficient),
    paste(re$resp, parameter, re$coef)
  )
  block <- unique(re$id[rows])
  stopif(
    anyNA(rows) || length(block) != 1 || any(re$cor[rows]),
    "Internal error: the random effects named in the latent structure were \\
    not found in one uncorrelated block. Please report this."
  )
  re <- re[rows, , drop = FALSE]
  prefix <- ifelse(
    nzchar(re$nlpar), paste0(re$resp, "_", re$nlpar),
    ifelse(nzchar(re$dpar) & re$dpar != "mu", paste0(re$dpar, "_", re$resp), re$resp)
  )
  ind$stan_r <- paste0("r_", block, "_", prefix, "_", re$cn)
  ind$cn <- re$cn
  ind$brms_sd <- paste0("sd_", structure$group, "__", prefix, "_", re$coef)
  structure$indicators <- ind
  structure$block <- block
  structure$levels <- attr(re, "levels")[[structure$group]]

  person_vars <- c(structure$covariates, structure$obs_ind)
  if (length(person_vars) > 0) {
    person_data <- structure$person_data
    idx <- match(structure$levels, as.character(person_data[[structure$group]]))
    stopif(
      anyNA(idx),
      "{sum(is.na(idx))} levels of '{structure$group}' in the components have \\
      no row in the person-level data of bmm_latent()."
    )
    values <- as.matrix(person_data[idx, person_vars, drop = FALSE])
    stopif(
      anyNA(values),
      "The person-level data of bmm_latent() contain missing values, which \\
      are not supported yet."
    )
    structure$scaling <- list(center = colMeans(values), scale = apply(values, 2, stats::sd))
    structure$person_values <- scale(values)
  }
  structure
}

latent_stanvars <- function(structure) {
  pt <- structure$partable
  N <- paste0("N_", structure$block)
  latents <- structure$latents
  lv_idx <- function(v) match(v, latents)
  x_idx <- function(v) match(v, structure$covariates)
  y_idx <- function(v) match(v, structure$obs_ind)
  term <- function(rows) {
    vapply(rows, function(i) {
      if (pt$op[i] == "=~") {
        glue("{pt$stan[i]} * lv_eta[, {lv_idx(pt$lhs[i])}]")
      } else if (pt$rhs_role[i] == "latent") {
        glue("{pt$stan[i]} * lv_eta[, {lv_idx(pt$rhs[i])}]")
      } else {
        glue("{pt$stan[i]} * lv_X[, {x_idx(pt$rhs[i])}]")
      }
    }, character(1))
  }
  n_of <- function(name) sum(grepl(paste0("^", name, "\\["), unique(pt$stan)))
  sd_of <- function(v) {
    s <- pt$stan[pt$kind == "sd" & pt$lhs == v]
    if (length(s) == 0) "1" else s
  }

  sv <- list()
  add <- function(scode, block, position = "start", ...) {
    sv[[length(sv) + 1]] <<- brms::stanvar(scode = scode, block = block, position = position, ...)
  }

  Q <- length(structure$covariates)
  P <- length(structure$obs_ind)
  if (Q > 0) {
    sv[[length(sv) + 1]] <- brms::stanvar(
      x = structure$person_values[, structure$covariates, drop = FALSE], name = "lv_X",
      scode = glue("  matrix[{N}, {Q}] lv_X;  // standardized person-level covariates"),
      block = "data"
    )
  }
  if (P > 0) {
    sv[[length(sv) + 1]] <- brms::stanvar(
      x = structure$person_values[, structure$obs_ind, drop = FALSE], name = "lv_Y",
      scode = glue("  matrix[{N}, {P}] lv_Y;  // standardized person-level indicators"),
      block = "data"
    )
  }

  F_all <- length(latents)
  F_exo <- length(structure$exo)
  par_lines <- c(
    glue("  matrix[{N}, {F_all}] z_lv;  // standardized latent variables"),
    if (n_of("lv_lambda") > 0) glue("  vector[{n_of('lv_lambda')}] lv_lambda;  // loadings"),
    if (n_of("lv_beta") > 0) glue("  vector[{n_of('lv_beta')}] lv_beta;  // latent regressions"),
    if (n_of("lv_sd") > 0) glue("  vector<lower=0>[{n_of('lv_sd')}] lv_sd;  // latent (residual) SDs"),
    if (structure$exo_cor) glue("  cholesky_factor_corr[{F_exo}] L_lv;  // exogenous latent correlations"),
    if (n_of("lv_rescor") > 0) glue("  vector<lower=-1, upper=1>[{n_of('lv_rescor')}] lv_rescor;  // residual correlations"),
    if (P > 0) glue("  vector[{P}] lv_nu;  // intercepts of observed indicators"),
    if (P > 0) glue("  vector<lower=0>[{P}] lv_sigma;  // residual SDs of observed indicators")
  )
  add(paste(par_lines, collapse = "\n"), "parameters")

  exo_sd <- vapply(structure$exo, sd_of, character(1))
  eta_lines <- glue("  matrix[{N}, {F_all}] lv_eta;  // latent variables")
  if (structure$exo_cor) {
    eta_lines <- c(eta_lines, glue(
      "  lv_eta[, 1:{F_exo}] = z_lv[, 1:{F_exo}] * diag_pre_multiply([{paste(exo_sd, collapse = ', ')}]', L_lv)';"
    ))
  } else {
    for (f in structure$exo) {
      eta_lines <- c(eta_lines, glue("  lv_eta[, {lv_idx(f)}] = {sd_of(f)} * z_lv[, {lv_idx(f)}];"))
    }
  }
  for (f in structure$endo) {
    rows <- which(pt$kind == "beta" & pt$lhs == f)
    eta_lines <- c(eta_lines, glue(
      "  lv_eta[, {lv_idx(f)}] = {paste(c(term(rows), paste0(sd_of(f), ' * z_lv[, ', lv_idx(f), ']')), collapse = ' + ')};"
    ))
  }
  add(paste(eta_lines, collapse = "\n"), "tparameters", "start")

  ind <- structure$indicators
  r_lines <- character()
  for (k in seq_len(nrow(ind))) {
    rows <- which((pt$kind == "lambda" & pt$rhs == ind$name[k]) |
      (pt$kind == "beta" & pt$lhs == ind$name[k]))
    if (length(rows) > 0) {
      r_lines <- c(r_lines, glue("  {ind$stan_r[k]} += {paste(term(rows), collapse = ' + ')};"))
    }
  }
  sd_b <- paste0("sd_", structure$block)
  z_b <- paste0("z_", structure$block)
  for (i in which(pt$kind == "rescor")) {
    a <- ind$cn[ind$name == pt$lhs[i]]
    b <- match(pt$rhs[i], ind$name)
    r_lines <- c(r_lines, glue(
      "  {ind$stan_r[b]} += {sd_b}[{ind$cn[b]}] * ({pt$stan[i]} * {z_b}[{a}] + (sqrt(1 - square({pt$stan[i]})) - 1) * {z_b}[{ind$cn[b]}]);"
    ))
  }
  if (length(r_lines) > 0) {
    add(paste(r_lines, collapse = "\n"), "tparameters", "end")
  }

  model_lines <- c(
    "  target += std_normal_lpdf(to_vector(z_lv));",
    if (n_of("lv_lambda") > 0) "  target += normal_lpdf(lv_lambda | 0, 1);",
    if (n_of("lv_beta") > 0) "  target += normal_lpdf(lv_beta | 0, 1);",
    if (n_of("lv_sd") > 0) "  target += exponential_lpdf(lv_sd | 1);",
    if (structure$exo_cor) "  target += lkj_corr_cholesky_lpdf(L_lv | 2);",
    if (P > 0) "  target += normal_lpdf(lv_nu | 0, 1);",
    if (P > 0) "  target += exponential_lpdf(lv_sigma | 1);"
  )
  for (p in structure$obs_ind) {
    rows <- which(pt$kind == "lambda" & pt$rhs == p)
    model_lines <- c(model_lines, glue(
      "  target += normal_lpdf(lv_Y[, {y_idx(p)}] | lv_nu[{y_idx(p)}] + {paste(term(rows), collapse = ' + ')}, lv_sigma[{y_idx(p)}]);"
    ))
  }
  add(paste(model_lines, collapse = "\n"), "model")

  if (structure$exo_cor) {
    add(
      glue("  corr_matrix[{F_exo}] Cor_lv = multiply_lower_tri_self_transpose(L_lv);"),
      "genquant"
    )
  }
  Reduce(`+`, sv)
}

############################################################################# !
# POSTPROCESSING                                                         ####
############################################################################# !

# The likelihood is invariant to flipping the sign of a factor whose loadings
# are all free, so each draw is reflected such that the factor's reference
# loading -- the one largest in absolute value on average -- is positive.
# Returns a draws x factors matrix of +1/-1 (+1 for factors with a fixed loading)
latent_signs <- function(structure, draws) {
  pt <- structure$partable
  vapply(structure$latents, function(f) {
    rows <- pt$kind == "lambda" & pt$lhs == f
    if (!any(rows) || any(pt$free[rows] == 0)) {
      return(rep(1, nrow(draws)))
    }
    lambda <- draws[, unique(pt$stan[rows]), drop = FALSE]
    ref <- which.max(colMeans(abs(lambda)))
    ifelse(lambda[, ref] < 0, -1, 1)
  }, numeric(nrow(draws)))
}

# per-row sign multiplier of a draw: loadings and latent regressions change
# sign with every latent variable they involve
latent_row_signs <- function(rows, signs) {
  vapply(seq_len(nrow(rows)), function(i) {
    s <- rep(1, nrow(signs))
    for (v in c(rows$lhs[i], rows$rhs[i])) {
      if (v %in% colnames(signs) && !(rows$op[i] == "~~" && rows$lhs[i] == rows$rhs[i])) {
        s <- s * signs[, v]
      }
    }
    s
  }, numeric(nrow(signs)))
}

#' @title Estimates of the latent structure of a multivariate bmm fit
#' @description Summarises the posterior of the parameters of a
#'   [bmm_latent()] structure in the layout of a `lavaan` parameter table.
#' @param fit A multivariate `bmmfit` with a latent structure.
#' @param probs The lower and upper quantiles of the credible intervals.
#' @return A data.frame with one row per parameter: `lhs`, `op`, `rhs`, the
#'   `label`, whether it was `free`, and the posterior mean, SD, quantiles,
#'   Rhat and bulk ESS. Loadings and regressions are on the link scale of the
#'   model parameters and on the standardized scale of person-level variables.
#'   Rows with `op = "~~"` and the same `lhs` and `rhs` are standard
#'   deviations (latent, unique or residual), not variances. The sign of each
#'   factor is aligned across draws so that its largest loading is positive.
#' @export
latent_summary <- function(fit, probs = c(0.025, 0.975)) {
  structure <- fit$bmm$latent
  stopif(is.null(structure), "This fit has no latent structure (see bmm_latent()).")
  pt <- structure$partable
  rows <- pt[pt$kind %in% c("lambda", "beta", "sd", "rescor"), c("lhs", "op", "rhs", "label", "free", "stan")]
  ind <- structure$indicators
  rows <- rbind(
    rows,
    data.frame(lhs = ind$name, op = "~~", rhs = ind$name, label = "", free = 1, stan = ind$brms_sd)
  )
  if (structure$exo_cor) {
    pairs <- utils::combn(seq_along(structure$exo), 2)
    rows <- rbind(rows, data.frame(
      lhs = structure$exo[pairs[1, ]], op = "~~", rhs = structure$exo[pairs[2, ]],
      label = "", free = 1, stan = paste0("Cor_lv[", pairs[1, ], ",", pairs[2, ], "]")
    ))
  }
  if (length(structure$obs_ind) > 0) {
    rows <- rbind(rows, data.frame(
      lhs = structure$obs_ind, op = "~~", rhs = structure$obs_ind, label = "", free = 1,
      stan = paste0("lv_sigma[", seq_along(structure$obs_ind), "]")
    ))
  }
  rows <- rows[grepl("[A-Za-z]", rows$stan), , drop = FALSE]
  variables <- unique(c(rows$stan, pt$stan[pt$kind == "lambda" & grepl("[A-Za-z]", pt$stan)]))
  draws <- posterior::as_draws_df(posterior::subset_draws(brms::as_draws_array(fit), variable = variables))
  values <- as.matrix(as.data.frame(draws)[, variables, drop = FALSE])
  aligned <- values[, rows$stan, drop = FALSE] *
    latent_row_signs(rows, latent_signs(structure, values))
  colnames(aligned) <- paste0("row", seq_len(nrow(rows)))
  aligned <- posterior::as_draws_df(cbind(
    as.data.frame(aligned),
    .chain = draws$.chain, .iteration = draws$.iteration
  ))
  est <- posterior::summarise_draws(
    aligned, mean, sd, ~ stats::quantile(.x, probs = probs),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
  )
  out <- cbind(rows[, c("lhs", "op", "rhs", "label")], free = TRUE, as.data.frame(est[, -1]))
  rownames(out) <- NULL
  out
}

#' @title Person-level scores on the latent variables
#' @description Posterior draws of each person's value on each latent variable
#'   of a [bmm_latent()] structure, with factor signs aligned as in
#'   [latent_summary()].
#' @param fit A multivariate `bmmfit` with a latent structure.
#' @return An array of draws x persons x latent variables, with the levels of
#'   the grouping variable and the latent variable names as dimnames.
#' @export
latent_scores <- function(fit) {
  structure <- fit$bmm$latent
  stopif(is.null(structure), "This fit has no latent structure (see bmm_latent()).")
  pt <- structure$partable
  lambda <- unique(pt$stan[pt$kind == "lambda" & grepl("[A-Za-z]", pt$stan)])
  draws <- unclass(posterior::as_draws_matrix(posterior::subset_draws(
    brms::as_draws_array(fit), variable = c("lv_eta", lambda)
  )))
  signs <- latent_signs(structure, draws)
  N <- length(structure$levels)
  F_all <- length(structure$latents)
  eta <- array(
    draws[, grepl("^lv_eta", colnames(draws))], dim = c(nrow(draws), N, F_all),
    dimnames = list(NULL, structure$levels, structure$latents)
  )
  for (f in seq_len(F_all)) {
    eta[, , f] <- eta[, , f] * signs[, f]
  }
  eta
}
