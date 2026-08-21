#' Transform kappa of the von Mises distribution to the circular standard
#' deviation
#' @description This function transforms the precision parameter kappa of the
#'   von Mises distribution to the circular standard deviation. Adapted from
#'   Matlab code by Paul Bays (https://www.paulbays.com/code.php)
#'
#' @param K numeric. A vector of kappa values.
#' @return A vector of sd values.
#' @keywords transform
#' @export
#' @examples
#' kappas <- runif(1000, 0.01, 100)
#'
#' # calcualte SD (in radians)
#' SDs <- k2sd(kappas)
#'
#' # transform SDs from radians to degrees
#' SDs_degress <- SDs * 180 / pi
#'
#' # plot the relationship between kappa and circular SD
#' plot(kappas, SDs)
#' plot(kappas, SDs_degress)
k2sd <- function(K) {
  log_bessel_ratio <- log(besselI(K, 1, expon.scaled = T)) - log(besselI(K, 0, expon.scaled = T))
  S <- sqrt(-2 * log_bessel_ratio)
  S[K == 0] <- Inf
  S[is.infinite(K)] <- 0
  S
}


#' Convert between parametrizations of the c parameter of the SDM distribution
#'
#' @name c_parametrizations
#' @inheritParams SDMdist
#' @return A numeric vector of the same length as `c` and `kappa`.
#' @details
#' `c_bessel2sqrtexp` converts the memory strength parameter (c)
#'   from the bessel parametrization to the sqrtexp parametrization,
#'   `c_sqrtexp2bessel` converts from the sqrtexp parametrization to the
#'   bessel parametrization.
#' @keywords transform
#' @details See [the online article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for details on the
#'   parameterization. The sqrtexp parametrization is the default in the
#'   `bmm` package.
#' @export
#'
#' @examples
#' c_bessel <- c_sqrtexp2bessel(c = 4, kappa = 3)
#' c_sqrtexp <- c_bessel2sqrtexp(c = c_bessel, kappa = 3)
#'
c_sqrtexp2bessel <- function(c, kappa) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(c < 0)), "c must be non-negative")
  c * besselI(kappa, 0, expon.scaled = TRUE) * sqrt(2 * pi * kappa)
}

#' @rdname c_parametrizations
#' @keywords transform
#' @export
c_bessel2sqrtexp <- function(c, kappa) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(c < 0)), "c must be non-negative")
  c / (besselI(kappa, 0, expon.scaled = TRUE) * sqrt(2 * pi * kappa))
}

#' @title Transform values from the native to parameter scale or vice versa according to a link function
#'
#' @description
#' This function transforms a vector of values from the native scale to the parameter scale,
#' according to the specified link function. The function is mainly used internally, to ensure proper
#' initial values
#'
#' @param values A numerical vector of values the transformation should be applied to
#' @param link A character specifying the link to be applied.
#' Available options are: "identity", "log", "softplus", "log1p", "logm1", "inverse", "sqrt", "logit", "probit", "tan_half", and "cloglog".
#' @param inverse A Boolean value indicating if values should be transformed from the native to
#' the parameter scale (FALSE), or from the parameter scale to the native scale (TRUE)
#'
#' @noRd
link_transform <- function(values, link, inverse = FALSE) {
  # Handle NULL or missing link as identity (no transformation)
  if (is.null(link)) link <- "identity"

  stopifnot(is.numeric(values), is.character(link), length(link) == 1L, is.logical(inverse), length(inverse) == 1L)
  if(inverse) {
    switch(
      link,
      identity = values,
      log = exp(values),
      softplus = log1p(exp(values)),
      log1p = expm1(values),
      logm1 = brms::expp1(values),
      inverse = 1 / values,
      sqrt = values^2,
      logit = plogis(values),
      probit = pnorm(values),
      tan_half = 2 * atan(values),
      loglog = exp(-exp(values)),
      cloglog = 1 - exp(-exp(values)),
      stop2("Link '{link}' not recognized.")
    )
  } else {
    switch(
      link,
      identity = values,
      log = log(values),
      softplus = log(expm1(values)),
      log1p = log1p(values),
      logm1 = brms::logm1(values),
      inverse = 1 / values,
      sqrt = sqrt(values),
      logit = qlogis(values),
      probit = qnorm(values),
      tan_half = tan(values / 2),
      loglog = log(-log(values)),
      cloglog = log(-log1p(-values)),
      stop2("Link '{link}' not recognized.")
    )
  }
}

#' Get parameter information for a bmm model
#'
#' @description Returns a data frame with information about the model
#'   parameters, including their descriptions, whether they are fixed,
#'   their link functions, and optionally their default priors.
#'
#' @param x A \code{bmmodel} object (e.g., \code{sdm(resp_error = "y")}) or a
#'   \code{bmmfit} object (a fitted model returned by \code{\link{bmm}}).
#' @param formula An optional \code{bmmformula} object. Only relevant for
#'   M3 custom models, where additional parameters are discovered from
#'   the formula. Ignored for all other models.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{bmm_parameters} with one row per
#'   parameter and columns: \code{parameter}, \code{description},
#'   \code{fixed}, \code{value}, and \code{link}.
#'
#' @export
#' @examples
#' # For an unfitted model
#' parameters(sdm(resp_error = "y"))
#'
#' # For an M3 model
#' parameters(m3(
#'   resp_cats = c("corr", "other", "npl"),
#'   num_options = c(1, 4, 5),
#'   version = "ss"
#' ))
parameters <- function(x, ...) {
  UseMethod("parameters")
}


#' @rdname parameters
#' @export
parameters.bmmodel <- function(x, formula = NULL, ...) {
  model <- x

  if (inherits(model, "m3_custom") && !is.null(formula)) {
    user_pars <- rhs_vars(formula[is_nl(formula)])
    user_pars <- setdiff(user_pars, names(formula[is_nl(formula)]))
    user_pars <- setdiff(user_pars, names(model$parameters))
    model$parameters <- c(model$parameters, setNames(
      as.list(user_pars), user_pars
    ))
  }

  pars <- names(model$parameters)
  if (length(pars) == 0) {
    message2("This model has no parameters defined.")
    return(invisible(data.frame()))
  }

  fixed <- pars %in% names(model$fixed_parameters)
  values <- rep(NA_character_, length(pars))
  values[fixed] <- as.character(model$fixed_parameters[pars[fixed]])

  links <- vapply(pars, function(p) {
    model$links[[p]] %||% "identity"
  }, character(1))

  out <- data.frame(
    parameter = pars,
    description = as.character(model$parameters),
    fixed = fixed,
    value = values,
    link = links,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  m3_note <- NULL
  if (inherits(model, "m3_custom") && is.null(formula)) {
    m3_note <- paste(
      "Note: This is a custom M3 model. Pass a formula to",
      "discover additional parameters."
    )
  }

  class(out) <- c("bmm_parameters", "data.frame")
  attr(out, "model_name") <- model$name
  attr(out, "m3_note") <- m3_note
  out
}


#' @rdname parameters
#' @export
parameters.bmmfit <- function(x, ...) {
  x <- restructure(x)
  parameters(x$bmm$model, formula = x$bmm$user_formula, ...)
}


#' @export
print.bmm_parameters <- function(x, max_desc_width = 50, ...) {
  model_name <- attr(x, "model_name")
  m3_note <- attr(x, "m3_note")

  if (!is.null(model_name) && nzchar(model_name)) {
    cat(style("purple1")("Model: "), model_name, "\n\n")
  }

  print_df <- x
  if ("description" %in% names(print_df)) {
    print_df$description <- vapply(print_df$description, function(d) {
      d <- gsub("\\s+", " ", d)
      if (nchar(d) > max_desc_width) {
        paste0(substr(d, 1, max_desc_width - 3), "...")
      } else {
        d
      }
    }, character(1))
  }

  for (col in names(print_df)) {
    if (is.character(print_df[[col]])) {
      print_df[[col]][is.na(print_df[[col]])] <- "--"
    }
  }

  if ("fixed" %in% names(print_df)) {
    print_df$fixed <- ifelse(print_df$fixed, "yes", "no")
  }

  print.data.frame(print_df, right = FALSE, row.names = FALSE)

  if (!is.null(m3_note)) {
    cat("\n", m3_note, "\n")
  }

  invisible(x)
}

#' Get parameter classification info for a bmmfit parameter
#'
#' @description
#' Determines whether a named model parameter is a distributional (`dpar`) or
#' non-linear (`nlpar`) parameter in the underlying brms model, retrieves its
#' link function, and checks for softmax transformation.
#'
#' @param bmmfit A bmmfit object
#' @param par Character string. Parameter name (e.g., `"kappa"`, `"c"`)
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`type`}{Character: `"dpar"` or `"nlpar"`}
#'     \item{`model_name`}{Character: the parameter name as specified in the bmmodel}
#'     \item{`brms_name`}{Character: the parameter name as used in brms}
#'     \item{`link`}{Character: the link function (e.g., `"log"`, `"identity"`)}
#'     \item{`softmax`}{Logical: whether the parameter uses softmax transformation}
#'   }
#'
#' @keywords internal
#' @noRd
.get_parameter_info <- function(bmmfit, par) {
  model <- bmmfit$bmm$model
  model_pars <- names(model$parameters)
  bterms <- brms::brmsterms(bmmfit$formula)

  if (!par %in% model_pars) {
    stop2(
      "Parameter '{par}' not found in model.\n",
      "Available parameters: {paste(model_pars, collapse = ', ')}"
    )
  }

  link <- model$links[[par]] %||% "identity"
  softmax <- .is_softmax_param(par, model)

  if (!is.null(bterms$dpars) && par %in% names(bterms$dpars)) {
    type <- "dpar"
  } else if (!is.null(bterms$nlpars) && par %in% names(bterms$nlpars)) {
    type <- "nlpar"
  } else {
    # parameters not in bterms (e.g. fixed or auxiliary params) default to nlpar
    type <- "nlpar"
  }

  list(
    type = type,
    model_name = par,
    brms_name = par,
    link = link,
    softmax = softmax
  )
}


#' Check if parameter uses softmax transformation
#'
#' @param par Character string. Parameter name
#' @param model A bmmodel object
#' @return Logical. TRUE if parameter uses softmax transformation
#'
#' @keywords internal
#' @noRd
.is_softmax_param <- function(par, model) {
  "mixture3p" %in% class(model) && par %in% c("thetat", "thetant")
}


#' @title Posterior draws of model parameters on the native scale
#'
#' @description
#' `bmm` samples all parameters on their link scale (`log` for precision and
#' boundary parameters, `logit` or `softmax` for mixture weights, `tan_half` for
#' circular locations). `native_parameters()` returns posterior draws of those
#' parameters on the scale they are interpreted and reported on, evaluated over a
#' grid of predictor values.
#'
#' Because draws are returned rather than summaries, any contrast is ordinary
#' arithmetic on the draws, with the baseline made explicit by the grid.
#'
#' @param x A `bmmfit` object returned by [bmm()].
#' @param newdata Data to evaluate the parameters on. If `NULL` (the default),
#'   the unique observed combinations of the parameters' predictors are used,
#'   taken as rows of the model's own data. Columns you do not supply are filled
#'   from the first row of the model data.
#' @param pars Character vector of parameters to return. If `NULL` (the default),
#'   all parameters of the model are returned, including those fixed to a
#'   constant. This filters the output only; parameters transformed jointly are
#'   always computed together.
#' @param re_formula Which group-level effects to include, as in
#'   [brms::posterior_linpred()]. See Details.
#' @param scale Either `"native"` (the default) to apply the model's inverse link
#'   functions, or `"sampling"` to return the untransformed linear predictor.
#' @param ndraws Number of posterior draws to use. If `NULL` (the default), all
#'   draws are used. The same draws are used for every parameter.
#' @param draw_ids Indices of the draws to use. Overrides `ndraws` and makes the
#'   result reproducible.
#' @param summary Logical. If `FALSE` (the default), return draws. If `TRUE`,
#'   return posterior summaries computed *after* the transformation.
#' @param prob Probability mass of the credible interval when `summary = TRUE`.
#' @param robust Logical. If `TRUE`, `summary = TRUE` reports the median and
#'   median absolute deviation instead of the mean and standard deviation.
#' @param ... Further arguments passed to [brms::posterior_linpred()], such as
#'   `allow_new_levels` and `sample_new_levels`.
#'
#' @return If `summary = FALSE`, a `data.frame` with one row per draw, grid cell
#'   and parameter, with columns `.chain`, `.iteration`, `.draw`, the grid
#'   variables, `parameter` and `value`. The draw indices refer to the draws of
#'   `x`, so they can be joined with the output of `posterior::as_draws_df(x)`.
#'   If `summary = TRUE`, a `data.frame` with one row per grid cell and
#'   parameter, with columns `Estimate`, `Est.Error` and the interval bounds.
#'
#' @details
#' # Transform first, then summarise
#'
#' Quantiles are preserved by monotone inverse links but means are not, so the
#' median and the credible interval of a native-scale parameter are exact while
#' its mean is not the inverse link of the linear predictor's mean. This function
#' transforms the draws and only then summarises, which is why `summary = TRUE`
#' is not the same as transforming the output of [summary.bmmfit()].
#'
#' # Group-level effects
#'
#' `re_formula = NULL` (the default) returns subject-specific parameters and adds
#' the grouping variables to the grid. `re_formula = NA` sets all group-level
#' effects to zero, which under a non-identity link gives the *median* subject,
#' not the population mean. The population mean requires marginalising over the
#' distribution of group-level effects, which is done by predicting for a new
#' level:
#'
#' ```r
#' native_parameters(fit, newdata = transform(nd, id = "new"),
#'                   allow_new_levels = TRUE, sample_new_levels = "gaussian")
#' ```
#'
#' These three quantities differ substantially at between-subject standard
#' deviations typical for working memory data.
#'
#' # What is returned
#'
#' Parameters fixed to a constant are returned at that constant, transformed to
#' the native scale. Because `bmm` fixes parameters on the *link* scale, the
#' native value can differ from the value shown by [parameters()]: the `ddm`
#' relative starting point `zr` is fixed at `0` under a `logit` link and is
#' therefore reported as `0.5`, and the `ezdm` and `cswald` diffusion constant `s`
#' is fixed at `0` under a `log` link and is reported as `1`.
#'
#' Trial-level derived quantities, such as the response probabilities of the
#' mixture models or the category probabilities of `m3`, are not returned; use
#' [brms::posterior_epred()] for those.
#'
#' The default grid contains the observed combinations of the predictors, not
#' their full crossing, so cells that were never presented do not appear. A
#' continuous predictor therefore produces one grid cell per observed value;
#' supply `newdata` for such models.
#'
#' @seealso [parameters()], [native_transform()], [conditional_effects.bmmfit()]
#' @keywords extract_info
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' fit <- bmm(
#'   bmf(c ~ 0 + set_size, kappa ~ 1),
#'   data = oberauer_lin_2017,
#'   model = sdm(resp_error = "dev_rad")
#' )
#'
#' # draws of c and kappa for every set size, on the native scale
#' np <- native_parameters(fit, re_formula = NA)
#' head(np)
#'
#' # a contrast is plain arithmetic on the draws
#' c_draws <- subset(np, parameter == "c")
#' quantile(
#'   c_draws$value[c_draws$set_size == 1] - c_draws$value[c_draws$set_size == 4],
#'   probs = c(0.025, 0.5, 0.975)
#' )
#'
#' # posterior summaries of the transformed draws
#' native_parameters(fit, re_formula = NA, summary = TRUE)
native_parameters <- function(x, newdata = NULL, pars = NULL, re_formula = NULL,
                              scale = c("native", "sampling"), ndraws = NULL,
                              draw_ids = NULL, summary = FALSE, prob = 0.95,
                              robust = FALSE, ...) {
  stopif(
    !is_bmmfit(x),
    "native_parameters() requires a bmmfit object, not an object of class \\
    {collapse_comma(class(x))}."
  )
  x <- restructure(x)
  scale <- match.arg(scale)

  model_pars <- names(x$bmm$model$parameters)
  stopif(length(model_pars) == 0, "This model has no parameters defined.")
  pars <- pars %||% model_pars
  unknown <- pars[not_in(pars, model_pars)]
  stopif(
    length(unknown) > 0,
    "Unknown parameter(s) {collapse_comma(unknown)}. \\
    The parameters of this model are {collapse_comma(model_pars)}."
  )

  if (is.null(draw_ids)) {
    draw_ids <- if (is.null(ndraws)) {
      seq_len(brms::ndraws(x))
    } else {
      sort(sample.int(brms::ndraws(x), ndraws))
    }
  }

  grid_vars <- .np_grid_vars(x, model_pars, re_formula)
  newdata <- .np_newdata(x, grid_vars, newdata)

  linpred <- .np_linpred(x, model_pars, pars, newdata, re_formula, draw_ids, list(...))
  if (scale == "native") {
    transformed <- native_transform(x$bmm$model, linpred, newdata)
    stopif(
      !setequal(names(transformed), names(linpred)),
      "The native_transform() method for model '{x$bmm$model$name}' must return \\
      one element per parameter it was given."
    )
    transformed <- transformed[names(linpred)]
    stopif(
      !identical(lapply(transformed, dim), lapply(linpred, dim)),
      "The native_transform() method for model '{x$bmm$model$name}' must preserve \\
      the dimensions of the draws it was given."
    )
    linpred <- transformed
  }
  linpred <- linpred[names(linpred) %in% pars]

  grid <- newdata[, grid_vars, drop = FALSE]
  row.names(grid) <- NULL
  if (summary) {
    .np_summary(linpred, grid, prob, robust)
  } else {
    .np_long(linpred, grid, .np_draw_index(x, draw_ids))
  }
}


#' @title Transform parameter draws from the sampling scale to the native scale
#'
#' @description
#' Applies a `bmmodel`'s inverse link functions to posterior draws of its
#' parameters. This is the extension point used by [native_parameters()]: the
#' transformation is a property of the *model*, not of an individual parameter,
#' because some models map several parameters jointly (e.g. mixture weights
#' through a softmax).
#'
#' The default method covers every transformation that can be expressed through
#' a `links` declaration in a `.model_*()` constructor, so **new models normally
#' need no method at all** — declaring `links` is sufficient. Write a method only
#' when a model's transformation cannot be written as an elementwise inverse link
#' or as a single softmax group.
#'
#' @param model A `bmmodel` object.
#' @param linpred A named list of matrices of linear predictor draws, one per
#'   model parameter, each with draws in rows and prediction grid cells in
#'   columns.
#' @param data The prediction grid the draws were computed on. Available for
#'   transformations that depend on the design.
#' @param ... Currently unused.
#'
#' @return A named list of matrices with the same names and dimensions as
#'   `linpred`, on the native scale.
#'
#' @details
#' Methods for specific models follow the `bmm` S3 chain, which runs from the
#' general classes to the specific ones. A method should therefore transform the
#' parameters it owns, remove them from `linpred`, and call `NextMethod()` so the
#' remaining parameters reach the default method:
#'
#' ```r
#' native_transform.mymodel <- function(model, linpred, data, ...) {
#'   own <- .my_joint_transform(linpred[c("p1", "p2")])
#'   linpred <- linpred[not_in(names(linpred), c("p1", "p2"))]
#'   c(NextMethod(), own)
#' }
#' ```
#'
#' @seealso [native_parameters()]
#' @keywords developer
#' @export
native_transform <- function(model, linpred, data, ...) {
  UseMethod("native_transform")
}

#' @rdname native_transform
#' @export
native_transform.default <- function(model, linpred, data, ...) {
  stopif(
    !is_bmmodel(model),
    "native_transform() requires a bmmodel object, not an object of class \\
    {collapse_comma(class(model))}."
  )

  softmax_pars <- .np_softmax_pars(model, names(linpred))
  if (length(softmax_pars) > 0) {
    linpred[softmax_pars] <- .np_softmax(linpred[softmax_pars])
  }

  for (par in names(linpred)[not_in(names(linpred), softmax_pars)]) {
    link <- model$links[[par]] %||% "identity"
    linpred[[par]] <- tryCatch(
      link_transform(linpred[[par]], link, inverse = TRUE),
      error = function(e) {
        stop2(
          "Cannot transform parameter '{par}' of model '{model$name}' to the \\
          native scale: {conditionMessage(e)}
          A model whose transformation is not an elementwise inverse link must \\
          provide a native_transform() method. See ?native_transform."
        )
      }
    )
  }

  linpred
}


#' Parameters transformed jointly through a softmax
#'
#' @param model A bmmodel object
#' @param pars Character vector of parameter names to consider
#' @return Character vector of the parameters declaring a `"softmax"` link
#'
#' @keywords internal
#' @noRd
.np_softmax_pars <- function(model, pars) {
  links <- unlist(model$links)
  pars[pars %in% names(links)[links == "softmax"]]
}


#' Joint softmax over draws matrices with an implicit reference category
#'
#' @description
#' The draws-matrix form of [softmax()] applied to `c(eta, 0)`. `brms` fixes the
#' linear predictor of the unpredicted mixture component to zero
#' (`brms:::stan_mixture`), so the reference category is structurally 0 and a
#' group of one parameter reduces to `plogis()`.
#'
#' @param mats A named list of matrices of linear predictor draws
#' @return A list of matrices of the same shape, on the probability scale
#'
#' @keywords internal
#' @noRd
.np_softmax <- function(mats) {
  if (length(mats) == 0) {
    return(mats)
  }

  shift <- Reduce(pmax, mats)
  shift[shift < 0] <- 0
  exps <- lapply(mats, function(m) exp(m - shift))
  denom <- Reduce(`+`, exps) + exp(-shift)
  lapply(exps, function(e) e / denom)
}


#' Variables spanning the prediction grid for a set of parameters
#'
#' @description
#' Collects the predictors of each parameter's user formula. Random-effects
#' grouping variables are retained according to `re_formula`, so that
#' `re_formula = NA` collapses the grid across grouping levels. A variable that
#' is both a grouping variable and a population-level predictor is always kept.
#'
#' @param x A bmmfit object
#' @param pars Character vector of parameter names
#' @param re_formula `NULL`, `NA`, or a formula, as in [brms::posterior_linpred()]
#' @return Character vector of variable names
#'
#' @keywords internal
#' @noRd
.np_grid_vars <- function(x, pars, re_formula) {
  keep_all_groups <- is.null(re_formula)
  keep_groups <- if (is_formula(re_formula)) {
    .extract_re_grouping_vars(re_formula)
  } else {
    character(0)
  }

  vars <- lapply(pars, function(par) {
    f <- x$bmm$user_formula[[par]]
    if (is.null(f) || is_constant(f)) {
      return(character(0))
    }
    rhs <- stats::formula(f)[-2]
    labels <- attr(stats::terms(rhs), "term.labels")
    fe_vars <- lapply(
      labels[!grepl("|", labels, fixed = TRUE)],
      function(label) all.vars(str2lang(label))
    )
    groups <- .extract_re_grouping_vars(rhs)
    c(
      unlist(fe_vars),
      if (keep_all_groups) groups else groups[groups %in% keep_groups]
    )
  })

  unique(unlist(vars))
}


#' Resolve the prediction grid for native_parameters()
#'
#' @description
#' By default the grid is the set of unique observed combinations of `grid_vars`,
#' taken as actual rows of the model's data. Real rows already carry the addition
#' terms, matrix columns and bmm-internal columns that `brms` requires, so this
#' generalises to any model without per-model code. User supplied `newdata` is
#' completed the same way: columns the user did not provide are filled from the
#' first row of the model's data.
#'
#' @param x A bmmfit object
#' @param grid_vars Character vector of variables spanning the grid
#' @param newdata User supplied prediction data, or `NULL`
#' @return A data frame suitable for [brms::posterior_linpred()]
#'
#' @keywords internal
#' @noRd
.np_newdata <- function(x, grid_vars, newdata) {
  data <- x$data
  stopif(nrow(data) == 0, "The fitted model contains no data to predict from.")

  unknown <- grid_vars[not_in(grid_vars, names(data))]
  stopif(
    length(unknown) > 0,
    "Cannot build a prediction grid: variable(s) {collapse_comma(unknown)} \\
    are not columns of the model data."
  )

  if (!is.null(newdata)) {
    return(.np_complete_newdata(newdata, data, grid_vars))
  }

  if (length(grid_vars) == 0) {
    return(data[1L, , drop = FALSE])
  }

  cells <- data[, grid_vars, drop = FALSE]
  rows <- which(!duplicated(cells))
  data[rows[do.call(order, unname(as.list(cells[rows, , drop = FALSE])))], , drop = FALSE]
}


.np_complete_newdata <- function(newdata, data, grid_vars) {
  newdata <- as.data.frame(newdata)
  missing <- names(data)[not_in(names(data), names(newdata))]
  stopif(
    any(missing %in% grid_vars),
    "'newdata' is missing the predictor(s) \\
    {collapse_comma(missing[missing %in% grid_vars])}."
  )

  reference <- data[rep(1L, nrow(newdata)), missing, drop = FALSE]
  for (var in missing) {
    newdata[[var]] <- reference[[var]]
  }
  newdata
}


#' Linear predictor draws for each model parameter
#'
#' @description
#' Resolves whether each parameter entered `brms` as a distributional or a
#' non-linear parameter and returns its linear predictor draws.
#' [brms::posterior_linpred()] returns the *link* scale for both kinds, so no
#' distinction is needed beyond choosing the argument name.
#'
#' A parameter declared by the model but absent from the fitted formula cannot be
#' predicted. This errors when the user asked for it, and warns and drops it
#' otherwise, so that one unresolvable parameter does not make the whole model
#' unusable.
#'
#' @param x A bmmfit object
#' @param pars Character vector of all parameters to compute
#' @param requested Character vector of the parameters the user asked for
#' @param newdata The prediction grid
#' @param re_formula Passed to [brms::posterior_linpred()]
#' @param draw_ids Indices of the draws to use, identical across parameters
#' @param dots Further arguments for [brms::posterior_linpred()]
#' @return A named list of draws x grid cell matrices on the link scale
#'
#' @keywords internal
#' @noRd
.np_linpred <- function(x, pars, requested, newdata, re_formula, draw_ids, dots) {
  bterms <- brms::brmsterms(x$formula)
  types <- ifelse(
    pars %in% names(bterms$dpars), "dpar",
    ifelse(pars %in% names(bterms$nlpars), "nlpar", NA_character_)
  )
  names(types) <- pars

  unresolved <- pars[is.na(types)]
  stopif(
    any(unresolved %in% requested),
    "Parameter(s) {collapse_comma(unresolved[unresolved %in% requested])} of model \\
    '{x$bmm$model$name}' are declared by the model but do not appear in the fitted \\
    formula, so they cannot be predicted. This usually means a parameter name \\
    collides with a data column or with another parameter's formula. Please report \\
    this at https://github.com/venpopov/bmm/issues"
  )
  if (length(unresolved) > 0) {
    warning2(
      "Skipping parameter(s) {collapse_comma(unresolved)}, which are declared by \\
      the model but do not appear in the fitted formula."
    )
    pars <- pars[not_in(pars, unresolved)]
  }

  linpred <- lapply(pars, function(par) {
    args <- c(
      nlist(object = x, newdata, re_formula, draw_ids),
      stats::setNames(list(par), types[[par]]),
      dots
    )
    brms::do_call(brms::posterior_linpred, args)
  })
  stats::setNames(linpred, pars)
}


#' Chain and iteration indices of the draws used
#'
#' @description
#' Draws returned by [brms::posterior_linpred()] are ordered chain-major and in
#' the same order as [as.matrix.brmsfit()], so the chain and iteration of each
#' draw can be recovered from its index.
#'
#' @param x A bmmfit object
#' @param draw_ids Indices of the draws used
#' @return A data frame with columns `.chain`, `.iteration` and `.draw`
#'
#' @keywords internal
#' @noRd
.np_draw_index <- function(x, draw_ids) {
  n_chains <- brms::nchains(x)
  n_iter <- brms::niterations(x)
  data.frame(
    .chain = rep(seq_len(n_chains), each = n_iter)[draw_ids],
    .iteration = rep(seq_len(n_iter), times = n_chains)[draw_ids],
    .draw = draw_ids
  )
}


#' Reshape parameter draws into long format
#'
#' @param linpred Named list of draws x grid cell matrices
#' @param grid Data frame of grid cells, one row per column of the matrices
#' @param draws Data frame of `.chain`, `.iteration` and `.draw`
#' @return A long data frame, one row per draw, cell and parameter
#'
#' @keywords internal
#' @noRd
.np_long <- function(linpred, grid, draws) {
  n_draws <- nrow(draws)
  n_cells <- nrow(grid)
  n_pars <- length(linpred)

  data.frame(
    draws[rep(seq_len(n_draws), times = n_cells * n_pars), , drop = FALSE],
    grid[rep(rep(seq_len(n_cells), each = n_draws), times = n_pars), , drop = FALSE],
    parameter = rep(names(linpred), each = n_draws * n_cells),
    value = unlist(lapply(linpred, as.vector), use.names = FALSE),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


#' Summarise parameter draws per grid cell
#'
#' @inheritParams .np_long
#' @param prob Probability mass of the credible interval
#' @param robust Whether to use median and MAD instead of mean and SD
#' @return A data frame, one row per cell and parameter
#'
#' @keywords internal
#' @noRd
.np_summary <- function(linpred, grid, prob, robust) {
  summaries <- lapply(linpred, .ce_summarize_draws, prob = prob, robust = robust)
  pull <- function(element) {
    unlist(lapply(summaries, `[[`, element), use.names = FALSE)
  }

  out <- data.frame(
    grid[rep(seq_len(nrow(grid)), times = length(linpred)), , drop = FALSE],
    parameter = rep(names(linpred), each = nrow(grid)),
    Estimate = pull("estimate"),
    Est.Error = pull("se"),
    lower = pull("lower"),
    upper = pull("upper"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  names(out)[names(out) %in% c("lower", "upper")] <- paste0("Q", probs * 100)
  out
}
