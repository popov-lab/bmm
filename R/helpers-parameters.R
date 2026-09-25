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
      softplus = .log1p_exp(values),
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
      softplus = .log_expm1(values),
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

#' The softplus link and its inverse, without overflow
#'
#' @description
#' `log1p(exp(x))` is `Inf` from x ~ 710 and `log(expm1(x))` from x ~ 709, where
#' both are within double precision of `x` itself. Same guard as
#' `brms:::log1p_exp()`.
#'
#' @param values A numeric vector
#' @return A numeric vector
#'
#' @keywords internal
#' @noRd
.log1p_exp <- function(values) {
  out <- log1p(exp(values))
  ifelse(out < Inf, out, values)
}

#' @noRd
.log_expm1 <- function(values) {
  out <- log(expm1(values))
  ifelse(out < Inf, out, values)
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
  bterms <- brms::brmsterms(bmmfit$formula, family = bmmfit$formula$family)

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
  length(.np_softmax_pars(model, par)) > 0
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
#' The same order is what makes `summary = TRUE` correct for the mixture weights,
#' whose softmax is not an elementwise map at all: summarising the transformed
#' draws reports quantiles of the weight's own marginal posterior, whereas
#' transforming the summaries would report the softmax of three separate
#' quantiles, which is not a quantile of anything.
#'
#' # Transforming coefficients is not the same thing
#'
#' The inverse link applies to the **linear predictor** of a grid cell, not to an
#' individual regression coefficient. For a model with `kappa ~ condition`,
#' `exp(b_kappa_conditionB)` is a multiplicative factor, not `kappa` in condition
#' B; `kappa` in condition B is `exp(b_kappa_Intercept + b_kappa_conditionB)`.
#' Building the grid is exactly what removes this step, which is why contrasts
#' are taken between rows of the output rather than read off the coefficients.
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
#' # Circular location parameters
#'
#' The circular models sample their location parameter (`mu` for `sdm`, `mu1` for
#' `mixture2p`, `mixture3p` and `imm`) through a `tan_half` link, and it is
#' returned in radians in `(-pi, pi)` — a response bias relative to the target,
#' since the response variable is the angular error. It is fixed to `0` unless
#' the `bmmformula` predicts it explicitly, so an all-zero `mu1` means the model
#' never estimated one.
#'
#' The inverse link is `2 * atan()`, and the caveat above applies to it with
#' particular force: `2 * atan(b_mu1_conditionB)` is neither the bias in
#' condition B nor the difference between conditions. The bias in condition B is
#' `2 * atan(b_mu1_Intercept + b_mu1_conditionB)`, and the difference between the
#' conditions is that value minus `2 * atan(b_mu1_Intercept)`, taken draw by
#' draw.
#'
#' # Mixture weights
#'
#' For `mixture3p` the returned `thetat` and `thetant` are probabilities but they
#' **do not sum to 1**: `brms` holds the linear predictor of one mixture
#' component at zero, and that component is the guessing distribution, so the
#' remaining mass `1 - thetat - thetant` is the probability of a guess. `thetant`
#' is the *total* probability of a non-target response, summed over the lures,
#' not the probability per lure.
#'
#' At set sizes where no non-target was presented the model switches its
#' non-target components off, so `thetant` is reported as exactly `0` and
#' `thetat` as `plogis()` of its linear predictor, rather than as the value the
#' set-size regression extrapolates to. This correction reads the lure indicators
#' out of the prediction grid; if you supply `newdata`, columns you do not supply
#' are filled from the first row of the model data rather than recomputed, and
#' `native_parameters()` warns when the two disagree.
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

  stopif(
    !is.numeric(prob) || length(prob) != 1L || is.na(prob) || prob <= 0 || prob >= 1,
    "'prob' must be a single number between 0 and 1, not {prob}."
  )
  stopif(
    !is.null(ndraws) &&
      (!is.numeric(ndraws) || length(ndraws) != 1L || is.na(ndraws) ||
        ndraws < 1 || ndraws > brms::ndraws(x)),
    "'ndraws' must be a single number between 1 and {brms::ndraws(x)}."
  )

  if (is.null(draw_ids)) {
    draw_ids <- if (is.null(ndraws)) {
      seq_len(brms::ndraws(x))
    } else {
      sample.int(brms::ndraws(x), ndraws)
    }
  }
  # brms subsets draws with posterior::subset_draws(), which returns them in
  # ascending id order whatever order it was given, so the index columns have to
  # be built from the same ordering or every row is labelled with another draw
  draw_ids <- sort(draw_ids)

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

  clash <- intersect(names(grid), .np_reserved_names(summary, prob))
  stopif(
    length(clash) > 0,
    "Predictor(s) {collapse_comma(clash)} share a name with a column of the \\
    output. Rename them in the data, or pass 'newdata' with different names."
  )

  if (summary) {
    .np_summary(linpred, grid, prob, robust)
  } else {
    warnif(
      length(linpred) * nrow(grid) * length(draw_ids) > 1e6,
      "This will return \\
      {length(linpred) * nrow(grid) * length(draw_ids)} rows. Supply 'newdata' \\
      with the grid cells you need, or set 'ndraws', to return fewer."
    )
    .np_long(linpred, grid, .np_draw_index(x, draw_ids))
  }
}


#' Output column names that a grid variable must not collide with
#'
#' @description
#' `data.frame()` silently disambiguates duplicate names by suffixing, and the
#' quantile columns are renamed after construction, so a predictor named `value`
#' or `Q2.5` would shadow the output rather than clash with it.
#'
#' @param summary Whether the summary or the draws format is returned
#' @param prob Probability mass of the credible interval
#' @return Character vector of reserved column names
#'
#' @keywords internal
#' @noRd
.np_reserved_names <- function(summary, prob) {
  if (!summary) {
    return(c(".chain", ".iteration", ".draw", "parameter", "value"))
  }
  c("parameter", "Estimate", "Est.Error", .np_interval_names(prob))
}


#' Names of the quantile columns of a summary
#'
#' @param prob Probability mass of the credible interval
#' @return Character vector of two column names
#'
#' @keywords internal
#' @noRd
.np_interval_names <- function(prob) {
  paste0("Q", c((1 - prob) / 2, 1 - (1 - prob) / 2) * 100)
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
#' **A model whose parameters transform elementwise, or through one softmax group
#' that is active on every row of the data, needs no method** — declaring `links`
#' in its `.model_*()` constructor is sufficient, and that covers every model
#' `bmm` currently ships. Write a method when either condition fails. Two cases
#' that are known to fail it: a model with several independent multinomial
#' branches, because a `"softmax"` link carries no group identity and all such
#' parameters are softmaxed together as one group; and a transformation that
#' depends on the design rather than only on the parameter's own value.
#'
#' @param model A `bmmodel` object.
#' @param linpred A named list of matrices of linear predictor draws, one per
#'   model parameter, each with draws in rows and prediction grid cells in
#'   columns. These are [brms::posterior_linpred()] output, so they are on the
#'   *link* scale for distributional parameters and are the raw value for
#'   non-linear parameters. For `bmm`'s `nlpar`-based models the link is a
#'   fiction maintained by the model's own `nlf()` expression — `imm` declares
#'   `c = "log"` because `configure_model.imm_abc()` writes `exp(c)` into the
#'   formula, not because `brms` applies a link — so a model that declares a link
#'   it does not actually apply in its `nlf()` will get silently wrong output
#'   here.
#' @param data The prediction grid the draws were computed on, with one row per
#'   column of the `linpred` matrices. These are rows of the fitted model's own
#'   data, so they carry the `bmm`-internal columns a design-dependent
#'   transformation keys on (`LureIdx*`, `inv_ss`, `Idx_*`, matrix columns).
#'   `native_transform.non_targets()` is the worked example.
#' @param ... Currently unused.
#'
#' @return A named list of matrices on the native scale. A method **must** return
#'   one element for every name it was given and preserve each matrix's
#'   dimensions; [native_parameters()] checks both and errors, naming the model,
#'   if either is broken. The order of the list is free.
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
#' Removing them before `NextMethod()` is what keeps a method idempotent with
#' respect to the parameters it does not own: a parameter left in `linpred` is
#' transformed a second time by the next method in the chain, which does not
#' error and does not warn.
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
    declared <- .np_softmax_pars(model, names(model$links))
    stopif(
      !setequal(softmax_pars, declared),
      "Parameter(s) {collapse_comma(setdiff(declared, softmax_pars))} of model \\
      '{model$name}' are part of a softmax group but their draws are missing, so \\
      the remaining weights cannot be computed. Transforming the group without \\
      them would silently return a different quantity."
    )
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

#' @rdname native_transform
#' @export
native_transform.non_targets <- function(model, linpred, data, ...) {
  # configure_model.mixture3p() gates each non-target component on its LureIdx, so
  # where no lure was presented the likelihood holds it at -100 and the weight is
  # zero rather than whatever the set-size regression extrapolates to
  non_target <- intersect(.np_softmax_pars(model, names(linpred)), "thetant")
  if (length(non_target) == 0) {
    return(NextMethod())
  }
  linpred[[non_target]][, .np_lure_free_rows(model, data)] <- -Inf
  NextMethod()
}


#' Grid cells in which no non-target was presented
#'
#' @description
#' The non-target mixture components are switched off in the likelihood wherever
#' every `LureIdx` is zero (`check_data.non_targets()` sets them from the set
#' size), so on those rows the model's non-target weight is exactly zero rather
#' than the value its regression coefficient extrapolates to.
#'
#' @param model A bmmodel object
#' @param data The prediction grid
#' @return Logical vector, one element per row of `data`
#'
#' @keywords internal
#' @noRd
.np_lure_free_rows <- function(model, data) {
  set_size_var <- model$other_vars$set_size
  if (is_data_var(set_size_var, data)) {
    return(as_numeric_vector(data[[set_size_var]]) == 1)
  }
  if ("ss_numeric" %in% names(data)) {
    return(data$ss_numeric == 1)
  }
  lure_cols <- grep("^LureIdx", names(data), value = TRUE)
  if (length(lure_cols) == 0) {
    return(rep(FALSE, nrow(data)))
  }
  rowSums(data[, lure_cols, drop = FALSE]) == 0
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
  names(out)[names(out) %in% c("lower", "upper")] <- .np_interval_names(prob)
  out
}
