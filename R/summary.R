#' Create a summary of a fitted model represented by a \code{bmmfit} object
#'
#' @inheritParams brms::summary.brmsfit
#' @param backend Choose whether to display the *bmm* summary method (default),
#'   or to display the *brms* summary method.
#' @seealso \code{\link[brms]{summary.brmsfit}}
#' @note You can turn off the color output by setting the option
#' options(bmm.color_summary = FALSE) or bmm_options(color_summary = FALSE)
#'
#' @return A list of class \code{bmmsummary} containing the summary of the model's
#'  parameters, the model formula, the model, and the data used to fit the model.
#'
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # generate artificial data from the Signal Discrimination Model
#' dat <- data.frame(y = rsdm(2000))
#'
#' # define formula
#' ff <- bmmformula(c ~ 1, kappa ~ 1)
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = sdm(resp_error = "y"),
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # summary of the model
#' summary(fit)
summary.bmmfit <- function(object, priors = FALSE, prob = 0.95, robust = FALSE, mc_se = FALSE, ..., backend = "bmm") {
  object <- restructure(object)
  backend <- match.arg(backend, c("bmm", "brms"))

  # get summary object from brms, since it contains a lot of necessary information:
  out <- NextMethod("summary")
  if (backend == "brms") {
    return(out)
  }

  out <- rename_mu_smry(out, get_mu_pars(object))

  # add bmm-specific info to smry object
  out$formula <- object$bmm$user_formula
  out$model <- object$bmm$model
  out$data <- object$data

  class(out) <- "bmmsummary"
  out
}


#' @export
print.bmmsummary <- function(x, digits = 2, color = getOption("bmm.color_summary", TRUE), ...) {
  options(bmm.color_summary = color)
  print_model_header(x$model, x$formula, attr(x$data, "data_name"), nrow(x$data))
  print_summary_draws(x)
  print_summary_random(x, digits)
  if (nrow(x$fixed)) {
    print_summary_fixed(.summary_fixed_rows(x$fixed, select_pars(x)), digits)
  }
  print_summary_footer(x)
  invisible(x)
}

# The blocks of a printed summary, shared by the univariate and the
# multivariate summary methods (the latter prints one header and one
# coefficient block per component)
print_model_header <- function(model, formula, data_name, nobs, links = TRUE) {
  cat(style("purple1")("  Model: "))
  cat(summarise_model(model, newline = TRUE, wsp = 9), "\n")
  if (links) {
    cat(style("purple1")("  Links: "))
    cat(summarise_links(model$links), "\n")
  }
  cat(style("purple1")("Formula: "))
  cat(summarise_formula.bmmformula(formula, newline = TRUE, wsp = 9, model = model), "\n")
  cat(style("purple1")("   Data:"), data_name, "(Number of observations:", paste0(nobs, ")\n"))
}

print_summary_draws <- function(x) {
  total_ndraws <- ceiling((x$iter - x$warmup) / x$thin * x$chains)
  cat(paste0(
    style("purple1")("  Draws: "), x$chains, " chains, each with iter = ", x$iter,
    "; warmup = ", x$warmup, "; thin = ", x$thin, ";\n",
    "         total post-warmup draws = ", total_ndraws, "\n\n"
  ))
}

print_summary_random <- function(x, digits) {
  if (length(x$random)) {
    cat(style("green")("Multilevel Hyperparameters:\n"))
    for (g in names(x$random)) {
      cat(paste0("~", g, " (Number of levels: ", x$ngrps[[g]], ") \n"))
      re <- x$random[[g]]
      print_format(re[!is.na(re$Rhat), , drop = FALSE], digits)
      cat("\n")
    }
  }
}

# constant parameters (fixed via the formula) have no Rhat and are listed
# separately with their value only
print_summary_fixed <- function(rows, digits, label = "") {
  is_constant <- is.na(rows$Rhat)
  if (any(!is_constant)) {
    cat(style("green")(paste0("Regression Coefficients", label, ":\n")))
    print_format(rows[!is_constant, , drop = FALSE], digits)
    cat("\n")
  }
  if (any(is_constant)) {
    cat(style("green")(paste0("Constant Parameters", label, ":\n")))
    constants <- rows[is_constant, , drop = FALSE]
    print_format(
      data.frame(Value = constants[, 1], row.names = paste0(rownames(constants), "    ")),
      digits
    )
    cat("\n")
  }
}

print_summary_footer <- function(x) {
  cat(paste0("Draws were sampled using ", x$sampler, ". "))
  if (x$algorithm == "sampling") {
    cat(
      glue(
        "For each parameter, Bulk_ESS
        and Tail_ESS are effective sample size measures, and Rhat is the potential
        scale reduction factor on split chains (at convergence, Rhat = 1)."
      )
    )
  }
  cat("\n")
}

rename_mu_smry <- function(x, mu_pars) {
  for (i in seq_along(x)) {
    if (is.data.frame(x[[i]])) {
      rownames(x[[i]])[rownames(x[[i]]) %in% mu_pars] <- paste0("mu_", mu_pars)
    }
  }
  x
}

select_pars <- function(x) {
  model_pars <- names(x$model$parameters)
  provided_dpars <- names(x$formula)[!is_nl(x$formula)]
  union(model_pars, provided_dpars)
}

# Fixed-effect rows whose brms prefix (the parameter name before the first
# underscore) is among the printed parameters. The prefix is tokenized and
# compared exactly — not via an unanchored grepl() — so "a_" does not also
# select rows of another parameter such as "kappa_Intercept" (#379). Plain
# data.frame subsetting with drop = FALSE also keeps a single fixed-effect row
# intact, which broke the old sapply+apply approach for single-coefficient
# models such as gumbel-min sdt_ranking (dprime ~ 1) (#369).
.summary_fixed_rows <- function(fixed, pars) {
  fixed[sub("_.*$", "", rownames(fixed)) %in% pars, , drop = FALSE]
}

summarise_links <- function(links) {
  out <- paste0(names(links), " = ", links)
  paste(out, sep = "", collapse = "; ")
}

summarise_formula.bmmformula <- function(formula, newline = TRUE, wsp = 0, model = NULL) {
  fixpars <- NULL
  if (!is.null(model)) {
    formula <- suppressMessages(add_missing_parameters(model, formula))
    fixpars <- model$fixed_parameters
    fixpars <- fixpars[names(fixpars) %in% names(model$parameters)]
    formula[names(fixpars)] <- fixpars
  }
  print(formula, newline = newline, wsp = wsp)
}

#' @export
print.bmmformula <- function(x, newline = TRUE, wsp = 0, ...) {
  wspace <- collapse(rep(" ", wsp))
  sep <- paste0(ifelse(newline, "\n", ","), wspace)
  for (i in seq_along(x)) {
    if (is.numeric(x[[i]])) {
      x[[i]] <- paste0(names(x)[i], " = ", x[[i]])
    }
  }
  cat(paste0(x, collapse = sep))
}

summarise_model <- function(model, ...) {
  UseMethod("summarise_model")
}

# TODO: build this up
#' @export
summarise_model.bmmodel <- function(model, ...) {
  construct_model_call(model, ...)
}

# creates a string representation of the call to the model with the user variables
# will likely also depend on the theme users have
construct_model_call <- function(model, newline = FALSE, wsp = 1, ...) {
  call <- attr(model, "call")
  if (is.null(call)) {
    classes <- class(model)
    model_name <- classes[length(classes)]

    # construct the inner part of the call
    allvars <- c(model$resp_vars, model$other_vars)
  } else {
    model_name <- deparse(call[[1]])
    allvars <- as.list(call)[-1]
  }
  vnames <- names(allvars)

  # create necessary padding
  wspace <- collapse(rep(" ", wsp + nchar(model_name) + 1))
  sep <- paste0(ifelse(newline, ",\n", ","), wspace)

  args <- sapply(vnames, function(var) {
    paste0(var, " = ", paste0(deparse(allvars[[var]]), collapse = ", "))
  })

  # combine with the padding and the function name
  args <- paste0(args, collapse = sep)
  paste0(style("coral2")(model_name), "(", args, ")")
}

# helper function to print summary matrices in nice format
# also displays -0.00 as a result of round negative values to zero (#263)
# @param x object to be printed; coerced to matrix
# @param digits number of digits to show
# @param no_digits names of columns for which no digits should be shown
# exported from brms
print_format <- function(x, digits = 2, no_digits = c("Bulk_ESS", "Tail_ESS")) {
  x <- as.matrix(x)
  digits <- as.numeric(digits)
  if (length(digits) != 1L) {
    stop2("'digits' should be a single numeric value.")
  }
  out <- x
  fmt <- paste0("%.", digits, "f")
  for (i in seq_cols(x)) {
    if (isTRUE(colnames(x)[i] %in% no_digits)) {
      out[, i] <- sprintf("%.0f", x[, i])
    } else {
      out[, i] <- sprintf(fmt, x[, i])
    }
  }
  print(out, quote = FALSE, right = TRUE)
  invisible(x)
}

style <- function(...) {
  if (getOption("bmm.color_summary", TRUE)) {
    crayon::make_style(...)
  } else {
    function(x) x
  }
}
