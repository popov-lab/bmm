############################################################################# !
# EXPRESSION HELPERS                                                    ####
############################################################################# !
# Branch expressions are stored as parsed R calls on the mpt_tree object. The
# helpers here operate on those calls and never on strings; strings are only
# produced for printing via deparse1().

# brms deparses the formula constants when emitting Stan code and spaces out
# operators, which turns scientific notation like 6.7e-05 into the subtraction
# '6.7e - 05' and breaks the Stan compilation; constants that R deparses in
# scientific notation are therefore collected so mpt_tree() can reject them
.mpt_scientific_constants <- function(expr) {
  if (is.numeric(expr) && length(expr) == 1L) {
    scientific <- withr::with_options(
      list(scipen = 0), grepl("e", deparse(expr), fixed = TRUE)
    )
    return(if (scientific) expr else numeric(0))
  }
  if (!is.call(expr)) {
    return(numeric(0))
  }
  unlist(lapply(as.list(expr)[-1], .mpt_scientific_constants)) %||% numeric(0)
}

# Stan compiles a bare numeric fraction like 1/4 as integer division (= 0),
# so constant divisions are folded into decimal literals before emission
.mpt_fold_numeric_division <- function(expr) {
  if (!is.call(expr)) {
    return(expr)
  }
  for (i in seq_along(expr)[-1]) {
    expr[[i]] <- .mpt_fold_numeric_division(expr[[i]])
  }
  if (identical(expr[[1]], quote(`(`)) && is.numeric(expr[[2]])) {
    return(expr[[2]])
  }
  if (identical(expr[[1]], quote(`/`)) && length(expr) == 3L &&
        is.numeric(expr[[2]]) && is.numeric(expr[[3]])) {
    return(expr[[2]] / expr[[3]])
  }
  expr
}

.mpt_expr_vars <- function(tree) {
  unique(unlist(lapply(tree$branches, all.vars)))
}

.mpt_substitute_symbols <- function(expr, values) {
  .mpt_fold_numeric_division(do.call(substitute, list(expr, values)))
}

.mpt_apply_restrictions <- function(tree, restrictions) {
  tree$branches[] <- lapply(tree$branches, .mpt_substitute_symbols, restrictions)
  tree
}

# restrictions in the MPTinR/TreeBUGS string syntax ("Dn = Do", "g = 0.5",
# "G1 = G2 = G3") or as a named list (list(Dn = "Do", g = 0.5)) become one
# named list mapping each restricted parameter to a symbol or a number
.mpt_parse_restrictions <- function(restrictions) {
  if (length(restrictions) == 0L) {
    return(list())
  }
  if (is.null(names(restrictions))) {
    stopif(
      !is.character(restrictions),
      "The restrictions argument must be a character vector such as \\
      c('Dn = Do', 'g = 0.5') or a named list such as list(Dn = 'Do', g = 0.5)."
    )
    return(do.call(c, lapply(restrictions, .mpt_parse_restriction_string)))
  }
  stopif(
    any(!nzchar(names(restrictions))),
    "Every element of a named restrictions list must be named after the \\
    parameter it restricts."
  )
  values <- lapply(names(restrictions), function(par) {
    value <- restrictions[[par]]
    if (is.character(value)) {
      parsed <- try(str2lang(value), silent = TRUE)
      stopif(
        is_try_error(parsed),
        "Cannot parse the restriction '{par} = {value}'."
      )
      value <- parsed
    }
    .mpt_restriction_value(value, glue("{par} = {deparse1(value)}"))
  })
  setNames(values, names(restrictions))
}

.mpt_parse_restriction_string <- function(text) {
  expr <- try(str2lang(text), silent = TRUE)
  stopif(is_try_error(expr), "Cannot parse the restriction '{text}'.")
  lhs <- character(0)
  while (is.call(expr) && identical(expr[[1]], quote(`=`))) {
    stopif(
      !is.symbol(expr[[2]]),
      "The left-hand side of the restriction '{text}' must be a parameter name."
    )
    lhs <- c(lhs, as.character(expr[[2]]))
    expr <- expr[[3]]
  }
  stopif(
    is.call(expr) && as.character(expr[[1]]) %in% c("<", ">", "<=", ">="),
    "Order constraints such as '{text}' are not supported by the \\
    restrictions argument. Reparameterize the larger parameter instead, e.g. \\
    Do ~ Dn + (1 - Dn) * inv_logit(phi) in the model formula; see the section \\
    'Ordered parameter constraints' of the MPT article."
  )
  stopif(
    length(lhs) == 0L,
    "Each restriction must have the form 'parameter = parameter' or \\
    'parameter = constant', not '{text}'."
  )
  setNames(rep(list(.mpt_restriction_value(expr, text)), length(lhs)), lhs)
}

# a right-hand side without symbols is a constant (1/4, 1 - 0.75); a single
# symbol equates two parameters; anything else has no MPT interpretation
.mpt_restriction_value <- function(expr, text) {
  if (is.symbol(expr)) {
    return(expr)
  }
  if (length(all.vars(expr)) == 0L) {
    value <- try(eval(expr, envir = baseenv()), silent = TRUE)
    stopif(
      is_try_error(value) || !is.numeric(value) || length(value) != 1L ||
        is.na(value),
      "The restriction '{text}' does not evaluate to a single number."
    )
    return(value)
  }
  stop2(
    "Restrictions can equate a parameter with another parameter or fix it to \\
    a numeric constant. '{text}' does neither."
  )
}

# chains such as A = B, B = C resolve to the final target; a cycle leaves a
# parameter pointing at a restricted name, which mpt() reports
.mpt_resolve_restrictions <- function(restrictions) {
  for (i in seq_along(restrictions)) {
    restrictions <- lapply(restrictions, function(value) {
      if (is.symbol(value) && as.character(value) %in% names(restrictions)) {
        restrictions[[as.character(value)]]
      } else {
        value
      }
    })
  }
  restrictions
}

.mpt_eval_branches <- function(tree, values) {
  vapply(tree$branches, function(branch) eval(branch, envir = values), numeric(1))
}

# branches of each tree must sum to 1 for any parameter values; evaluating at
# several distinct test points catches swapped-complement errors that a single
# symmetric point (e.g. all 0.5) would miss. Returns the first deviating
# branch sum per tree, NA where every test point sums to 1.
.mpt_tree_sum_deviations <- function(trees, parameters, covariates, simplex,
                                     tolerance = 1e-6) {
  symbols <- c(parameters, covariates)
  test_vals <- c(0.137, 0.421, 0.683, 0.852)
  vapply(trees, function(tree) {
    for (shift in seq_along(test_vals)) {
      vals <- setNames(
        test_vals[(seq_along(symbols) + shift - 2L) %% length(test_vals) + 1L],
        symbols
      )
      for (grp in simplex) {
        vals[grp] <- vals[grp] / sum(vals[grp])
      }
      total <- sum(.mpt_eval_branches(tree, as.list(vals)))
      if (abs(total - 1) > tolerance) {
        return(total)
      }
    }
    NA_real_
  }, numeric(1))
}
