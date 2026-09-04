############################################################################# !
# EXPRESSION HELPERS                                                    ####
############################################################################# !

.mpt_canonical_expr <- function(expr, tree_name) {
  stopif(
    !is.character(expr) || length(expr) != 1L || !nzchar(expr),
    "Branch probability expressions must be single character strings \\
    (tree '{tree_name}')."
  )
  parsed <- try(str2lang(expr), silent = TRUE)
  stopif(
    is_try_error(parsed),
    "Cannot parse the branch expression '{expr}' in tree '{tree_name}'."
  )
  folded <- .mpt_fold_numeric_division(parsed)
  .mpt_check_stan_constants(folded, tree_name)
  paste(deparse(folded), collapse = " ")
}

# brms deparses the formula constants when emitting Stan code and spaces out
# operators, which turns scientific notation like 6.7e-05 into the subtraction
# '6.7e - 05' and breaks the Stan compilation; constants that R deparses in
# scientific notation are therefore rejected at construction
.mpt_check_stan_constants <- function(node, tree_name) {
  if (is.numeric(node) && length(node) == 1L) {
    scientific <- withr::with_options(
      list(scipen = 0), grepl("e", deparse(node), fixed = TRUE)
    )
    stopif(
      scientific,
      "The numeric constant {node} in tree '{tree_name}' is too extreme to \\
      be written into the generated Stan code (brms emits it in scientific \\
      notation, which breaks the Stan syntax). Please provide such values as \\
      a data column declared in the covariates argument, or use a larger \\
      constant."
    )
  }
  if (is.call(node)) {
    for (i in seq_along(node)[-1]) {
      .mpt_check_stan_constants(node[[i]], tree_name)
    }
  }
  invisible(NULL)
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

.mpt_tree_vars <- function(tree) {
  unique(unlist(lapply(tree$branches, function(branch) all.vars(str2lang(branch)))))
}

# branches of each tree must sum to 1 for any parameter values; evaluating at
# several distinct test points catches swapped-complement errors that a single
# symmetric point (e.g. all 0.5) would miss
.mpt_validate_tree_sums <- function(trees, parameters, covariates, simplex,
                                    tolerance = 1e-6) {
  symbols <- c(parameters, covariates)
  if (length(symbols) == 0L) {
    return(invisible(NULL))
  }
  test_vals <- c(0.137, 0.421, 0.683, 0.852)
  for (tree in trees) {
    for (shift in seq_along(test_vals)) {
      vals <- setNames(
        test_vals[(seq_along(symbols) + shift - 2L) %% length(test_vals) + 1L],
        symbols
      )
      for (grp in simplex) {
        vals[grp] <- vals[grp] / sum(vals[grp])
      }
      total <- sum(vapply(
        tree$branches,
        function(branch) eval(str2lang(branch), envir = as.list(vals)),
        numeric(1)
      ))
      if (abs(total - 1) > tolerance) {
        warning2(
          "The branch probabilities of tree '{tree$name}' sum to \\
          {signif(total, 6)} instead of 1 when evaluated at numeric test \\
          values. Please check the branch expressions."
        )
        break
      }
    }
  }
  invisible(NULL)
}
