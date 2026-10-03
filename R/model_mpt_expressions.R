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

.mpt_eval_branches <- function(tree, values) {
  vapply(tree$branches, function(branch) eval(branch, envir = values), numeric(1))
}

# branches of each tree must sum to 1 for any parameter values; evaluating at
# several distinct test points catches swapped-complement errors that a single
# symmetric point (e.g. all 0.5) would miss. Returns the first deviating
# branch sum per tree, NA where every test point sums to 1.
.mpt_tree_sum_deviations <- function(trees, parameters, tolerance = 1e-6) {
  test_vals <- c(0.137, 0.421, 0.683, 0.852)
  vapply(trees, function(tree) {
    for (shift in seq_along(test_vals)) {
      vals <- setNames(
        test_vals[(seq_along(parameters) + shift - 2L) %% length(test_vals) + 1L],
        parameters
      )
      total <- sum(.mpt_eval_branches(tree, as.list(vals)))
      if (abs(total - 1) > tolerance) {
        return(total)
      }
    }
    NA_real_
  }, numeric(1))
}
