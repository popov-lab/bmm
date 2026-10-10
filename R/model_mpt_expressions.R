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

# the constants as the user wrote them in the branch string, each followed by
# the form brms writes when the two differ; a value folded from constant
# arithmetic (1/10000) has no written form and is shown as brms writes it
.mpt_written_constants <- function(expr, values) {
  as_brms_writes <- function(value) {
    withr::with_options(list(scipen = 0), deparse(value))
  }
  tokens <- utils::getParseData(parse(text = expr, keep.source = TRUE))
  written <- tokens$text[tokens$token == "NUM_CONST"]
  written <- unique(written[suppressWarnings(as.numeric(written)) %in% values])
  written <- c(
    written,
    vapply(setdiff(values, as.numeric(written)), as_brms_writes, character(1))
  )
  emitted <- vapply(as.numeric(written), as_brms_writes, character(1))
  paste0("'", written, "'", ifelse(written == emitted, "", paste0(" (", emitted, ")")))
}

# Stan compiles a bare integer fraction like 1/4 or 1/(2*2) as integer
# division (= 0), so every variable-free arithmetic subexpression is folded
# into its numeric value before emission
.mpt_fold_numeric_division <- function(expr) {
  if (!is.call(expr)) {
    return(expr)
  }
  for (i in seq_along(expr)[-1]) {
    expr[[i]] <- .mpt_fold_numeric_division(expr[[i]])
  }
  operands <- as.list(expr)[-1]
  if (is.symbol(expr[[1]]) &&
        as.character(expr[[1]]) %in% c("(", "+", "-", "*", "/", "^") &&
        length(operands) %in% 1:2 &&
        all(vapply(operands, is.numeric, logical(1)))) {
    return(do.call(as.character(expr[[1]]), operands))
  }
  expr
}

.mpt_expr_vars <- function(tree) {
  unique(unlist(lapply(tree$branches, all.vars)))
}

.mpt_tree_parameters <- function(trees) {
  unique(unlist(lapply(trees, .mpt_expr_vars)))
}

.mpt_eval_branches <- function(tree, values) {
  vapply(tree$branches, function(branch) eval(branch, envir = values), numeric(1))
}

# branches of each tree must sum to 1 for any parameter values; evaluating at
# several distinct test points catches swapped-complement errors that a single
# symmetric point (e.g. all 0.5) would miss. A golden-ratio sequence gives
# every symbol its own value at every point, however many symbols there are
# (a short cycle of values would equate symbols a fixed number of places apart).
# The first four points are linear in the symbol index, so v_i + v_j equals
# v_k + v_l at all four whenever i + j = k + l; the fifth raises the index to
# the irrational power sqrt(2), which no such linear relation survives, so the
# rank check cannot meet the same coincidence at every point
.mpt_test_points <- function(symbols) {
  index <- seq_along(symbols)
  lapply(1:5, function(point) {
    position <- if (point < 5) {
      index * 0.6180339887 + point * 0.2718281828
    } else {
      index^1.4142135624 * 0.6180339887
    }
    setNames(0.05 + 0.9 * (position %% 1), symbols)
  })
}

# Local identifiability: the rank of the Jacobian of all category probabilities
# (all trees stacked) with respect to the free parameters, at the interior test
# points; fixed parameters enter at their values. The maximum rank over the
# points is reported, so one point on a singular set cannot raise a false alarm,
# and it is capped at the degrees of freedom, which no Jacobian can exceed.
# Derivatives are symbolic (stats::D() on the stored calls) rather than finite
# differences: MPT branches are nearly always polynomials, which D()
# differentiates exactly, so a null direction has a singular value at rounding
# level (about 1e-16 relative) instead of the 1e-13 to 1e-11 left by central
# differences. A branch that D() cannot differentiate (a function outside its
# table) makes the check unavailable, never a false result.
# Columns are scaled to unit norm, so a parameter that moves the probabilities
# little at a test point is not mistaken for a redundant one. A parameter that
# cancels from every branch can leave a rounding residue instead of an exact
# zero (D(a * b * c * q + c * b * a * (1 - q), "q") = a * b * c - c * b * a),
# which scaling would blow up into a full column, so a column at most 1e-12 of
# the largest column is a zero column. The largest column is taken over the
# fixed parameters too, so a residue whose partners are all fixed still has a
# genuine column to be measured against. In every case measured, residues stay
# at about 1e-16 of the largest column and the genuine columns of a 32-deep
# chain at about 3e-12 or more. Among the points of maximum rank, the one with
# the fewest zero columns names the parameters.
# The relative singular-value tolerance 1e-8 sits eight orders above the null
# singular values and six below the smallest singular value of the identified
# imports (hybrid.eqn and unitization.eqn with their design constants fixed:
# 0.015 and 0.43 relative). Parameters whose unit vector is not orthogonal to
# the null space take part in a non-identified combination; zero columns are
# reported separately as parameters that appear not to affect any category
# probability.
.mpt_jacobian_rank <- function(trees, free, fixed = list(), tolerance = 1e-8) {
  # the counts are reported even when the rank cannot be computed
  counts <- list(
    n_free = length(free), free = free,
    df = sum(lengths(lapply(trees, `[[`, "branches")) - 1L)
  )
  if (length(free) == 0L) {
    return(c(counts, list(
      rank = 0L, involved = character(0), absent = character(0)
    )))
  }
  branches <- unlist(lapply(unname(trees), `[[`, "branches"), use.names = FALSE)
  parameters <- c(free, names(fixed))
  derivs <- try(unlist(lapply(branches, function(branch) {
    lapply(parameters, function(par) {
      if (par %in% all.vars(branch)) stats::D(branch, par) else 0
    })
  }), recursive = FALSE), silent = TRUE)
  if (is_try_error(derivs)) {
    return(c(counts, list(error = conditionMessage(attr(derivs, "condition")))))
  }
  symbols <- .mpt_tree_parameters(trees)
  points <- .mpt_test_points(symbols)
  # one vector per symbol holding its value at every test point, so each
  # derivative is evaluated once for all points
  vals <- lapply(setNames(nm = symbols), function(symbol) {
    vapply(points, `[[`, numeric(1), symbol)
  })
  vals[names(fixed)] <- fixed
  values <- vapply(derivs, function(deriv) {
    rep_len(eval(deriv, vals), length(points))
  }, numeric(length(points)))
  decompositions <- lapply(seq_along(points), function(point) {
    jacobian <- matrix(values[point, ], ncol = length(parameters), byrow = TRUE)
    norms <- sqrt(colSums(jacobian^2))
    zero <- norms[seq_along(free)] <= 1e-12 * max(norms)
    jacobian <- jacobian[, seq_along(free), drop = FALSE]
    norms <- norms[seq_along(free)]
    jacobian[, zero] <- 0
    norms[zero] <- 1
    decomposition <- svd(
      sweep(jacobian, 2, norms, "/"), nu = 0, nv = length(free)
    )
    decomposition$rank <- min(
      sum(decomposition$d > tolerance * decomposition$d[1]), counts$df
    )
    decomposition$absent <- free[zero]
    decomposition
  })
  best <- decompositions[[order(
    -vapply(decompositions, `[[`, integer(1), "rank"),
    lengths(lapply(decompositions, `[[`, "absent"))
  )[1]]]
  null_space <- best$v[, seq_along(free) > best$rank, drop = FALSE]
  c(counts, list(
    rank = best$rank,
    involved = free[rowSums(abs(null_space) > 1e-6) > 0],
    absent = best$absent
  ))
}

# the first branch-sum or branch-range violation per tree, NA where every test
# point gives branches in (0, 1] that sum to 1. Branches can sum to 1 for every
# value and still leave (0, 1] (2 * a and 1 - 2 * a); a branch below or at 0 is
# log(p) of a non-positive number in Stan. All test points are interior, so an
# exact 0 is a branch that is 0 for every value, e.g. (1 - a) * 0, or one that
# underflows at a test value, e.g. (1 - D)^392 at D = 0.851; refusing the
# latter too is accepted, as no realistic tree has such a power
.mpt_tree_branch_errors <- function(trees, parameters, tolerance = 1e-6) {
  points <- .mpt_test_points(parameters)
  vapply(trees, function(tree) {
    for (vals in points) {
      probs <- .mpt_eval_branches(tree, as.list(vals))
      if (abs(sum(probs) - 1) > tolerance) {
        return(glue(
          "The branch probabilities of tree '{tree$name}' sum to \\
          {signif(sum(probs), 6)} instead of 1 when evaluated at numeric test \\
          values. Please check the branch expressions. To equate parameters, \\
          give them the same name in the branch expressions or tie them in the \\
          formula (e.g. Dn ~ Do)."
        ))
      }
      outside <- names(probs)[probs <= 0 | probs > 1 + tolerance][1]
      if (!is.na(outside)) {
        symbols <- all.vars(tree$branches[[outside]])
        zero_hint <- if (probs[[outside]] == 0) {
          " A zero probability makes the likelihood undefined."
        } else {
          ""
        }
        return(glue(
          "The branch probability of category '{outside}' in tree \\
          '{tree$name}' is {signif(probs[[outside]], 6)} at the test values \\
          {paste(symbols, '=', signif(vals[symbols], 3), collapse = ', ')}, \\
          outside (0, 1]. Please check the branch expressions.{zero_hint}"
        ))
      }
    }
    NA_character_
  }, character(1))
}
