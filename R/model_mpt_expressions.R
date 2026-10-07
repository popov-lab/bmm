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

# the constants as the user wrote them in the branch string; a value folded
# from constant arithmetic (1/10000) has no written form and is deparsed
.mpt_written_constants <- function(expr, values) {
  tokens <- utils::getParseData(parse(text = expr, keep.source = TRUE))
  written <- tokens$text[tokens$token == "NUM_CONST"]
  written <- unique(written[suppressWarnings(as.numeric(written)) %in% values])
  c(written, vapply(setdiff(values, as.numeric(written)), deparse, character(1)))
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

.mpt_tree_parameters <- function(trees, covariates = NULL) {
  setdiff(unique(unlist(lapply(trees, .mpt_expr_vars))), covariates)
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
    # MPTinR and TreeBUGS pass inline restrictions as an unnamed list of strings
    if (is.list(restrictions) && all(vapply(restrictions, is.character, logical(1)))) {
      restrictions <- unlist(restrictions)
    }
    stopif(
      !is.character(restrictions),
      "The restrictions argument must be a character vector such as \\
      c('Dn = Do', 'g = 0.5') or a named list such as list(Dn = 'Do', g = 0.5)."
    )
    # readLines() returns blank and comment lines of a restriction file as is
    restrictions <- restrictions[!grepl("^\\s*(#|$)", restrictions)]
    return(Reduce(
      c, lapply(restrictions, .mpt_parse_restriction_string), list()
    ))
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
  # a restriction has an `=`, so a line without one that names a file is a path
  # (the extension is judged without a trailing comment)
  stopif(
    !grepl("=", text, fixed = TRUE) &&
      (file.exists(text) ||
        grepl(
          "\\.(restr|txt)\\s*$", sub("#.*$", "", text), ignore.case = TRUE
        )),
    "A restriction must not be a file path. bmm does not read restriction \\
    files; pass readLines({encodeString(text, quote = '\"')}) instead (blank \\
    and comment lines are skipped)."
  )
  expr <- try(str2lang(text), silent = TRUE)
  # R cannot parse a chain such as G1 < G2 < G3
  if (is_try_error(expr) && grepl("[<>]", text)) {
    .mpt_stop_order_constraint(text)
  }
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
  if (is.call(expr) && as.character(expr[[1]]) %in% c("<", ">", "<=", ">=")) {
    .mpt_stop_order_constraint(text)
  }
  stopif(
    length(lhs) == 0L,
    "Each restriction must have the form 'parameter = parameter' or \\
    'parameter = constant', not '{text}'."
  )
  setNames(rep(list(.mpt_restriction_value(expr, text)), length(lhs)), lhs)
}

.mpt_stop_order_constraint <- function(text) {
  stop2(
    "Order constraints such as '{text}' are not supported by the \\
    restrictions argument. Reparameterize the larger parameter instead, e.g. \\
    Do ~ Dn + (1 - Dn) * inv_logit(phi) in the model formula; see the section \\
    'Ordered parameter constraints' of the MPT article."
  )
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
    stopif(
      value < 0 || value > 1,
      "The restriction '{text}' fixes a parameter to {value}. A constant must \\
      lie strictly between 0 and 1."
    )
    stopif(
      value == 0 || value == 1,
      "The restriction '{text}' fixes a parameter to {value}. Constants must \\
      be strictly between 0 and 1: a constant 0 or 1 can turn a branch \\
      into 0, which has no log probability. To model a parameter at 0 or 1, \\
      remove it from the branch expressions (write the reduced tree) and \\
      declare the categories that no branch reaches with \\
      mpt_tree(impossible = )."
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
# symmetric point (e.g. all 0.5) would miss. A golden-ratio sequence gives
# every symbol its own value at every point, however many symbols there are
# (a short cycle of values would equate symbols a fixed number of places apart).
# The first four points are linear in the symbol index, so v_i + v_j equals
# v_k + v_l at all four whenever i + j = k + l; the fifth raises the index to
# the irrational power sqrt(2), which no such linear relation survives, so the
# rank check cannot meet the same coincidence at every point. The members of a
# simplex group are rescaled to sum to 1.
.mpt_test_points <- function(symbols, simplex) {
  index <- seq_along(symbols)
  lapply(1:5, function(point) {
    position <- if (point < 5) {
      index * 0.6180339887 + point * 0.2718281828
    } else {
      index^1.4142135624 * 0.6180339887
    }
    vals <- setNames(0.05 + 0.9 * (position %% 1), symbols)
    for (grp in simplex) {
      vals[grp] <- vals[grp] / sum(vals[grp])
    }
    vals
  })
}

# all parameters near 0 and all near 1, for the range checks in mpt() and,
# with the observed covariate values, in check_data(). These two corners catch
# a branch (or a covariate) that leaves (0, 1] near the edge of the parameter
# space when the parameters enter the branch with one orientation; branches
# that mix a parameter with its complement can leave (0, 1] at other vertices,
# which are not evaluated. A simplex group sits at a corner (one member takes the rest of the
# mass), so its members stay positive and sum to 1.
.mpt_boundary_points <- function(symbols, simplex, eps = 0.001) {
  lapply(c(eps, 1 - eps), function(value) {
    vals <- setNames(rep(value, length(symbols)), symbols)
    for (grp in simplex) {
      big <- if (value < 0.5) 1L else length(grp)
      vals[grp] <- eps / length(grp)
      vals[grp[big]] <- 1 - sum(vals[grp[-big]])
    }
    vals
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
# table) makes the check unavailable, never a false result. A test point where
# a derivative is not finite (a cusp such as ((c - k)^2)^0.25 at c = k, here or
# at one covariate setting) is left out, as a point on a singular set is
# outvoted, so a deficit at the other points is still reported. Outvoting needs
# a second point, so with fewer than two points left the check is unavailable;
# a derivative that is not finite at a covariate value in the data for every
# point mostly comes with branches that are undefined there, which check_data()
# reports by row.
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
# A simplex group is free through its stick-breaking components (`sticks`, named
# by member), not its members: the member columns are multiplied by the
# derivative of the members with respect to the sticks, at test points that lie
# on the simplex. A fixed stick keeps its column out of the free set.
# Covariates identify parameters through rows with different covariate values,
# so `settings` (one data frame of covariate values per tree) stacks one block
# of rows per tree and setting; the degrees of freedom grow with the number of
# blocks. The blocks enter `chunk_size` settings at a time through the R factor
# of a QR decomposition, which leaves the rank, the zero columns and the null
# space as they are for the full stack (the stick derivative is a
# right-multiplication, so it applies to R as it would to the stack).
.mpt_jacobian_rank <- function(trees, free, fixed = list(), simplex = list(),
                               sticks = character(0), settings = NULL,
                               tolerance = 1e-8, chunk_size = 2000L) {
  # the counts are reported even when the rank cannot be computed; df is that
  # of one design cell, as print() states it
  counts <- list(
    n_free = length(free), free = free,
    df = sum(lengths(lapply(trees, `[[`, "branches")) - 1L)
  )
  if (length(free) == 0L) {
    return(c(counts, list(
      rank = 0L, involved = character(0), absent = character(0), left_out = 0L
    )))
  }
  trees <- unname(trees)
  parameters <- c(setdiff(free, sticks), unlist(simplex), names(fixed))
  derivs <- try(lapply(trees, function(tree) {
    unlist(lapply(tree$branches, function(branch) {
      lapply(parameters, function(par) {
        if (par %in% all.vars(branch)) stats::D(branch, par) else 0
      })
    }), recursive = FALSE)
  }), silent = TRUE)
  if (is_try_error(derivs)) {
    return(c(counts, list(error = glue(
      "stats::D() cannot differentiate the branch expressions \\
      ({conditionMessage(attr(derivs, 'condition'))})"
    ))))
  }
  settings <- settings %||% lapply(trees, function(tree) data.frame(row.names = 1L))
  # no rank exceeds the rows stacked over every setting
  max_rank <- sum(
    (lengths(lapply(trees, `[[`, "branches")) - 1L) * vapply(settings, nrow, integer(1))
  )
  symbols <- .mpt_tree_parameters(trees)
  points <- .mpt_test_points(symbols, simplex)
  # one vector per symbol holding its value at every test point, so each
  # derivative is evaluated once for all points
  vals <- lapply(setNames(nm = symbols), function(symbol) {
    vapply(points, `[[`, numeric(1), symbol)
  })
  vals[names(fixed)] <- fixed
  # per test point, the R factor of a QR decomposition of the rows stacked so
  # far: R'R = J'J, so column norms and the singular values of the
  # column-scaled matrix are those of the full stack, while memory stays bounded
  # by one chunk of settings. tol = 0 turns off LINPACK's column pivoting, so
  # the columns stay in parameter order
  n_points <- length(points)
  factors <- vector("list", n_points)
  finite <- rep(TRUE, n_points)
  for (tree in seq_along(trees)) {
    n_rows <- nrow(settings[[tree]])
    for (rows in split(seq_len(n_rows), ceiling(seq_len(n_rows) / chunk_size))) {
      # each derivative evaluated once over the chunk's settings and every test
      # point, the points varying fastest
      row_vals <- lapply(vals, rep, times = length(rows))
      row_vals[names(settings[[tree]])] <- lapply(
        settings[[tree]], function(column) rep(column[rows], each = n_points)
      )
      block <- vapply(derivs[[tree]], function(deriv) {
        rep_len(eval(deriv, row_vals), length(rows) * n_points)
      }, numeric(length(rows) * n_points))
      for (point in which(finite)) {
        point_block <- block[seq(point, nrow(block), by = n_points), , drop = FALSE]
        # qr() stops on a non-finite entry with a bare foreign-call error
        if (!all(is.finite(point_block))) {
          finite[point] <- FALSE
          next
        }
        factors[[point]] <- qr.R(qr(rbind(factors[[point]], matrix(
          t(point_block), ncol = length(parameters), byrow = TRUE
        )), tol = 0))
      }
    }
  }
  if (sum(finite) < 2L) {
    return(c(counts, list(error = glue(
      "at {length(finite) - sum(finite)} of the {length(finite)} interior \\
      test values some derivative of the branch expressions is not finite"
    ))))
  }
  decompositions <- lapply(which(finite), function(point) {
    jacobian <- factors[[point]]
    colnames(jacobian) <- parameters
    for (grp in simplex) {
      jacobian <- cbind(
        jacobian[, setdiff(colnames(jacobian), grp), drop = FALSE],
        .mpt_stick_jacobian(
          jacobian[, grp, drop = FALSE], points[[point]][grp],
          sticks[grp[-length(grp)]]
        )
      )
    }
    norms <- sqrt(colSums(jacobian^2))
    zero <- norms[free] <= 1e-12 * max(norms)
    jacobian <- jacobian[, free, drop = FALSE]
    norms <- norms[free]
    jacobian[, zero] <- 0
    norms[zero] <- 1
    decomposition <- svd(
      sweep(jacobian, 2, norms, "/"), nu = 0, nv = length(free)
    )
    decomposition$rank <- min(
      sum(decomposition$d > tolerance * decomposition$d[1]), max_rank
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
    absent = best$absent,
    left_out = sum(!finite), n_points = n_points
  ))
}

# chain rule through the stick-breaking map: member k is
# s_k * prod_{j < k} (1 - s_j) and the last member takes the remainder, so
# d member_k / d s_k is the remainder before k, and every later member falls in
# proportion to its own value, d member_k / d s_m = -member_k / (1 - s_m)
.mpt_stick_jacobian <- function(member_jacobian, members, sticks) {
  n_sticks <- length(sticks)
  remainder <- 1 - c(0, cumsum(members[seq_len(n_sticks)]))
  stick_values <- members[seq_len(n_sticks)] / remainder[seq_len(n_sticks)]
  map <- matrix(0, length(members), n_sticks, dimnames = list(NULL, sticks))
  for (m in seq_len(n_sticks)) {
    map[m, m] <- remainder[m]
    later <- seq.int(m + 1L, length(members))
    map[later, m] <- -members[later] / (1 - stick_values[m])
  }
  member_jacobian %*% map
}

# the first undefined, branch-sum or branch-range violation per tree, NA where
# every test point gives branches in (0, 1] that sum to 1. A branch that is not
# a number (0 / 0) is NaN in Stan too. Branches can sum to 1 for every value and
# still leave (0, 1] (2 * a and 1 - 2 * a); a branch below or at 0 is log(p) of
# a non-positive number in Stan. At an interior point an exact 0 is a branch
# that is 0 for every value, e.g. (1 - a) * 0; at the two boundary corners a
# valid branch may underflow to 0 or, written as 1 minus the others, cancel to
# just below it (-8.5e-20 for seven stages), so only a branch below -tolerance
# counts there. The corners catch a branch that leaves (0, 1] only near the
# edge of the parameter space (1.2 * a - 0.2)
.mpt_tree_branch_errors <- function(trees, parameters, simplex,
                                    tolerance = 1e-6) {
  interior <- .mpt_test_points(parameters, simplex)
  points <- c(interior, .mpt_boundary_points(parameters, simplex))
  vapply(trees, function(tree) {
    for (point in seq_along(points)) {
      vals <- points[[point]]
      probs <- .mpt_eval_branches(tree, as.list(vals))
      at_values <- function(category) {
        symbols <- all.vars(tree$branches[[category]])
        paste(symbols, "=", signif(vals[symbols], 3), collapse = ", ")
      }
      undefined <- names(probs)[is.na(probs)][1]
      if (!is.na(undefined)) {
        return(glue(
          "The branch probability of category '{undefined}' in tree \\
          '{tree$name}' is not a number at the test values \\
          {at_values(undefined)}, so the likelihood is undefined there. \\
          Please check the branch expressions."
        ))
      }
      if (abs(sum(probs) - 1) > tolerance) {
        return(glue(
          "The branch probabilities of tree '{tree$name}' sum to \\
          {signif(sum(probs), 6)} instead of 1 when evaluated at numeric test \\
          values. Please check the branch expressions. To equate parameters or \\
          fix one to a constant, use the restrictions argument (e.g. \\
          restrictions = 'Dn = Do'), or tie them in the formula (e.g. Dn ~ Do). \\
          Parameters that must sum to 1 across branches (e.g. guessing over \\
          three options) belong in the simplex argument."
        ))
      }
      too_low <- if (point <= length(interior)) probs <= 0 else probs < -tolerance
      outside <- names(probs)[too_low | probs > 1 + tolerance][1]
      if (!is.na(outside)) {
        zero_hint <- if (probs[[outside]] == 0) {
          " A zero probability makes the likelihood undefined."
        } else {
          ""
        }
        return(glue(
          "The branch probability of category '{outside}' in tree \\
          '{tree$name}' is {signif(probs[[outside]], 6)} at the test values \\
          {at_values(outside)}, outside (0, 1]. Please check the branch \\
          expressions.{zero_hint}"
        ))
      }
    }
    NA_character_
  }, character(1))
}
