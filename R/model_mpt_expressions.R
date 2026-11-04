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

# mpt_tree() formulas, category ~ expression, become the same named list of
# language objects that the string form produces; the right-hand sides are
# never evaluated and the formula environments are dropped
.mpt_formula_branches <- function(formulas, tree_name) {
  is_formula <- vapply(formulas, inherits, logical(1), "formula")
  stopif(
    !all(is_formula),
    "Cannot mix formulas and other branch expressions in tree '{tree_name}'. \\
    Write every branch as a formula, category ~ expression, or pass a named \\
    list of character strings."
  )
  one_sided <- lengths(formulas) != 3L
  stopif(
    any(one_sided),
    "Branch formulas must be two-sided, category ~ expression (tree \\
    '{tree_name}'). One-sided: {collapse_comma(vapply(formulas[one_sided], deparse1, ''))}"
  )
  lhs_is_symbol <- vapply(formulas, function(f) is.symbol(f[[2]]), logical(1))
  stopif(
    !all(lhs_is_symbol),
    "The left-hand side of a branch formula must be a single response category \\
    name (tree '{tree_name}'). Offending: \\
    {collapse_comma(vapply(formulas[!lhs_is_symbol], deparse1, ''))}"
  )
  setNames(
    lapply(formulas, function(f) f[[3]]),
    vapply(formulas, function(f) as.character(f[[2]]), character(1))
  )
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
    'Order constraints' of the MPT article."
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
# which .mpt_branch_vertices() supplies, in mpt() and in check_data(). A
# simplex group sits at a corner (one member takes the rest of the mass), so
# its members stay positive and sum to 1.
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

# The parameter values at which one branch is checked to stay in (0, 1], one
# row per point and one column per parameter of the branch.
# A branch that is multilinear in its parameters (each enters with degree <= 1,
# also when it appears in several terms) takes its extremes over the parameter
# box at the vertices, so the 2^k vertices of its k parameters make the range
# check exact. A parameter of higher degree can peak between the vertices, as u
# does in the pair-clustering term u * (1 - u), so it takes 0.5 as well; that
# covers the common u * (1 - u) terms, other higher-degree terms are checked
# approximately only. Degree is read from stats::D(): a parameter is of higher
# degree when its derivative still contains it (or another member of its
# simplex group). A derivative that D() cannot form counts as higher degree.
# A simplex group is one dimension whose vertices put the mass on one member,
# as the box vertices would not sum to 1; with a higher-degree term it also
# takes the midpoints between members and the centroid.
# The vertices sit eps inside the box, where a branch such as a / (a + b) is
# still defined; eps is far below the tolerance of the range check, so a branch
# that passes them stays within the tolerance at the true vertices. At 0.001 a
# multilinear branch could exceed 1 by 0.1% per parameter unnoticed.
# A branch whose grid has more than max_points rows (k = 10 multilinear
# parameters) falls back to the two corners and n_fallback further points.
# Those come from a Weyl sequence (steps sqrt(prime)) instead of a random
# number generator, so they are the same on every call and leave the global
# random number state alone; the attribute "exact" says whether the grid is
# complete. The attributes "higher_free" and "higher_groups" name the
# parameters and simplex groups of higher degree, which .mpt_tree_range()
# searches between the vertices; "differentiable" says whether D() could read
# the degree of every parameter.
.mpt_branch_vertices <- function(branch, parameters, simplex, eps = 1e-9,
                                 max_points = 1024L, n_fallback = 64L) {
  used <- intersect(all.vars(branch), parameters)
  if (length(used) == 0L) {
    return(structure(matrix(numeric(0), 1L, 0L), exact = TRUE,
                     higher_free = character(0), higher_groups = list(),
                     differentiable = TRUE))
  }
  differentiable <- !any(vapply(used, function(par) {
    is_try_error(try(stats::D(branch, par), silent = TRUE))
  }, logical(1)))
  groups <- Filter(function(grp) any(grp %in% used), simplex)
  free <- setdiff(used, unlist(groups))
  higher_free <- free[vapply(free, function(par) .mpt_higher_degree(branch, par), logical(1))]
  higher_groups <- Filter(function(grp) .mpt_higher_degree(branch, intersect(grp, used)), groups)
  candidates <- c(
    lapply(free, function(par) {
      values <- if (par %in% higher_free) c(eps, 0.5, 1 - eps) else c(eps, 1 - eps)
      matrix(values, ncol = 1L, dimnames = list(NULL, par))
    }),
    lapply(groups, function(grp) {
      members <- intersect(grp, used)
      higher <- any(vapply(higher_groups, identical, logical(1), grp))
      points <- .mpt_simplex_points(grp, eps, higher)
      unique(points[, members, drop = FALSE])
    })
  )

  n_points <- prod(vapply(candidates, nrow, integer(1)))
  if (n_points <= max_points) {
    index <- expand.grid(lapply(candidates, function(points) seq_len(nrow(points))))
    return(structure(
      do.call(cbind, Map(function(points, i) points[i, , drop = FALSE], candidates, index)),
      exact = TRUE, higher_free = higher_free, higher_groups = higher_groups,
      differentiable = differentiable
    ))
  }

  unit <- outer(seq_len(n_fallback), sqrt(.mpt_first_primes(length(candidates)))) %% 1
  sampled <- do.call(cbind, lapply(seq_along(candidates), function(dim) {
    candidates[[dim]][1 + floor(nrow(candidates[[dim]]) * unit[, dim]), , drop = FALSE]
  }))
  corners <- do.call(rbind, .mpt_boundary_points(parameters, simplex, eps))
  structure(
    rbind(corners[, colnames(sampled), drop = FALSE], sampled),
    exact = FALSE, higher_free = higher_free, higher_groups = higher_groups,
    differentiable = differentiable
  )
}

# the first vertex at which `branch` is not a number or leaves (0, 1], with the
# covariate values given as a named list of equal-length vectors (empty for a
# branch without covariates); NULL when no vertex does. At the vertices a branch
# may underflow to 0 or cancel to just below it, so only a value below
# -tolerance counts.
.mpt_first_range_violation <- function(branch, vertices, covariate_values = list(),
                                       tolerance = 1e-6) {
  n_rows <- if (length(covariate_values) > 0L) length(covariate_values[[1]]) else 1L
  for (vertex in seq_len(nrow(vertices))) {
    values <- setNames(vertices[vertex, ], colnames(vertices))
    branch_value <- rep(
      eval(branch, envir = c(as.list(values), covariate_values)),
      length.out = n_rows
    )
    outside <- is.na(branch_value) | branch_value < -tolerance |
      branch_value > 1 + tolerance
    if (any(outside)) {
      first <- which(outside)[1]
      at <- if (ncol(vertices) > 0L) {
        paste(colnames(vertices), "=", round(values, 4), collapse = ", ")
      }
      return(list(row = first, value = branch_value[first], at = at))
    }
  }
  NULL
}

# For the search, the names of `labels` say why the branch is searched rather
# than certified, and after "|" what its budget left out (.mpt_tree_range()),
# one message per reason; an unnamed label is outside the certified form.
.mpt_message_sampled <- function(labels, reason = c("vertices", "search")) {
  if (length(labels) == 0L) {
    return(invisible(NULL))
  }
  if (match.arg(reason) == "vertices") {
    labels <- paste(labels, collapse = ", ")
    message2(
      "The range check of the branches of {labels} uses selected points of \\
      the parameter box only, not every vertex, because they have too many \\
      vertices to evaluate. A covariate value or parameter combination that \\
      takes one of them outside (0, 1] elsewhere is not detected."
    )
    return(invisible(NULL))
  }
  keys <- names(labels) %||% rep("", length(labels))
  causes <- sub("\\|.*$", "", keys)
  causes[!causes %in% c("grid", "grid_covariates", "identity", "derivative")] <- "form"
  skipped <- ifelse(grepl("|", keys, fixed = TRUE), sub("^[^|]*\\|", "", keys), NA)
  trees <- sub("^.* in tree '(.*)'$", "\\1", labels)
  for (cause in unique(causes)) {
    in_cause <- causes == cause
    one <- length(unique(trees[in_cause])) == 1L
    tree_has <- if (one) "their tree has" else "their trees have"
    # glue() keeps the line continuations of an interpolated value, so the
    # explanation is pasted
    explanation <- switch(cause,
      form = paste(
        "a parameter enters them more than once in a product (or several",
        "members of a simplex group do), and",
        if (one) "their tree is" else "their trees are",
        "not written with parameters, complements (1 - p), non-negative",
        "constants, sums, products and whole-number powers alone, the form for",
        "which the check is exact"
      ),
      grid = paste(
        tree_has, "the form for which the check is exact, but the exact check",
        "would need more than 65,536 points (a grid over the parameters)"
      ),
      grid_covariates = paste(
        tree_has, "the form for which the check is exact, but the exact check",
        "would need more than 65,536 points (a grid over the parameters for each",
        "distinct covariate setting in", if (one) "the tree's" else "each tree's", "rows)"
      ),
      identity = paste(
        tree_has, "the form for which the check is exact, but the sum of",
        if (one) "its branches" else "the branches of each", "could not be",
        "confirmed as 1 within the numerical precision of the check, as happens",
        "with high powers such as (1 - g)^60 of a member of a simplex group"
      ),
      derivative = paste(
        "stats::D() cannot differentiate them, so it is unknown whether their",
        "extremes lie at the vertices"
      )
    )
    budget <- !is.na(skipped) & in_cause
    left_out <- if (any(budget)) {
      paste0(" The search ", paste(vapply(unique(skipped[budget]), function(skip) {
        paste(skip, "for", paste(labels[budget & skipped == skip], collapse = ", "))
      }, character(1)), collapse = "; it "), ".")
    }
    message2(
      "The range check of the branches of {paste(labels[in_cause], collapse = ', ')} \\
      is approximate: {explanation}. A local search between the \\
      vertices found no value outside (0, 1]; a parameter combination that it \\
      misses is not detected.{left_out %||% ''}"
    )
  }
}

# The range check of one tree at the rows of covariate_values (a named list of
# equal-length vectors, empty without covariates; `probabilities` and
# `exponents` say which covariates the certificate may read as such): every
# branch at its vertices (.mpt_branch_vertices()), which settle a branch
# multilinear in its parameters. Only a tree with a branch they do not settle
# (a higher-degree parameter or group, or a sampled grid) goes on: to the
# certificate (.mpt_certified_tree()), which settles the whole tree, and
# failing that to a local search between the vertices of each higher-degree
# branch. The result holds the category and violation of the first branch
# found outside (0, 1], or `approximate`, the labels of the branches whose
# check is not exact by reason: "vertices" when the grid was only sampled,
# "search" when the search stood in. A search label is named by the reason the
# certificate did not settle it (.mpt_certified_tree()), or "derivative" for a
# branch that stats::D() cannot differentiate, followed after "|" by what the
# search left out for its budget, if anything, which .mpt_message_sampled()
# states.
.mpt_tree_range <- function(tree, parameters, simplex, covariate_values = list(),
                            probabilities = parameters, exponents = numeric(0),
                            tolerance = 1e-6) {
  unsettled <- list()
  for (resp_cat in names(tree$branches)) {
    branch <- tree$branches[[resp_cat]]
    vertices <- .mpt_branch_vertices(branch, parameters, simplex)
    violation <- .mpt_first_range_violation(branch, vertices, covariate_values, tolerance)
    if (!is.null(violation)) {
      return(list(category = resp_cat, violation = violation))
    }
    higher <- length(attr(vertices, "higher_free")) + length(attr(vertices, "higher_groups"))
    if (!attr(vertices, "exact") || higher > 0L) {
      unsettled[[resp_cat]] <- vertices
    }
  }
  approximate <- list(vertices = character(0), search = character(0))
  if (length(unsettled) == 0L) {
    return(list(approximate = approximate))
  }
  certified <- .mpt_certified_tree(
    tree, parameters, simplex, probabilities, covariate_values, exponents
  )
  if (certified) {
    return(list(approximate = approximate))
  }
  for (resp_cat in names(unsettled)) {
    vertices <- unsettled[[resp_cat]]
    label <- glue("category '{resp_cat}' in tree '{tree$name}'")
    if (!attr(vertices, "exact")) {
      approximate$vertices <- c(approximate$vertices, label)
    }
    if (length(attr(vertices, "higher_free")) + length(attr(vertices, "higher_groups")) == 0L) {
      next
    }
    searched <- .mpt_search_range_violation(
      tree$branches[[resp_cat]], vertices, attr(vertices, "higher_free"),
      attr(vertices, "higher_groups"), covariate_values, tolerance
    )
    if (!is.null(searched$violation)) {
      return(list(category = resp_cat, violation = searched$violation))
    }
    approximate$search <- c(approximate$search, setNames(
      as.character(label),
      paste(c(
        if (attr(vertices, "differentiable")) attr(certified, "reason") else "derivative",
        searched$skipped
      ), collapse = "|")
    ))
  }
  list(approximate = approximate)
}

# The search between the vertices, over the higher-degree coordinates only.
# The other parameters stay at the start's vertex: for fixed values of the
# higher-degree coordinates the branch is affine in each of them, so its
# extremes lie at their vertices. The starts come per cell
# (.mpt_search_starts()), at most max_starts + max_further towards the
# largest values and as many towards the smallest; `skipped` in the result
# says how many cells got no start in either direction, or that the starts
# after a cell's first were capped, when the budget did not reach them all.
# Two local searches run from each start, as neither finds every case alone
# (measured on 1-D branches, #560 review): L-BFGS-B on the probability scale,
# within the box, whose first step runs to the edge of the box, finds the
# narrow dip of 0.06 - 3.7095 * u * (1 - u)^19 from the start u = 0.5, where
# BFGS on the logit scale stalls; BFGS on the logit scale finds the peak of
# 60 * u * (1 - u)^20 beside the broad one of 8 * u^2 * (1 - u)^2 from the
# grid point u = 1 - 1e-9, which to_z() moves to u = 1 - 1e-3 on the logit
# scale, while L-BFGS-B starts at 1 - 1e-9 and ends on the broad peak. Neither
# finds that peak from the grid points u = 1e-9 or u = 0.5, so a cell needs its
# third start there.
# Without a free coordinate both would move the same softmax coordinates
# without bounds, so only BFGS runs. Both climb log(value) towards the
# largest values and log(1 - value) towards the smallest: on the logit scale
# the tail of a peak such as 60 * u * (1 - u)^20 is convex, where BFGS resets
# its curvature estimate at every step and crawls, while the log of a product
# of parameters and complements is concave there.
# A free parameter stays in [eps, 1 - eps] (L-BFGS-B bounds its coordinate,
# plogis() maps the logit), and a simplex group moves through softmax() over
# all its members, clamped and renormalised, so every point visited lies in
# the box and on the simplex, and any value visited outside
# [-tolerance, 1 + tolerance] is a true witness; the most extreme one is
# reported, and a value that is not finite is recorded before it ends the
# search. Fixed starts and no random numbers make the search deterministic.
.mpt_search_range_violation <- function(branch, vertices, free, groups,
                                        covariate_values, tolerance,
                                        max_starts = 256L, max_further = 64L,
                                        eps = 1e-9) {
  n_rows <- if (length(covariate_values) > 0L) length(covariate_values[[1]]) else 1L
  grid_values <- matrix(vapply(seq_len(nrow(vertices)), function(vertex) {
    values <- setNames(vertices[vertex, ], colnames(vertices))
    rep(eval(branch, c(as.list(values), covariate_values)), length.out = n_rows)
  }, numeric(n_rows)), nrow = n_rows)
  moved <- c(free, unlist(groups))
  other_keys <- do.call(paste, c(
    list(rep("", nrow(vertices))),
    lapply(setdiff(colnames(vertices), moved), function(par) vertices[, par]),
    sep = "|"
  ))
  cell_vertices <- unname(split(seq_len(nrow(vertices)), match(other_keys, unique(other_keys))))
  on_grid <- moved %in% colnames(vertices)
  columns <- match(moved[on_grid], colnames(vertices))
  n_members <- length(moved) - length(free)
  # z holds the free coordinates first, then every member of each group; on
  # the probability scale the free coordinates are the values themselves.
  # pmin.int() and pmax.int() rather than pmin() and pmax(), which cost more
  # than the branch itself on these short vectors
  to_values <- function(z, scale) {
    values <- z[seq_along(free)]
    if (scale == "logit") values <- stats::plogis(values)
    values <- pmin.int(pmax.int(values, eps), 1 - eps)
    offset <- length(free)
    for (grp in groups) {
      share <- exp(z[offset + seq_along(grp)] - max(z[offset + seq_along(grp)]))
      share <- pmin.int(pmax.int(share / sum(share), eps), 1 - eps)
      values <- c(values, share / sum(share))
      offset <- offset + length(grp)
    }
    values[on_grid]
  }
  to_z <- function(start, scale) {
    z <- if (scale == "logit") {
      stats::qlogis(pmin(pmax(start[free], 1e-3), 1 - 1e-3))
    } else {
      pmin(pmax(start[free], eps), 1 - eps)
    }
    for (grp in groups) {
      members <- intersect(grp, names(start))
      shares <- setNames(numeric(length(grp)), grp)
      shares[members] <- start[members]
      rest <- setdiff(grp, members)
      shares[rest] <- (1 - sum(start[members])) / length(rest)
      z <- c(z, log(pmax(shares, 1e-3)))
    }
    z
  }
  local_search <- function(z, objective, scale) {
    if (scale == "logit") {
      return(stats::optim(z, objective, method = "BFGS", control = list(maxit = 100L)))
    }
    stats::optim(
      z, objective, method = "L-BFGS-B",
      lower = c(rep(eps, length(free)), rep(-Inf, n_members)),
      upper = c(rep(1 - eps, length(free)), rep(Inf, n_members)),
      control = list(maxit = 100L)
    )
  }
  directions <- c(1, -1)
  starts <- lapply(directions, function(direction) {
    .mpt_search_starts(direction * grid_values, cell_vertices, max_starts, max_further)
  })
  for (d in seq_along(directions)) {
    direction <- directions[d]
    for (start in seq_along(starts[[d]]$row)) {
      row <- starts[[d]]$row[start]
      vertex <- setNames(vertices[starts[[d]]$vertex[start], ], colnames(vertices))
      env <- c(as.list(vertex), lapply(covariate_values, `[`, row))
      extreme <- NULL
      for (scale in c(if (length(free) > 0L) "probability", "logit")) {
        objective <- function(z) {
          moved_values <- to_values(z, scale)
          env[columns] <- moved_values
          value <- eval(branch, env)
          if (is.null(extreme) || !is.finite(value) ||
                direction * value > direction * extreme$value) {
            values <- vertex
            values[columns] <- moved_values
            extreme <<- list(value = value, values = values)
          }
          if (!is.finite(value)) stop2("undefined")
          -log(max(if (direction > 0) value else 1 - value, .Machine$double.xmin))
        }
        try(local_search(to_z(vertex, scale), objective, scale), silent = TRUE)
        if (!is.finite(extreme$value) || extreme$value > 1 + tolerance ||
              extreme$value < -tolerance) {
          return(list(violation = list(
            row = row, value = extreme$value,
            at = paste(colnames(vertices), "=", round(extreme$values, 4), collapse = ", ")
          )))
        }
      }
    }
  }
  cell_of <- integer(nrow(vertices))
  cell_of[unlist(cell_vertices)] <- rep(seq_along(cell_vertices), lengths(cell_vertices))
  started <- lapply(starts, function(s) s$row + n_rows * (cell_of[s$vertex] - 1L))
  capped <- any(vapply(started, function(keys) {
    sum(lengths(cell_vertices)[(unique(keys) - 1L) %/% n_rows + 1L]) > length(keys)
  }, logical(1)))
  n_cells <- n_rows * length(cell_vertices)
  unsearched <- n_cells - length(unique(unlist(started)))
  if (unsearched == 0L && !capped) {
    return(list(violation = NULL))
  }
  big <- function(n) formatC(n, format = "d", big.mark = ",")
  cells <- if (length(cell_vertices) == 1L) {
    "distinct covariate settings"
  } else if (n_rows == 1L) {
    "vertices of the parameters it holds fixed"
  } else {
    "combinations of a distinct covariate setting with a vertex of the parameters it holds fixed"
  }
  list(violation = NULL, skipped = paste(c(
    if (unsearched > 0L) glue("left unsearched {big(unsearched)} of the {big(n_cells)} {cells}"),
    if (unsearched == 0L && n_cells > 1L) glue("searched all {big(n_cells)} {cells}"),
    if (capped) "capped its further starts"
  ), collapse = " and "))
}

# The starts of the search, as the indices `row` (covariate row) and `vertex`
# of grid points, picked by `scores` (the grid values times the direction of
# the search, one row per covariate row and one column per vertex). A cell is
# one covariate row with one vertex of the parameters that the search leaves
# alone (`cell_vertices` holds the vertices of each), and its grid points
# differ only in the higher-degree coordinates. Starting from the best scores
# overall would leave a cell whose grid values are all poor unsearched,
# however far its value rises between them, so every cell gets its best grid
# point before any cell gets its second best, and its second best before any
# gets its third (a peak can be reachable from one grid point only). The best
# and second best points take up to max_starts starts, the third and later
# ones up to max_further more, as a branch with several higher-degree
# parameters has many grid points per cell (243 for five) and each start
# costs more there. Among the best points, the best of each covariate
# row come first, so a budget below the number of cells still reaches every
# row it can. Within each of these stages the points are ordered by score.
# Each rank takes one max.col() per vertex of the other parameters, with no
# sort over all grid points.
.mpt_search_starts <- function(scores, cell_vertices, max_starts, max_further) {
  rows <- seq_len(nrow(scores))
  starts <- list(row = integer(0), vertex = integer(0))
  limit <- max_starts
  for (rank in seq_len(max(lengths(cell_vertices)))) {
    if (rank == 3L) {
      limit <- length(starts$row) + max_further
    }
    if (length(starts$row) >= limit) {
      break
    }
    cells <- cell_vertices[lengths(cell_vertices) >= rank]
    best <- score <- vector("list", length(cells))
    for (cell in seq_along(cells)) {
      best[[cell]] <- cells[[cell]][
        max.col(scores[, cells[[cell]], drop = FALSE], ties.method = "first")
      ]
      score[[cell]] <- scores[cbind(rows, best[[cell]])]
      scores[cbind(rows, best[[cell]])] <- -Inf
    }
    vertex <- unlist(best)
    row <- rep(rows, length(cells))
    by_score <- order(-unlist(score))
    if (rank == 1L) {
      by_score <- by_score[order(duplicated(row[by_score]))]
    }
    by_score <- utils::head(by_score, limit - length(starts$row))
    starts$row <- c(starts$row, row[by_score])
    starts$vertex <- c(starts$vertex, vertex[by_score])
  }
  starts
}

# A tree is certified when every branch is non-negative by its form alone
# (.mpt_form()) and the branches sum to 1 as polynomials: each branch is then
# also at most 1, as the others are non-negative, on the whole parameter box
# and simplex, so no vertex or search is needed. Covariates count as
# probabilities only where named in `probabilities`, which check_data() does
# when their values in the tree's rows lie in [0, 1], and as exponents (named
# by their largest value in `exponents`) when those values are whole numbers.
# An uncertified tree carries the reason in the attribute "reason": "form"
# (a branch outside the form), "grid" or, with covariate rows,
# "grid_covariates" (the sum is too costly to check) or "identity" (the sum is
# not 1 at every grid point within the tolerance).
.mpt_certified_tree <- function(tree, parameters, simplex, probabilities = parameters,
                                covariate_values = list(), exponents = numeric(0)) {
  used <- intersect(.mpt_expr_vars(tree), parameters)
  forms <- lapply(tree$branches, .mpt_form, probabilities, used, exponents)
  if (any(vapply(forms, `[[`, integer(1), "bound") < 1L)) {
    return(structure(FALSE, reason = "form"))
  }
  degree <- setNames(do.call(pmax, lapply(forms, `[[`, "degree")), used)
  sums_to_one <- .mpt_branches_sum_to_one(tree, degree, simplex, covariate_values)
  if (is.na(sums_to_one)) {
    return(structure(
      FALSE, reason = if (length(covariate_values) > 0L) "grid_covariates" else "grid"
    ))
  }
  if (!sums_to_one) {
    return(structure(FALSE, reason = "identity"))
  }
  TRUE
}

# The bound an expression keeps by its form alone, for all values in [0, 1] of
# `probabilities`: 2 for [0, 1] (a probability, a constant in [0, 1], a
# product or whole-number power of such terms, 1 minus one of them), 1 for
# non-negative (a non-negative constant, a sum, product or whole-number power
# of non-negative terms), 0 otherwise. A part without symbols is evaluated, so
# 1/5 is a constant, while a division or function of a parameter is outside
# the form. A covariate named in `exponents` is a whole-number exponent of at
# most that value, as in (1 - D)^n. Where the bound is 1 or 2 the expression
# is a polynomial, and `degree` bounds the degree of each of `symbols` in it.
# One pass over the expression.
.mpt_form <- function(expr, probabilities, symbols, exponents = numeric(0)) {
  zero <- numeric(length(symbols))
  constant <- function(value) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      return(list(bound = 0L, degree = zero))
    }
    list(bound = if (value < 0) 0L else if (value <= 1) 2L else 1L, degree = zero,
         value = value)
  }
  if (is.numeric(expr)) {
    return(constant(expr))
  }
  if (is.symbol(expr)) {
    name <- as.character(expr)
    return(list(bound = if (name %in% probabilities) 2L else 0L,
                degree = as.numeric(symbols == name)))
  }
  if (!is.call(expr)) {
    return(list(bound = 0L, degree = zero))
  }
  args <- lapply(as.list(expr)[-1L], .mpt_form, probabilities, symbols, exponents)
  if (all(vapply(args, function(arg) !is.null(arg$value), logical(1)))) {
    value <- try(eval(expr, envir = baseenv()), silent = TRUE)
    return(if (is_try_error(value)) list(bound = 0L, degree = zero) else constant(value))
  }
  bounds <- vapply(args, `[[`, integer(1), "bound")
  degrees <- lapply(args, `[[`, "degree")
  op <- as.character(expr[[1L]])
  exponent <- if (op == "^") args[[2L]]$value
  if (op == "^" && is.symbol(expr[[3L]]) && as.character(expr[[3L]]) %in% names(exponents)) {
    exponent <- exponents[[as.character(expr[[3L]])]]
  }
  whole <- !is.null(exponent) && exponent >= 0 && exponent == round(exponent)
  bound <- switch(op,
    "(" = bounds[1L],
    "*" = if (length(args) == 2L) min(bounds),
    "+" = if (length(args) == 2L) min(bounds, 1L),
    "-" = if (length(args) == 2L && isTRUE(args[[1L]]$value == 1) && bounds[2L] == 2L) 2L,
    "^" = if (whole) bounds[1L]
  )
  degree <- switch(op,
    "*" = Reduce(`+`, degrees),
    "^" = if (whole) degrees[[1L]] * exponent else zero,
    do.call(pmax, degrees)
  )
  list(bound = bound %||% 0L, degree = degree)
}

# The branches of a certified form are polynomials, and their sum is checked to
# be 1 identically, not only at the test points of .mpt_tree_branch_errors(): a
# polynomial of degree at most d_j in variable j that vanishes on a grid of
# d_j + 1 distinct values per variable is 0. `degree` holds d_j per parameter.
# A simplex group loses its last member to 1 minus the others, which adds that
# member's degree to theirs. Every other member of the group goes on the grid,
# with degree 0 when the tree does not use it: the sum must be 1 on the whole
# simplex, not on the face where an unused member is 0 (a + c with group
# {a, b, c} is 1 - b), and with the last member alone in the tree the others
# are what it is computed from. The identity is algebraic, so the grid need
# not lie on the simplex, and the last member may leave [0, 1] there.
# The grid values are Chebyshev-Lobatto nodes (1 - cos(pi j / d)) / 2,
# j = 0, ..., d, on [0, 1]. A sum that differs from 1 by at most r at every
# grid point then differs from 1 by at most Lambda * r on the whole box, where
# Lambda is the product of the one-dimensional Lebesgue constants: 1 at degree
# 1, 1.25 at 2, 1.67 at 3, 2.72 at 16, 4.59 at 300 and 8.02 at 65535, and at
# most 68.53 for any grid of up to max_points points (computed by a search over
# all degree combinations; equally spaced interior values give 1.3e5 at degree
# 16 and 2.2e12 at degree 40). With tolerance 1e-9 every branch of a
# certified tree is therefore at most 1 + 7e-8 on the parameter box and
# simplex, inside the tolerance of the range check. Off the simplex a valid
# branch can be large ((1 - gC)^k = (gA + gB)^k reaches 2^k in a group of
# three), and its rounding error then exceeds the tolerance (from k = 54 there;
# k = 53 passes): such a tree is not certified and goes to the search, which
# says so. With covariates, the sum
# must be 1 at each row of covariate_values. A grid of more than max_points
# points (16 parameters of degree 1) leaves the sum undecided (NA).
.mpt_branches_sum_to_one <- function(tree, degree, simplex, covariate_values = list(),
                                     max_points = 2^16, tolerance = 1e-9) {
  last <- character(0)
  for (grp in Filter(function(grp) any(grp %in% names(degree)), simplex)) {
    final <- grp[length(grp)]
    members <- grp[-length(grp)]
    degree[setdiff(members, names(degree))] <- 0
    if (final %in% names(degree)) {
      degree[members] <- degree[members] + degree[[final]]
      last <- c(last, final)
    }
  }
  variables <- setdiff(names(degree), last)
  n_rows <- if (length(covariate_values) > 0L) length(covariate_values[[1]]) else 1L
  n_points <- prod(degree[variables] + 1)
  if (n_points * n_rows > max_points) {
    return(NA)
  }
  grid <- expand.grid(lapply(setNames(nm = variables), function(par) {
    if (degree[[par]] == 0) 0.5 else (1 - cos(pi * seq.int(0, degree[[par]]) / degree[[par]])) / 2
  }))
  for (grp in simplex) {
    if (grp[length(grp)] %in% last) {
      grid[[grp[length(grp)]]] <- 1 - rowSums(grid[grp[-length(grp)]])
    }
  }
  env <- c(
    lapply(grid, rep, times = n_rows),
    lapply(covariate_values, rep, each = n_points)
  )
  total <- Reduce(`+`, lapply(tree$branches, function(branch) {
    rep_len(eval(branch, envir = env), n_points * n_rows)
  }))
  !anyNA(total) && all(abs(total - 1) <= tolerance)
}

.mpt_higher_degree <- function(branch, members) {
  any(vapply(members, function(par) {
    derivative <- try(stats::D(branch, par), silent = TRUE)
    is_try_error(derivative) || any(members %in% all.vars(derivative))
  }, logical(1)))
}

# one row per point, one column per member of the group; the vertices put eps
# spread over all members but one, which takes the rest
.mpt_simplex_points <- function(group, eps, with_interior = FALSE) {
  n <- length(group)
  points <- matrix(eps / n, n, n, dimnames = list(NULL, group))
  diag(points) <- 1 - eps * (n - 1) / n
  if (!with_interior) {
    return(points)
  }
  pairs <- which(upper.tri(diag(n)), arr.ind = TRUE)
  midpoints <- (points[pairs[, 1], , drop = FALSE] + points[pairs[, 2], , drop = FALSE]) / 2
  rbind(points, midpoints, colMeans(points))
}

.mpt_first_primes <- function(n) {
  candidates <- seq_len(20 * n + 20) + 1L
  is_prime <- vapply(candidates, function(x) {
    all(x %% seq_len(floor(sqrt(x)))[-1L] != 0L)
  }, logical(1))
  candidates[is_prime][seq_len(n)]
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

# the first branch that is not a number or leaves (0, 1] at a vertex of its own
# parameter box or at a point the local search finds, per tree (NA where none
# does); the attributes "sampled" and "searched" name the branches whose range
# check is approximate, by reason (.mpt_tree_range()). The points of
# .mpt_tree_branch_errors() visit two corners and five interior points, which
# miss a branch that mixes a parameter with another one's complement
# (2 * a * (1 - b) reaches 2 as a -> 1 and b -> 0).
.mpt_tree_vertex_errors <- function(trees, parameters, simplex, tolerance = 1e-6) {
  results <- lapply(trees, function(tree) {
    range <- .mpt_tree_range(tree, parameters, simplex, tolerance = tolerance)
    violation <- range$violation
    if (is.null(violation)) {
      return(list(error = NA_character_, approximate = range$approximate))
    }
    resp_cat <- range$category
    error <- if (is.na(violation$value)) {
      glue(
        "The branch probability of category '{resp_cat}' in tree \\
        '{tree$name}' is not a number at the test values {violation$at}, so \\
        the likelihood is undefined there. Please check the branch expressions."
      )
    } else {
      glue(
        "The branch probability of category '{resp_cat}' in tree \\
        '{tree$name}' is {signif(violation$value, 6)} at the test values \\
        {violation$at}, outside (0, 1]. Please check the branch expressions."
      )
    }
    list(error = error, approximate = list())
  })
  # unname() keeps the reason names of the search labels free of tree names
  approximate <- function(reason) {
    unlist(lapply(unname(results), function(result) result$approximate[[reason]])) %||%
      character(0)
  }
  structure(
    vapply(results, `[[`, character(1), "error"),
    sampled = approximate("vertices"), searched = approximate("search")
  )
}

# the first undefined, branch-sum or branch-range violation per tree, NA where
# every test point gives branches in (0, 1] that sum to 1. A branch that is not
# a number (0 / 0) is NaN in Stan too. Branches can sum to 1 for every value and
# still leave (0, 1] (2 * a and 1 - 2 * a); a branch below or at 0 is log(p) of
# a non-positive number in Stan. At an interior point an exact 0 is a branch
# that is 0 for every value, e.g. (1 - a) * 0, or one that underflows there,
# e.g. (1 - D)^392 at D = 0.851; refusing the latter too is accepted, as no
# realistic tree has such a power. At the two boundary corners a valid branch
# may underflow to 0 or, written as 1 minus the others, cancel to just below it
# (-8.5e-20 for seven stages), so only a branch below -tolerance counts there.
# The corners catch a branch that leaves (0, 1] only near the edge of the
# parameter space (1.2 * a - 0.2); .mpt_tree_vertex_errors() takes the rest
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
        # interpolated into the glue() below, which leaves the line
        # continuations of an interpolated value in place
        zero_hint <- if (probs[[outside]] == 0) {
          paste(
            " A zero probability makes the likelihood undefined. Declare a",
            "response category that a tree cannot produce with",
            "mpt_tree(impossible = ) instead of giving it a zero branch."
          )
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
