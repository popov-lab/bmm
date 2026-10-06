############################################################################# !
# PLOT METHODS                                                           ####
############################################################################# !

#' @title Plot the structure of MPT trees
#'
#' @description Draws the processing-tree diagram implied by the branch
#'   probability expressions: each branch expression is expanded into its
#'   root-to-leaf paths (a sum of products, with the factors in writing
#'   order), paths with a common prefix share the corresponding tree nodes,
#'   and the leaves show the response categories. Calling `plot()` on a model
#'   created with [mpt()] draws all trees of the model.
#'
#' @param x An `mpt_tree` object or a `bmmodel` object created with [mpt()].
#' @param cex Character expansion factor for the edge and leaf labels.
#' @param ... Ignored.
#'
#' @details Edges are drawn in the order the branch expressions first mention
#'   them, bare constants included. A complement edge `1 - X` is placed directly
#'   below its sibling `X`. A numeric constant factor above 1 is not a separate
#'   edge: the constant and all factors after it form one edge, as in
#'   Batchelder & Riefer (1999, Fig. 1). For example, `B * 2 * u * (1 - u)` draws
#'   one edge `B`, then one edge `2 * u * (1 - u)`. If the constant is at the start
#'   (`2 * u * (1 - u)`), the entire expression is one edge. If it is in the middle
#'   (`u * 2 * (1 - u)`), the constant and the factors after it join into one edge,
#'   and a trailing constant joins the factor before it. The plot warns when this
#'   merging hides a process step: a factor after the constant that has no
#'   complement among the merged factors (`D * 2 * g * (1 - g) * r`), a factor
#'   separated from its complement by the constant (`u * 2 * (1 - u)`), or a
#'   quotient that contains a complement (`(1 - Pi) * 1 / 15` is one edge; write
#'   `(1 - Pi) * (1/15)`).
#'   Expressions are expanded distributively, so `Pm * (Pb + (1 - Pb) * 0.25)`
#'   displays two paths sharing the `Pm` edge.
#'
#'   The plot warns when edges leaving an internal node do not sum to 1, but only
#'   for nodes whose edge labels are all constants, plain parameters, or
#'   `1 - parameter`; product edges, quotient edges, function calls, and covariate
#'   fans that lack complements are not checked, nor are nodes with a single edge
#'   or, when plotting a bare tree, fans of plain parameters (which a model may
#'   declare as a simplex). Write factors in the same order across all categories
#'   to avoid mispairing due to different factor orders.
#'
#'   For a model created with [mpt()], the plot draws the trees as written and
#'   marks restricted parameters: the edge of a fixed parameter reads `g = 0.5`,
#'   the edge of an equated parameter `Dn (= Do)`, and an edge that contains a
#'   restricted parameter shows the restriction in parentheses, as in
#'   `1 - g (g = 0.5)`. Impossible categories are omitted. A literal 0 in an
#'   expression (as in `D + 0 * x`) is drawn as an edge labelled 0. The `...`
#'   argument is ignored; the title cannot be changed.
#'   Leaf labels wider than about 90% of a panel are clipped, and long edge
#'   labels near the root can be clipped earlier, when the leaf labels leave
#'   little room.
#'
#' @return The object `x`, invisibly.
#'
#' @references Batchelder, W. H., & Riefer, D. M. (1999). Theoretical and
#'   empirical review of multinomial process tree modeling. Psychonomic
#'   Bulletin & Review, 6(1), 57-86. https://doi.org/10.3758/BF03210812
#'
#' @keywords transform
#'
#' @examples
#' tree_old <- mpt_tree("old", list(
#'   old = "D + (1 - D) * g",
#'   new = "(1 - D) * (1 - g)"
#' ))
#' plot(tree_old)
#'
#' model <- mpt(
#'   list(
#'     tree_old,
#'     mpt_tree("new", list(
#'       old = "(1 - D) * g",
#'       new = "D + (1 - D) * (1 - g)"
#'     ))
#'   ),
#'   tree_id = "item_type"
#' )
#' plot(model)
#' @export
plot.mpt_tree <- function(x, cex = 0.9, ...) {
  .mpt_plot_tree(x, cex, simplex = NULL)
  invisible(x)
}

#' @rdname plot.mpt_tree
#' @export
plot.mpt <- function(x, cex = 0.9, ...) {
  n_trees <- length(x$other_vars$trees)
  if (n_trees > 1 && prod(graphics::par("mfrow")) == 1L) {
    n_row <- floor(sqrt(n_trees))
    old_par <- graphics::par(mfrow = c(n_row, ceiling(n_trees / n_row)))
    on.exit(graphics::par(old_par))
  }
  for (tree in x$other_vars$unrestricted_trees %||% x$other_vars$trees) {
    .mpt_plot_tree(
      tree, cex,
      simplex = x$other_vars$simplex, covariates = x$other_vars$covariates,
      restrictions = x$other_vars$restrictions
    )
  }
  invisible(x)
}

# simplex = NULL means the caller does not know the model's simplex groups
.mpt_plot_tree <- function(x, cex, simplex, covariates = character(0),
                           restrictions = NULL) {
  graph <- .mpt_tree_graph(x)
  graph$edges$label <- .mpt_restricted_labels(graph$edges, restrictions)
  .mpt_warn_sibling_sums(graph, x$name, simplex, covariates, restrictions)
  .mpt_warn_merged_factors(x)
  old_par <- graphics::par(mar = c(0.5, 0.5, 2, 0.5))
  on.exit(graphics::par(old_par))

  is_leaf <- !is.na(graph$nodes$category)
  ylim <- c(min(graph$nodes$y) - 0.5, max(graph$nodes$y) + 0.5)
  graphics::plot.new()
  # the leaf labels hang right of the deepest nodes: the window grows by their
  # width in inches, which only the device can tell
  right_inches <- max(graphics::strwidth(
    graph$nodes$category[is_leaf], units = "inches", cex = cex, font = 2
  )) + 0.3 * graphics::strwidth("m", units = "inches", cex = cex)
  # edge labels near the root are pushed left by half their width. The room
  # they need depends on the window they are placed in, so it is found by a
  # secant step: the overhang falls by less than the room added, because the
  # tree is compressed and its edges steepen
  left_inches <- 0
  previous <- NULL
  for (pass in 1:3) {
    .mpt_set_window(max(graph$nodes$x), ylim, left_inches, right_inches)
    label_pos <- .mpt_edge_label_positions(graph$edges, cex)
    if (pass == 3L) break
    overhang <- .mpt_edge_label_overhang(label_pos, graph$edges$label, cex)
    gain <- if (is.null(previous)) {
      1
    } else {
      (previous$overhang - overhang) / (left_inches - previous$room)
    }
    # no gain: the tree already sits at its minimum width
    if (overhang <= 0 || gain < 0.05) break
    previous <- list(room = left_inches, overhang = overhang)
    left_inches <- left_inches + overhang / gain
  }
  graphics::segments(
    graph$edges$x0, graph$edges$y0, graph$edges$x1, graph$edges$y1,
    col = "grey40"
  )
  graphics::text(label_pos$x, label_pos$y, labels = graph$edges$label, cex = cex)
  graphics::points(
    graph$nodes$x[!is_leaf], graph$nodes$y[!is_leaf],
    pch = 16, cex = 0.6, col = "grey40"
  )
  graphics::text(
    graph$nodes$x[is_leaf], graph$nodes$y[is_leaf],
    labels = graph$nodes$category[is_leaf],
    pos = 4, offset = 0.3, cex = cex, font = 2
  )
  graphics::title(main = glue("MPT tree '{x$name}'"), cex.main = 1)
}

# the tree keeps a tenth of the plot width at least, so absurdly long labels
# cannot invert the window; the leaf labels get their room first and the edge
# labels what remains. xaxs = "i" because the default 4% padding would eat into
# the room reserved for the labels
.mpt_set_window <- function(x_max, ylim, left_inches, right_inches) {
  plot_inches <- graphics::par("pin")[1]
  left_inches <- min(left_inches, max(0.9 * plot_inches - right_inches, 0))
  tree_inches <- max(plot_inches - left_inches - right_inches, 0.1 * plot_inches)
  x_per_inch <- (x_max + 0.2) / tree_inches
  graphics::plot.window(
    xlim = c(-0.1 - left_inches * x_per_inch, x_max + 0.1 + right_inches * x_per_inch),
    ylim = ylim, xaxs = "i"
  )
}

# how far, in inches, the farthest edge label reaches past the left edge of the
# plot region
.mpt_edge_label_overhang <- function(position, labels, cex) {
  usr <- graphics::par("usr")
  x_per_inch <- (usr[2] - usr[1]) / graphics::par("pin")[1]
  left <- position$x -
    graphics::strwidth(labels, units = "inches", cex = cex) / 2 * x_per_inch
  max(0, (usr[1] - min(left)) / x_per_inch)
}

# centers each edge label beside its edge, at a fixed gap from the line: the
# label box is pushed along the edge's normal by the gap plus the box's own
# extent in that direction, so it clears flat and steep edges alike. Distances
# are measured in inches because user units differ between the axes. The label
# sits 60% along the edge rather than at the midpoint, where the edges of a fan
# are still close together.
.mpt_edge_label_positions <- function(edges, cex) {
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  x_per_inch <- (usr[2] - usr[1]) / pin[1]
  y_per_inch <- (usr[4] - usr[3]) / pin[2]
  dx <- (edges$x1 - edges$x0) / x_per_inch
  dy <- (edges$y1 - edges$y0) / y_per_inch
  length_inches <- sqrt(dx^2 + dy^2)
  # edges rising from their parent carry the label above, falling edges below,
  # so the labels of a fan point away from each other instead of into the wedge
  side <- ifelse(dy < 0, -1, 1)
  normal_x <- -side * dy / length_inches
  normal_y <- side * dx / length_inches
  width <- graphics::strwidth(edges$label, units = "inches", cex = cex)
  height <- graphics::strheight(edges$label, units = "inches", cex = cex)
  distance <- 0.4 * height + abs(normal_x) * width / 2 + abs(normal_y) * height / 2
  list(
    x = edges$x0 + 0.6 * (edges$x1 - edges$x0) + normal_x * distance * x_per_inch,
    y = edges$y0 + 0.6 * (edges$y1 - edges$y0) + normal_y * distance * y_per_inch
  )
}

# the edges leaving a node are the branches of one process step, so their
# probabilities sum to 1 for every parameter value. The drawing follows the
# writing order of the factors, which can pair edges that do not belong
# together when categories write their factors in different orders. Edges are
# evaluated from their original expressions (the labels are rounded). Only
# nodes whose labels are all simple probabilities are checked, see
# .mpt_node_is_checkable()
.mpt_warn_sibling_sums <- function(graph, tree_name, simplex, covariates = character(0),
                                   restrictions = NULL, tolerance = 1e-6) {
  edges <- graph$edges
  symbols <- unique(c(
    unlist(lapply(edges$expr, function(e) all.vars(str2lang(e)))),
    unlist(simplex), unlist(lapply(restrictions, all.vars))
  ))
  # the edges come from the trees before the restrictions, which mpt() checked
  # only with the restrictions in place (a, 1 - b sums to 1 once b = a); chains
  # (b = c, c = 0.4) are resolved first so the order of binding does not matter
  restrictions <- .mpt_resolve_restrictions(restrictions)
  points <- lapply(.mpt_test_points(symbols, simplex %||% list()), function(vals) {
    for (par in intersect(names(restrictions), symbols)) {
      vals[[par]] <- eval(restrictions[[par]], as.list(vals))
    }
    vals
  })
  for (node in unique(edges$from)) {
    out <- edges[edges$from == node, ]
    if (!.mpt_node_is_checkable(out$expr, !is.null(simplex), covariates)) next
    sums <- vapply(points, function(vals) {
      sum(vapply(out$expr, function(e) eval(str2lang(e), as.list(vals)), numeric(1)))
    }, numeric(1))
    if (all(abs(sums - 1) <= tolerance)) next
    where <- .mpt_node_description(edges, node)
    siblings <- collapse_comma(out$label)
    warning2(
      "In the tree '{tree_name}', the edges leaving {where} do not sum to 1: \\
      {siblings}. The diagram follows the order in which the factors are \\
      written, so categories that write their factors in different orders \\
      end up under different nodes. Write the factors of all categories in the \\
      same order to get a tree in which every node's edges sum to 1."
    )
  }
  invisible(NULL)
}

# A fan of edges is a probability distribution only when its labels are simple
# probabilities (a constant, a parameter or covariate, or 1 minus one of them).
# Products, divisions, calls and multiplicity edges are not: they mix several
# process steps, so a deviating sum does not show a factor-order mismatch.
.mpt_node_is_checkable <- function(exprs, simplex_known, covariates) {
  # a lone edge continues a path in series (rp * rp) instead of branching
  if (length(exprs) < 2L) {
    return(FALSE)
  }
  parts <- lapply(exprs, .mpt_simple_label)
  if (any(vapply(parts, is.null, logical(1)))) {
    return(FALSE)
  }
  kinds <- vapply(parts, `[[`, character(1), "kind")
  symbols <- vapply(parts, `[[`, character(1), "symbol")
  # members of a simplex sum to 1 as a group, which a bare tree cannot know
  if (!simplex_known && all(kinds == "symbol")) {
    return(FALSE)
  }
  # covariates may sum to 1 in the data only, unless a covariate comes with its complement
  unpaired <- vapply(seq_along(kinds), function(i) {
    kinds[i] == "symbol" && symbols[i] %in% covariates &&
      !any(symbols[-i] == symbols[i] & kinds[-i] == "complement")
  }, logical(1))
  !any(unpaired)
}

# classifies an edge expression as a simple probability: a symbol, a numeric
# constant, or `1 - symbol`, with redundant parentheses ignored; NULL otherwise
.mpt_simple_label <- function(expr) {
  parsed <- .mpt_strip_parens(str2lang(expr))
  if (is.symbol(parsed)) {
    return(list(kind = "symbol", symbol = as.character(parsed)))
  }
  if (length(all.vars(parsed)) == 0L) {
    value <- tryCatch(eval(parsed), error = function(e) NULL)
    return(if (is.numeric(value) && length(value) == 1L) {
      list(kind = "constant", symbol = NA_character_)
    })
  }
  base <- .mpt_complement_base(parsed)
  if (is.symbol(base)) list(kind = "complement", symbol = as.character(base))
}

.mpt_strip_parens <- function(parsed) {
  while (is.call(parsed) && identical(parsed[[1]], quote(`(`))) {
    parsed <- parsed[[2]]
  }
  parsed
}

# X for an expression of the form `1 - X`, NULL for anything else
.mpt_complement_base <- function(parsed) {
  parsed <- .mpt_strip_parens(parsed)
  is_complement <- is.call(parsed) && identical(parsed[[1]], quote(`-`)) &&
    length(parsed) == 3L && identical(.mpt_strip_parens(parsed[[2]]), 1)
  if (is_complement) .mpt_strip_parens(parsed[[3]])
}

# the diagram keeps a restricted parameter's name, so the restriction is
# written beside it: "g = 0.5" and "Dn (= Do)" for the parameter itself, and in
# parentheses after an edge that only contains it ("1 - g (g = 0.5)")
.mpt_restricted_labels <- function(edges, restrictions) {
  if (length(restrictions) == 0L) {
    return(edges$label)
  }
  values <- vapply(restrictions, function(value) {
    if (is.numeric(value)) as.character(signif(value, 3)) else deparse1(value)
  }, character(1))
  vapply(seq_len(nrow(edges)), function(i) {
    parsed <- .mpt_strip_parens(str2lang(edges$expr[i]))
    restricted <- intersect(all.vars(parsed), names(restrictions))
    if (length(restricted) == 0L) {
      edges$label[i]
    } else if (is.symbol(parsed) && is.numeric(restrictions[[restricted]])) {
      glue("{restricted} = {values[[restricted]]}")
    } else if (is.symbol(parsed)) {
      glue("{restricted} (= {values[[restricted]]})")
    } else {
      glue("{edges$label[i]} ({paste(restricted, '=', values[restricted], collapse = ', ')})")
    }
  }, character(1))
}

# names a node by the labels of the edges from the root to it
.mpt_node_description <- function(edges, node) {
  path <- character(0)
  while (node %in% edges$to) {
    edge <- edges[edges$to == node, ]
    path <- c(edge$label, path)
    node <- edge$from
  }
  if (length(path) == 0L) "the root" else glue("the node '{paste(path, collapse = ' > ')}'")
}

# expands a branch probability expression into its root-to-leaf paths: sums
# split into separate paths, products append their factors in writing order,
# and any other subexpression (parameters, complements, constants, covariates)
# is one atomic edge label
.mpt_branch_paths <- function(expr) {
  lapply(.mpt_raw_paths(expr), .mpt_collapse_multiplicity)
}

.mpt_raw_paths <- function(expr) {
  expand <- function(node) {
    if (is.call(node)) {
      op <- if (is.symbol(node[[1]])) as.character(node[[1]]) else ""
      if (op == "+" && length(node) == 3L) {
        return(c(expand(node[[2]]), expand(node[[3]])))
      }
      if (op == "*" && length(node) == 3L) {
        left <- expand(node[[2]])
        right <- expand(node[[3]])
        return(unlist(
          lapply(left, function(l) lapply(right, function(r) c(l, r))),
          recursive = FALSE
        ))
      }
      if (op == "(") {
        return(expand(node[[2]]))
      }
    }
    list(deparse1(node))
  }
  expand(expr)
}

# a constant factor above 1 (the 2 in 2 * u * (1 - u)) counts outcomes rather
# than branching, so it is not an edge of its own, as in the pair-clustering
# tree of Batchelder & Riefer (1999). Written first, it makes the whole path
# one edge; written last, it joins the factor before it; otherwise it and
# every factor after it form one edge
.mpt_collapse_multiplicity <- function(path) {
  first <- .mpt_multiplicity_position(path)
  if (is.na(first)) {
    return(path)
  }
  join_from <- if (first == length(path)) first - 1L else first
  edge_factors <- vapply(path[join_from:length(path)], .mpt_parenthesise_sum, character(1))
  c(path[seq_len(join_from - 1L)], paste(edge_factors, collapse = " * "))
}

# position of the first multiplicity constant, NA if the path has none (a path
# of one factor is an edge of its own, whatever it is)
.mpt_multiplicity_position <- function(path) {
  if (length(path) == 1L) {
    return(NA_integer_)
  }
  which(vapply(path, .mpt_is_multiplicity, logical(1)))[1]
}

.mpt_is_multiplicity <- function(factor_label) {
  parsed <- str2lang(factor_label)
  length(all.vars(parsed)) == 0L &&
    isTRUE(tryCatch(eval(parsed) > 1, error = function(e) FALSE))
}

# only sums and differences bind more loosely than the product they join
.mpt_parenthesise_sum <- function(factor_label) {
  parsed <- str2lang(factor_label)
  is_sum <- is.call(parsed) && is.symbol(parsed[[1]]) &&
    as.character(parsed[[1]]) %in% c("+", "-")
  if (is_sum) paste0("(", factor_label, ")") else factor_label
}

# the drawing merges the factors after a multiplicity constant, and a quotient
# is one edge, so some expressions are drawn without a node the writer meant:
# a factor swallowed by the constant, a factor split from its complement by it,
# or a complement inside a quotient. One warning per tree lists the paths
.mpt_warn_merged_factors <- function(tree) {
  notes <- unlist(lapply(names(tree$branches), function(resp_cat) {
    paths <- .mpt_raw_paths(tree$branches[[resp_cat]])
    lapply(paths, function(path) {
      note <- .mpt_merged_factors_note(path)
      if (!is.null(note)) {
        glue("branch '{resp_cat}': {note}")
      }
    })
  }))
  if (length(notes) == 0L) {
    return(invisible(NULL))
  }
  details <- paste(unique(notes), collapse = "; ")
  warning2(
    "In the tree '{tree$name}', the drawing merges factors into one edge where \\
    the expression may mean separate process steps: {details}. A constant above \\
    1 takes the factors after it into its edge, and a quotient is a single \\
    edge. Write the constant directly before the factors it counts, as in \\
    (1 - c) * 2 * u * (1 - u), and fractions of a path as (1/15)."
  )
  invisible(NULL)
}

# NULL when the path draws as written
.mpt_merged_factors_note <- function(path) {
  factors <- lapply(path, function(factor_label) .mpt_strip_parens(str2lang(factor_label)))
  quotient <- vapply(factors, .mpt_is_quotient_with_complement, logical(1))
  path_text <- paste(vapply(path, .mpt_parenthesise_sum, character(1)), collapse = " * ")
  if (any(quotient)) {
    return(glue("{path_text} is one edge, so the complement it contains is not a node"))
  }
  constant <- .mpt_multiplicity_position(path)
  # a leading constant makes the whole path one edge, as documented
  if (is.na(constant) || constant == 1L) {
    return(NULL)
  }
  own <- vapply(factors, deparse1, character(1))
  complement_of <- vapply(factors, function(f) {
    base <- .mpt_complement_base(f)
    if (is.null(base)) NA_character_ else deparse1(base)
  }, character(1))
  has_partner <- function(i, pool) {
    any(complement_of[pool] == own[i] | own[pool] == complement_of[i], na.rm = TRUE)
  }
  merged <- constant:length(path)
  merged_steps <- merged[lengths(lapply(factors[merged], all.vars)) > 0L]
  split_pair <- any(vapply(merged_steps, has_partner, logical(1), pool = seq_len(constant - 1L)))
  swallowed <- length(merged_steps) >= 2L &&
    !all(vapply(merged_steps, has_partner, logical(1), pool = merged_steps))
  if (split_pair || swallowed) {
    glue("{path_text}, where the drawing merges the factors after the constant {path[constant]} into one edge")
  }
}

.mpt_is_quotient_with_complement <- function(parsed) {
  is_quotient <- is.call(parsed) && identical(parsed[[1]], quote(`/`))
  if (!is_quotient) {
    return(FALSE)
  }
  numerator_factors <- function(node) {
    node <- .mpt_strip_parens(node)
    if (is.call(node) && identical(node[[1]], quote(`*`))) {
      c(numerator_factors(node[[2]]), numerator_factors(node[[3]]))
    } else {
      list(node)
    }
  }
  any(vapply(numerator_factors(parsed[[2]]), function(f) !is.null(.mpt_complement_base(f)), logical(1)))
}

# merges the paths of all branch expressions into a prefix trie and computes
# plot coordinates. Node ids count in row order, and every edge records its
# parent and child ids and the original, unrounded expression of its label. At
# every node, leaves and child subtrees are stacked top to bottom in the order
# the expressions first mention them, with each complement `1 - X` placed below
# its sibling `X`; internal nodes are centered on their children. Only the
# shared process prefixes of the paths merge: the final factor of each path is
# the category-specific weight and stays a separate leaf edge, so parallel
# branches with equal weights (e.g. uniform guessing fans) do not collapse
.mpt_tree_graph <- function(tree) {
  trie <- list(children = list(), leaves = list())
  for (resp_cat in names(tree$branches)) {
    for (path in .mpt_branch_paths(tree$branches[[resp_cat]])) {
      trie <- .mpt_trie_insert(trie, path, resp_cat)
    }
  }
  layout <- .mpt_layout_subtree(trie, depth = 0, next_y = 0, next_id = 1L)
  nodes <- do.call(rbind, layout$nodes)
  edges <- do.call(rbind, layout$edges)
  rownames(nodes) <- NULL
  rownames(edges) <- NULL
  nlist(nodes, edges)
}

.mpt_trie_insert <- function(node, factors, category) {
  if (length(factors) == 1L) {
    node$leaves <- c(node$leaves, list(list(
      label = factors[[1]], category = category,
      rank = length(node$leaves) + length(node$children)
    )))
    return(node)
  }
  key <- factors[[1]]
  child <- node$children[[key]] %||% list(
    children = list(), leaves = list(),
    rank = length(node$leaves) + length(node$children)
  )
  node$children[[key]] <- .mpt_trie_insert(child, factors[-1], category)
  node
}

# a complement edge `1 - X` is drawn right below its sibling `X` (the success
# edge on top, as in published tree diagrams); every other edge keeps the order
# in which the expressions first mention it
.mpt_sibling_order <- function(labels, ranks) {
  parsed <- lapply(labels, str2lang)
  texts <- vapply(parsed, deparse1, character(1))
  keys <- vapply(seq_along(labels), function(i) {
    base <- .mpt_complement_base(parsed[[i]])
    anchors <- if (!is.null(base)) ranks[texts == deparse1(base)]
    if (length(anchors) > 0L) max(anchors) + 0.5 else ranks[i]
  }, numeric(1))
  order(keys, ranks)
}

# lays out one subtree: leaves and child subtrees take the next free rows from
# the top in the order given by .mpt_sibling_order(), and the node itself is
# centered on its children. Node ids count up in the order the node rows are
# created, so an id is also the row of its node in the final table. Returns the
# node's id and y, the next free row and id, and the node and edge rows of the
# subtree.
.mpt_layout_subtree <- function(node, depth, next_y, next_id) {
  nodes <- list()
  edges <- list()
  children <- list()
  entries <- c(
    lapply(node$leaves, function(leaf) c(leaf, kind = "leaf")),
    lapply(names(node$children), function(key) {
      list(kind = "child", label = key, rank = node$children[[key]]$rank)
    })
  )
  sibling_order <- .mpt_sibling_order(
    vapply(entries, `[[`, character(1), "label"),
    vapply(entries, `[[`, numeric(1), "rank")
  )
  for (entry in entries[sibling_order]) {
    if (entry$kind == "leaf") {
      nodes <- c(nodes, list(
        data.frame(id = next_id, x = depth + 1, y = next_y, category = entry$category)
      ))
      children <- c(children, list(
        list(id = next_id, x = depth + 1, y = next_y, label = entry$label)
      ))
      next_id <- next_id + 1L
      next_y <- next_y - 1
    } else {
      subtree <- .mpt_layout_subtree(
        node$children[[entry$label]], depth + 1, next_y, next_id
      )
      next_y <- subtree$next_y
      next_id <- subtree$next_id
      nodes <- c(nodes, subtree$nodes)
      edges <- c(edges, subtree$edges)
      children <- c(children, list(
        list(id = subtree$id, x = depth + 1, y = subtree$y, label = entry$label)
      ))
    }
  }
  y <- mean(vapply(children, `[[`, numeric(1), "y"))
  nodes <- c(nodes, list(
    data.frame(id = next_id, x = depth, y = y, category = NA_character_)
  ))
  edges <- c(edges, lapply(children, function(child) {
    data.frame(
      from = next_id, to = child$id,
      x0 = depth, y0 = y, x1 = child$x, y1 = child$y,
      label = .mpt_edge_label(child$label), expr = child$label
    )
  }))
  nlist(id = next_id, y, next_y, next_id = next_id + 1L, nodes, edges)
}

# constant factors are displayed rounded (1/15 folds to 0.0666666666666667,
# which would be unreadable as an edge label); a variable-free call that
# cannot be evaluated falls back to the unevaluated label
.mpt_edge_label <- function(factor_label) {
  parsed <- str2lang(factor_label)
  if (length(all.vars(parsed)) == 0L) {
    tryCatch(
      as.character(signif(eval(parsed), 3)),
      error = function(e) factor_label
    )
  } else {
    factor_label
  }
}
