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
#' @details The diagram reflects the branch expressions as written: the first
#'   factor of each product is the edge leaving the root, the second factor
#'   the next edge, and so on. Expressions are expanded distributively, so a
#'   factored subtree like `Pm * (Pb + (1 - Pb) * 0.25)` is displayed as two
#'   paths sharing the `Pm` edge.
#'
#' @return The object `x`, invisibly.
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
  graph <- .mpt_tree_graph(x)
  old_par <- graphics::par(mar = c(0.5, 0.5, 2, 0.5))
  on.exit(graphics::par(old_par))

  x_max <- max(graph$nodes$x)
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(-0.1, x_max + 0.9),
    ylim = c(min(graph$nodes$y) - 0.5, max(graph$nodes$y) + 0.5)
  )
  graphics::segments(
    graph$edges$x0, graph$edges$y0, graph$edges$x1, graph$edges$y1,
    col = "grey40"
  )
  graphics::text(
    (graph$edges$x0 + graph$edges$x1) / 2,
    (graph$edges$y0 + graph$edges$y1) / 2,
    labels = graph$edges$label, pos = 3, offset = 0.2, cex = cex
  )
  is_leaf <- !is.na(graph$nodes$category)
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
  invisible(x)
}

#' @rdname plot.mpt_tree
#' @export
plot.mpt <- function(x, cex = 0.9, ...) {
  n_trees <- length(x$trees)
  if (n_trees > 1) {
    n_row <- floor(sqrt(n_trees))
    old_par <- graphics::par(mfrow = c(n_row, ceiling(n_trees / n_row)))
    on.exit(graphics::par(old_par))
  }
  for (tree in x$trees) {
    plot(tree, cex = cex)
  }
  invisible(x)
}

# expands a branch probability expression into its root-to-leaf paths: sums
# split into separate paths, products append their factors in writing order,
# and any other subexpression (parameters, complements, constants, covariates)
# is one atomic edge label
.mpt_branch_paths <- function(expr) {
  expand <- function(node) {
    if (is.call(node)) {
      op <- as.character(node[[1]])
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

# merges the paths of all branch expressions into a prefix trie and computes
# plot coordinates: leaves are stacked top to bottom in specification order,
# internal nodes are centered on their children. Only the shared process
# prefixes of the paths merge; the final factor of each path is the
# category-specific weight and stays a separate leaf edge, so parallel
# branches with equal weights (e.g. uniform guessing fans) do not collapse
.mpt_tree_graph <- function(tree) {
  trie <- list(children = list(), leaves = list())
  for (resp_cat in names(tree$branches)) {
    for (path in .mpt_branch_paths(tree$branches[[resp_cat]])) {
      trie <- .mpt_trie_insert(trie, path, resp_cat)
    }
  }
  layout <- .mpt_layout_subtree(trie, depth = 0, next_y = 0)
  nodes <- do.call(rbind, layout$nodes)
  edges <- do.call(rbind, layout$edges)
  rownames(nodes) <- NULL
  rownames(edges) <- NULL
  nlist(nodes, edges)
}

.mpt_trie_insert <- function(node, factors, category) {
  if (length(factors) == 1L) {
    node$leaves <- c(
      node$leaves, list(list(label = factors[[1]], category = category))
    )
    return(node)
  }
  key <- factors[[1]]
  child <- node$children[[key]] %||% list(children = list(), leaves = list())
  node$children[[key]] <- .mpt_trie_insert(child, factors[-1], category)
  node
}

# lays out one subtree: leaves take the next free rows from the top, child
# subtrees follow, and the node itself is centered on its children. Returns
# the node's y, the next free row, and the node and edge rows of the subtree.
.mpt_layout_subtree <- function(node, depth, next_y) {
  nodes <- list()
  edges <- list()
  children <- list()
  for (leaf in node$leaves) {
    nodes <- c(nodes, list(
      data.frame(x = depth + 1, y = next_y, category = leaf$category)
    ))
    children <- c(children, list(
      list(x = depth + 1, y = next_y, label = .mpt_edge_label(leaf$label))
    ))
    next_y <- next_y - 1
  }
  for (key in names(node$children)) {
    subtree <- .mpt_layout_subtree(node$children[[key]], depth + 1, next_y)
    next_y <- subtree$next_y
    nodes <- c(nodes, subtree$nodes)
    edges <- c(edges, subtree$edges)
    children <- c(children, list(
      list(x = depth + 1, y = subtree$y, label = .mpt_edge_label(key))
    ))
  }
  y <- mean(vapply(children, `[[`, numeric(1), "y"))
  nodes <- c(nodes, list(data.frame(x = depth, y = y, category = NA_character_)))
  edges <- c(edges, lapply(children, function(child) {
    data.frame(x0 = depth, y0 = y, x1 = child$x, y1 = child$y, label = child$label)
  }))
  nlist(y, next_y, nodes, edges)
}

# constant factors are displayed rounded (1/15 folds to 0.0666666666666667,
# which would be unreadable as an edge label)
.mpt_edge_label <- function(factor_label) {
  parsed <- str2lang(factor_label)
  if (length(all.vars(parsed)) == 0L) {
    as.character(signif(eval(parsed), 3))
  } else {
    factor_label
  }
}
