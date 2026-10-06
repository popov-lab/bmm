pairs_tree <- function(...) {
  mpt_tree("pairs", list(
    E1 = "c * r",
    E2 = "(1 - c) * u^2",
    E3 = "(1 - c) * 2 * u * (1 - u)",
    E4 = "c * (1 - r) + (1 - c) * (1 - u)^2"
  ), ...)
}

source_trees <- function() {
  list(
    mpt_tree("sourceA", list(
      A = "D1 * d1 + D1 * (1 - d1) * a + (1 - D1) * b * g",
      B = "D1 * (1 - d1) * (1 - a) + (1 - D1) * b * (1 - g)",
      N = "(1 - D1) * (1 - b)"
    )),
    mpt_tree("sourceB", list(
      A = "D2 * (1 - d2) * a + (1 - D2) * b * g",
      B = "D2 * d2 + D2 * (1 - d2) * (1 - a) + (1 - D2) * b * (1 - g)",
      N = "(1 - D2) * (1 - b)"
    )),
    mpt_tree("new_items", list(A = "b * g", B = "b * (1 - g)", N = "1 - b"))
  )
}

oberauer_2019_tree <- function() {
  mpt_tree("ss", list(
    correct = "Pb + (1 - Pb) * Pi * GcorrPi + (1 - Pb) * (1 - Pi) * GcorrNoPi",
    other = paste0(
      "(1 - Pb) * Pi * (1 - GcorrPi) + ",
      "(1 - Pb) * (1 - Pi) * (1 - GcorrNoPi) * GotherNoPi"
    ),
    npl = "(1 - Pb) * (1 - Pi) * (1 - GcorrNoPi) * (1 - GotherNoPi)"
  ))
}

oberauer_lin_trees <- function() {
  list(
    mpt_tree("newdist", list(
      corr = "Pi * Pb + Pi * (1 - Pb) * (1 - Pd) * (1/5) + (1 - Pi) * (1/15)",
      other = "Pi * (1 - Pb) * (1 - Pd) * (4/5) + (1 - Pi) * (4/15)",
      dist = "Pi * (1 - Pb) * Pd + (1 - Pi) * (5/15)",
      npl = "(1 - Pi) * (5/15)"
    )),
    mpt_tree("nodist", list(
      corr = "Pi * Pb + Pi * (1 - Pb) * (1/5) + (1 - Pi) * (1/15)",
      other = "Pi * (1 - Pb) * (4/5) + (1 - Pi) * (4/15)",
      npl = "(1 - Pi) * (10/15)"
    ), impossible = "dist")
  )
}

simplex_trees <- function() {
  list(
    mpt_tree("sourceA", list(
      A = "dA + (1 - dA) * gA", B = "(1 - dA) * gB", New = "(1 - dA) * gNew"
    )),
    mpt_tree("sourceB", list(
      A = "(1 - dB) * gA", B = "dB + (1 - dB) * gB", New = "(1 - dB) * gNew"
    )),
    mpt_tree("new", list(A = "gA", B = "gB", New = "gNew"))
  )
}

# one wide, one narrow leaf label and long edge labels, which near the root
# stick out to the left of steep edges
long_label_tree <- function() {
  mpt_tree("long", list(
    Correct_Correct_Long_Category_Name =
      "Very_Long_Parameter_Name_One * Very_Long_Parameter_Name_Two",
    Other_Response_Category_With_A_Long_Name =
      "Very_Long_Parameter_Name_One * (1 - Very_Long_Parameter_Name_Two)",
    Rest = "1 - Very_Long_Parameter_Name_One"
  ))
}

# draws the tree on a pdf device and returns where every label starts and ends,
# in inches from the left edge of the plot region, with the plot width
label_extents <- function(tree, width, height, cex = 0.9, pointsize = 12) {
  extents <- data.frame()
  record <- function(x, y = NULL, labels = y, ...) {
    args <- list(...)
    usr <- graphics::par("usr")
    start <- (x - usr[1]) / (usr[2] - usr[1]) * graphics::par("pin")[1]
    size <- graphics::strwidth(labels, units = "inches", cex = args$cex, font = args$font %||% 1)
    shift <- if (is.null(args$pos)) {
      -size / 2
    } else {
      args$offset * graphics::strwidth("m", units = "inches", cex = args$cex)
    }
    extents <<- rbind(extents, data.frame(
      kind = if (is.null(args$pos)) "edge" else "leaf",
      left = start + shift, right = start + shift + size, plot_width = graphics::par("pin")[1]
    ))
  }
  withr::with_pdf(
    NULL,
    testthat::with_mocked_bindings(plot(tree, cex = cex), text = record, .package = "graphics"),
    width = width, height = height, pointsize = pointsize
  )
  extents
}

# draws the model on a pdf device and returns the pairs of edge labels whose
# boxes overlap, measured in inches and 30% taller than the text height (which
# leaves out descenders and parentheses)
overlapping_edge_labels <- function(model, width, height, cex = 0.9) {
  boxes <- data.frame()
  record <- function(x, y = NULL, labels = y, ...) {
    args <- list(...)
    if (!is.null(args$pos)) {
      return(invisible())
    }
    usr <- graphics::par("usr")
    boxes <<- rbind(boxes, data.frame(
      x = (x - usr[1]) / (usr[2] - usr[1]) * graphics::par("pin")[1],
      y = (y - usr[3]) / (usr[4] - usr[3]) * graphics::par("pin")[2],
      width = graphics::strwidth(labels, units = "inches", cex = args$cex),
      height = 1.3 * graphics::strheight(labels, units = "inches", cex = args$cex),
      label = labels
    ))
  }
  withr::with_pdf(
    NULL,
    testthat::with_mocked_bindings(
      suppressWarnings(plot(model, cex = cex)), text = record, .package = "graphics"
    ),
    width = width, height = height
  )
  overlap_x <- outer(boxes$width, boxes$width, `+`) / 2 - abs(outer(boxes$x, boxes$x, `-`))
  overlap_y <- outer(boxes$height, boxes$height, `+`) / 2 - abs(outer(boxes$y, boxes$y, `-`))
  pairs <- which(overlap_x > 0 & overlap_y > 0 & upper.tri(overlap_x), arr.ind = TRUE)
  as.character(paste(boxes$label[pairs[, 1]], "<>", boxes$label[pairs[, 2]])[nrow(pairs) > 0])
}

test_that("branch expressions expand into root-to-leaf paths", {
  paths <- .mpt_branch_paths(quote(D + (1 - D) * g))
  expect_equal(paths, list("D", c("1 - D", "g")))

  # factored subtrees are expanded distributively
  paths_factored <- .mpt_branch_paths(quote(Pm * (Pb + (1 - Pb) * 0.25)))
  expect_equal(paths_factored, list(c("Pm", "Pb"), c("Pm", "1 - Pb", "0.25")))
})

test_that("a constant factor above 1 and the factors after it form one edge", {
  expect_equal(
    .mpt_branch_paths(quote((1 - c) * 2 * u * (1 - u))),
    list(c("1 - c", "2 * u * (1 - u)"))
  )
  expect_equal(
    .mpt_branch_paths(quote(k * 2 * (1 - u) * v)),
    list(c("k", "2 * (1 - u) * v"))
  )
  # a folded constant counts by its value, not by how it is written
  expect_equal(
    .mpt_branch_paths(quote(k * (6 / 3) * u)),
    list(c("k", "6/3 * u"))
  )
})

test_that("a leading constant above 1 makes the whole path one edge", {
  expect_equal(
    .mpt_branch_paths(quote(2 * (1 - c) * u * (1 - u))),
    list("2 * (1 - c) * u * (1 - u)")
  )
  expect_equal(.mpt_branch_paths(quote(2 * u)), list("2 * u"))
  expect_equal(.mpt_branch_paths(quote(2)), list("2"))
})

test_that("a trailing constant above 1 joins the factor before it", {
  expect_equal(
    .mpt_branch_paths(quote((1 - c) * u * (1 - u) * 2)),
    list(c("1 - c", "u", "(1 - u) * 2"))
  )
  expect_equal(.mpt_branch_paths(quote((1 - c) * 3)), list("(1 - c) * 3"))
})

test_that("a constant between factors takes everything after it", {
  expect_equal(
    .mpt_branch_paths(quote(u * 2 * (1 - u))),
    list(c("u", "2 * (1 - u)"))
  )
})

test_that("constants of 1 or less and a variable-free call stay ordinary edges", {
  expect_equal(
    .mpt_branch_paths(quote(0.2 * x + (1 - x) * (1 / 3))),
    list(c("0.2", "x"), c("1 - x", "1/3"))
  )
  expect_equal(.mpt_branch_paths(quote(1 * x)), list(c("1", "x")))
  expect_equal(.mpt_branch_paths(quote(foo() * D)), list(c("foo()", "D")))
})

test_that("pair clustering draws the multiplicity as one edge below 1 - c", {
  graph <- .mpt_tree_graph(pairs_tree())
  expect_equal(
    graph$edges$label,
    c("r", "1 - r", "u^2", "2 * u * (1 - u)", "(1 - u)^2", "c", "1 - c")
  )
  expect_equal(graph$nodes$category, c("E1", "E4", NA, "E2", "E3", "E4", NA, NA))
  expect_equal(graph$edges$from, c(3, 3, 7, 7, 7, 8, 8))
  expect_equal(graph$edges$to, c(1, 2, 4, 5, 6, 3, 7))
})

test_that("a package tree with the constant written first keeps it in one root edge", {
  study <- mpt_tree("study", list(
    C = "cp + (1 - cp) * rp * rp",
    E = "2 * (1 - cp) * rp * (1 - rp)",
    U = "(1 - cp) * (1 - rp) * (1 - rp)"
  ))
  graph <- .mpt_tree_graph(study)
  root_edges <- graph$edges[graph$edges$x0 == 0, ]
  expect_setequal(root_edges$label, c("cp", "1 - cp", "2 * (1 - cp) * rp * (1 - rp)"))
  expect_false("2" %in% graph$edges$label)
})

test_that("plot() warns once per tree when a constant above 1 hides a later process step", {
  withr::local_pdf(NULL)
  swallowed <- mpt_tree("swallowed", list(
    A = "D * 2 * g * (1 - g) * r",
    B = "D * 2 * g * (1 - g) * (1 - r)",
    C = "D * (g^2 + (1 - g)^2)",
    E = "1 - D"
  ))
  messages <- character(0)
  withCallingHandlers(
    plot(swallowed),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(messages, 1)
  expect_match(messages, "In the tree 'swallowed'", fixed = TRUE)
  expect_match(messages, "branch 'A': D * 2 * g * (1 - g) * r", fixed = TRUE)
  expect_match(messages, "branch 'B': D * 2 * g * (1 - g) * (1 - r)", fixed = TRUE)
  expect_match(messages, "merges the factors after the constant 2", fixed = TRUE)

  # the constant splits a factor from its complement
  expect_warning(
    plot(mpt_tree("split", list(a = "u * 2 * (1 - u)", b = "1 - u"))),
    "branch 'a': u * 2 * (1 - u)", fixed = TRUE
  )
})

test_that("plot() warns when a quotient hides a complement", {
  withr::local_pdf(NULL)
  quotient <- mpt_tree("quotient", list(
    a = "(1 - Pi) * 1 / 15", b = "Pi + (1 - Pi) * 14 / 15"
  ))
  expect_warning(
    plot(quotient),
    "branch 'a': (1 - Pi) * 1/15 is one edge", fixed = TRUE
  )
})

test_that("plot() stays silent about constants that merge cleanly", {
  withr::local_pdf(NULL)
  silent <- list(
    pairs_tree(),
    mpt_tree("leading", list(a = "2 * (1 - c) * u * (1 - u)", b = "c")),
    mpt_tree("trailing", list(a = "(1 - c) * u * (1 - u) * 2", b = "c")),
    mpt_tree("single", list(a = "(1 - c) * 3", b = "c")),
    mpt_tree("one_step", list(a = "D * 2 * r", b = "1 - D")),
    mpt_tree("constants", list(a = "(1 - Pi) * 10 * (1/15)", b = "Pi")),
    mpt_tree("folded", list(a = "(1 - Pi) * (1/15)", b = "Pi + (1 - Pi) * (14/15)"))
  )
  for (tree in silent) {
    expect_silent(plot(tree))
  }
})

test_that("the Oberauer (2019) tree keeps its structure", {
  graph <- .mpt_tree_graph(oberauer_2019_tree())
  expect_equal(nrow(graph$nodes), 11)
  expect_equal(nrow(graph$edges), 10)
  leaves <- graph$nodes[!is.na(graph$nodes$category), ]
  expect_equal(
    leaves$category,
    c("correct", "correct", "other", "correct", "other", "npl")
  )
  expect_equal(
    graph$edges$label[match(leaves$id, graph$edges$to)],
    c("Pb", "GcorrPi", "1 - GcorrPi", "GcorrNoPi", "GotherNoPi", "1 - GotherNoPi")
  )
})

test_that("node ids identify the nodes and the edges refer to them", {
  trees <- c(
    mpt_2htm_trees()[1], list(pairs_tree(), mpt_tree("one", list(only = "1")))
  )
  for (tree in trees) {
    graph <- .mpt_tree_graph(tree)
    nodes <- graph$nodes
    edges <- graph$edges
    expect_type(nodes$id, "integer")
    expect_equal(nodes$id, seq_len(nrow(nodes)))
    expect_true(all(c(edges$from, edges$to) %in% nodes$id))
    # a tree: every node but the root has exactly one incoming edge
    expect_equal(sort(edges$to), sort(setdiff(nodes$id, setdiff(edges$from, edges$to))))
    expect_false(anyDuplicated(edges$to) > 0)
    from <- nodes[match(edges$from, nodes$id), ]
    to <- nodes[match(edges$to, nodes$id), ]
    expect_equal(c(edges$x0, edges$y0), c(from$x, from$y))
    expect_equal(c(edges$x1, edges$y1), c(to$x, to$y))
  }
})

test_that("the tree graph merges shared path prefixes only", {
  graph <- .mpt_tree_graph(mpt_2htm_trees()[[1]])
  leaves <- graph$nodes[!is.na(graph$nodes$category), ]
  expect_equal(leaves$category, c("old", "old", "new"))
  # internal nodes: the root and the shared (1 - D) node
  expect_equal(sum(is.na(graph$nodes$category)), 2)
  expect_setequal(graph$edges$label, c("D", "1 - D", "g", "1 - g"))
})

test_that("equal terminal weights stay separate leaf edges", {
  fan <- mpt_tree("t", list(
    a = "m + (1 - m) * (1/3)",
    b = "(1 - m) * (1/3)",
    c = "(1 - m) * (1/3)"
  ))
  graph <- .mpt_tree_graph(fan)
  leaves <- graph$nodes[!is.na(graph$nodes$category), ]
  expect_equal(leaves$category, c("a", "a", "b", "c"))
  # constant edge labels are rounded for display
  expect_equal(sum(graph$edges$label == "0.333"), 3)
})

test_that("a three-level tree has one node per process state and centers its internal nodes", {
  tree <- mpt_tree("t", list(
    a = "D * b * c",
    b = "D * b * (1 - c)",
    c = "D * (1 - b)",
    d = "1 - D"
  ))
  graph <- .mpt_tree_graph(tree)
  # root, the D node, the D > b node, and four leaves
  expect_equal(nrow(graph$nodes), 7)
  expect_equal(nrow(graph$edges), 6)
  leaves <- graph$nodes[!is.na(graph$nodes$category), ]
  expect_equal(leaves$category, c("a", "b", "c", "d"))
  expect_equal(leaves$y, c(0, -1, -2, -3))
  expect_equal(leaves$x, c(3, 3, 2, 1))
  expect_true(all(graph$edges$x1 == graph$edges$x0 + 1))
  expect_equal(sum(graph$nodes$x == 0), 1)
  # the D > b node has leaves at 0 and -1, the D node has that node and a
  # leaf at -2, the root has the D node and a leaf at -3
  internal <- graph$nodes[is.na(graph$nodes$category), ]
  expect_equal(internal$y[order(internal$x)], c(-2.125, -1.25, -0.5))
})

test_that("a node lists its process children in the order they are first mentioned", {
  tree <- mpt_tree("t", list(
    a = "D * m * x",
    b = "D * m * (1 - x)",
    c = "(1 - D) * n * y",
    d = "(1 - D) * n * (1 - y)"
  ))
  graph <- .mpt_tree_graph(tree)
  root_edges <- graph$edges[graph$edges$x0 == 0, ]
  expect_equal(root_edges$label[order(-root_edges$y1)], c("D", "1 - D"))
})

test_that("leaves and subtrees keep the order of first mention", {
  # source monitoring, new items: the success edge b comes before 1 - b
  new_items <- source_trees()[[3]]
  graph <- .mpt_tree_graph(new_items)
  root_edges <- graph$edges[graph$edges$x0 == 0, ]
  expect_equal(root_edges$label[order(-root_edges$y1)], c("b", "1 - b"))
  leaves <- graph$nodes[!is.na(graph$nodes$category), ]
  expect_equal(leaves$category, c("A", "B", "N"))

  leaf_categories <- function(branches) {
    nodes <- .mpt_tree_graph(mpt_tree("t", branches))$nodes
    nodes$category[!is.na(nodes$category)]
  }
  # a leaf listed first is drawn first, unless it is the complement of a sibling
  expect_equal(leaf_categories(list(z = "g", a = "D * b", b = "D * (1 - b)")), c("z", "a", "b"))
  # a leaf mentioned after a subtree is drawn below it
  expect_equal(leaf_categories(list(a = "D * b", b = "D * (1 - b)", n = "k")), c("a", "b", "n"))
  expect_equal(
    leaf_categories(list(d = "1 - D", a = "D * b", b = "D * (1 - b)")),
    c("a", "b", "d")
  )
})

test_that("a complement edge follows its sibling and other edges keep first-mention order", {
  expect_equal(.mpt_sibling_order(c("1 - D", "g", "D"), c(0, 1, 2)), c(2, 3, 1))
  # a complement without its sibling, plain symbols and constants keep their rank
  expect_equal(.mpt_sibling_order(c("1 - b", "gB", "gA", "0.2"), c(0, 1, 2, 3)), 1:4)
  expect_equal(.mpt_sibling_order(c("gA", "gB", "gNew"), c(0, 1, 2)), 1:3)
  # repeated siblings: the complement goes below the last X, ties by first mention
  expect_equal(.mpt_sibling_order(c("g", "g", "1 - g"), c(0, 1, 2)), 1:3)
  expect_equal(.mpt_sibling_order(c("1 - g", "1 - g", "g"), c(2, 1, 0)), c(3, 2, 1))
  expect_equal(.mpt_sibling_order(c("2 - D", "D"), c(0, 1)), 1:2)
  # redundant parentheses do not hide the pair
  expect_equal(.mpt_sibling_order(c("1 - (D)", "D"), c(0, 1)), c(2, 1))
})

test_that("the success edge X is drawn above its complement 1 - X in published trees", {
  trees <- c(
    mpt_2htm_trees(), source_trees(), oberauer_lin_trees()[1], list(
      pairs_tree(),
      mpt_tree("three_levels", list(
        a = "D * b * c", b = "D * b * (1 - c)", c = "D * (1 - b)", d = "1 - D"
      )),
      mpt_tree("complement_listed_first", list(
        d = "1 - D", a = "D * b", b = "D * (1 - b)"
      ))
    )
  )
  complement_pairs_checked <- 0
  for (tree in trees) {
    edges <- .mpt_tree_graph(tree)$edges
    for (siblings in split(edges, edges$from)) {
      for (label in siblings$label[startsWith(siblings$label, "1 - ")]) {
        success <- siblings$y1[siblings$label == substring(label, 5)]
        if (length(success) == 0L) next
        complement_pairs_checked <- complement_pairs_checked + 1
        expect_true(
          all(success > siblings$y1[siblings$label == label]),
          label = glue("{tree$name}: {substring(label, 5)} above {label}")
        )
      }
    }
  }
  expect_gt(complement_pairs_checked, 15)
})

test_that("plot methods leave par() as it was, on success and on error", {
  tree <- mpt_tree("t", list(a = "D", b = "1 - D"))
  withr::local_pdf(NULL)
  before <- graphics::par(c("mfrow", "mar"))
  plot(tree)
  expect_equal(graphics::par(c("mfrow", "mar")), before)

  local_mocked_bindings(.mpt_edge_label_positions = function(edges, cex) stop("boom"))
  expect_error(plot(tree), "boom")
  expect_equal(graphics::par(c("mfrow", "mar")), before)
  expect_error(plot(mpt(mpt_2htm_trees(), tree_id = "item_type")), "boom")
  expect_equal(graphics::par(c("mfrow", "mar")), before)
})

test_that("plot methods return x invisibly", {
  tree <- mpt_tree("t", list(a = "D", b = "1 - D"))
  model <- mpt(list(tree), tree_id = "item")

  withr::local_pdf(NULL)
  tree_result <- withVisible(plot(tree))
  expect_identical(tree_result$value, tree)
  expect_false(tree_result$visible)

  model_result <- withVisible(plot(model))
  expect_identical(model_result$value, model)
  expect_false(model_result$visible)
})

test_that("namespaced function calls in edge labels plot correctly", {
  tree_ns <- mpt_tree("ns_call", list(
    A = "stats::plogis(a) * g",
    B = "stats::plogis(a) * (1 - g)",
    C = "1 - stats::plogis(a)"
  ))
  expect_true(any(grepl("stats::plogis", .mpt_tree_graph(tree_ns)$edges$label)))
  withr::local_pdf(NULL)
  expect_silent(plot(tree_ns))
})

test_that("variable-free calls that cannot be evaluated fall back to label", {
  tree_no_eval <- mpt_tree("no_eval", list(D = "foo() * D", not_D = "1 - D"))
  expect_silent(graph <- .mpt_tree_graph(tree_no_eval))
  expect_true(any(grepl("foo()", graph$edges$label, fixed = TRUE)))
  withr::local_pdf(NULL)
  expect_silent(plot(tree_no_eval))
})

test_that("models whose sibling edges sum to 1 plot without a warning", {
  withr::local_pdf(NULL)
  old_new <- mpt_2htm_trees()
  pair_clustering <- list(
    pairs_tree(impossible = c("F1", "F2")),
    mpt_tree("singletons", list(F1 = "u", F2 = "1 - u"), impossible = c("E1", "E2", "E3", "E4"))
  )
  models <- list(
    two_htm = mpt(old_new, tree_id = "item_type"),
    pair_clustering = mpt(pair_clustering, tree_id = "item_type"),
    source_monitoring = mpt(source_trees(), tree_id = "source"),
    oberauer_2019 = mpt(
      oberauer_2019_tree(),
      covariates = c("GcorrPi", "GcorrNoPi", "GotherNoPi")
    ),
    oberauer_lin = mpt(oberauer_lin_trees(), tree_id = "tree"),
    simplex = mpt(simplex_trees(), tree_id = "source", simplex = list(c("gA", "gB", "gNew"))),
    constants_summing_to_1 = mpt(old_new, tree_id = "item_type", restrictions = "g = 0.5"),
    equated = mpt(list(
      mpt_tree("old", list(old = "Do + (1 - Do) * g", new = "(1 - Do) * (1 - g)")),
      mpt_tree("new", list(old = "(1 - Dn) * g", new = "Dn + (1 - Dn) * (1 - g)"))
    ), tree_id = "item_type", restrictions = "Dn = Do")
  )
  for (name in names(models)) {
    expect_silent(plot(models[[name]]))
  }
  # a bare tree cannot know the simplex: fans of plain symbols are not checked
  expect_silent(plot(simplex_trees()[[1]]))
  expect_silent(plot(simplex_trees()[[3]]))
  expect_silent(plot(pair_clustering[[1]]))
})

test_that("valid models whose edges are not simple probabilities plot without a warning", {
  withr::local_pdf(NULL)
  trees <- list(
    complement_of_a_product = mpt_tree("t", list(hit = "D * G", miss = "1 - D * G")),
    factor_with_a_constant = mpt_tree("t", list(
      yes = "(Do) + ((1 - Do) * guess * 0.25)", no = "(1 - Do) * (1 - guess * 0.25)"
    )),
    multiplicity_first = mpt_tree("study", list(
      C = "cp + (1 - cp) * rp * rp",
      E = "2 * (1 - cp) * rp * (1 - rp)",
      U = "(1 - cp) * (1 - rp) * (1 - rp)"
    )),
    multiplicity_last = mpt_tree("t", list(
      a = "(1 - c) * u^2", b = "(1 - c) * u * (1 - u) * 2", d = "(1 - c) * (1 - u)^2", e = "c"
    )),
    zero_product = mpt_tree("t", list(x = "0 * D", y = "1 - 0 * D"))
  )
  for (name in names(trees)) {
    expect_silent(plot(trees[[name]]))
    expect_silent(plot(mpt(trees[[name]])))
  }
  in_a_node <- mpt_tree("t", list(
    hit = "G * D + (1 - G * D) * gA", other = "(1 - G * D) * gB", miss = "(1 - G * D) * gC"
  ))
  expect_silent(plot(mpt(in_a_node, simplex = list(c("gA", "gB", "gC")))))
  # quotient edges are not checked for their sum; the quotient warning is
  # tested apart
  quotients <- list(
    mpt(mpt_tree("t", list(
      a = "D + (1 - D) * g / 2", b = "(1 - D) * g / 2", c = "(1 - D) * (1 - g)"
    ))),
    mpt(mpt_tree("t", list(
      hit = "D + (1 - D)/ss", miss = "(1 - D) * (1 - 1/ss)"
    )), covariates = "ss")
  )
  for (model in quotients) {
    messages <- character(0)
    withCallingHandlers(
      plot(model),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_false(any(grepl("do not sum to 1", messages, fixed = TRUE)))
  }
})

test_that("covariates that sum to 1 only in the data are not checked", {
  withr::local_pdf(NULL)
  fan <- mpt_tree("cued", list(correct = "D + (1 - D) * Gcorr", incorrect = "(1 - D) * Gother"))
  model <- suppressWarnings(mpt(fan, covariates = c("Gcorr", "Gother")))
  expect_silent(plot(model))
  # a covariate with its complement does sum to 1, and is checked with the rest of the node
  wrong <- mpt_tree("w", list(
    a = "D * Gc", b = "D * (1 - Gc)", c = "D * 0.5", d = "1 - D"
  ))
  expect_warning(
    plot(suppressWarnings(mpt(wrong, covariates = "Gc"))),
    "the node 'D' do not sum to 1: 'Gc', '1 - Gc', '0.5'",
    fixed = TRUE
  )
  swapped <- mpt_tree("s", list(a = "D * G", b = "(1 - G) * D", c = "1 - D"))
  expect_warning(
    plot(mpt(swapped, covariates = "G")),
    "the root do not sum to 1: 'D', '1 - D', '1 - G'",
    fixed = TRUE
  )
})

test_that("a node whose edges do not sum to 1 warns once, names the node, and still draws", {
  withr::local_pdf(NULL)
  mismatch <- mpt_tree("mismatch", list(a = "x * y", b = "(1 - y) * x", c = "1 - x"))
  expect_warning(
    result <- plot(mismatch),
    "In the tree 'mismatch', the edges leaving the root do not sum to 1: 'x', '1 - x', '1 - y'",
    fixed = TRUE
  )
  expect_identical(result, mismatch)

  # one warning per offending node, each naming its node and its edges
  bad <- mpt_tree("bad", list(
    a = "u * v", b = "u * (1 - w)", c = "1 - u", d = "(1 - v) * 0.5", e = "(1 - v) * 0.25"
  ))
  messages <- character(0)
  withCallingHandlers(
    plot(bad),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(messages, 3)
  expect_match(messages[1], "the node 'u' do not sum to 1: 'v', '1 - w'", fixed = TRUE)
  expect_match(messages[2], "the node '1 - v' do not sum to 1: '0.5', '0.25'", fixed = TRUE)
  expect_match(messages[3], "the root do not sum to 1: 'u', '1 - u', '1 - v'", fixed = TRUE)

  # a deeper node is named by its whole path of edges, from the root down
  deep <- mpt_tree("deep", list(
    a = "x * y * z", b = "x * y * (1 - w)", c = "x * (1 - y)", d = "1 - x"
  ))
  expect_warning(plot(deep), "the node 'x > y' do not sum to 1: 'z', '1 - w'", fixed = TRUE)
})

test_that("a swapped complement in a published tree is reported", {
  withr::local_pdf(NULL)
  swapped <- mpt_tree("swapped", list(old = "D + (1 - D) * (1 - g)", new = "(1 - D) * (1 - g)"))
  expect_warning(
    plot(swapped),
    "the node '1 - D' do not sum to 1: '1 - g', '1 - g'",
    fixed = TRUE
  )
  swapped_root <- mpt_tree("swapped_root", list(old = "D + (1 - D) * g", new = "D * (1 - g)"))
  expect_warning(plot(swapped_root), "the root do not sum to 1: 'D', 'D', '1 - D'", fixed = TRUE)
})

test_that("a parenthesised complement is classified like the plain one", {
  withr::local_pdf(NULL)
  expect_silent(plot(mpt_tree("p", list(a = "D", b = "1 - (D)"))))
  expect_warning(
    plot(mpt_tree("p", list(a = "D", b = "1 - (D)", c = "D"))),
    "'D', 'D', '1 - (D)'",
    fixed = TRUE
  )
})

test_that("sibling sums use the original expressions and a tolerance of 1e-6", {
  withr::local_pdf(NULL)
  # the edge labels show 0.333, the sum of the expressions is 1 within 1e-6
  within_tolerance <- mpt_tree("t", list(
    a = "m", b = "(1 - m) * 0.3333333", c = "(1 - m) * 0.3333333", d = "(1 - m) * 0.3333333"
  ))
  expect_silent(plot(within_tolerance))
  fraction <- mpt_tree("t", list(
    a = "m", b = "(1 - m) * (1/3)", c = "(1 - m) * (1/3)", d = "(1 - m) * (1/3)"
  ))
  expect_silent(plot(fraction))
  rounded <- mpt_tree("t", list(
    a = "m", b = "(1 - m) * 0.333", c = "(1 - m) * 0.333", d = "(1 - m) * 0.333"
  ))
  expect_warning(plot(rounded), "the node '1 - m' do not sum to 1", fixed = TRUE)
})

test_that("a simplex fan is checked only when the simplex is known", {
  tree <- simplex_trees()[[1]]
  graph <- .mpt_tree_graph(tree)
  expect_silent(.mpt_warn_sibling_sums(graph, "sA", list(c("gA", "gB", "gNew"))))
  expect_warning(
    .mpt_warn_sibling_sums(graph, "sA", list()),
    "'gA', 'gB', 'gNew'",
    fixed = TRUE
  )
  expect_silent(.mpt_warn_sibling_sums(graph, "sA", NULL))
  # a group member missing from the tree's edges still gets its test value
  partial <- .mpt_tree_graph(mpt_tree("t", list(A = "dA * gA", B = "dA * gB", C = "1 - dA")))
  expect_warning(
    .mpt_warn_sibling_sums(partial, "t", list(c("gA", "gB", "gNew"))),
    "'gA', 'gB'",
    fixed = TRUE
  )
})

test_that("plot.mpt() checks plain-symbol fans against the model's simplex", {
  withr::local_pdf(NULL)
  tree <- mpt_tree("s", list(A = "gA", B = "gB", C = "x * gNew", D = "gNew * (1 - x)"))
  model <- mpt(tree, simplex = list(c("gA", "gB", "gNew")))
  expect_warning(
    plot(model),
    "the root do not sum to 1: 'gA', 'gB', 'x', 'gNew'",
    fixed = TRUE
  )
  expect_silent(plot(tree))
})

test_that("leaf labels fit inside the plot region when they take less than 90% of its width", {
  tree <- long_label_tree()
  fits <- logical(0)
  overflow <- numeric(0)
  settings <- list(
    list(7, 7, 12), list(12, 4, 12), list(6, 4, 12), list(5, 3, 12),
    list(4, 3, 12), list(2.5, 2.5, 12), list(4, 3, 16)
  )
  for (setting in settings) {
    extents <- label_extents(
      tree, setting[[1]], setting[[2]], pointsize = setting[[3]]
    )
    leaves <- extents[extents$kind == "leaf", ]
    widest <- max(leaves$right - leaves$left)
    fits <- c(fits, widest <= 0.9 * leaves$plot_width[1])
    overflow <- c(overflow, max(leaves$right) / leaves$plot_width[1] - 1)
  }
  expect_gte(sum(fits), 4)
  expect_true(any(!fits))
  expect_lte(max(overflow[fits]), 0.02)
})

test_that("edge labels and leaf labels stay inside the plot region", {
  extents <- label_extents(long_label_tree(), 10, 7)
  expect_gte(min(extents$left), -0.02 * extents$plot_width[1])
  expect_lte(max(extents$right), 1.02 * extents$plot_width[1])
  expect_equal(sum(extents$kind == "edge"), 4)
})

test_that("labels of neighbouring edges do not overlap in a dense tree", {
  # four binary process steps with restricted parameters: the long restriction
  # labels of neighbouring edges used to overprint at these sizes
  steps <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
  branches <- apply(steps, 1, function(row) {
    paste(ifelse(row == 1, names(row), paste0("(1 - ", names(row), ")")), collapse = " * ")
  })
  dense <- mpt(
    mpt_tree("dense", as.list(setNames(branches, paste0("c", seq_along(branches))))),
    restrictions = c("B = 0.571", "C = 0.2", "D = 0.125")
  )
  for (size in list(c(4, 3), c(5, 3), c(6, 4), c(8, 6))) {
    expect_identical(overlapping_edge_labels(dense, size[1], size[2]), character(0))
  }
})

test_that("the plot window makes room for the farthest edge label and no more", {
  # the room a label needs depends on the window it is placed in, so the window
  # is found in a few passes; without them the label leaves the plot region
  root_long <- mpt_tree("root_long", list(
    a = "Very_Long_Parameter_Name_One_abc * g",
    b = "Very_Long_Parameter_Name_One_abc * (1 - g)",
    c = "1 - Very_Long_Parameter_Name_One_abc"
  ))
  for (size in list(c(5, 4), c(4.5, 3.5), c(4, 3))) {
    extents <- label_extents(root_long, size[1], size[2])
    leftmost <- min(extents$left[extents$kind == "edge"])
    expect_gte(leftmost, -0.02)
    expect_lt(leftmost, 0.05)
  }
})

test_that("edge labels sit beside their edge on the outer side of a fan", {
  withr::local_pdf(NULL)
  edges <- .mpt_tree_graph(mpt_tree("t", list(
    a = "D", b = "(1 - D) * g", c = "(1 - D) * (1 - g)", d = "(1 - D) * 0.5"
  )))$edges
  graphics::plot.new()
  # unequal user units per inch on the two axes
  graphics::plot.window(xlim = c(-0.1, 40), ylim = c(-4, 1))
  cex <- 0.9
  position <- .mpt_edge_label_positions(edges, cex)
  pin <- graphics::par("pin")
  usr <- graphics::par("usr")
  x_per_inch <- (usr[2] - usr[1]) / pin[1]
  y_per_inch <- (usr[4] - usr[3]) / pin[2]
  dx <- (edges$x1 - edges$x0) / x_per_inch
  dy <- (edges$y1 - edges$y0) / y_per_inch
  # distance of the label centre from the (infinite) edge line, in inches
  distance <- abs(
    (position$x - edges$x0) / x_per_inch * dy - (position$y - edges$y0) / y_per_inch * dx
  ) / sqrt(dx^2 + dy^2)
  width <- graphics::strwidth(edges$label, units = "inches", cex = cex)
  height <- graphics::strheight(edges$label, units = "inches", cex = cex)
  # the label box extends from its centre by this much towards the line
  # (|nx| * w / 2 + |ny| * h / 2), and the gap beyond that is 0.4 label heights
  extent <- (abs(dy) * width / 2 + abs(dx) * height / 2) / sqrt(dx^2 + dy^2)
  expect_equal(distance, 0.4 * height + extent)
  # rising edges carry the label above the line, falling edges below
  line_y <- edges$y0 + (position$x - edges$x0) * (edges$y1 - edges$y0) / (edges$x1 - edges$x0)
  expect_equal(sign(position$y - line_y), ifelse(edges$y1 >= edges$y0, 1, -1))
})

test_that("plot.mpt() lays the trees out on the expected grid", {
  withr::local_pdf(NULL)
  grids <- list()
  local_mocked_bindings(.mpt_plot_tree = function(x, cex, simplex, covariates, restrictions) {
    grids[[length(grids) + 1]] <<- graphics::par("mfrow")
  })
  for (n_trees in 2:6) {
    categories <- paste0(c("a", "b"), rep(seq_len(n_trees), each = 2))
    trees <- lapply(seq_len(n_trees), function(i) {
      own <- categories[(2 * i - 1):(2 * i)]
      mpt_tree(
        paste0("t", i), setNames(list("p", "1 - p"), own),
        impossible = setdiff(categories, own)
      )
    })
    grids <- list()
    plot(mpt(trees, tree_id = "item"))
    expect_length(grids, n_trees)
    expect_true(all(vapply(grids, identical, logical(1), grids[[1]])), label = n_trees)
    expect_equal(
      grids[[1]],
      list(c(1, 2), c(1, 3), c(2, 2), c(2, 3), c(2, 3))[[n_trees - 1]],
      label = glue("grid for {n_trees} trees")
    )
  }
})

test_that("plot.mpt() leaves a multi-panel layout of the user alone", {
  withr::local_pdf(NULL)
  grids <- list()
  local_mocked_bindings(.mpt_plot_tree = function(x, cex, simplex, covariates, restrictions) {
    grids[[length(grids) + 1]] <<- graphics::par("mfrow")
  })
  trees <- list(
    mpt_tree("t1", list(a = "p", b = "1 - p"), impossible = c("c", "d")),
    mpt_tree("t2", list(c = "q", d = "1 - q"), impossible = c("a", "b"))
  )
  graphics::par(mfrow = c(2, 2))
  plot(mpt(trees, tree_id = "item"))
  expect_equal(grids, list(c(2, 2), c(2, 2)))
  expect_equal(graphics::par("mfrow"), c(2, 2))
})

test_that("cex reaches every edge and leaf label and the title shows the tree name", {
  withr::local_pdf(NULL)
  tree <- mpt_2htm_trees()[[1]]
  graph <- .mpt_tree_graph(tree)
  texts <- list()
  titles <- list()
  local_mocked_bindings(
    text = function(x, y = NULL, labels = y, ...) {
      texts[[length(texts) + 1]] <<- list(labels = labels, cex = list(...)$cex)
    },
    title = function(...) titles[[length(titles) + 1]] <<- list(...),
    .package = "graphics"
  )
  plot(tree, cex = 1.7)
  expect_length(texts, 2)
  expect_true(all(vapply(texts, function(call) identical(call$cex, 1.7), logical(1))))
  expect_setequal(unlist(lapply(texts, `[[`, "labels")), c(graph$edges$label, "old", "new"))
  expect_equal(sum(lengths(lapply(texts, `[[`, "labels"))), nrow(graph$edges) + 3)
  expect_length(titles, 1)
  expect_equal(titles[[1]]$main, "MPT tree 'old'")
})

test_that("a leaf label wider than the plot still leaves a valid window", {
  tree <- mpt_tree("long", list(
    Correct_Correct_Long_Category_Name = "p", Other = "1 - p"
  ))
  usr <- NULL
  label_positions <- .mpt_edge_label_positions
  local_mocked_bindings(.mpt_edge_label_positions = function(edges, cex) {
    usr <<- graphics::par("usr")
    label_positions(edges, cex)
  })
  withr::with_pdf(NULL, plot(tree), width = 1.5, height = 1.5)
  expect_true(all(is.finite(usr)))
  expect_gt(usr[2], usr[1])
})

drawn_labels <- function(x) {
  texts <- character(0)
  local_mocked_bindings(
    text = function(x, y = NULL, labels = y, ...) texts <<- c(texts, labels),
    .package = "graphics"
  )
  withr::with_pdf(NULL, expect_no_warning(plot(x)))
  texts
}

alias_trees <- function() {
  list(
    mpt_tree("old", list(old = "Do + (1 - Do) * g", new = "(1 - Do) * (1 - g)")),
    mpt_tree("new", list(old = "(1 - Dn) * g", new = "Dn + (1 - Dn) * (1 - g)"))
  )
}

test_that("plot() of a model labels the edges of a fixed parameter with the restriction", {
  fixed <- mpt(mpt_2htm_trees(), tree_id = "item_type", restrictions = "g = 0.5")
  labels <- drawn_labels(fixed)
  expect_contains(labels, c("D", "1 - D", "g = 0.5", "1 - g (g = 0.5)"))
  expect_false("0.5" %in% labels)
})

test_that("plot() of a model labels the edges of an equated parameter with the restriction", {
  aliased <- mpt(alias_trees(), tree_id = "item_type", restrictions = "Dn = Do")
  labels <- drawn_labels(aliased)
  expect_contains(labels, c("Do", "1 - Do", "Dn (= Do)", "1 - Dn (Dn = Do)", "g", "1 - g"))
  expect_contains(drawn_labels(alias_trees()[[2]]), c("Dn", "1 - Dn"))
})

test_that("a restricted constant keeps the multiplicity edge and the sibling sums intact", {
  restricted <- mpt(pairs_tree(), restrictions = "u = 0.5")
  expect_contains(drawn_labels(restricted), c(
    "u^2 (u = 0.5)", "2 * u * (1 - u) (u = 0.5)", "(1 - u)^2 (u = 0.5)", "c", "1 - c"
  ))
  # the fan a, 1 - b sums to 1 only once b = a
  equated <- mpt(mpt_tree("t", list(A = "a", B = "1 - b")), restrictions = "b = a")
  expect_contains(drawn_labels(equated), c("a", "1 - b (b = a)"))

  # a chain resolves in either order: b = c with c = 0.4 makes the fan c, 1 - b sum to 1
  chain_tree <- mpt_tree("t", list(A = "a * c", B = "a * (1 - b)", C = "1 - a"))
  for (order in list(c("b = c", "c = 0.4"), c("c = 0.4", "b = c"))) {
    expect_contains(drawn_labels(mpt(chain_tree, restrictions = order)), "c = 0.4")
  }
})
