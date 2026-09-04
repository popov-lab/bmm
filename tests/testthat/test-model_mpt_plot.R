test_that("branch expressions expand into root-to-leaf paths", {
  paths <- .mpt_branch_paths(quote(D + (1 - D) * g))
  expect_equal(paths, list("D", c("1 - D", "g")))

  # factored subtrees are expanded distributively
  paths_factored <- .mpt_branch_paths(quote(Pm * (Pb + (1 - Pb) * 0.25)))
  expect_equal(paths_factored, list(c("Pm", "Pb"), c("Pm", "1 - Pb", "0.25")))
})

test_that("the tree graph merges shared path prefixes only", {
  tree <- mpt_tree("old", list(
    old = "D + (1 - D) * g",
    new = "(1 - D) * (1 - g)"
  ))
  graph <- .mpt_tree_graph(tree)
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

test_that("plot methods run for trees and models", {
  tree_old <- mpt_tree("old", list(
    old = "D + (1 - D) * g",
    new = "(1 - D) * (1 - g)"
  ))
  model <- mpt(
    list(tree_old, mpt_tree("new", list(
      old = "(1 - D) * g",
      new = "D + (1 - D) * (1 - g)"
    ))),
    tree_id = "item_type"
  )
  grDevices::pdf(NULL)
  expect_silent(plot(tree_old))
  expect_silent(plot(model))
  grDevices::dev.off()
})
