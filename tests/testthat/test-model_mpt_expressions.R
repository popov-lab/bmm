test_that("mpt_tree folds integer fractions into decimal literals", {
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1/4)",
    b = "(1 - p) * (3/4)"
  ))
  expect_equal(deparse1(tree$branches$a), "p + (1 - p) * 0.25")
  expect_equal(deparse1(tree$branches$b), "(1 - p) * 0.75")
})

test_that("mpt_tree stores branch expressions as parsed calls", {
  tree <- mpt_tree("t", list(a = "D + (1 - D) * g", b = "(1 - D) * (1 - g)"))
  expect_identical(tree$branches$a, quote(D + (1 - D) * g))
  expect_equal(.mpt_expr_vars(tree), c("D", "g"))
  expect_equal(
    .mpt_eval_branches(tree, list(D = 0.7, g = 0.5)),
    c(a = 0.85, b = 0.15)
  )
})

test_that("mpt warns when branch probabilities do not sum to 1", {
  bad_tree <- mpt_tree("t", list(a = "D * g", b = "(1 - D) * g"))
  expect_warning(mpt(bad_tree), "sum to")
})

test_that("constants that Stan would receive in scientific notation error", {
  # brms deparses formula constants into the Stan code and spaces out
  # operators, so 6.7e-05 would become the subtraction '6.7e - 05'
  expect_error(
    mpt_tree("t", list(
      a = "p + (1 - p) * (1 - 0.001/15.001)",
      b = "(1 - p) * (0.001/15.001)"
    )),
    "scientific"
  )
  # constants that deparse in fixed notation are fine
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1 - 0.001)",
    b = "(1 - p) * 0.001"
  ))
  expect_equal(deparse1(tree$branches$b), "(1 - p) * 0.001")
})
