test_that("mpt_tree folds integer fractions into decimal literals", {
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1/4)",
    b = "(1 - p) * (3/4)"
  ))
  expect_false(grepl("/", tree$branches$a, fixed = TRUE))
  expect_true(grepl("0.25", tree$branches$a, fixed = TRUE))
  expect_true(grepl("0.75", tree$branches$b, fixed = TRUE))
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
  expect_equal(tree$branches$b, "(1 - p) * 0.001")
})
