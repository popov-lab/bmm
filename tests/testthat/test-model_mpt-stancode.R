test_that("decimal constants are emitted literally, never in scientific notation", {
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1 - 0.001)",
    b = "(1 - p) * 0.001"
  ))
  dat <- data.frame(a = c(30L, 28L), b = c(0L, 2L))
  code <- stancode(bmf(p ~ 1), data = dat, model = mpt(tree))
  expect_match(code, "0.001", fixed = TRUE)
  expect_false(grepl("e-", code, fixed = TRUE))
})
