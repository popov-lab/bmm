test_that("impossible categories reach the Stan code as a -100 predictor", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  formula <- bmf(Pm ~ 1, Pb ~ 1)

  # brms replaces the indicator columns by C_<parameter>_<i>[n] in the code
  code <- stancode(formula, data = dat, model = model)
  expect_match(
    code,
    "mudist[n] = (C_mudist_1[n] * log(nlp_dist[n]) + (1 - C_mudist_1[n]) * ( - 100))",
    fixed = TRUE
  )
  # the placeholder keeps log() defined for the tree without the category
  expect_match(code, "nlp_dist[n] = (C_dist_1[n] * ((1 - inv_logit(nlp_Pm[n])) * 0.2) + C_dist_2[n] * (1))", fixed = TRUE)
  expect_match(code, "mucorr[n] = (log(nlp_corr[n]))", fixed = TRUE)
  expect_match(code, "multinomial_logit2_lpmf(Y[n] | mu[n])", fixed = TRUE)

  # brms passes the generated indicator columns as non-linear covariates,
  # named per category formula in writing order: Idx_withdist, Idx_nodist
  # for each category and Poss_dist for the guarded linear predictor
  sdata <- standata(formula, data = dat, model = model)
  expect_equal(sdata$N, nrow(dat))
  expect_equal(as.integer(sdata$C_corr_1), as.integer(dat$tree == "withdist"))
  expect_equal(as.integer(sdata$C_corr_2), as.integer(dat$tree == "nodist"))
  expect_equal(as.integer(sdata$C_mudist_1), as.integer(dat$cond == "withdist"))
})

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
