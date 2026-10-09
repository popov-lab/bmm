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

test_that("simplex sticks get the beta priors of a symmetric Dirichlet(1)", {
  tree <- mpt_tree("t", list(A = "gA", B = "gB", C = "gC", D = "gD"))
  dat <- data.frame(cond = factor(c("x", "y")), A = 5L, B = 3L, C = 2L, D = 1L)
  formula <- bmf(gA ~ 1, gB ~ 1, gC ~ 1)
  for (link in c("logit", "probit")) {
    model <- mpt(tree, simplex = c("gA", "gB", "gC", "gD"), links = link)
    code <- suppressMessages(stancode(formula, dat, model))
    expect_match(code, paste0("real ", link, "beta_lpdf(real x"), fixed = TRUE)
    expect_match(code, paste0("real ", link, "beta_rng("), fixed = TRUE)
    expect_match(code, paste0(link, "beta_lpdf(b_gAraw[1] | 1, 3)"), fixed = TRUE)
    expect_match(code, paste0(link, "beta_lpdf(b_gBraw[1] | 1, 2)"), fixed = TRUE)
    expect_match(code, "_lpdf(b_gCraw[1] | 0, 1)", fixed = TRUE)
  }
  # under cell-means coding every cell of a stick gets its beta prior
  code <- suppressMessages(stancode(
    bmf(gA ~ 0 + cond, gB ~ 0 + cond, gC ~ 1), dat,
    mpt(tree, simplex = c("gA", "gB", "gC", "gD"))
  ))
  expect_match(code, "logitbeta_lpdf(b_gAraw | 1, 3)", fixed = TRUE)
})

test_that("models without a simplex group carry no simplex prior functions", {
  pcm_tree <- mpt_tree("study", list(
    C = "cp + (1 - cp) * rp * rp",
    E = "2 * (1 - cp) * rp * (1 - rp)",
    U = "(1 - cp) * (1 - rp) * (1 - rp)"
  ))
  pcm_dat <- data.frame(id = factor(1:2), C = 20L, E = 12L, U = 8L)
  codes <- list(
    suppressMessages(stancode(
      bmf(D ~ 1 + (1 | id), g ~ 1), mpt_2htm_data(n_id = 2),
      mpt(mpt_2htm_trees(), tree_id = "item_type")
    )),
    suppressMessages(stancode(bmf(cp ~ 1 + (1 | id), rp ~ 1), pcm_dat, mpt(pcm_tree)))
  )
  for (code in codes) {
    expect_no_match(code, "beta_lpdf|beta_rng")
  }
})
