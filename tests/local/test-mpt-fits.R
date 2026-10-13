test_that("conditional_effects() shows mpt parameters on the native scale", {
  # unequal trials per row, with a mean that is not an integer
  dat <- data.frame(
    item_type = rep(c("old", "new"), 2), cond = rep(c("x", "y"), each = 2),
    old = c(70, 20, 45, 30), new = c(30, 61, 15, 70)
  )
  fit <- bmm(
    bmf(D ~ cond, g ~ 1), dat, mpt(mpt_2htm_trees(), "item_type"),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, silent = 2
  )

  ce <- conditional_effects(fit, par = "D", robust = TRUE)
  np <- native_parameters(fit, pars = "D")
  medians <- c(tapply(np$value, np$cond, stats::median))
  expect_equal(ce$cond$estimate__, unname(medians[as.character(ce$cond$cond)]), tolerance = 1e-3)
  expect_named(conditional_effects(fit), "D.cond")
})

test_that("the default simplex prior gives every member prior mean 1/K under both links", {
  withr::local_options("bmm.silent" = 2)
  tree <- mpt_tree("t", list(A = "gA", B = "gB", C = "gC"))
  dat <- data.frame(A = 4L, B = 3L, C = 3L)
  inv_links <- list(logit = stats::plogis, probit = stats::pnorm)
  for (link in names(inv_links)) {
    fit <- bmm(
      bmf(gA ~ 1, gB ~ 1), dat, mpt(tree, simplex = c("gA", "gB", "gC"), links = link),
      sample_prior = "only", chains = 2, iter = 3000, warmup = 500, refresh = 0,
      seed = 1, rename = FALSE
    )
    draws <- brms::as_draws_df(fit)
    stick_a <- inv_links[[link]](draws[[grep("^b_gAraw", names(draws))]])
    stick_b <- inv_links[[link]](draws[[grep("^b_gBraw", names(draws))]])
    # the per-stick uniform prior this replaces gives 0.50/0.25/0.25; the
    # Monte Carlo error of each mean is below 0.01 at these draws
    expect_equal(
      c(mean(stick_a), mean((1 - stick_a) * stick_b), mean((1 - stick_a) * (1 - stick_b))),
      rep(1 / 3, 3),
      tolerance = 0.03
    )
  }
})

test_that("conditional_effects() shows simplex members on the probability scale", {
  trees <- list(
    mpt_tree("x", list(A = "D + (1 - D) * a", B = "(1 - D) * b", C = "(1 - D) * c")),
    mpt_tree("y", list(A = "(1 - D) * a", B = "D + (1 - D) * b", C = "(1 - D) * c"))
  )
  dat <- data.frame(
    tt = rep(c("x", "y"), 2), cond = rep(c("p", "q"), each = 2),
    A = c(60, 15, 50, 30), B = c(25, 55, 20, 40), C = c(15, 21, 30, 35)
  )
  fit <- bmm(
    bmf(D ~ 1, a ~ cond, b ~ cond), dat, mpt(trees, "tt", simplex = c("a", "b", "c")),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, silent = 2
  )

  np <- native_parameters(fit, pars = c("a", "b"))
  # b = (1 - a) * stick, so its panel depends on both stick-breaking components
  for (p in c("a", "b")) {
    ce <- conditional_effects(fit, par = p, robust = TRUE)$cond
    medians <- c(tapply(np$value[np$parameter == p], np$cond[np$parameter == p], stats::median))
    expect_true(all(ce$lower__ > 0 & ce$upper__ < 1))
    expect_equal(ce$estimate__, unname(medians[as.character(ce$cond)]), tolerance = 1e-3)
  }
})
