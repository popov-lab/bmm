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
