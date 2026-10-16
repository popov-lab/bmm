test_that("ppp() does not flag a well-specified mpt model and flags a misspecified one", {
  withr::local_seed(2)
  withr::local_options("bmm.silent" = 2)
  n_id <- 15
  dat <- expand.grid(id = factor(seq_len(n_id)), item_type = c("old", "new"),
                     stringsAsFactors = FALSE)
  D <- stats::plogis(stats::qlogis(0.6) + stats::rnorm(n_id, 0, 0.5))[dat$id]
  p_old <- ifelse(dat$item_type == "old", D + (1 - D) * 0.5, (1 - D) * 0.5)
  dat$old <- stats::rbinom(nrow(dat), 60, p_old)
  dat$new <- 60L - dat$old

  fit <- bmm(
    bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id)), dat, mpt(mpt_2htm_trees(), "item_type"),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, seed = 2
  )
  res <- ppp(fit)
  expect_equal(attr(res, "by"), "item_type")
  expect_true(all(res$ppp > 0.05))
  expect_equal(nrow(ppp(fit, level = "participant")), n_id)

  misfit <- bmm(
    bmf(D ~ 1, g = 0.9), dat, mpt(mpt_2htm_trees(), "item_type"),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, seed = 2
  )
  withr::local_options("bmm.silent" = 1)
  expect_message(res_misfit <- ppp(misfit), "computes T1 only")
  expect_lt(res_misfit$ppp, 0.01)
})

test_that("ppp() handles the impossible categories of a fitted mpt model", {
  withr::local_seed(3)
  withr::local_options("bmm.silent" = 2)
  dat <- mpt_impossible_data()
  dat$corr <- dat$corr + stats::rbinom(nrow(dat), 5, 0.5)

  fit <- bmm(
    bmf(Pm ~ 1 + (1 | id), Pb ~ 1), dat, mpt(mpt_impossible_trees(), "tree"),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, seed = 3
  )
  res <- ppp(fit)
  expect_equal(attr(res, "by"), "tree")
  expect_true(all(is.finite(c(res$observed, res$replicated, res$ppp))))
})
