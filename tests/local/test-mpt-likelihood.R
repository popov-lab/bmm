test_that("an mpt fit's likelihood and expected counts equal dmpt() at the same draws", {
  withr::local_seed(21)
  withr::local_options(bmm.sort_data = FALSE)
  # an impossible category, so that the -100 linear predictor of bmf2bf.mpt()
  # is part of what is compared
  dat <- mpt_impossible_data(n_id = 2)
  dat$corr <- dat$corr + stats::rbinom(nrow(dat), 10, 0.5)
  model <- mpt(mpt_impossible_trees(), "tree")
  # sig_figs = 18: lp__ is compared at full precision, not at the CSV default
  fit <- bmm(
    bmf(Pm ~ 1, Pb ~ 1), dat, model, backend = "cmdstanr", chains = 2,
    iter = 1000, refresh = 0, seed = 21, sig_figs = 18
  )
  counts <- as.matrix(dat[c("corr", "dist", "npl")])
  expect_equal(unname(fit$data$Y), unname(counts))

  draws <- posterior::as_draws_df(fit)
  ids <- round(seq(1, nrow(draws), length.out = 20))
  inv_links <- list(logit = stats::plogis, probit = stats::pnorm)
  native <- function(par, i) inv_links[[model$links[[par]]]](draws[[paste0("b_", par, "_Intercept")]][i])
  rows <- seq_len(nrow(dat))
  r_loglik <- t(vapply(ids, function(i) {
    pars <- c(Pm = native("Pm", i), Pb = native("Pb", i))
    vapply(rows, function(n) dmpt(counts[n, ], pars, model, tree = dat$tree[n]), numeric(1))
  }, numeric(nrow(dat))))
  expect_lt(max(abs(brms::log_lik(fit, draw_ids = ids) - r_loglik)), 1e-8)

  # brms evaluates log_lik() in R from the formula; lp__ is Stan's own sum:
  # the multinomial log-likelihood with its normalizing constant plus the
  # logistic(0, 1) priors of the two intercepts, with no Jacobian
  r_lp <- rowSums(r_loglik) +
    stats::dlogis(draws$b_Pm_Intercept[ids], log = TRUE) +
    stats::dlogis(draws$b_Pb_Intercept[ids], log = TRUE)
  expect_lt(max(abs(draws$lp__[ids] - r_lp)), 1e-8)

  # one response in a category has the probability of that category
  unit <- diag(3)
  dimnames(unit) <- list(NULL, colnames(counts))
  r_epred <- aperm(vapply(ids, function(i) {
    pars <- c(Pm = native("Pm", i), Pb = native("Pb", i))
    t(vapply(rows, function(n) {
      rowSums(counts)[n] * exp(apply(unit, 1, dmpt, pars = pars, mpt_model = model, tree = dat$tree[n]))
    }, numeric(3)))
  }, matrix(0, nrow(dat), 3)), c(3, 1, 2))
  epred <- brms::posterior_epred(fit, draw_ids = ids)
  expect_equal(apply(epred, c(1, 2), sum), matrix(rowSums(counts), length(ids), nrow(dat), byrow = TRUE), ignore_attr = TRUE)
  expect_lt(max(abs(unname(epred) - r_epred)), 1e-8)
})
