# conditional_effects() on real fits (moved from tests/testthat in #422)

skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")

test_that("conditional_effects works on sdt_rating fits (#512)", {
  withr::local_seed(512)

  dat <- expand.grid(id = 1:6, stimulus = c(0L, 1L), cond = c("A", "B"))
  dat <- cbind(dat, rsdt_rating(nrow(dat), 150, dat$stimulus,
                                d = ifelse(dat$cond == "A", 1, 2),
                                thresholds = c(-0.5, 0, 0.5)))
  dat$x <- stats::rnorm(nrow(dat))
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", links = list(d = "log"))
  fit <- bmm(bmf(d ~ 1 + cond + x, criterion ~ 1, spacing ~ 1), dat, model,
             backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  expect_no_message(ce <- conditional_effects(fit))
  expect_setequal(names(ce), c("d.cond", "d.x"))

  ce <- conditional_effects(
    fit, par = "d", effects = "x:cond", int_conditions = list(cond = "B")
  )
  expect_true(all(ce[[1]]$cond == "B"))

  draw_ids <- 1:11
  ce <- conditional_effects(fit, par = "d", effects = "x", spaghetti = TRUE,
                            draw_ids = draw_ids, scale = "sampling")
  draws <- brms::posterior_linpred(fit, newdata = ce$x, nlpar = "d",
                                   re_formula = NA, draw_ids = draw_ids)
  expect_equal(ce$x$estimate__, apply(draws, 2, median))
  expect_equal(attr(ce$x, "spaghetti")$estimate__, as.numeric(t(draws)))

  ce <- conditional_effects(fit, par = "d", effects = "x", spaghetti = TRUE,
                            draw_ids = draw_ids)
  expect_equal(ce$x$estimate__, exp(apply(draws, 2, median)))
  expect_equal(attr(ce$x, "spaghetti")$estimate__, exp(as.numeric(t(draws))))
})

test_that("m3 conditional_effects does not evaluate the category dpar on other formulas' effects (#512)", {
  # n_dist is 5, 0, 0 across conditions and enters the first category dpar
  # through log(), so brms's grid of Idx_dist:n_dist (n_dist at mean - sd < 0)
  # makes it NaN
  formula <- bmf(corr ~ b + a + c, dist ~ b + d, other ~ b + a, npl ~ b,
                 c ~ 0 + cond, a ~ 1, d ~ 1, b = 0)
  model <- m3(resp_cats = c("dist", "corr", "other", "npl"),
              num_options = c("n_dist", "n_corr", "n_other", "n_npl"),
              version = "custom", links = list(a = "log", c = "log", d = "log"))
  fit <- suppressWarnings(suppressMessages(bmm(
    formula, oberauer_lewandowsky_2019_e1, model, backend = "cmdstanr",
    chains = 1, iter = 300, warmup = 150, refresh = 0, silent = 2, seed = 1
  )))

  expect_no_warning(ce_c <- conditional_effects(fit, par = "c"))
  expect_no_warning(ce_a <- conditional_effects(fit, par = "a"))
  expect_no_warning(ce_all <- conditional_effects(fit))
  expect_equal(ce_c, conditional_effects(fit, par = "c", effects = "cond"))
  expect_length(ce_a, 0)
  expect_named(ce_all, "c.cond")
  expect_equal(
    .ce_brms_default_effects(fit, NA),
    unname(lapply(suppressWarnings(suppressMessages(.brms_conditional_effects(
      fit, dpar = names(brms::brmsterms(fit$formula)$dpars)[1], draw_ids = 1
    ))), attr, "effects"))
  )
})
