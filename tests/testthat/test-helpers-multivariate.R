test_that("add_subset() attaches subset() with the correct separator", {
  bare <- add_subset(brms::bf(y ~ x), "s1")
  expect_equal(deparse1(bare$formula), "y | subset(s1) ~ x")

  with_aterms <- add_subset(brms::bf(rt | dec(resp) ~ x), "s2")
  expect_equal(deparse1(with_aterms$formula), "rt | dec(resp) + subset(s2) ~ x")
})

test_that("add_subset() preserves formula attributes and environment", {
  env <- new.env()
  f <- brms::bf(y ~ exp(a) + b)
  environment(f$formula) <- env
  attr(f$formula, "nl") <- TRUE

  out <- add_subset(f, "s1")
  expect_true(attr(out$formula, "nl"))
  expect_identical(environment(out$formula), env)
})

test_that("sanitize_resp_name() mirrors the brms resp suffixes", {
  expect_equal(sanitize_resp_name("error"), "error")
  expect_equal(sanitize_resp_name("mean_rt"), "meanrt")
  expect_equal(sanitize_resp_name("y.2"), "y2")
})

test_that("stack_component_data() fills missing columns with typed NAs", {
  d1 <- data.frame(id = factor(c("a", "b")), error = c(0.1, 0.2))
  d2 <- data.frame(id = factor(c("b", "c")), rt = c(1.5, 2.5))

  stacked <- stack_component_data(list(d1, d2), c(".subset_error", ".subset_rt"))

  expect_equal(nrow(stacked), 4)
  expect_equal(stacked$.subset_error, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(stacked$.subset_rt, c(FALSE, FALSE, TRUE, TRUE))
  expect_false(anyNA(stacked$.subset_error))
  expect_equal(stacked$error, c(0.1, 0.2, NA, NA))
  expect_equal(stacked$rt, c(NA, NA, 1.5, 2.5))
  expect_equal(levels(stacked$id), c("a", "b", "c"))
})

test_that("stack_component_data() preserves matrix response columns", {
  d1 <- data.frame(id = 1:2)
  d1$Y <- matrix(1:6, nrow = 2, dimnames = list(NULL, c("corr", "other", "npl")))
  d2 <- data.frame(id = 3:4, rt = c(1, 2))

  stacked <- stack_component_data(list(d1, d2), c(".subset_Y", ".subset_rt"))

  expect_true(is.matrix(stacked$Y))
  expect_equal(colnames(stacked$Y), c("corr", "other", "npl"))
  expect_equal(stacked$Y[1:2, ], d1$Y)
  expect_true(all(is.na(stacked$Y[3:4, ])))
})

test_that("stack_component_data() preserves factor columns absent in one dataset", {
  d1 <- data.frame(id = 1:2, cond = factor(c("x", "y")))
  d2 <- data.frame(id = 3:4)

  stacked <- stack_component_data(list(d1, d2), c(".subset_a", ".subset_b"))

  expect_s3_class(stacked$cond, "factor")
  expect_equal(levels(stacked$cond), c("x", "y"))
  expect_true(all(is.na(stacked$cond[3:4])))
})

mv_dat_vwm <- data.frame(
  id = factor(rep(1:8, each = 10)),
  error = runif(80, -pi, pi)
)
mv_dat_rt <- data.frame(
  id = factor(rep(1:8, each = 10)),
  rt = rlnorm(80, meanlog = -0.5, sdlog = 0.3),
  resp = rbinom(80, 1, 0.7)
)
mv_f_mixture <- bmf(thetat ~ 1 + (1 | p | id), kappa ~ 1 + (1 | p | id))
mv_f_ddm <- bmf(
  drift ~ 1 + (1 | p | id), bound ~ 1 + (1 | p | id),
  ndt ~ 1, zr ~ 1
)
mv_comp_mixture <- bmm_component(mv_f_mixture,
  model = mixture2p(resp_error = "error"), data = mv_dat_vwm
)
mv_comp_lognormal <- bmm_component(bmf(rt ~ 1 + (1 | p | id)),
  family = brms::lognormal(), data = mv_dat_rt
)
mv_comp_ddm <- bmm_component(mv_f_ddm,
  model = ddm(rt = "rt", response = "resp"), data = mv_dat_rt
)

test_that("mvbmm_config() requires at least two components", {
  expect_error(
    bmm:::mvbmm_config(bmm:::new_mvbmmformula(list(mv_comp_mixture))),
    "at least two components"
  )
})

test_that("models with Stan code outside the functions block are rejected", {
  withr::local_options(bmm.sort_data = FALSE)
  sdm_component <- list(
    formula = bmf(c ~ 1, kappa ~ 1),
    model = sdm(resp_error = "error"),
    data = mv_dat_vwm,
    prior = NULL,
    data_name = "mv_dat_vwm"
  )
  expect_error(
    bmm:::mv_configure_component(sdm_component),
    "injects Stan code into"
  )
})

test_that("response variables must be unique after brms name sanitization", {
  d1 <- data.frame(id = factor(1:10), meanrt = rnorm(10))
  d2 <- data.frame(id = factor(1:10), mean_rt = rnorm(10))
  joint <- bmm_component(bmf(meanrt ~ 1 + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(mean_rt ~ 1 + (1 | p | id)), family = gaussian(), data = d2)
  expect_error(bmm:::mvbmm_config(joint), "must be unique across the components")
})

test_that("data columns clashing with subset indicators are rejected", {
  d1 <- data.frame(id = factor(1:10), y1 = rnorm(10), .subset_y2 = TRUE)
  d2 <- data.frame(id = factor(1:10), y2 = rnorm(10))
  joint <- bmm_component(bmf(y1 ~ 1 + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | p | id)), family = gaussian(), data = d2)
  expect_error(bmm:::mvbmm_config(joint), "conflict with the subset indicators")
})

test_that("a missing shared random-effects ID triggers a warning", {
  d1 <- data.frame(id = factor(1:10), y1 = rnorm(10))
  d2 <- data.frame(id = factor(1:10), y2 = rnorm(10))
  joint <- bmm_component(bmf(y1 ~ 1 + (1 | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | id)), family = gaussian(), data = d2)
  expect_warning(bmm:::mvbmm_config(joint), "No random-effects ID is shared")

  expect_silent(bmm:::mvbmm_config(mv_comp_mixture + mv_comp_lognormal))
})

test_that("formula_re_pairs() extracts ID labels with their grouping variables", {
  pairs <- bmm:::formula_re_pairs(
    bmf(kappa ~ 1 + (1 | p | id), thetat ~ cond + (cond | q | session))
  )
  expect_equal(pairs$id, c("p", "q"))
  expect_equal(pairs$group, c("id", "session"))

  expect_null(bmm:::formula_re_pairs(bmf(kappa ~ 1 + (1 | id))))
})

test_that("components with no shared grouping values trigger a warning", {
  d1 <- data.frame(id = factor(1:10), y1 = rnorm(10))
  d2 <- data.frame(id = factor(paste0("P", 1:10)), y2 = rnorm(10))
  joint <- bmm_component(bmf(y1 ~ 1 + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | p | id)), family = gaussian(), data = d2)
  expect_warning(
    bmm:::mvbmm_config(joint),
    "values of the grouping variable 'id' in common"
  )
})

test_that("partially overlapping grouping values are reported", {
  d1 <- data.frame(id = factor(1:6), y1 = rnorm(6))
  d2 <- data.frame(id = factor(4:9), y2 = rnorm(6))
  joint <- bmm_component(bmf(y1 ~ 1 + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | p | id)), family = gaussian(), data = d2)
  expect_message(bmm:::mvbmm_config(joint), "3 of 9 unique values of 'id'")
})

test_that("factor predictors with different levels across components warn", {
  d1 <- data.frame(
    id = factor(rep(1:5, 2)), y1 = rnorm(10),
    cond = factor(rep(c("a", "b"), 5))
  )
  d2 <- data.frame(
    id = factor(rep(1:5, 2)), y2 = rnorm(10),
    cond = factor(rep(c("a", "c"), 5))
  )
  joint <- bmm_component(bmf(y1 ~ cond + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ cond + (1 | p | id)), family = gaussian(), data = d2)
  expect_warning(bmm:::mvbmm_config(joint), "different levels across")
})

test_that("component priors are tagged with the sanitized response name", {
  user_prior <- brms::set_prior("normal(1, 0.5)",
    class = "b", coef = "Intercept", nlpar = "kappa"
  )
  joint <- bmm_component(mv_f_mixture,
    model = mixture2p(resp_error = "error"), data = mv_dat_vwm, prior = user_prior
  ) + mv_comp_lognormal
  cfg <- bmm:::mvbmm_config(joint)

  kappa_row <- cfg$prior[cfg$prior$nlpar == "kappa" & cfg$prior$coef == "Intercept", ]
  expect_equal(kappa_row$prior, "normal(1, 0.5)")
  expect_equal(kappa_row$resp, "error")
  expect_true(all(cfg$prior$resp %in% c("error", "rt")))
})

test_that("a global prior is passed through without a resp tag", {
  cfg <- bmm:::mvbmm_config(
    mv_comp_mixture + mv_comp_lognormal,
    prior = brms::set_prior("lkj(2)", class = "cor")
  )
  cor_row <- cfg$prior[cfg$prior$class == "cor", ]
  expect_equal(cor_row$prior, "lkj(2)")
  expect_equal(cor_row$resp, "")
})

test_that("multivariate standata matches the univariate standata per component", {
  sd_mv <- standata(mv_comp_mixture + mv_comp_lognormal)
  sd_uni <- standata(mv_f_mixture, data = mv_dat_vwm, model = mixture2p(resp_error = "error"))

  expect_equal(sd_mv$N, nrow(mv_dat_vwm) + nrow(mv_dat_rt))
  expect_equal(sd_mv$N_error, sd_uni$N)
  expect_equal(sd_mv$Y_error, sd_uni$Y)
  expect_equal(sd_mv$X_error_kappa, sd_uni$X_kappa)
  expect_equal(sd_mv$X_kappa2_error, sd_uni$X_kappa2)
  expect_equal(sd_mv$Z_1_error_kappa_1, sd_uni$Z_1_kappa_1)
  expect_equal(sd_mv$N_rt, nrow(mv_dat_rt))
  expect_equal(as.vector(sd_mv$Y_rt), mv_dat_rt$rt)
  expect_equal(sd_mv$M_1, 3)
})

test_that("custom family variables are partitioned correctly (ddm + gaussian)", {
  joint <- suppressMessages(
    mv_comp_ddm +
      bmm_component(bmf(error ~ 1 + (1 | p | id)), family = gaussian(), data = mv_dat_vwm)
  )
  sd_mv <- suppressMessages(standata(joint))
  sd_uni <- suppressMessages(
    standata(mv_f_ddm, data = mv_dat_rt, model = ddm(rt = "rt", response = "resp"))
  )

  expect_equal(sd_mv$N_rt, sd_uni$N)
  expect_equal(sd_mv$Y_rt, sd_uni$Y)
  expect_equal(sd_mv$dec_rt, sd_uni$dec)
  expect_equal(sd_mv$X_drift_rt, sd_uni$X_drift)
  expect_equal(sd_mv$Z_1_drift_rt_1, sd_uni$Z_1_drift_1)
})

test_that("two measurement models share one correlation matrix", {
  sd_mv <- standata(mv_comp_mixture + mv_comp_ddm)
  expect_equal(sd_mv$N_error, nrow(mv_dat_vwm))
  expect_equal(sd_mv$N_rt, nrow(mv_dat_rt))
  expect_equal(sd_mv$M_1, 4)

  scode <- stancode(mv_comp_mixture + mv_comp_ddm)
  expect_match(scode, "ddm_lpdf", fixed = TRUE)
  expect_match(scode, "von_mises", fixed = TRUE)
})

test_that("the m3 matrix response survives multivariate composition", {
  m3dat <- oberauer_lewandowsky_2019_e1
  rt_dat <- data.frame(
    ID = unique(m3dat$ID),
    rt = rlnorm(length(unique(m3dat$ID)))
  )
  m3_model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "simple", version = "ss"
  )
  m3_formula <- bmf(c ~ 1 + (1 | p | ID), a ~ 1 + (1 | p | ID))
  joint <- suppressMessages(
    bmm_component(m3_formula, model = m3_model, data = m3dat) +
      bmm_component(bmf(rt ~ 1 + (1 | p | ID)), family = brms::lognormal(), data = rt_dat)
  )
  sd_mv <- suppressMessages(standata(joint))
  sd_uni <- suppressMessages(standata(m3_formula, data = m3dat, model = m3_model))

  expect_equal(sd_mv$N_Y, sd_uni$N)
  expect_equal(sd_mv$Y_Y, sd_uni$Y)
  expect_equal(sd_mv$trials_Y, sd_uni$trials)
  expect_equal(sd_mv$X_Y_c, sd_uni$X_c)
})

test_that("the same model on two tasks deduplicates its Stan functions", {
  dat_two_tasks <- data.frame(
    id = factor(rep(1:8, each = 10)),
    rt1 = rlnorm(80, -0.5, 0.3), resp1 = rbinom(80, 1, 0.7),
    rt2 = rlnorm(80, -0.3, 0.3), resp2 = rbinom(80, 1, 0.8)
  )
  joint <- bmm_component(
    bmf(drift ~ 1 + (1 | p | id), bound ~ 1 + (1 | p | id), ndt ~ 1, zr ~ 1),
    model = ddm(rt = "rt1", response = "resp1"), data = dat_two_tasks
  ) +
    bmm_component(
      bmf(drift ~ 1 + (1 | p | id), bound ~ 1 + (1 | p | id), ndt ~ 1, zr ~ 1),
      model = ddm(rt = "rt2", response = "resp2"), data = dat_two_tasks
    )

  scode <- stancode(joint)
  expect_equal(
    lengths(regmatches(scode, gregexpr("real ddm_lpdf", scode, fixed = TRUE))), 1
  )

  sd_mv <- standata(joint)
  sd_uni <- standata(
    bmf(drift ~ 1 + (1 | p | id), bound ~ 1 + (1 | p | id), ndt ~ 1, zr ~ 1),
    data = dat_two_tasks, model = ddm(rt = "rt1", response = "resp1")
  )
  expect_equal(sd_mv$Y_rt1, sd_uni$Y)
  expect_equal(sd_mv$dec_rt1, sd_uni$dec)
})

test_that("components without init_ranges produce the default init", {
  cfg <- bmm:::mvbmm_config(mv_comp_mixture + mv_comp_lognormal)
  expect_equal(cfg$config_args$init, 1)
})

test_that("components with init_ranges produce a valid init function", {
  joint <- suppressMessages(
    mv_comp_ddm +
      bmm_component(bmf(error ~ 1 + (1 | p | id)), family = gaussian(), data = mv_dat_vwm)
  )
  cfg <- suppressMessages(bmm:::mvbmm_config(joint))
  initfun <- cfg$config_args$init
  expect_true(is.function(initfun))

  inits <- initfun()
  scode <- suppressMessages(stancode(joint))
  spars <- names(bmm:::extract_parameter_dimensions(
    bmm:::extract_stan_blocks(scode)$parameters
  ))
  expect_true(all(names(inits) %in% spars))

  ranges <- bmm:::.model_ddm()$init_ranges
  expect_gte(inits$Intercept_ndt_rt, log(ranges$ndt[1]))
  expect_lte(inits$Intercept_ndt_rt, log(ranges$ndt[2]))
  expect_gte(inits$Intercept_drift_rt, ranges$drift[1])
  expect_lte(inits$Intercept_drift_rt, ranges$drift[2])
  expect_gte(inits$Intercept_bound_rt, log(ranges$bound[1]))
  expect_lte(inits$Intercept_bound_rt, log(ranges$bound[2]))

  expect_equal(dim(inits$z_1), c(3, length(unique(mv_dat_rt$id))))
  expect_equal(inits$L_1, diag(3))
  expect_length(inits$sd_1, 3)

  expect_false("Intercept_error" %in% names(inits))
  expect_false("sigma_error" %in% names(inits))
})

test_that("init matching distinguishes the same model on two tasks", {
  dat_two_tasks <- data.frame(
    id = factor(rep(1:8, each = 10)),
    rt1 = rlnorm(80, -0.5, 0.3), resp1 = rbinom(80, 1, 0.7),
    rt2 = rlnorm(80, -0.3, 0.3), resp2 = rbinom(80, 1, 0.8)
  )
  f_ddm <- bmf(drift ~ 1 + (1 | p | id), bound ~ 1 + (1 | p | id), ndt ~ 1, zr ~ 1)
  joint <- bmm_component(f_ddm, model = ddm(rt = "rt1", response = "resp1"), data = dat_two_tasks) +
    bmm_component(f_ddm, model = ddm(rt = "rt2", response = "resp2"), data = dat_two_tasks)
  cfg <- bmm:::mvbmm_config(joint)

  inits <- cfg$config_args$init()
  expect_true(all(
    c("Intercept_ndt_rt1", "Intercept_ndt_rt2", "Intercept_drift_rt1", "Intercept_drift_rt2") %in%
      names(inits)
  ))
})

test_that("bmm() fits multivariate specifications and postprocesses the fit", {
  joint <- mv_comp_mixture + mv_comp_lognormal
  fit <- bmm(joint, backend = "mock", mock_fit = 1, rename = FALSE)

  expect_s3_class(fit, c("mvbmmfit", "bmmfit", "brmsfit"), exact = TRUE)
  expect_length(fit$bmm$components, 2)
  expect_named(
    fit$bmm$components[[1]],
    c("model", "user_formula", "resp", "resp_name", "subset_var", "data_name")
  )
  expect_equal(
    vapply(fit$bmm$components, function(x) x$resp_name, character(1)),
    c("error", "rt")
  )
  expect_s3_class(fit$bmm$components[[1]]$model, "mixture2p")
  expect_s3_class(fit$bmm$components[[1]]$user_formula, "bmmformula")
  expect_equal(attr(fit$data, "data_name"), "mv_dat_vwm + mv_dat_rt")
  expect_equal(nrow(fit$data), nrow(mv_dat_vwm) + nrow(mv_dat_rt))
  expect_equal(fit$version$bmm, utils::packageVersion("bmm"))
})

test_that("bmm() rejects invalid multivariate calls", {
  joint <- mv_comp_mixture + mv_comp_lognormal
  expect_error(
    bmm(joint, data = mv_dat_vwm, backend = "mock", mock_fit = 1),
    "Omit the 'data' and 'model' arguments"
  )
  expect_error(
    bmm(joint, model = mixture2p(resp_error = "error"), backend = "mock", mock_fit = 1),
    "Omit the 'data' and 'model' arguments"
  )
  expect_error(
    bmm(mv_comp_mixture, backend = "mock", mock_fit = 1),
    "single bmm_component"
  )
})

test_that("multivariate fits are cached with the file argument", {
  file <- tempfile(fileext = ".rds")
  joint <- mv_comp_mixture + mv_comp_lognormal
  fit <- bmm(joint, backend = "mock", mock_fit = 1, rename = FALSE, file = file)
  expect_true(file.exists(file))

  cached <- bmm(joint, backend = "mock", mock_fit = 2, rename = FALSE, file = file)
  expect_equal(cached$fit, fit$fit)
  unlink(file)
})

test_that("merge_stanvars() deduplicates identical Stan code", {
  sv <- brms::stanvar(scode = "real foo_lpdf(real y) { return y; }", block = "functions")
  sv_other <- brms::stanvar(scode = "real bar_lpdf(real y) { return y; }", block = "functions")

  expect_null(merge_stanvars(list(NULL, NULL)))
  expect_length(merge_stanvars(list(sv, NULL)), 1)
  expect_length(merge_stanvars(list(sv, sv)), 1)

  combined <- merge_stanvars(list(sv, sv_other, sv))
  expect_length(combined, 2)
  expect_s3_class(combined, "stanvars")
})
