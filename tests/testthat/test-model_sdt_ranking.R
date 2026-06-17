# Test Ranking SDT model specification and integration (native multinomial)

ranks4 <- paste0("rank", 1:4)
ranks5 <- paste0("rank", 1:5)

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_ranking model can be created", {
  model <- sdt_ranking(ranks4, m = 4)
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_ranking")
})

test_that("sdt_ranking has correct parameters (only dprime, no criterion, no mu)", {
  model <- sdt_ranking(ranks4, m = 4)
  expect_true("dprime" %in% names(model$parameters))
  expect_false("mu" %in% names(model$parameters))
  expect_false("criterion" %in% names(model$parameters))
  expect_length(model$fixed_parameters, 0)
})

test_that("sdt_ranking stores response and m correctly", {
  model <- sdt_ranking(ranks4, m = 4)
  expect_equal(model$resp_vars$response, ranks4)
  expect_equal(model$other_vars$m, 4L)
})

test_that("sdt_ranking accepts m as a constant or a column name", {
  m_const <- sdt_ranking(ranks4, m = 4)
  expect_equal(m_const$other_vars$m, 4L)

  m_col <- sdt_ranking(ranks5, m = "set_size")
  expect_identical(m_col$other_vars$m, "set_size")
})

test_that("sdt_ranking has correct links and accepts custom links", {
  model <- sdt_ranking(ranks4, m = 4)
  expect_equal(model$links$dprime, "identity")

  custom <- sdt_ranking(ranks4, m = 4, links = list(dprime = "log"))
  expect_equal(custom$links$dprime, "log")
})

test_that("sdt_ranking requires a valid m", {
  expect_error(sdt_ranking(ranks4, m = 1), "m must be")
  expect_error(sdt_ranking(ranks4, m = c(3, 4)), "m must be")
})

test_that("sdt_ranking requires response and m", {
  expect_error(sdt_ranking(m = 4))
  expect_error(sdt_ranking(ranks4))
})

test_that("constant m must match the number of response columns", {
  expect_error(sdt_ranking(ranks4, m = 5), "exactly m columns")
  expect_error(sdt_ranking("rank1", m = 4), "at least 2")
})

test_that("sdt_ranking only supports gumbel_min and normal", {
  expect_silent(sdt_ranking(ranks4, m = 4, dist = "gumbel_min"))
  expect_silent(sdt_ranking(ranks4, m = 4, dist = "normal"))
  expect_error(sdt_ranking(ranks4, m = 4, dist = "logistic"))
})

test_that("dist='normal' adds sdratio (fixed to 0); gumbel_min does not", {
  normal <- sdt_ranking(ranks4, m = 4, dist = "normal")
  expect_true("sdratio" %in% names(normal$parameters))
  expect_equal(normal$fixed_parameters$sdratio, 0)
  expect_equal(normal$init_ranges$sdratio, c(-0.3, 0.3))

  gumbel <- sdt_ranking(ranks4, m = 4, dist = "gumbel_min")
  expect_false("sdratio" %in% names(gumbel$parameters))
})

test_that("sdt_ranking has init_ranges with dprime (no mu)", {
  model <- sdt_ranking(ranks4, m = 4)
  expect_equal(model$init_ranges$dprime, c(0.5, 1.5))
  expect_false("mu" %in% names(model$init_ranges))
})


############################################################################# !
# R-SIDE PROBABILITY TESTS                                               ####
############################################################################# !

test_that(".ranking_all_probs_r sums to 1 and is monotone (Gumbel)", {
  probs <- bmm:::.ranking_all_probs_r(1.5, 4L, dist = "gumbel_min")
  expect_equal(sum(probs), 1.0, tolerance = 1e-10)

  p_low <- bmm:::.ranking_all_probs_r(0.5, 4L, dist = "gumbel_min")
  p_high <- bmm:::.ranking_all_probs_r(2.0, 4L, dist = "gumbel_min")
  expect_true(p_high[1] > p_low[1])
})

test_that(".ranking_all_probs_r at d'=0 gives uniform for Gumbel", {
  probs <- bmm:::.ranking_all_probs_r(0, 4L, dist = "gumbel_min")
  expect_equal(probs, rep(0.25, 4), tolerance = 1e-6)
})

test_that(".ranking_all_probs_r works for m=3 and for Gaussian UV-SDT", {
  probs3 <- bmm:::.ranking_all_probs_r(1.0, 3L, dist = "gumbel_min")
  expect_length(probs3, 3)
  expect_equal(sum(probs3), 1.0, tolerance = 1e-10)

  probs_n <- bmm:::.ranking_all_probs_r(1.5, 4L, dist = "normal", sdratio = 0)
  expect_equal(sum(probs_n), 1.0, tolerance = 1e-6)
  expect_true(probs_n[1] > probs_n[4])
})

test_that("sdt_ranking_logmu R companion matches .ranking_all_probs_r and masks", {
  for (info in list(list(d = "gumbel_min", id = 2L, s = 0),
                    list(d = "normal", id = 1L, s = 0.2))) {
    internal <- bmm:::.ranking_all_probs_r(1.3, 4L, info$d, info$s)
    comp <- vapply(1:4, function(k) {
      exp(sdt_ranking_logmu(k, 4, matrix(1.3), sdratio = info$s, dist = info$id))
    }, numeric(1))
    expect_equal(comp, internal, tolerance = 1e-10)
    # ranks beyond the set size switch off to the -100 sentinel
    expect_equal(as.numeric(sdt_ranking_logmu(5L, 4, matrix(1.3),
                            sdratio = info$s, dist = info$id)), -100)
  }
})

test_that("sdt_ranking_logmu preserves the draws-by-observation shape", {
  dprime <- matrix(c(1, 1.5, 0.5, 2), nrow = 2)
  out <- sdt_ranking_logmu(1L, c(4, 4), dprime, dist = 2L)
  expect_equal(dim(out), dim(dprime))
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_ranking generates valid wide ranking data", {
  dat <- rsdt_ranking(n_per_cell = 100, n_subjects = 5, dprime = 1.5, m = 4)

  expect_equal(colnames(dat), c("id", "set_size", paste0("rank", 1:4)))
  expect_equal(nrow(dat), 5)
  expect_true(all(dat$set_size == 4))
  expect_true(all(rowSums(dat[paste0("rank", 1:4)]) == 100))
})

test_that("rsdt_ranking validates inputs", {
  expect_error(rsdt_ranking(n_per_cell = 50, n_subjects = 1, dprime = 1, m = NULL))
  expect_error(rsdt_ranking(n_per_cell = 50, n_subjects = 1, dprime = 1, m = 1))
})

test_that("dsdt_ranking computes valid ranking densities", {
  dens <- dsdt_ranking(counts = c(40, 30, 20, 10), dprime = 1.5, m = 4)
  expect_true(dens > 0)

  log_dens <- dsdt_ranking(counts = c(40, 30, 20, 10), dprime = 1.5, m = 4,
                           log = TRUE)
  expect_equal(log(dens), log_dens, tolerance = 1e-10)
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("check_data builds the multinomial response matrix and covariates", {
  model <- sdt_ranking(ranks4, m = 4)
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 4)
  result <- check_data(model, dat, bmf(dprime ~ 1))

  expect_true(all(c("Y", "nTrials", "max_rank") %in% colnames(result)))
  expect_equal(dim(result$Y), c(3L, 4L))
  expect_equal(result$nTrials, rowSums(result$Y))
  expect_equal(unique(result$max_rank), 4)
  expect_false(any(ranks4 %in% colnames(result)))
})

test_that("check_data errors on missing response columns", {
  model <- sdt_ranking(ranks4, m = 4)
  dat <- data.frame(id = 1:5, rank1 = 10, rank2 = 10)
  expect_error(check_data(model, dat, bmf(dprime ~ 1)), "missing")
})

test_that("check_data resolves per-row set size from a data column", {
  model <- sdt_ranking(ranks5, m = "set_size")
  dat <- rbind(
    data.frame(id = 1, set_size = 3L,
               rank1 = 30, rank2 = 15, rank3 = 5, rank4 = 0, rank5 = 0),
    data.frame(id = 1, set_size = 4L,
               rank1 = 30, rank2 = 12, rank3 = 6, rank4 = 2, rank5 = 0)
  )
  result <- check_data(model, dat, bmf(dprime ~ 1))
  expect_equal(result$max_rank, c(3, 4))
})

test_that("check_data errors when the set-size column is missing", {
  model <- sdt_ranking(ranks5, m = "set_size")
  dat <- data.frame(id = 1, rank1 = 30, rank2 = 12, rank3 = 6, rank4 = 2, rank5 = 0)
  expect_error(check_data(model, dat, bmf(dprime ~ 1)), "Set-size column")
})

test_that("check_data errors on counts beyond the per-row set size", {
  model <- sdt_ranking(ranks5, m = "set_size")
  # set size 3 but a nonzero count in rank4 (a structural-zero column)
  dat <- data.frame(id = 1, set_size = 3L,
                    rank1 = 30, rank2 = 15, rank3 = 5, rank4 = 6, rank5 = 0)
  expect_error(check_data(model, dat, bmf(dprime ~ 1)),
               "beyond the row's set size")
})

test_that("check_data errors when set size exceeds the number of rank columns", {
  model <- sdt_ranking(ranks4, m = "set_size")
  dat <- data.frame(id = 1, set_size = 5L,
                    rank1 = 10, rank2 = 10, rank3 = 10, rank4 = 10)
  expect_error(check_data(model, dat, bmf(dprime ~ 1)),
               "number of rank columns")
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("bmf2bf.sdt_ranking builds a multinomial non-linear formula", {
  model <- sdt_ranking(ranks4, m = 4)
  bf <- bmf2bf(model, bmf(dprime ~ 1))

  expect_s3_class(bf, "brmsformula")
  bf_str <- deparse(bf$formula, width.cutoff = 500)
  expect_true(grepl("Y | trials(nTrials)", bf_str, fixed = TRUE))
  expect_true(grepl("sdt_ranking_logmu", bf_str))
  # one non-linear mu formula per non-reference category
  expect_true(all(vapply(paste0("mu", ranks4[-1]), function(p) {
    any(grepl(p, names(bf$pforms), fixed = TRUE))
  }, logical(1))))
})


############################################################################# !
# STAN CODE GENERATION TESTS                                              ####
############################################################################# !

test_that("sdt_ranking produces multinomial stancode with the ranking functions", {
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 4)
  model <- sdt_ranking(ranks4, m = 4)
  code <- stancode(bmf(dprime ~ 1), data = dat, model = model)
  expect_true(grepl("multinomial", code))
  expect_true(grepl("sdt_ranking_logmu", code))
  expect_true(grepl("sdt_ranking_logp", code))
  expect_false(grepl("sdt_ranking_lpmf", code))
})

test_that("sdt_ranking stancode includes the UV function for the normal dist", {
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 3, dprime = 1.2, m = 4,
                      dist = "normal")
  normal <- sdt_ranking(ranks4, m = 4, dist = "normal")
  code <- stancode(bmf(dprime ~ 1), data = dat, model = normal)
  expect_true(grepl("sdt_ranking_uv_logp", code))
})

test_that("Gaussian ranking fixes sdratio for EV and frees it for UV", {
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 3, dprime = 1.2, m = 4,
                      dist = "normal")
  normal <- sdt_ranking(ranks4, m = 4, dist = "normal")

  ev_sdratio <- subset(default_prior(bmf(dprime ~ 1), data = dat, model = normal),
                       nlpar == "sdratio")
  expect_true(any(grepl("constant", ev_sdratio$prior)))

  uv_sdratio <- subset(default_prior(bmf(dprime ~ 1, sdratio ~ 1), data = dat,
                                     model = normal),
                       nlpar == "sdratio")
  expect_false(all(grepl("constant", uv_sdratio$prior)))
})

test_that("sdt_ranking stancode works with predictors and random effects", {
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 5, dprime = 1.5, m = 4)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt_ranking(ranks4, m = 4)
  expect_true(nchar(stancode(bmf(dprime ~ condition), data = dat, model = model)) > 0)
  expect_true(nchar(stancode(bmf(dprime ~ 1 + (1 | id)), data = dat, model = model)) > 0)
})


############################################################################# !
# PIPELINE INTEGRATION TESTS                                             ####
############################################################################# !

test_that("sdt_ranking integrates with the bmm pipeline via mock backend", {
  dat <- rsdt_ranking(n_per_cell = 80, n_subjects = 8, dprime = 1.2, m = 4)
  for (d in c("gumbel_min", "normal")) {
    model <- sdt_ranking(ranks4, m = 4, dist = d)
    expect_silent(
      bmm(bmf(dprime ~ 1), dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
    )
  }
})

test_that("sdt_ranking fits mixed set sizes through the pipeline (mock)", {
  pad <- function(d, cols = ranks5) {
    for (p in setdiff(cols, names(d))) d[[p]] <- 0L
    d[c("id", "set_size", cols)]
  }
  dat <- rbind(
    pad(rsdt_ranking(n_per_cell = 80, n_subjects = 4, dprime = 1.2, m = 3)),
    pad(rsdt_ranking(n_per_cell = 80, n_subjects = 4, dprime = 1.2, m = 5))
  )
  model <- sdt_ranking(ranks5, m = "set_size")
  expect_silent(
    bmm(bmf(dprime ~ 1), dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_ranking unequal-variance ranking runs through the pipeline (mock)", {
  dat <- rsdt_ranking(n_per_cell = 80, n_subjects = 8, dprime = 1.2, m = 4,
                      dist = "normal", sdratio = 1.3)
  model <- sdt_ranking(ranks4, m = 4, dist = "normal")
  expect_silent(
    bmm(bmf(dprime ~ 1, sdratio ~ 1), dat, model,
        backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_ranking default_prior returns a valid prior object", {
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 4)
  model <- sdt_ranking(ranks4, m = 4)
  expect_s3_class(default_prior(bmf(dprime ~ 1), data = dat, model = model),
                  "brmsprior")
})

test_that("sdt_ranking recovers d' and returns a multinomial count array", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("cmdstanr")

  dat <- rsdt_ranking(n_per_cell = 200, n_subjects = 10, dprime = 1.4, m = 4)
  model <- sdt_ranking(ranks4, m = 4)
  fit <- bmm(bmf(dprime ~ 1 + (1 | id)), dat, model, backend = "cmdstanr",
             chains = 2, iter = 600, warmup = 300, init = 0.5,
             refresh = 0, silent = 2)

  est <- brms::fixef(fit)["dprime_Intercept", "Estimate"]
  expect_equal(est, 1.4, tolerance = 0.3)

  pp <- brms::posterior_predict(fit, ndraws = 25)
  expect_length(dim(pp), 3)
  expect_equal(dim(pp)[3], 4L)
  expect_false(any(is.nan(pp)))
})
