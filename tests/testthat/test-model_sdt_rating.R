# Test Rating SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_rating model can be created with vector response", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_rating")
  expect_equal(model$other_vars$n_ratings, 4L)
})

test_that("sdt_rating model infers n_ratings from response length", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  expect_equal(model$other_vars$n_ratings, 6L)
})

test_that("sdt_rating model rejects mismatched n_ratings", {
  expect_error(
    sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", n_ratings = 6),
    "must match"
  )
})

test_that("sdt_rating model rejects n_ratings <= 2", {
  expect_error(sdt_rating(c("r1", "r2"), "stimulus"), "n_ratings > 2")
})

test_that("sdt_rating model has parsimonious threshold params by default", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("spacing" %in% names(model$parameters))
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_false("mu" %in% names(model$parameters))
})

test_that("sdt_rating model has log_distance threshold params", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  expect_true("delta1" %in% names(model$parameters))
  expect_true("delta3" %in% names(model$parameters))
  expect_false("delta2" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
})

test_that("sdt_rating model does not estimate a base mu", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_false("mu" %in% names(model$parameters))
})

test_that("sdt_rating stores threshold_type correctly", {
  model_pa <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_equal(model_pa$other_vars$threshold_type, "parsimonious")

  model_ld <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                         threshold_type = "log_distance")
  expect_equal(model_ld$other_vars$threshold_type, "log_distance")
})

test_that("sdt_rating model stores all distribution options", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", dist = di)
    expect_equal(model$other_vars$dist, di)
  }
})

test_that("sdt_rating model has log_ratio threshold params", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_ratio")
  expect_true("delta1" %in% names(model$parameters))
  expect_true("delta3" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
  expect_false("delta2" %in% names(model$parameters))
})

test_that("sdt_rating model has softmax threshold params", {
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "softmax")
  expect_true("spacing" %in% names(model$parameters))
  # K = 6 -> K - 3 = 3 allocation deltas
  expect_true(all(c("delta1", "delta2", "delta3") %in% names(model$parameters)))
})

test_that("sdt_rating always has sdratio as parameter (fixed to 0 by default)", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("sdratio" %in% names(model$parameters))
  expect_true("sdratio" %in% names(model$links))
  expect_true("sdratio" %in% names(model$default_priors))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_rating model accepts custom links", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      links = list(dprime = "log"))
  expect_equal(model$links$dprime, "log")
  expect_equal(model$links$criterion, "identity")
})

test_that("sdt_rating supplies init_ranges for every estimated parameter", {
  # create_initfun() looks up init_ranges per parameter; a missing entry yields
  # NA inits, and the flexible threshold types reject brms' default random init.
  for (tt in c("parsimonious", "equidistant", "log_distance",
               "log_ratio", "softmax")) {
    model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = tt)
    est_pars <- setdiff(names(model$parameters), "sdratio")
    expect_true(all(est_pars %in% names(model$init_ranges)),
                info = paste("threshold_type:", tt))
    expect_true("sdratio" %in% names(model$init_ranges))
    for (r in model$init_ranges) {
      expect_length(r, 2)
      expect_lte(r[1], r[2])
    }
  }
})

test_that(".sdt_make_thresholds mid matches the Stan builders for odd and even K", {
  # Regression: the R companion's middle-threshold index must equal the Stan
  # builders' ((K-1) %/% 2 + 1) so prediction reproduces the likelihood. The two
  # diverged for odd K (R used K %/% 2), shifting equidistant/anchored thresholds.
  for (K in 4:7) {
    mid <- (K - 1L) %/% 2L + 1L
    for (tt in c("log_distance", "log_ratio", "softmax")) {
      nd <- if (tt == "softmax") max(0L, K - 3L) else K - 2L
      deltas <- if (nd > 0) seq(0.1, 0.3, length.out = nd) else NULL
      thr <- .sdt_make_thresholds(0.3, K, tt, spacing = 0.2, deltas = deltas)
      expect_equal(thr[mid], 0.3, info = paste(tt, "K =", K))    # criterion sits at mid
      expect_false(is.unsorted(thr), info = paste(tt, "K =", K))
    }
  }
  # equidistant closed form criterion + (k - mid) * exp(spacing); odd K was the bug
  expect_equal(.sdt_make_thresholds(0, 5L, "equidistant", spacing = 0), c(-2, -1, 0, 1))
  expect_equal(.sdt_make_thresholds(0, 7L, "equidistant", spacing = 0), c(-3, -2, -1, 0, 1, 2))
  expect_equal(.sdt_make_thresholds(0, 6L, "equidistant", spacing = 0), c(-2, -1, 0, 1, 2))
})

test_that("odd-K delta labels skip the same middle threshold as the builders", {
  model <- sdt_rating(paste0("r", 1:5), "stimulus",
                      threshold_type = "log_distance")
  expect_true(all(c("delta1", "delta2", "delta4") %in% names(model$parameters)))
  expect_false("delta3" %in% names(model$parameters))
})

############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

sim_rating <- function(n_subjects, n_trials, dprime, criterion, n_ratings,
                       spacing = NULL, deltas = NULL,
                       threshold_type = "parsimonious",
                       sdratio = 1, dist = "normal") {
  thr <- bmm:::.sdt_make_thresholds(criterion, n_ratings, threshold_type,
                                    spacing, deltas)
  dat <- expand.grid(id = seq_len(n_subjects), stimulus = c(0L, 1L))
  cbind(dat, as.data.frame(rsdt_rating(nrow(dat), n_trials, dat$stimulus,
                                       dprime, thr, sdratio = sdratio,
                                       dist = dist)))
}

test_that("sdt_rating check_data validates response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  valid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 1L)
  )
  result <- check_data(model, valid_data, formula)
  # multinomial encoding: one row per observation with a count matrix Y
  expect_equal(nrow(result), nrow(valid_data))
  expect_true(all(c("Y", "nTrials", "stimulus") %in% colnames(result)))
  expect_equal(ncol(result$Y), 4)
  expect_equal(result$nTrials, c(50, 100))
})

test_that("sdt_rating check_data rejects missing response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10),
    stimulus = c(0L, 1L)
  )
  expect_error(check_data(model, invalid_data, formula), "missing in the data")
})

test_that("sdt_rating check_data rejects negative counts", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(-1, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 1L)
  )
  expect_error(check_data(model, invalid_data, formula), "non-negative")
})

test_that("sdt_rating check_data validates stimulus coding", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 2L)
  )
  expect_error(check_data(model, invalid_data, formula), "must be coded as 0")
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_rating generates counts matching the design", {
  counts <- rsdt_rating(10, 50, rep(c(0L, 1L), 5), dprime = 1.5,
                        thresholds = c(-0.5, 0, 0.5))
  expect_true(is.matrix(counts))
  expect_equal(dim(counts), c(10L, 4L))
  expect_equal(colnames(counts), paste0("r", 1:4))
  expect_true(all(rowSums(counts) == 50))
})

test_that("rsdt_rating recycles vectorized parameters per observation", {
  counts <- rsdt_rating(6, c(50, 100), c(0L, 1L), dprime = rnorm(6, 1.5, 0.2),
                        thresholds = c(-0.5, 0, 0.5), sdratio = c(1, 1.3))
  expect_equal(dim(counts), c(6L, 4L))
  expect_true(all(rowSums(counts) == c(50, 100)))

  thr <- rbind(c(-0.5, 0, 0.5), c(-1, 0, 1))
  counts2 <- rsdt_rating(2, 100, c(0L, 1L), dprime = 1.5, thresholds = thr)
  expect_equal(dim(counts2), c(2L, 4L))
})

test_that("rsdt_rating validates inputs", {
  expect_error(rsdt_rating(c(2, 3), 50, 1L, dprime = 1,
                           thresholds = c(-0.5, 0, 0.5)),
               "single positive integer")
  expect_error(rsdt_rating(2, 50, c(0L, 2L), dprime = 1,
                           thresholds = c(-0.5, 0, 0.5)),
               "0 \\(noise\\) or 1 \\(signal\\)")
  expect_error(rsdt_rating(2, 50, c(0L, 1L), dprime = 1,
                           thresholds = c(-0.5, 0, 0.5), dist = "cauchy"),
               "should be one of")
})

test_that("sim helper round-trips every threshold parameterization", {
  configs <- list(
    list(n_ratings = 4L, spacing = 0.5, threshold_type = "parsimonious"),
    list(n_ratings = 4L, spacing = 0.5, threshold_type = "equidistant"),
    list(n_ratings = 4L, deltas = c(0.5, 0.5), threshold_type = "log_distance"),
    list(n_ratings = 4L, deltas = c(0, 0), threshold_type = "log_ratio"),
    list(n_ratings = 6L, spacing = 0.3, deltas = c(0, 0, 0),
         threshold_type = "softmax")
  )
  for (cfg in configs) {
    dat <- do.call(sim_rating, c(list(3, 50, dprime = 1.5, criterion = 0), cfg))
    rcols <- paste0("r", seq_len(cfg$n_ratings))
    expect_equal(nrow(dat), 6, info = cfg$threshold_type)
    expect_true(all(rowSums(dat[, rcols]) == 50), info = cfg$threshold_type)
  }
})

test_that("dsdt_rating returns valid density", {
  d <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                   dprime = 1.5, thresholds = c(-0.5, 0, 0.5))
  expect_true(d > 0)
  expect_true(d <= 1)
})

test_that("dsdt_rating log matches log of density", {
  d <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                   dprime = 1.5, thresholds = c(-0.5, 0, 0.5))
  ld <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                    dprime = 1.5, thresholds = c(-0.5, 0, 0.5), log = TRUE)
  expect_equal(log(d), ld, tolerance = 1e-10)
})

test_that("dsdt_rating validates inputs", {
  expect_error(
    dsdt_rating(counts = c(-1, 20, 30, 40), stimulus = 1,
                dprime = 1.5, thresholds = c(-0.5, 0, 0.5)),
    "non-negative"
  )
  expect_error(
    dsdt_rating(counts = c(10, 20, 30), stimulus = 1,
                dprime = 1.5, thresholds = c(-0.5, 0, 0.5)),
    "K - 1"
  )
})

test_that("dsdt_rating works for all distributions", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    d <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                     dprime = 1.5, thresholds = c(-0.5, 0, 0.5), dist = di)
    expect_true(d > 0, info = paste("dist:", di))
  }
})

test_that("dsdt_rating with sdratio != 1 changes density", {
  d_ev <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                      dprime = 1.5, thresholds = c(-0.5, 0, 0.5), sdratio = 1)
  d_uv <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                      dprime = 1.5, thresholds = c(-0.5, 0, 0.5), sdratio = 1.3)
  expect_false(d_ev == d_uv)
})

test_that("dsdt_rating matches dmultinom and is vectorized over rows", {
  thr <- c(-0.5, 0, 0.5)
  probs <- bmm:::.sdt_category_probs(thr, 1.5, 1, 1, "normal")
  expect_equal(dsdt_rating(c(10, 20, 30, 40), 1, 1.5, thr),
               dmultinom(c(10, 20, 30, 40), prob = probs), tolerance = 1e-12)

  counts <- rbind(c(10, 20, 30, 40), c(40, 30, 20, 10))
  d <- dsdt_rating(counts, stimulus = c(1L, 0L), dprime = 1.5,
                   thresholds = thr, log = TRUE)
  expect_length(d, 2)
  expect_true(all(is.finite(d)))
})

test_that(".sdt_make_thresholds vectorized path matches per-draw evaluation", {
  cr <- c(-0.2, 0.1, 0.4)
  sp <- c(0.1, 0.3, -0.1)
  for (tt in c("parsimonious", "equidistant", "log_distance",
               "log_ratio", "softmax")) {
    K <- 5L
    nd <- if (tt == "softmax") K - 3L else K - 2L
    dl <- if (tt %in% c("log_distance", "log_ratio", "softmax")) {
      matrix(seq(-0.3, 0.3, length.out = 3 * nd), 3, nd)
    }
    vec <- bmm:::.sdt_make_thresholds(cr, K, tt, sp, dl)
    ref <- t(vapply(1:3, function(i) {
      bmm:::.sdt_make_thresholds(cr[i], K, tt, sp[i],
                                 if (is.null(dl)) NULL else dl[i, ])
    }, numeric(K - 1L)))
    expect_equal(vec, ref, tolerance = 1e-12, info = tt)
  }
})

test_that(".sdt_category_probs vectorized path matches per-row evaluation", {
  thr <- rbind(c(-0.5, 0, 0.5), c(-1, 0.2, 0.9))
  probs <- bmm:::.sdt_category_probs(thr, c(1.2, 0.8), c(1.3, 1), c(1L, 0L),
                                     "normal")
  ref <- t(vapply(1:2, function(i) {
    bmm:::.sdt_category_probs(thr[i, ], c(1.2, 0.8)[i], c(1.3, 1)[i],
                              c(1L, 0L)[i], "normal")
  }, numeric(4)))
  expect_equal(probs, ref, tolerance = 1e-12)
})

test_that("sdt_rating_logmu is vectorized and preserves the draw shape", {
  dprime <- matrix(c(1, 1.5, 0.5, 2), nrow = 2)
  out <- sdt_rating_logmu(2L, 4L, 1L, 1L, dprime, criterion = 0.1,
                          spacing = 0.3, sdratio = 0.2, stimulus = c(0, 1))
  expect_equal(dim(out), dim(dprime))

  thr <- bmm:::.sdt_make_thresholds(0.1, 4L, "parsimonious", 0.3)
  stim <- rep_len(c(0, 1), 4)
  ref <- vapply(seq_along(dprime), function(j) {
    log(bmm:::.sdt_category_probs(thr, as.vector(dprime)[j], exp(0.2),
                                  stim[j], "normal")[2])
  }, numeric(1))
  expect_equal(as.vector(out), ref, tolerance = 1e-12)
})

test_that("log_ratio requires n_ratings >= 4", {
  expect_error(sdt_rating(c("r1", "r2", "r3"), "stimulus",
                          threshold_type = "log_ratio"),
               "n_ratings >= 4")
})

test_that("K = 3 softmax reduces to a single symmetric interval", {
  thr <- bmm:::.sdt_make_thresholds(0.1, 3L, "softmax", spacing = 0.2)
  expect_equal(thr, c(0.1 - exp(0.2), 0.1), tolerance = 1e-12)

  dat <- sim_rating(3, 50, dprime = 1.2, criterion = 0, n_ratings = 3,
                    spacing = 0.2, threshold_type = "softmax")
  model <- sdt_rating(paste0("r", 1:3), "stimulus", threshold_type = "softmax")
  code <- stancode(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1),
                   data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("sdt_rating produces valid stancode with parsimonious thresholds", {
  dat <- sim_rating(3, 50, dprime = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
  expect_true(grepl("sdt_thresholds_parsimonious_rating", code, fixed = TRUE))
})

test_that("sdt_rating stancode loads the shared distribution dispatch", {
  # the noise-distribution dispatch lives in sdt_dist_funs.stan; the rating
  # likelihood relies on the log-scale dispatchers from that shared chunk.
  dat <- sim_rating(3, 50, dprime = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(grepl("sdt_log_cumprob", code, fixed = TRUE))
  expect_true(grepl("sdt_log_one_minus_cumprob", code, fixed = TRUE))
  expect_true(grepl("log_diff_exp", code, fixed = TRUE))
})

test_that("sdt_rating produces valid stancode with log_distance thresholds", {
  dat <- sim_rating(3, 50, dprime = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0.5, 0.5), threshold_type = "log_distance")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
  expect_true(grepl("sdt_thresholds_log_distance_rating", code, fixed = TRUE))
})

test_that("sdt_rating produces valid stancode with log_ratio thresholds", {
  dat <- sim_rating(3, 50, dprime = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0, 0), threshold_type = "log_ratio")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_ratio")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_thresholds_log_ratio_rating", code, fixed = TRUE))
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
})

test_that("sdt_rating produces valid stancode with softmax thresholds", {
  dat <- sim_rating(3, 50, dprime = 1.5, criterion = 0, n_ratings = 6,
                    spacing = 0.3, deltas = c(0, 0, 0),
                    threshold_type = "softmax")
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "softmax")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1,
                 delta1 ~ 1, delta2 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_thresholds_softmax_rating", code, fixed = TRUE))
})

test_that("sdt_rating handles K=6 parsimonious", {
  dat <- sim_rating(2, 50, dprime = 1.5, criterion = 0, n_ratings = 6,
                    spacing = 0.3)
  model <- sdt_rating(paste0("r", 1:6), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_rating handles predictors on dprime", {
  dat <- sim_rating(3, 50, dprime = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ condition, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# UV-SDT TESTS (sdratio overridable fixed parameter)                     ####
############################################################################# !

test_that("sdt_rating with sdratio ~ 1 produces valid stancode", {
  dat <- sim_rating(3, 50, dprime = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5, sdratio = 1.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
  expect_true(grepl("sdratio", code))
})

test_that("sdt_rating UV-SDT with log_distance thresholds produces valid stancode", {
  dat <- sim_rating(3, 50, dprime = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0.5, 0.5), threshold_type = "log_distance",
                    sdratio = 1.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1,
                 sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdratio", code))
})


############################################################################# !
# PIPELINE INTEGRATION TESTS (mock backend)                              ####
############################################################################# !

test_that("sdt_rating integrates with the bmm pipeline via mock backend", {
  dat <- sim_rating(6, 80, dprime = 1.2, criterion = 0, n_ratings = 4,
                    spacing = 0.4)
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", dist = di)
    expect_silent(
      bmm(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1),
          dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
    )
  }
})

test_that("sdt_rating UV-SDT integrates with the bmm pipeline via mock backend", {
  dat <- sim_rating(6, 80, dprime = 1.2, criterion = 0, n_ratings = 4,
                    spacing = 0.4, sdratio = 1.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_silent(
    bmm(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
        dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_rating log_distance integrates with the bmm pipeline via mock backend", {
  # log_distance has variable-arity delta parameters, so this exercises the
  # generated multinomial formula and Stan function for that case.
  dat <- sim_rating(5, 60, dprime = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0.5, 0.5), threshold_type = "log_distance")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  expect_silent(
    bmm(bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1),
        dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_rating posterior_predict draws joint multinomial counts (real fit)", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)))

  dat <- sim_rating(6, 120, dprime = 1.4, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  fit <- bmm(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1), dat, model,
             backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  pp <- brms::posterior_predict(fit)
  expect_length(dim(pp), 3)
  row_sums <- apply(pp, c(1, 2), sum)
  expect_true(all(row_sums == 120))
})

test_that("sdt_rating log_distance fits with finite likelihood at K=6 (real fit)", {
  # Regression: at K=6 the log_distance threshold builder has mid=3, so its
  # lower thresholds were filled by a descending Stan for-loop. Stan for-loops
  # only count up, so the loop was skipped and thresholds[1:2] stayed NaN,
  # failing every chain. Earlier tests only used K=4 (mid=2), where the loop ran.
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)))

  dat <- sim_rating(6, 150, dprime = 1.4, criterion = 0, n_ratings = 6,
                    deltas = c(-0.4, -0.4, -0.4, -0.4),
                    threshold_type = "log_distance")
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "log_distance")
  fit <- bmm(bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1,
                 delta4 ~ 1, delta5 ~ 1),
             dat, model, backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  expect_true(all(is.finite(brms::log_lik(fit))))
})
