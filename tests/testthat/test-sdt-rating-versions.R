# Tests for the dual-process (dpsdt) and meta-d' (metad) versions of sdt_rating:
# constructor wiring, the fixed-parameter freeing of recollection, the reduction
# identities to the standard model, distribution/simulation functions, and the
# version-aware post-processing (roc_sdt / latent_sdt) via the mock fits in
# helper-sdt-analysis.R.

############################################################################# !
# CONSTRUCTOR & VERSION WIRING                                           ####
############################################################################# !

test_that("dpsdt version adds Ro/Rn fixed off by default", {
  m <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
  expect_s3_class(m, "sdt_rating_dpsdt")
  expect_equal(m$version, "dpsdt")
  expect_true(all(c("Ro", "Rn") %in% names(m$parameters)))
  expect_equal(m$fixed_parameters, list(sdratio = 0, Ro = -100, Rn = -100))
  expect_equal(bmm:::.sdt_rating_variant(m)$logmu_fun, "sdt_dpsdt_logmu")
  expect_equal(bmm:::.sdt_rating_variant(m)$extra_params, c("Ro", "Rn"))
})

test_that("metad version adds logmratio (log M-ratio) estimated by default", {
  m <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", version = "metad")
  expect_s3_class(m, "sdt_rating_metad")
  expect_equal(m$version, "metad")
  expect_true("logmratio" %in% names(m$parameters))
  expect_false("metad" %in% names(m$parameters))
  expect_equal(m$links$logmratio, "identity")
  expect_equal(m$fixed_parameters, list(sdratio = 0))
  expect_equal(bmm:::.sdt_rating_variant(m)$logmu_fun, "sdt_metad_logmu")
  expect_equal(bmm:::.sdt_rating_variant(m)$extra_params, "logmratio")
})

test_that("legacy model objects without a variant resolve to the standard one", {
  m <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  m$version <- "rating"
  expect_equal(bmm:::.sdt_rating_variant(m)$logmu_fun, "sdt_rating_logmu")
})

test_that("standard version is the default and keeps the base parameters", {
  m <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_equal(m$version, "standard")
  expect_s3_class(m, "sdt_rating_standard")
  expect_false(any(c("Ro", "Rn", "metad", "logmratio") %in% names(m$parameters)))
})

test_that("recollection is freed from fixed_parameters via the formula", {
  m <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
  f <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  expect_equal(names(update_model_fixed_parameters(m, f)$fixed_parameters),
               c("sdratio", "Ro", "Rn"))

  f1 <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1)
  expect_equal(names(update_model_fixed_parameters(m, f1)$fixed_parameters),
               c("sdratio", "Rn"))

  f2 <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1)
  expect_equal(names(update_model_fixed_parameters(m, f2)$fixed_parameters),
               "sdratio")
})

############################################################################# !
# CATEGORY PROBABILITIES: REDUCTION IDENTITIES & STRUCTURE               ####
############################################################################# !

test_that("dpsdt with no recollection reduces to the standard category probs", {
  thr <- bmm:::.sdt_make_thresholds(0, 5, "parsimonious", spacing = 0.3)
  for (stim in c(0L, 1L)) {
    base <- bmm:::.sdt_category_probs(thr, 1.4, 1.2, stim, "normal")
    off  <- bmm:::.sdt_dpsdt_category_probs(thr, 1.4, 1.2, stim, "normal",
                                            plogis(-100), plogis(-100))
    expect_equal(off, base, tolerance = 1e-12)
  }
})

test_that("metad equal to dprime reduces to the standard category probs", {
  thr <- bmm:::.sdt_make_thresholds(0, 5, "parsimonious", spacing = 0.3)
  for (stim in c(0L, 1L)) {
    base <- bmm:::.sdt_category_probs(thr, 1.4, 1.2, stim, "normal")
    md   <- bmm:::.sdt_metad_category_probs(thr, 1.4, 1.4, stim, 1.2, "normal")
    expect_equal(md, base, tolerance = 1e-12)
  }
})

test_that("dpsdt category probs match an independent hand computation (K=4)", {
  thr <- c(-0.6, 0.1, 0.7)
  d <- 1.5
  Ro <- 0.35
  base_s <- diff(c(0, pnorm(thr - d / 2), 1))            # signal, sr = 1
  hand_s <- (1 - Ro) * base_s
  hand_s[4] <- hand_s[4] + Ro
  got_s <- bmm:::.sdt_dpsdt_category_probs(thr, d, 1, 1L, "normal", Ro, 0)
  expect_equal(got_s, hand_s / sum(hand_s), tolerance = 1e-12)
})

test_that("recollection moves mass to the most-confident category and sums to 1", {
  thr <- bmm:::.sdt_make_thresholds(0, 5, "parsimonious", spacing = 0.3)
  base_s <- bmm:::.sdt_category_probs(thr, 1.4, 1, 1L, "normal")
  ps     <- bmm:::.sdt_dpsdt_category_probs(thr, 1.4, 1, 1L, "normal", 0.4, 0)
  expect_gt(ps[5], base_s[5])
  expect_equal(sum(ps), 1, tolerance = 1e-12)

  base_n <- bmm:::.sdt_category_probs(thr, 1.4, 1, 0L, "normal")
  pn     <- bmm:::.sdt_dpsdt_category_probs(thr, 1.4, 1, 0L, "normal", 0, 0.4)
  expect_gt(pn[1], base_n[1])
})

############################################################################# !
# DISTRIBUTION & SIMULATION FUNCTIONS                                    ####
############################################################################# !

test_that("rsdt_dpsdt / rsdt_metad return rating count matrices", {
  thr <- bmm:::.sdt_make_thresholds(0, 4L, "parsimonious", 0.3)
  d <- rsdt_dpsdt(8, 50, rep(c(0L, 1L), 4), dprime = 1.4, thresholds = thr,
                  Ro = 0.3, Rn = 0.1)
  expect_true(is.matrix(d))
  expect_equal(colnames(d), paste0("r", 1:4))
  expect_true(all(rowSums(d) == 50))

  m <- rsdt_metad(8, 50, rep(c(0L, 1L), 4), dprime = 1.4, thresholds = thr,
                  metad = 1.0)
  expect_equal(dim(m), c(8L, 4L))
})

test_that("dsdt_dpsdt / dsdt_metad return finite densities and validate inputs", {
  expect_true(is.finite(dsdt_dpsdt(c(2, 8, 20, 70), 1L, 1.5,
                                   c(-0.5, 0, 0.5), Ro = 0.3, Rn = 0)))
  expect_true(is.finite(dsdt_metad(c(5, 15, 25, 55), 1L, 1.5,
                                   c(-0.5, 0, 0.5), metad = 1.0)))
  expect_error(dsdt_dpsdt(c(1, 2, 3, 4), 1L, 1, c(0, 0.5, 1), Ro = 1.5, Rn = 0),
               "Ro must be a probability")
  expect_error(rsdt_metad(c(2, 3), 50, 1L, 1.4, c(-0.5, 0, 0.5), metad = 1),
               "single positive integer")
})

test_that("dsdt version densities match dmultinom and vectorize over rows", {
  thr <- c(-0.5, 0, 0.5)
  p_dp <- bmm:::.sdt_dpsdt_category_probs(thr, 1.5, 1, 1L, "normal", 0.3, 0)
  expect_equal(dsdt_dpsdt(c(2, 8, 20, 70), 1L, 1.5, thr, Ro = 0.3, Rn = 0),
               dmultinom(c(2, 8, 20, 70), prob = p_dp), tolerance = 1e-12)

  p_md <- bmm:::.sdt_metad_category_probs(thr, 1.5, 1.0, 1L, 1, "normal")
  expect_equal(dsdt_metad(c(5, 15, 25, 55), 1L, 1.5, thr, metad = 1.0),
               dmultinom(c(5, 15, 25, 55), prob = p_md), tolerance = 1e-12)

  counts <- rbind(c(2, 8, 20, 70), c(70, 20, 8, 2))
  d <- dsdt_dpsdt(counts, c(1L, 0L), 1.5, thr, Ro = c(0.3, 0.2),
                  Rn = c(0, 0.1), log = TRUE)
  expect_length(d, 2)
  expect_true(all(is.finite(d)))
})

test_that("version category probs vectorized path matches per-draw evaluation", {
  thr <- rbind(c(-0.5, 0, 0.5), c(-1, 0.2, 0.9), c(-0.8, -0.1, 0.4))
  dp <- c(1.2, 0.8, 1.6)
  sr <- c(1.3, 1, 1.1)
  ro <- c(0.3, 0.1, 0.5)
  rn <- c(0.05, 0.2, 0)
  md <- c(1.0, 0.6, 1.4)
  for (stim in c(0L, 1L)) {
    vec_dp <- bmm:::.sdt_dpsdt_category_probs(thr, dp, sr, stim, "normal",
                                              ro, rn)
    ref_dp <- t(vapply(1:3, function(i) {
      bmm:::.sdt_dpsdt_category_probs(thr[i, ], dp[i], sr[i], stim, "normal",
                                      ro[i], rn[i])
    }, numeric(4)))
    expect_equal(vec_dp, ref_dp, tolerance = 1e-12, info = paste("dpsdt", stim))

    vec_md <- bmm:::.sdt_metad_category_probs(thr, dp, md, stim, sr, "normal")
    ref_md <- t(vapply(1:3, function(i) {
      bmm:::.sdt_metad_category_probs(thr[i, ], dp[i], md[i], stim, sr[i],
                                      "normal")
    }, numeric(4)))
    expect_equal(vec_md, ref_md, tolerance = 1e-12, info = paste("metad", stim))
  }
})

############################################################################# !
# PIPELINE INTEGRATION (mock backend)                                    ####
############################################################################# !

test_that("all rating versions configure through the bmm pipeline", {
  thr <- bmm:::.sdt_make_thresholds(0, 4L, "parsimonious", 0.3)
  dat <- expand.grid(id = 1:4, stimulus = c(0L, 1L))
  dat <- cbind(dat, as.data.frame(rsdt_dpsdt(nrow(dat), 40, dat$stimulus,
                                             dprime = 1.4, thresholds = thr,
                                             Ro = 0.3, Rn = 0.1)))
  md <- sdt_rating(paste0("r", 1:4), "stimulus", version = "dpsdt")
  mm <- sdt_rating(paste0("r", 1:4), "stimulus", version = "metad")
  expect_silent(bmm(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1),
                    dat, md, backend = "mock", mock_fit = 1, rename = FALSE))
  expect_silent(bmm(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, logmratio ~ 1),
                    dat, mm, backend = "mock", mock_fit = 1, rename = FALSE))
})

test_that("the generated Stan wrapper carries the version-specific call", {
  md <- sdt_rating(paste0("r", 1:4), "stimulus", version = "dpsdt")
  code <- bmm:::.sdt_rating_logmu_stan(md)
  expect_match(code, "real sdt_dpsdt_logmu\\(")
  expect_match(code, "real Ro, real Rn, real stimulus")
  expect_match(code, "sdt_dpsdt_logmu_cat\\(cat, thr, dprime, sdratio, stimulus, dist_type, Ro, Rn\\)")

  mm <- sdt_rating(paste0("r", 1:4), "stimulus", version = "metad")
  code_md <- bmm:::.sdt_rating_logmu_stan(mm)
  # the estimated parameter is the log M-ratio; meta-d' is derived in the call
  expect_match(code_md, "real logmratio, real stimulus")
  expect_match(code_md, "sdt_metad_logmu_cat\\(cat, thr, dprime, exp\\(logmratio\\) \\* dprime")
})

test_that("metad R companion with logmratio = 0 reduces to the standard model", {
  mm  <- sdt_rating(paste0("r", 1:5), "stimulus", version = "metad")
  std <- sdt_rating(paste0("r", 1:5), "stimulus")
  args <- list(K = 5L, dist = 1L, thresh = 1L, dprime = matrix(1.4),
               criterion = matrix(0), spacing = matrix(0.3),
               sdratio = matrix(0), stimulus = matrix(1))
  for (k in 1:5) {
    md_k <- do.call(sdt_metad_logmu, c(list(cat = k), args, list(logmratio = matrix(0))))
    st_k <- do.call(sdt_rating_logmu, c(list(cat = k), args))
    expect_equal(as.numeric(md_k), as.numeric(st_k), tolerance = 1e-12)
  }
})

############################################################################# !
# VERSION-AWARE POST-PROCESSING (mock fits)                              ####
############################################################################# !

test_that("roc_sdt reflects dual-process recollection (higher AUC, lifted curve)", {
  fit_std <- fake_rating_fit(n_ratings = 6L)
  fit_dp  <- fake_rating_fit(n_ratings = 6L, version = "dpsdt")
  draws_std <- list(dprime = 1.2, criterion = 0, spacing = 0)
  draws_dp  <- c(draws_std, list(Ro = qlogis(0.4), Rn = qlogis(0.2)))

  local_mocked_bindings(ranef = function(...) list(), .package = "brms")
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(draws_std), .package = "brms")
  roc_std <- roc_sdt(fit_std, n_points = 50)

  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(draws_dp), .package = "brms")
  roc_dp <- roc_sdt(fit_dp, n_points = 50)

  # recollection raises the smooth curve toward the top-left -> larger area
  expect_gt(mean(attr(roc_dp, "summary")$Hit_mean),
            mean(attr(roc_std, "summary")$Hit_mean))
  # the K-1 threshold operating points lie on the lifted dual-process curve
  expect_true(all(attr(roc_dp, "points")$Hit_mean >= -1e-9))
})

test_that("default dpsdt roc_sdt (recollection off) matches the standard roc", {
  fit_std <- fake_rating_fit(n_ratings = 6L)
  fit_dp  <- fake_rating_fit(n_ratings = 6L, version = "dpsdt")
  draws_std <- list(dprime = 1.2, criterion = 0, spacing = 0)
  draws_off <- c(draws_std, list(Ro = -100, Rn = -100))

  local_mocked_bindings(ranef = function(...) list(), .package = "brms")
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(draws_std), .package = "brms")
  roc_std <- roc_sdt(fit_std, n_points = 40)

  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(draws_off), .package = "brms")
  roc_dp <- roc_sdt(fit_dp, n_points = 40)

  expect_equal(attr(roc_dp, "summary")$Hit_mean,
               attr(roc_std, "summary")$Hit_mean, tolerance = 1e-8)
})

test_that("latent_sdt reports the response-process parameters as an attribute", {
  fit_dp <- fake_rating_fit(n_ratings = 6L, version = "dpsdt")
  local_mocked_bindings(ranef = function(...) list(), .package = "brms")
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(
      list(dprime = 1.2, criterion = 0, spacing = 0,
           Ro = qlogis(0.4), Rn = qlogis(0.2))),
    .package = "brms")
  lat <- latent_sdt(fit_dp)
  extra <- attr(lat, "extra")
  expect_false(is.null(extra))
  expect_setequal(extra$parameter, c("Ro", "Rn"))
  expect_equal(extra$mean[extra$parameter == "Ro"], 0.4, tolerance = 1e-6)
})
