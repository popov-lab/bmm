# Tests for roc_sdt(), roc_observed(), auc_sdt() in helpers-analysis.R.
# Mock/structural only: no Stan fit. Posterior draws are supplied by mocking
# brms::posterior_linpred (and brms::ranef / brms::variables) so the ROC/AUC
# orchestration and the binary multi-criteria logic are exercised against
# controlled, deterministic inputs. Fake-fit builders live in
# helper-sdt-analysis.R.


############################################################################# !
# INPUT VALIDATION                                                       ####
############################################################################# !

test_that("roc_sdt / auc_sdt / roc_observed reject non-bmmfit input", {
  expect_error(roc_sdt(list()), "bmmfit")
  expect_error(auc_sdt(list()), "bmmfit")
  expect_error(roc_observed(list()), "bmmfit")
})

test_that("roc_sdt / auc_sdt error for criterion-free models", {
  fake <- function(cls) {
    structure(list(bmm = list(model = structure(list(),
              class = c("bmmodel", "sdt", cls)))), class = c("bmmfit", "brmsfit"))
  }
  expect_error(roc_sdt(fake("sdt_mafc")), "not defined")
  expect_error(roc_sdt(fake("sdt_ranking")), "not defined")
  expect_error(auc_sdt(fake("sdt_mafc")), "not defined")
  expect_error(auc_sdt(fake("sdt_ranking")), "not defined")
})

test_that("roc_sdt errors for a non-SDT model", {
  fake <- structure(list(bmm = list(model = structure(list(),
            class = c("bmmodel", "sdm")))), class = c("bmmfit", "brmsfit"))
  expect_error(roc_sdt(fake), "only available for SDT")
})


############################################################################# !
# RATING ROC MATH (pure helpers, no fit)                                 ####
############################################################################# !

test_that("rating category probs yield valid, monotone ROC points (all dists)", {
  thr <- .sdt_make_thresholds(0, 6L, "parsimonious", spacing = 0)
  for (dist in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    pn <- .sdt_category_probs(thr, 1.5, 1, 0L, dist)
    ps <- .sdt_category_probs(thr, 1.5, 1, 1L, dist)
    expect_equal(sum(pn), 1, tolerance = 1e-8)
    expect_equal(sum(ps), 1, tolerance = 1e-8)

    fa  <- c(1, 1 - cumsum(pn)[1:5], 0)
    hit <- c(1, 1 - cumsum(ps)[1:5], 0)
    expect_true(all(fa >= 0 & fa <= 1 & hit >= 0 & hit <= 1))
    expect_true(all(diff(fa) <= 1e-9), info = dist)
    expect_true(all(diff(hit) <= 1e-9), info = dist)
    expect_true(all(hit >= fa - 1e-9), info = dist)  # d' > 0 -> ROC above diagonal
  }
})

test_that("threshold reconstruction works for all 5 types incl softmax boundaries", {
  specs <- list(
    list(tt = "parsimonious", K = 6L, spacing = 0,    deltas = NULL),
    list(tt = "equidistant",  K = 6L, spacing = 0,    deltas = NULL),
    list(tt = "log_distance", K = 6L, spacing = NULL, deltas = rep(0, 4)),
    list(tt = "log_ratio",    K = 6L, spacing = NULL, deltas = rep(0, 4)),
    list(tt = "softmax",      K = 6L, spacing = 0,    deltas = rep(0, 3)),
    list(tt = "softmax",      K = 4L, spacing = 0,    deltas = 0)
  )
  for (s in specs) {
    thr <- .sdt_make_thresholds(0, s$K, s$tt, spacing = s$spacing, deltas = s$deltas)
    expect_length(thr, s$K - 1L)
    expect_true(all(diff(thr) > 0), info = paste(s$tt, s$K))
    p <- .sdt_category_probs(thr, 1.2, 1, 1L, "normal")
    expect_equal(sum(p), 1, tolerance = 1e-8)
  }
})


############################################################################# !
# CRITERION-POINT DETECTION                                              ####
############################################################################# !

test_that(".sdt_criterion_point_dims classifies criterion-only predictors", {
  fit <- fake_binary_fit(uv = TRUE, multi = TRUE)
  local_mocked_bindings(ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
                        .package = "brms")
  dims <- .sdt_criterion_point_dims(fit)
  expect_equal(dims$points, "condition")
  expect_length(dims$curves, 0L)
})

test_that(".sdt_criterion_point_dims treats dprime predictors as curves", {
  fit <- fake_binary_fit()
  fit$bmm$user_formula <- bmf(dprime ~ 0 + condition, criterion ~ 0 + condition,
                              sdratio ~ 1)
  local_mocked_bindings(ranef = function(...) list(), .package = "brms")
  dims <- .sdt_criterion_point_dims(fit)
  expect_true("condition" %in% dims$curves)
  expect_length(dims$points, 0L)
})

test_that(".sdt_criterion_point_dims is empty for intercept-only fits", {
  fit <- fake_binary_fit()
  local_mocked_bindings(ranef = function(...) list(), .package = "brms")
  dims <- .sdt_criterion_point_dims(fit)
  expect_length(dims$points, 0L)
  expect_length(dims$curves, 0L)
})


############################################################################# !
# sdratio DETECTION                                                      ####
############################################################################# !

test_that(".sdt_has_estimated_sdratio uses the model object", {
  m_ev <- sdt_binary(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  expect_false(.sdt_has_estimated_sdratio(m_ev))
  m_uv <- m_ev
  m_uv$fixed_parameters$sdratio <- NULL
  expect_true(.sdt_has_estimated_sdratio(m_uv))
})

test_that(".sdt_has_estimated_sdratio cross-checks brms::variables", {
  m_ev <- sdt_binary(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  fit <- structure(list(), class = c("bmmfit", "brmsfit"))
  local_mocked_bindings(
    variables = function(...) c("b_dprime_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  expect_true(.sdt_has_estimated_sdratio(m_ev, fit))
})


############################################################################# !
# ROC — BINARY                                                           ####
############################################################################# !

test_that("roc_sdt() binary single criterion returns a smooth curve", {
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_dprime_Intercept"),
    .package = "brms"
  )
  roc <- roc_sdt(fit, n_points = 20)
  expect_s3_class(roc, "bmm_sdt_roc")
  expect_true(all(c("FA", "Hit", ".draw") %in% names(roc)))
  expect_true(all(roc$FA >= 0 & roc$FA <= 1 & roc$Hit >= 0 & roc$Hit <= 1))
  expect_equal(nrow(roc), length(unique(roc$.draw)) * (20L + 2L))
  expect_null(attr(roc, "points"))

  d1 <- roc[roc$.draw == 1L, ]
  expect_true(any(d1$FA == 0 & d1$Hit == 0))
  expect_true(any(d1$FA == 1 & d1$Hit == 1))
})

test_that("roc_sdt() binary auto-detects multi-criteria operating points", {
  fit <- fake_binary_fit(uv = TRUE, multi = TRUE)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(
      dprime = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8), sdratio = log(1.3))),
    ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
    variables = function(...) c("b_dprime_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  roc <- roc_sdt(fit)
  pts <- attr(roc, "points")
  expect_false(is.null(pts))
  expect_equal(nrow(pts), 5L)
  expect_true(all(c("FA_mean", "Hit_mean", "condition") %in% names(pts)))
  expect_true(all(pts$FA_mean >= 0 & pts$FA_mean <= 1))
})

test_that("roc_sdt() binary multi-criteria points fall on the curve (UV + spread)", {
  fit <- fake_binary_fit(uv = TRUE, multi = TRUE)
  z      <- stats::qnorm(stats::ppoints(n_draws_mock))
  spread <- list(dprime = 0.18, criterion = 0.10, sdratio = 0.12)
  base   <- list(dprime = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8),
                 sdratio = log(1.3))
  local_mocked_bindings(
    posterior_linpred = function(object, dpar = NULL, nlpar = NULL,
                                 newdata = NULL, ...) {
      par <- if (!is.null(dpar)) dpar else nlpar
      vals <- rep_len(base[[par]], if (is.null(newdata)) 1L else nrow(newdata))
      matrix(rep(vals, each = n_draws_mock) + spread[[par]] * z,
             nrow = n_draws_mock)
    },
    ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
    variables = function(...) c("b_dprime_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  roc  <- roc_sdt(fit, n_points = 100)
  summ <- attr(roc, "summary")
  pts  <- attr(roc, "points")
  hit_on_curve <- stats::approx(summ$FA, summ$Hit_mean, xout = pts$FA_mean)$y
  expect_lt(max(abs(pts$Hit_mean - hit_on_curve)), 0.01)
})

test_that("roc_sdt() criterion_points = FALSE disables multi-criteria points", {
  fit <- fake_binary_fit(uv = TRUE, multi = TRUE)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(
      dprime = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8), sdratio = log(1.3))),
    ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
    variables = function(...) c("bsp_sdratio"),
    .package = "brms"
  )
  roc <- roc_sdt(fit, n_points = 10, criterion_points = FALSE)
  expect_null(attr(roc, "points"))
  expect_true("condition" %in% names(roc))   # one curve per condition instead
})


############################################################################# !
# ROC — RATING                                                           ####
############################################################################# !

test_that("roc_sdt() rating returns K+1 points per draw with endpoints", {
  fit <- fake_rating_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  roc <- roc_sdt(fit)
  expect_s3_class(roc, "bmm_sdt_roc")
  expect_true(isTRUE(attr(roc, "is_rating")))
  expect_equal(nrow(roc), length(unique(roc$.draw)) * (6L + 1L))
  expect_true(all(roc$FA >= 0 & roc$FA <= 1 & roc$Hit >= 0 & roc$Hit <= 1))

  d1 <- roc[roc$.draw == 1L, ]
  expect_true(any(abs(d1$FA - 1) < 1e-8 & abs(d1$Hit - 1) < 1e-8))
  expect_true(any(abs(d1$FA) < 1e-8 & abs(d1$Hit) < 1e-8))
})


############################################################################# !
# AUC                                                                    ####
############################################################################# !

test_that("auc_sdt() binary normal EV uses the analytical Phi(d/sqrt(2))", {
  fit <- fake_binary_fit()
  dpr <- 1.5
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = dpr, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_dprime_Intercept"),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  expect_s3_class(auc, "bmm_sdt_auc")
  expect_equal(mean(auc$AUC), stats::pnorm(dpr / sqrt(2)), tolerance = 1e-10)
})

test_that("auc_sdt() rating uses the numerical path and stays in (0.5, 1)", {
  fit <- fake_rating_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  expect_true(all(auc$AUC > 0.5 & auc$AUC < 1))
})


############################################################################# !
# OBSERVED ROC                                                           ####
############################################################################# !

test_that("roc_observed() rating returns K+1 rows with endpoints", {
  fit <- fake_rating_fit()
  obs <- roc_observed(fit)
  expect_s3_class(obs, "bmm_sdt_roc_observed")
  expect_equal(nrow(obs), 6L + 1L)
  expect_true(all(obs$FA >= 0 & obs$FA <= 1 & obs$Hit >= 0 & obs$Hit <= 1))
})

test_that("roc_observed() binary gives one operating point per criterion", {
  fit <- fake_binary_fit(multi = TRUE)
  obs <- roc_observed(fit, conditions = "condition")
  expect_s3_class(obs, "bmm_sdt_roc_observed")
  expect_equal(nrow(obs), 5L)
  expect_true("condition" %in% names(obs))
  expect_true(all(obs$Hit >= obs$FA))   # positive sensitivity
})


############################################################################# !
# LATENT DECISION VARIABLE                                               ####
############################################################################# !

test_that("latent_sdt() errors for criterion-free models", {
  fake <- function(cls) {
    structure(list(bmm = list(model = structure(list(),
              class = c("bmmodel", "sdt", cls)))), class = c("bmmfit", "brmsfit"))
  }
  expect_error(latent_sdt(list()), "bmmfit")
  expect_error(latent_sdt(fake("sdt_mafc")), "not defined")
  expect_error(latent_sdt(fake("sdt_ranking")), "not defined")
})

test_that("latent_sdt() binary returns noise/signal densities + criterion line", {
  fit <- fake_binary_fit()
  dp  <- 1.5
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = dp, criterion = 0.3)),
    ranef = function(...) list(),
    variables = function(...) "b_dprime_Intercept",
    .package = "brms"
  )
  lat <- latent_sdt(fit, n_grid = 200L)
  expect_s3_class(lat, "bmm_sdt_latent")
  expect_setequal(unique(lat$distribution), c("noise", "signal"))
  expect_equal(nrow(lat), 2L * 200L)
  expect_true(all(lat$density >= 0))

  # densities are the standard normal at -d/2 (noise) and +d/2 (signal, SD 1)
  noise  <- lat[lat$distribution == "noise", ]
  signal <- lat[lat$distribution == "signal", ]
  expect_equal(noise$density,  dnorm(noise$x, -dp / 2), tolerance = 1e-10)
  expect_equal(signal$density, dnorm(signal$x, dp / 2), tolerance = 1e-10)

  ln <- attr(lat, "lines")
  expect_equal(nrow(ln), 1L)
  expect_equal(ln$position, 0.3)
  expect_equal(ln$marker, "criterion")

  # area beyond the criterion reproduces the model false-alarm rate
  keep <- noise$x > ln$position
  area <- sum(diff(noise$x[keep]) *
              (utils::head(noise$density[keep], -1L) +
               utils::tail(noise$density[keep], -1L)) / 2)
  expect_equal(area, pnorm(-dp / 2 - ln$position), tolerance = 0.02)
})

test_that("latent_sdt() rating returns ordered threshold lines per condition", {
  fit <- fake_rating_fit(n_ratings = 6L)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  lat <- latent_sdt(fit, n_grid = 100L)
  ln  <- attr(lat, "lines")
  expect_equal(nrow(ln), 5L)
  expect_equal(ln$marker, paste0("t", 1:5))
  expect_true(all(diff(ln$position) > 0))
})

test_that("latent_sdt() UV widens the signal distribution by exp(sdratio)", {
  fit  <- fake_binary_fit(uv = TRUE)
  sdr  <- log(1.4)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.2, criterion = 0,
                                                  sdratio = sdr)),
    ranef = function(...) list(),
    variables = function(...) c("b_dprime_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  lat    <- latent_sdt(fit, n_grid = 200L)
  signal <- lat[lat$distribution == "signal", ]
  expect_equal(signal$density, dnorm(signal$x, 1.2 / 2, exp(sdr)),
               tolerance = 1e-10)
})


############################################################################# !
# SUMMARIES & PRINT                                                      ####
############################################################################# !

test_that(".roc_sdt_summary collapses draws to a quantile band per FA", {
  df <- data.frame(
    FA    = rep(c(0, 0.5, 1), 4L),
    Hit   = c(0, 0.40, 1, 0, 0.60, 1, 0, 0.50, 1, 0, 0.55, 1),
    .draw = rep(1:4, each = 3L)
  )
  s <- .roc_sdt_summary(df, c(0.025, 0.975))
  expect_equal(nrow(s), 3L)
  expect_true(all(c("FA", "Hit_mean", "Hit_lower", "Hit_upper") %in% names(s)))
  expect_true(all(s$Hit_lower <= s$Hit_mean & s$Hit_mean <= s$Hit_upper))
})

test_that(".auc_sdt_summary reports the posterior mean and band", {
  df <- data.frame(AUC = c(0.70, 0.75, 0.80, 0.72), .draw = 1:4)
  s <- .auc_sdt_summary(df, c(0.025, 0.975))
  expect_true(all(c("AUC_mean", "AUC_lower", "AUC_upper") %in% names(s)))
  expect_equal(s$AUC_mean, mean(df$AUC))
})

test_that("print methods emit a short header", {
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  expect_output(print(roc_sdt(fit, n_points = 10)), "SDT ROC curve")
  expect_output(print(auc_sdt(fit)), "SDT AUC")
  expect_output(print(latent_sdt(fit, n_grid = 20)), "SDT latent distributions")
})
