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

test_that("model-implied functions work when fit has no group-level effects", {
  fit_binary <- fake_binary_fit()
  fit_rating <- fake_rating_fit()
  local_mocked_bindings(
    ranef = function(...) stop("The model does not contain group-level effects."),
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0)),
    variables = function(...) c("b_d_Intercept", "b_criterion_Intercept"),
    .package = "brms"
  )
  expect_no_error(roc_sdt(fit_binary))
  expect_no_error(auc_sdt(fit_binary))
  expect_no_error(latent_sdt(fit_binary))
  expect_no_error(sdt_sensitivity(fit_binary))
  expect_no_error(sdt_thresholds(fit_rating))
})

test_that("nested interaction groupings (id:session) are excluded from conditions", {
  data <- data.frame(
    stimulus = c(0, 1), n_old = c(20L, 80L), n_trials = 100L, dist_type = 1L,
    id = c(1L, 2L), session = c(1L, 1L)
  )
  data$`id:session` <- c("1_1", "2_1")
  fit <- structure(
    list(
      data = data,
      bmm = list(
        model = sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials"),
        user_formula = bmf(d ~ 1, criterion ~ 1)
      )
    ),
    class = c("bmmfit", "brmsfit")
  )
  local_mocked_bindings(
    ranef = function(...) list(`id:session` = array(0, dim = c(2, 1, 1))),
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0)),
    variables = function(...) c("b_d_Intercept", "b_criterion_Intercept"),
    .package = "brms"
  )
  conditions <- .sdt_resolve_conditions(fit, NULL)
  expect_equal(nrow(conditions), 1L)
  expect_false(any(c("id", "session", "id:session") %in% names(conditions)))
})

sdt_entry_calls <- function(fit_binary, fit_rating) {
  list(
    roc_sdt         = function(...) roc_sdt(fit_binary, ...),
    auc_sdt         = function(...) auc_sdt(fit_binary, ...),
    latent_sdt      = function(...) latent_sdt(fit_binary, ...),
    sdt_sensitivity = function(...) sdt_sensitivity(fit_binary, ...),
    sdt_thresholds  = function(...) sdt_thresholds(fit_rating, ...)
  )
}

test_that("ndraws is refused before any posterior draw is taken", {
  n_calls <- 0L
  local_mocked_bindings(
    posterior_linpred = function(...) {
      n_calls <<- n_calls + 1L
      stop("posterior_linpred should not be reached")
    },
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  calls <- sdt_entry_calls(fake_binary_fit(uv = TRUE), fake_rating_fit(uv = TRUE))
  for (nm in names(calls)) {
    expect_error(calls[[nm]](ndraws = 10), "draw_ids", info = nm)
  }
  expect_identical(n_calls, 0L)
})

test_that("draw_ids reaches every posterior_linpred call unchanged", {
  seen <- list()
  base <- mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0,
                                    sdratio = 0.3))
  local_mocked_bindings(
    posterior_linpred = function(object, ..., draw_ids = NULL) {
      seen[[length(seen) + 1L]] <<- draw_ids
      base(object, ...)
    },
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  sdt_sensitivity(fake_binary_fit(uv = TRUE), draw_ids = 1:5)
  roc_sdt(fake_rating_fit(uv = TRUE), draw_ids = 1:5)
  expect_gt(length(seen), 4L)
  for (ids in seen) expect_identical(ids, 1:5)
})

test_that("conditions must be a data frame of columns in the data", {
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  calls <- sdt_entry_calls(fake_binary_fit(), fake_rating_fit())
  for (nm in names(calls)) {
    expect_error(calls[[nm]](conditions = "stimulus"), "must be a data frame",
                 info = nm)
    expect_error(calls[[nm]](conditions = data.frame(base_rate = "br1")),
                 "not in the data:\\s+'base_rate'", info = nm)
  }
})

test_that("roc_observed() takes conditions as names of data columns", {
  fit <- fake_binary_fit(multi = TRUE)
  expect_error(roc_observed(fit, conditions = data.frame(condition = "br1")),
               "character vector")
  expect_error(roc_observed(fit, conditions = "base_rate"),
               "not in the data:\\s+'base_rate'")
  expect_equal(nrow(roc_observed(fit, conditions = "condition")), 5L)
})

test_that("probs must be two increasing probabilities", {
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  calls <- sdt_entry_calls(fake_binary_fit(), fake_rating_fit())
  bad <- list(c(0.975, 0.025), 0.5, c(-0.1, 0.9), "a", c(0.1, NA))
  for (nm in names(calls)) {
    for (p in bad) {
      expect_error(calls[[nm]](probs = p), "probs must be",
                   info = paste(nm, deparse(p)))
    }
  }
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

test_that(".sdt_criterion_point_dims treats d predictors as curves", {
  fit <- fake_binary_fit()
  fit$bmm$user_formula <- bmf(d ~ 0 + condition, criterion ~ 0 + condition,
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
  m_ev <- sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  expect_false(.sdt_has_estimated_sdratio(m_ev))
  m_uv <- m_ev
  m_uv$fixed_parameters$sdratio <- NULL
  expect_true(.sdt_has_estimated_sdratio(m_uv))
})

test_that(".sdt_has_estimated_sdratio cross-checks brms::variables", {
  m_ev <- sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  fit <- structure(list(), class = c("bmmfit", "brmsfit"))
  local_mocked_bindings(
    variables = function(...) c("b_d_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  expect_true(.sdt_has_estimated_sdratio(m_ev, fit))
})

test_that(".sdt_unit_sdratio is model-only and keys off the natural-scale value", {
  m_ev <- sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  expect_true(.sdt_unit_sdratio(m_ev))

  m_fixed_int <- m_ev
  m_fixed_int$fixed_parameters$sdratio <- 0L
  expect_true(.sdt_unit_sdratio(m_fixed_int))

  m_fixed <- m_ev
  m_fixed$fixed_parameters$sdratio <- 0.3
  expect_false(.sdt_unit_sdratio(m_fixed))

  m_uv <- m_ev
  m_uv$fixed_parameters$sdratio <- NULL
  expect_false(.sdt_unit_sdratio(m_uv))

  expect_true(.sdt_unit_sdratio(fake_mafc_fit()$bmm$model))
  expect_true(.sdt_unit_sdratio(fake_ranking_fit()$bmm$model))
})


############################################################################# !
# ROC — BINARY                                                           ####
############################################################################# !

test_that("roc_sdt() binary single criterion returns a smooth curve", {
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept"),
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
      d = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8), sdratio = log(1.3))),
    ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
    variables = function(...) c("b_d_Intercept", "bsp_sdratio"),
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
  spread <- list(d = 0.18, criterion = 0.10, sdratio = 0.12)
  base   <- list(d = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8),
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
    variables = function(...) c("b_d_Intercept", "bsp_sdratio"),
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
      d = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8), sdratio = log(1.3))),
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
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0)),
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

test_that("roc_sdt() rating attaches a smooth implied curve + threshold points", {
  fit <- fake_rating_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  roc  <- roc_sdt(fit, n_points = 50)
  summ <- attr(roc, "summary")
  pts  <- attr(roc, "points")

  expect_true(all(c("FA", "Hit_mean", "Hit_lower", "Hit_upper") %in% names(summ)))
  expect_equal(nrow(summ), 50L + 2L)
  expect_equal(summ$FA[1L], 0)
  expect_equal(utils::tail(summ$FA, 1L), 1)
  expect_false(is.unsorted(summ$FA))
  expect_true(all(summ$Hit_mean >= 0 & summ$Hit_mean <= 1))

  expect_equal(nrow(pts), 5L)
  expect_s3_class(pts$threshold, "factor")
  expect_equal(levels(pts$threshold), paste0("t", 1:5))
  expect_true(all(c("FA_mean", "FA_lower", "FA_upper",
                    "Hit_mean", "Hit_lower", "Hit_upper") %in% names(pts)))
})

test_that("roc_sdt() rating threshold points fall on the model-implied ROC", {
  # Checked against the ROC evaluated at each point's own FA rather than against
  # a linear interpolation of `summary`: a threshold can land at an FA below the
  # curve's grid floor (0.001) -- for gumbel_min the ROC is already at Hit ~ 0.1
  # by FA = 1e-5 -- where interpolating from the (0, 0) endpoint is inaccurate
  # no matter how large n_points is.
  base <- list(d = 1.4, criterion = 0.1, spacing = 0, sdratio = log(1.2))
  z    <- stats::qnorm(stats::ppoints(n_draws_mock))

  roc_points <- function(dist, spread) {
    fit <- fake_rating_fit(uv = TRUE)
    fit$bmm$model$other_vars$dist <- dist
    local_mocked_bindings(
      posterior_linpred = function(object, dpar = NULL, nlpar = NULL,
                                   newdata = NULL, ...) {
        par <- if (!is.null(dpar)) dpar else nlpar
        n_cond <- if (is.null(newdata)) 1L else nrow(newdata)
        matrix(rep(rep_len(base[[par]], n_cond), each = n_draws_mock) +
                 spread[[par]] * z, nrow = n_draws_mock)
      },
      ranef = function(...) list(),
      variables = function(...) c("bsp_sdratio"),
      .package = "brms"
    )
    roc <- roc_sdt(fit, n_points = 200)
    list(points = attr(roc, "points"), summary = attr(roc, "summary"))
  }

  # FA = 1 - cdf(t + sep/2) and Hit = 1 - cdf((t - sep/2) / sigma), so
  # eliminating t gives Hit as a function of FA. `d` is d_a, so the separation
  # on the noise-standardized axis is d * sqrt((1 + sigma^2) / 2).
  hit_on_roc <- function(dist, fa) {
    sigma <- exp(base$sdratio)
    sep   <- base$d * sqrt((1 + sigma^2) / 2)
    1 - .sdt_dists[[dist]]$cdf(
      (.sdt_dists[[dist]]$qf(1 - fa) - sep) / sigma
    )
  }

  flat <- list(d = 0, criterion = 0, spacing = 0, sdratio = 0)
  wide <- list(d = 0.15, criterion = 0.10, spacing = 0.05, sdratio = 0.10)

  for (dist in names(.sdt_dists)) {
    # With no posterior spread the operating points must sit on the ROC exactly.
    # This is the sharp check: a mirrored distribution convention breaks it by
    # >0.05 even though every marginal probability still looks plausible.
    flat_res <- roc_points(dist, flat)
    expect_equal(flat_res$points$Hit_mean,
                 hit_on_roc(dist, flat_res$points$FA_mean),
                 tolerance = 1e-8, info = dist)

    # With spread the points are posterior means of a non-linear map, so they
    # sit slightly off the mean-parameter curve. The gap is largest for
    # gumbel_min, whose ROC is steepest near the origin (0.063 here; the other
    # three stay under 0.023). Reporting d_a widens it: the separation is
    # d * sqrt((1 + sigma^2)/2), so the draws pass through a second non-linear
    # step before reaching the ROC. This is a slack bound on that Jensen gap,
    # not an accuracy claim -- the sharp check is the flat case above.
    wide_res <- roc_points(dist, wide)
    expect_lt(max(abs(wide_res$points$Hit_mean -
                        hit_on_roc(dist, wide_res$points$FA_mean))), 0.07)

    expect_true(all(diff(wide_res$summary$FA) >= 0), info = dist)
    expect_equal(range(wide_res$summary$FA), c(0, 1), info = dist)
  }
})


############################################################################# !
# AUC                                                                    ####
############################################################################# !

test_that("auc_sdt() binary normal EV uses the analytical Phi(d/sqrt(2))", {
  fit <- fake_binary_fit()
  dpr <- 1.5
  local_mocked_bindings(
    # A real EV fit's sdratio is fixed at 0 via a constant prior, not omitted,
    # so brms::variables() lists it just like any other parameter.
    posterior_linpred = mock_linpred_factory(list(d = dpr, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_criterion_Intercept",
                                "b_sdratio_Intercept"),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  expect_s3_class(auc, "bmm_sdt_auc")
  expect_equal(mean(auc$AUC), stats::pnorm(dpr / sqrt(2)), tolerance = 1e-10)
})

test_that("auc_sdt() binary gumbel_max EV uses the analytical plogis(d)", {
  fit <- fake_binary_fit()
  fit$bmm$model$other_vars$dist <- "gumbel_max"
  dpr <- 0.8
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = dpr, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_criterion_Intercept",
                                "b_sdratio_Intercept"),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  expect_equal(mean(auc$AUC), stats::plogis(dpr), tolerance = 1e-10)
})

test_that("auc_sdt() does not take the closed form when sdratio is fixed away from 0", {
  fit_fixed <- fake_binary_fit()
  fit_fixed$bmm$model$fixed_parameters$sdratio <- 0.3
  fit_uv <- fake_binary_fit(uv = TRUE)
  dpr <- 1.5
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = dpr, criterion = 0,
                                                  sdratio = 0.3)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_criterion_Intercept",
                                "b_sdratio_Intercept"),
    .package = "brms"
  )
  auc_fixed <- auc_sdt(fit_fixed)
  auc_uv    <- auc_sdt(fit_uv)
  closed    <- stats::pnorm(dpr / sqrt(2))

  expect_true(all(abs(auc_fixed$AUC - closed) > 1e-6))
  expect_equal(auc_fixed$AUC, auc_uv$AUC, tolerance = 1e-12)
})

test_that("auc_sdt() rating uses the numerical path and stays in (0.5, 1)", {
  fit <- fake_rating_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  expect_true(all(auc$AUC > 0.5 & auc$AUC < 1))
})

test_that("auc_sdt() rating EV matches the closed-form Phi(d_a/sqrt(2)) oracle", {
  fit <- fake_rating_fit()
  d_true <- 1.5
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = d_true, sdratio = 0,
                                                  criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  expect_equal(mean(auc$AUC), stats::pnorm(d_true / sqrt(2)), tolerance = 1e-3)
})

test_that("auc_sdt() rating UV (positive sdratio) matches the Phi(d_a/sqrt(2)) oracle", {
  fit <- fake_rating_fit(uv = TRUE)
  d_true <- 1.5
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = d_true, sdratio = log(1.35),
                                                  criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  expect_equal(mean(auc$AUC), stats::pnorm(d_true / sqrt(2)), tolerance = 1e-3)
})

test_that("auc_sdt() rating UV (negative sdratio) matches the Phi(d_a/sqrt(2)) oracle", {
  fit <- fake_rating_fit(uv = TRUE)
  d_true <- 1.5
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = d_true, sdratio = -0.5,
                                                  criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  expect_equal(mean(auc$AUC), stats::pnorm(d_true / sqrt(2)), tolerance = 1e-3)
})

test_that("auc_sdt() rating gumbel_max UV matches an independent stats::integrate() oracle", {
  fit <- fake_rating_fit(uv = TRUE)
  fit$bmm$model$other_vars$dist <- "gumbel_max"
  d_true <- 1.2
  sdratio_log <- 0.2
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = d_true, sdratio = sdratio_log,
                                                  criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  auc <- auc_sdt(fit)

  sdratio <- exp(sdratio_log)
  sep <- d_true * sqrt((1 + sdratio^2) / 2)
  gumbel_max_cdf <- function(x) exp(-exp(-x))
  gumbel_max_qf  <- function(p) -log(-log(p))
  hit_of_fa <- function(fa) {
    t <- gumbel_max_qf(1 - fa) - sep / 2
    1 - gumbel_max_cdf((t - sep / 2) / sdratio)
  }
  oracle <- stats::integrate(hit_of_fa, lower = 1e-8, upper = 1 - 1e-8,
                             rel.tol = 1e-10, subdivisions = 1000L)$value

  expect_equal(mean(auc$AUC), oracle, tolerance = 1e-3)
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

test_that("latent_sdt() errors for non-bmmfit input", {
  expect_error(latent_sdt(list()), "bmmfit")
})

test_that("latent_sdt() draws densities without boundary lines (mafc, ranking)", {
  for (maker in list(fake_mafc_fit, fake_ranking_fit)) {
    fit <- maker()
    local_mocked_bindings(
      posterior_linpred = mock_linpred_factory(list(d = 1.4)),
      ranef = function(...) list(),
      variables = function(...) character(0),
      .package = "brms"
    )
    lat <- latent_sdt(fit, n_grid = 120L)
    expect_s3_class(lat, "bmm_sdt_latent")
    expect_setequal(unique(lat$distribution), c("noise", "signal"))
    expect_equal(nrow(lat), 2L * 120L)
    expect_true(all(lat$density >= 0))
    expect_null(attr(lat, "lines"))            # no criterion -> no boundary lines
  }
})

test_that("latent_sdt() binary returns noise/signal densities + criterion line", {
  fit <- fake_binary_fit()
  dp  <- 1.5
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = dp, criterion = 0.3)),
    ranef = function(...) list(),
    variables = function(...) "b_d_Intercept",
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
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  lat <- latent_sdt(fit, n_grid = 100L)
  ln  <- attr(lat, "lines")
  expect_equal(nrow(ln), 5L)
  expect_equal(ln$marker, paste0("t", 1:5))
  expect_true(all(diff(ln$position) > 0))
  expect_equal(as.character(ln$level), paste0("t", 1:5))   # colour key = threshold
})

test_that("latent_sdt() collapses criterion-only dimensions into one panel", {
  fit <- fake_binary_fit(uv = TRUE, multi = TRUE)   # criterion ~ 0 + condition
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(
      d = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8), sdratio = log(1.3))),
    ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
    variables = function(...) "bsp_sdratio",
    .package = "brms"
  )
  lat <- latent_sdt(fit, n_grid = 60L)
  # one density panel, five criterion lines colour-coded by condition
  expect_equal(nrow(attr(lat, "conditions")), 1L)
  expect_equal(nrow(lat), 2L * 60L)
  ln <- attr(lat, "lines")
  expect_equal(nrow(ln), 5L)
  expect_setequal(as.character(ln$level), paste0("br", 1:5))

  # collapse = FALSE restores one panel per criterion level
  lat2 <- latent_sdt(fit, n_grid = 60L, collapse = FALSE)
  expect_equal(nrow(attr(lat2, "conditions")), 5L)
})

test_that("latent_sdt() show_competitors overlays max-of-distractors per set size", {
  fit <- fake_ranking_fit()                          # constant m = 3
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.4)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  lat <- latent_sdt(fit, n_grid = 80L, show_competitors = TRUE)
  comp <- attr(lat, "competitors")
  expect_false(is.null(comp))
  expect_setequal(levels(comp$set_size), "3")
  expect_true(all(comp$density >= 0))
  # competitor density integrates to ~1 (it is a proper order-statistic density)
  cd <- comp[order(comp$x), ]
  area <- sum(diff(cd$x) * (utils::head(cd$density, -1L) + utils::tail(cd$density, -1L)) / 2)
  expect_equal(area, 1, tolerance = 0.02)

  # off by default and ignored for binary
  expect_null(attr(latent_sdt(fit, n_grid = 40L), "competitors"))
})

test_that("latent_sdt() UV widens the signal distribution by exp(sdratio)", {
  fit  <- fake_binary_fit(uv = TRUE)
  sdr  <- log(1.4)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.2, criterion = 0,
                                                  sdratio = sdr)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  lat    <- latent_sdt(fit, n_grid = 200L)
  signal <- lat[lat$distribution == "signal", ]
  noise  <- lat[lat$distribution == "noise", ]
  # d is d_a, so the densities sit at +/- sep/2 with sep the separation in
  # noise-SD units -- not at +/- d/2, which would understate it whenever
  # sdratio > 0.
  sep <- 1.2 * sqrt((1 + exp(sdr)^2) / 2)
  expect_equal(signal$density, dnorm(signal$x, sep / 2, exp(sdr)),
               tolerance = 1e-10)
  expect_equal(noise$density, dnorm(noise$x, -sep / 2, 1), tolerance = 1e-10)
})

test_that("latent_sdt() density separation is d_a in RMS-SD units", {
  # the defining property of d_a: the gap between the two density modes,
  # divided by the root-mean-square of the two SDs, returns the d parameter
  fit <- fake_binary_fit(uv = TRUE)
  sdr <- log(1.7)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.2, criterion = 0,
                                                  sdratio = sdr)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  lat  <- latent_sdt(fit, n_grid = 4001L)
  mode_of <- function(which) {
    sub <- lat[lat$distribution == which, ]
    sub$x[which.max(sub$density)]
  }
  sep <- mode_of("signal") - mode_of("noise")
  expect_equal(sep / sqrt((1 + exp(sdr)^2) / 2), 1.2, tolerance = 1e-2)
})


############################################################################# !
# THRESHOLDS                                                             ####
############################################################################# !

test_that("sdt_thresholds() returns per-draw samples + a K-1 summary", {
  fit <- fake_rating_fit(n_ratings = 6L)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  thr <- sdt_thresholds(fit)
  expect_s3_class(thr, "bmm_sdt_thresholds")
  expect_true(all(c("marker", "position", ".draw") %in% names(thr)))
  expect_equal(length(unique(thr$.draw)), n_draws_mock)
  expect_setequal(unique(thr$marker), paste0("t", 1:5))
  expect_equal(nrow(thr), n_draws_mock * 5L)

  s <- attr(thr, "summary")
  expect_equal(nrow(s), 5L)
  expect_true(all(c("marker", "position", "lower", "upper") %in% names(s)))

  d1 <- thr[thr$.draw == 1L, ]
  d1 <- d1[order(match(d1$marker, paste0("t", 1:5))), ]
  expect_true(all(diff(d1$position) > 0))                # ordered within a draw

  expect_output(print(thr), "SDT decision thresholds")
})

test_that("sdt_thresholds() errors for fits without rating thresholds", {
  expect_error(sdt_thresholds(list()), "bmmfit")
  expect_error(sdt_thresholds(fake_binary_fit()), "rating")
  expect_error(sdt_thresholds(fake_mafc_fit()), "rating")
  expect_error(sdt_thresholds(fake_ranking_fit()), "rating")
})

test_that("sdt_thresholds() summary matches latent_sdt() lines across parameterizations", {
  specs <- list(
    list(tt = "equidistant",  K = 6L),   # spacing-based, even K
    list(tt = "softmax",      K = 5L),   # spacing + deltas, odd K
    list(tt = "log_ratio",    K = 6L),   # delta-based, even K
    list(tt = "log_distance", K = 5L)    # delta-based, odd K
  )
  for (s in specs) {
    fit  <- fake_rating_fit(threshold_type = s$tt, n_ratings = s$K)
    pars <- names(fit$bmm$model$parameters)
    draws <- list(d = 1.5, criterion = 0.2)
    if ("spacing" %in% pars) draws$spacing <- 0
    for (d in grep("^delta", pars, value = TRUE)) draws[[d]] <- 0

    local_mocked_bindings(
      posterior_linpred = mock_linpred_factory(draws),
      ranef = function(...) list(),
      variables = function(...) character(0),
      .package = "brms"
    )
    summ  <- attr(sdt_thresholds(fit), "summary")
    lines <- attr(latent_sdt(fit, n_grid = 20L), "lines")
    info  <- paste(s$tt, s$K)
    expect_equal(nrow(summ), s$K - 1L, info = info)
    expect_equal(summ$position, lines$position, info = info)
    expect_true(all(diff(summ$position) > 0), info = info)
  }
})


test_that("sdt_thresholds() places criterion by the rating model's rule", {
  # even K: criterion is the middle threshold; odd K: the centre of the middle
  # category, whatever the threshold_type
  for (tt in c("parsimonious", "equidistant", "softmax", "log_ratio", "log_distance")) {
    for (K in 5:6) {
      fit  <- fake_rating_fit(threshold_type = tt, n_ratings = K)
      pars <- names(fit$bmm$model$parameters)
      draws <- list(d = 1.5, criterion = 0.2)
      if ("spacing" %in% pars) draws$spacing <- -0.4
      deltas <- grep("^delta", pars, value = TRUE)
      draws[deltas] <- as.list(seq(-0.3, 0.3, length.out = length(deltas)))

      local_mocked_bindings(
        posterior_linpred = mock_linpred_factory(draws),
        ranef = function(...) list(),
        variables = function(...) character(0),
        .package = "brms"
      )
      pos <- attr(sdt_thresholds(fit), "summary")$position
      centre <- if (K %% 2L == 0L) pos[K / 2] else mean(pos[(K - 1) / 2 + 0:1])
      expect_equal(centre, 0.2, info = paste(tt, K))
    }
  }
})


############################################################################# !
# LINKS                                                                  ####
############################################################################# !

# The kernels receive d and criterion through the model's inverse link, so the
# post-processing must read them on that scale too, not as linear predictors.

test_that("rating post-processing applies the links on d and criterion", {
  fit <- fake_rating_fit(uv = TRUE, links = list(d = "log", criterion = "softplus"))
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(
      d = log(1.5), criterion = log(expm1(0.2)), spacing = 0, sdratio = 0
    )),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  da <- attr(sdt_sensitivity(fit, measure = "da"), "summary")$mean
  expect_equal(da, 1.5)

  pos <- attr(sdt_thresholds(fit), "summary")$position
  expect_equal(pos, .sdt_make_thresholds(0.2, 6L, "parsimonious", spacing = 0))

  pts <- attr(roc_sdt(fit, n_points = 10), "points")
  fa  <- 1 - cumsum(.sdt_category_probs(pos, 1.5, 1, 0L, "normal"))[1:5]
  expect_equal(pts$FA_mean, fa)
})

test_that("binary post-processing applies the links on d and criterion", {
  fit <- fake_binary_fit(links = list(d = "log", criterion = "softplus"))
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(
      d = log(1.5), criterion = log(expm1(0.2))
    )),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  auc <- attr(auc_sdt(fit), "summary")$AUC_mean
  expect_equal(auc, stats::pnorm(1.5 / sqrt(2)))

  lines <- attr(latent_sdt(fit, n_grid = 20L), "lines")
  expect_equal(lines$position, 0.2)
})


############################################################################# !
# SUMMARIES & PRINT                                                      ####
############################################################################# !

test_that(".auc_sdt_summary reports the posterior mean and band", {
  df <- data.frame(AUC = c(0.70, 0.75, 0.80, 0.72), .draw = 1:4)
  s <- .auc_sdt_summary(df, c(0.025, 0.975))
  expect_true(all(c("AUC_mean", "AUC_lower", "AUC_upper") %in% names(s)))
  expect_equal(s$AUC_mean, mean(df$AUC))
})

test_that("print methods emit a short header", {
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  expect_output(print(roc_sdt(fit, n_points = 10)), "SDT ROC curve")
  expect_output(print(auc_sdt(fit)), "SDT AUC")
  expect_output(print(latent_sdt(fit, n_grid = 20)), "SDT latent distributions")
})


############################################################################# !
# SENSITIVITY                                                            ####
############################################################################# !

test_that("sdt_sensitivity() converts d_a to the noise and signal scales", {
  fit <- fake_binary_fit(uv = TRUE)
  sdr <- log(1.6)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.3, criterion = 0,
                                                  sdratio = sdr)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  res <- sdt_sensitivity(fit)
  val <- function(m) unique(res$value[res$measure == m])

  # first principles: separation = d_a * rms scale; divide by each reference SD
  sigma <- exp(sdr)
  sep   <- 1.3 * sqrt((1 + sigma^2) / 2)
  expect_equal(val("da"), 1.3, tolerance = 1e-12)
  expect_equal(val("dn"), sep, tolerance = 1e-12)
  expect_equal(val("ds"), sep / sigma, tolerance = 1e-12)

  # the ratio of the two reference scales is the SD ratio itself
  expect_equal(val("dn") / val("ds"), sigma, tolerance = 1e-12)
  # and d_a orders between them whenever the SDs differ
  expect_true(val("ds") < val("da") && val("da") < val("dn"))
})

test_that("sdt_sensitivity() collapses to one value under equal variance", {
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  res <- sdt_sensitivity(fit)
  expect_equal(unique(res$value), 1.5, tolerance = 1e-12)
  expect_setequal(unique(res$measure), c("da", "dn", "ds"))
})

test_that("sdt_sensitivity() converts per draw, not from point estimates", {
  # d and sdratio are correlated across draws here, so a conversion applied to
  # the posterior means would give a different mean than converting draw-wise
  fit <- fake_binary_fit(uv = TRUE)
  z   <- stats::qnorm(stats::ppoints(n_draws_mock))
  local_mocked_bindings(
    posterior_linpred = function(object, dpar = NULL, nlpar = NULL,
                                 newdata = NULL, ...) {
      par <- if (!is.null(dpar)) dpar else nlpar
      base <- c(d = 1.3, criterion = 0, sdratio = log(1.6))[[par]]
      spread <- c(d = 0.4, criterion = 0, sdratio = 0.35)[[par]]
      matrix(base + spread * z, nrow = n_draws_mock)
    },
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  res <- sdt_sensitivity(fit, measure = "dn")

  d_draws     <- 1.3 + 0.4 * z
  sigma_draws <- exp(log(1.6) + 0.35 * z)
  expect_equal(res$value, d_draws * sqrt((1 + sigma_draws^2) / 2),
               tolerance = 1e-12)
  # the plug-in shortcut would be biased -- confirm the two really differ
  plug_in <- mean(d_draws) * sqrt((1 + mean(sigma_draws)^2) / 2)
  expect_gt(abs(mean(res$value) - plug_in), 1e-3)
})

test_that("sdt_sensitivity() honours measure selection and summarises draws", {
  fit <- fake_binary_fit(uv = TRUE)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.3, criterion = 0,
                                                  sdratio = log(1.6))),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "bsp_sdratio"),
    .package = "brms"
  )
  res <- sdt_sensitivity(fit, measure = c("dn", "ds"))
  expect_setequal(unique(res$measure), c("dn", "ds"))
  expect_equal(nrow(res), 2L * n_draws_mock)

  summ <- attr(res, "summary")
  expect_true(all(c("measure", "mean", "lower", "upper") %in% names(summ)))
  expect_equal(summ$mean[summ$measure == "dn"],
               mean(res$value[res$measure == "dn"]), tolerance = 1e-12)
  expect_error(sdt_sensitivity(fit, measure = "dprime"), "arg")
})

test_that("sdt_sensitivity() rejects non-SDT and non-bmmfit input", {
  expect_error(sdt_sensitivity(list()), "bmmfit")
})

test_that("sdt_sensitivity() print method labels the three scales", {
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  expect_output(print(sdt_sensitivity(fit)), "SDT sensitivity")
})


############################################################################# !
# SUMMARY NOTES                                                          ####
############################################################################# !

test_that("summary_notes.sdt() fires whenever sdratio departs from 0, not just when estimated", {
  m_ev <- fake_binary_fit()$bmm$model
  expect_null(summary_notes(m_ev, NULL))

  m_fixed <- m_ev
  m_fixed$fixed_parameters$sdratio <- 0.3
  expect_match(summary_notes(m_fixed, NULL), "d_a")

  m_uv <- fake_binary_fit(uv = TRUE)$bmm$model
  expect_match(summary_notes(m_uv, NULL), "d_a")

  expect_null(summary_notes(fake_mafc_fit()$bmm$model, NULL))
})


############################################################################# !
# PRINT METHOD COUNTS                                                    ####
############################################################################# !

test_that("print.bmm_sdt_roc() counts FA points per curve, not per condition", {
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(
      d = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8), sdratio = log(1.3))),
    ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  one_curve <- capture.output(print(roc_sdt(fake_binary_fit(uv = TRUE, multi = TRUE))))
  expect_true(any(grepl("Smooth curve: 102 FA points per draw", one_curve, fixed = TRUE)))
  expect_false(any(grepl(" x 5 ", one_curve, fixed = TRUE)))

  fit <- fake_binary_fit(uv = TRUE, multi = TRUE)
  fit$bmm$user_formula <- bmf(d ~ 0 + condition, criterion ~ 1, sdratio ~ 1)
  five_curves <- capture.output(print(roc_sdt(fit)))
  expect_true(any(grepl("x 5 curves", five_curves, fixed = TRUE)))
  expect_true(any(grepl("Smooth curve: 102 FA points per draw", five_curves, fixed = TRUE)))
})

test_that("print.bmm_sdt_roc() counts the K-1 rating thresholds", {
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0,
                                                  sdratio = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  out <- capture.output(print(roc_sdt(fake_rating_fit(n_ratings = 6L))))
  expect_true(any(grepl("Rating model: 5 threshold ROC points per draw", out, fixed = TRUE)))
})

test_that("a column subset that lost its attributes prints as a data frame", {
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(d = 1.5, criterion = 0, spacing = 0,
                                                  sdratio = 0)),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  objects <- list(
    roc_sdt(fake_binary_fit(), n_points = 10),
    latent_sdt(fake_binary_fit(), n_grid = 20L),
    auc_sdt(fake_binary_fit()),
    sdt_thresholds(fake_rating_fit()),
    sdt_sensitivity(fake_binary_fit())
  )
  for (obj in objects) {
    sub <- obj[, 1:2]
    expect_null(attr(sub, "model_class"))
    expect_identical(capture.output(print(sub)),
                     capture.output(print.data.frame(sub)),
                     info = class(obj)[1])
  }
})
