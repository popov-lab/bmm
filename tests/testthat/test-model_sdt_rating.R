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

test_that("sdt_rating takes the number of categories from the response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  expect_equal(model$other_vars$n_ratings, 6L)
})

test_that("sdt_rating model rejects fewer than 3 response columns", {
  expect_error(sdt_rating(c("r1", "r2"), "stimulus"), "more than 2")
  expect_error(sdt_rating("rating", "stimulus"), "more than 2")
})

test_that("sdt_rating refuses column names brms cannot turn into parameters", {
  expect_error(sdt_rating(paste0("conf_", 1:4), "stimulus"),
               "must not contain '.' or '_'")
  expect_error(sdt_rating(c("r1", "r.2", "r3"), "stimulus"), "'r.2'")
})

test_that("sdt_rating model has parsimonious threshold params by default", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("spacing" %in% names(model$parameters))
  expect_true("d" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_false("mu" %in% names(model$parameters))
})

test_that("sdt_rating model has log_distance threshold params", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  expect_true(all(c("delta1", "delta2") %in% names(model$parameters)))
  expect_false("delta3" %in% names(model$parameters))
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
  expect_true(all(c("delta1", "delta2") %in% names(model$parameters)))
  expect_true("criterion" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
  expect_false("delta3" %in% names(model$parameters))
  # the roles differ, so the labels must say which delta is the spread
  expect_match(model$parameters$delta2, "spread")
  expect_match(model$parameters$delta1, "ratio")
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
                      links = list(d = "log"))
  expect_equal(model$links$d, "log")
  expect_equal(model$links$criterion, "identity")
})

test_that("sdt_rating refuses a link on sdratio or a threshold parameter", {
  # sdratio and the thresholds are read through exp() in Stan, and sdratio is
  # additionally fixed at 0, so neither survives a change of link
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_equal(settable_links(model), c("d", "criterion"))
  for (par in c("sdratio", "spacing")) {
    expect_error(
      sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                 links = stats::setNames(list("log"), par)),
      paste0("link of '", par, "' cannot be changed")
    )
  }
  expect_error(
    sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
               links = list(sensitivity = "log")),
    "Unrecognized link target"
  )
  model$links$spacing <- "log"
  expect_error(check_links(model), "link of 'spacing' cannot be changed")
})

test_that("sdt_rating offers only the links it can invert in its formula", {
  expect_error(
    sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", links = list(d = "sqrt")),
    "Unknown link function"
  )
})

test_that("sdt_rating gives every parameter an sd default prior", {
  # Ro and Rn are logit-scale probabilities, so they join d at rate 1;
  # everything else is a log-scale quantity and takes rate 2
  rate1 <- c("d", "Ro", "Rn")
  for (v in names(bmm:::.sdt_rating_variants)) {
    for (tt in c("parsimonious", "log_distance", "softmax")) {
      model <- sdt_rating(paste0("r", 1:6), "stimulus",
                          threshold_type = tt, version = v)
      sds <- vapply(model$default_priors, function(p) p$sd %||% NA_character_,
                    character(1))
      info <- paste(v, tt)
      expect_false(anyNA(sds), info = info)
      expect_true(all(sds[intersect(names(sds), rate1)] == "exponential(1)"), info = info)
      expect_true(all(sds[setdiff(names(sds), rate1)] == "exponential(2)"), info = info)
    }
  }
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

test_that("every threshold type places criterion the same way for odd and even K", {
  # even K: criterion is the middle threshold; odd K: the centre of the middle
  # category, with the two thresholds around it half an interval away. The R
  # builder must agree with sdt_place_thresholds_rating() in Stan (parity test
  # below), and the two closed forms (parsimonious, equidistant) with both.
  for (K in 3:8) {
    for (tt in c("parsimonious", "equidistant", "log_distance", "log_ratio", "softmax")) {
      if (tt == "log_ratio" && K < 4) next
      nd <- if (tt == "softmax") max(0L, K - 3L) else K - 2L
      deltas <- if (nd > 0) seq(0.1, 0.3, length.out = nd) else NULL
      thr <- .sdt_make_thresholds(0.3, K, tt, spacing = 0.2, deltas = deltas)
      expect_false(is.unsorted(thr), info = paste(tt, "K =", K))
      if (K %% 2 == 0) {
        expect_equal(thr[K / 2], 0.3, info = paste(tt, "K =", K))
      } else {
        g <- (K - 1) / 2
        expect_equal((thr[g] + thr[g + 1]) / 2, 0.3, info = paste(tt, "K =", K))
      }
    }
  }
  expect_equal(.sdt_make_thresholds(0, 5L, "equidistant", spacing = 0), c(-1.5, -0.5, 0.5, 1.5))
  expect_equal(.sdt_make_thresholds(0, 7L, "equidistant", spacing = 0), seq(-2.5, 2.5))
  expect_equal(.sdt_make_thresholds(0, 6L, "equidistant", spacing = 0), c(-2, -1, 0, 1, 2))
  # with equal intervals the anchored types coincide with equidistant
  for (K in 5:6) {
    expect_equal(.sdt_make_thresholds(0, K, "log_distance", deltas = rep(0, K - 2)),
                 .sdt_make_thresholds(0, K, "equidistant", spacing = 0))
    expect_equal(.sdt_make_thresholds(0, K, "log_ratio", deltas = rep(0, K - 2)),
                 .sdt_make_thresholds(0, K, "equidistant", spacing = 0))
    expect_equal(.sdt_make_thresholds(0, K, "softmax", spacing = 0, deltas = rep(0, K - 3)),
                 .sdt_make_thresholds(0, K, "equidistant", spacing = 0))
  }
})

test_that("log_ratio deltas act as spread and ratios, in interval order", {
  # even K = 6: delta3 is the spread (interval above the criterion threshold),
  # delta2 the ratio of the interval below to it, delta1 and delta4 ratios to
  # the first interval on their side
  gaps <- function(thr) diff(thr)
  expect_equal(gaps(.sdt_make_thresholds(0, 6L, "log_ratio", deltas = c(0, 0, 1, 0))),
               rep(exp(1), 4))
  expect_equal(gaps(.sdt_make_thresholds(0, 6L, "log_ratio", deltas = c(0, 1, 0, 0))),
               c(exp(1), exp(1), 1, 1))
  expect_equal(gaps(.sdt_make_thresholds(0, 6L, "log_ratio", deltas = c(1, 0, 0, 0))),
               c(exp(1), 1, 1, 1))
  # odd K = 5: delta2 is the middle category, delta1 and delta3 ratios to it
  expect_equal(gaps(.sdt_make_thresholds(0, 5L, "log_ratio", deltas = c(0, 1, 0))),
               rep(exp(1), 3))
  expect_equal(gaps(.sdt_make_thresholds(0, 5L, "log_ratio", deltas = c(0, 0, 1))),
               c(1, 1, exp(1)))
  # odd K = 7: delta3 is the middle category, delta2 and delta4 ratios to it,
  # delta1 and delta5 ratios to the first interval on their side
  expect_equal(gaps(.sdt_make_thresholds(0, 7L, "log_ratio", deltas = c(1, 0, 0, 0, 0))),
               c(exp(1), 1, 1, 1, 1))
  expect_equal(gaps(.sdt_make_thresholds(0, 7L, "log_ratio", deltas = c(0, 1, 0, 0, 0))),
               c(exp(1), exp(1), 1, 1, 1))
  expect_equal(gaps(.sdt_make_thresholds(0, 7L, "log_ratio", deltas = c(0, 0, 0, 1, 0))),
               c(1, 1, 1, exp(1), exp(1)))
})

test_that("deltas are named in interval order", {
  for (tt in c("log_distance", "log_ratio")) {
    model <- sdt_rating(paste0("r", 1:5), "stimulus", threshold_type = tt)
    expect_equal(grep("^delta", names(model$parameters), value = TRUE),
                 paste0("delta", 1:3), info = tt)
  }
})

############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

sim_rating <- function(n_subjects, n_trials, d, criterion, n_ratings,
                       spacing = NULL, deltas = NULL,
                       threshold_type = "parsimonious",
                       sdratio = 1, dist = "normal") {
  thr <- bmm:::.sdt_make_thresholds(criterion, n_ratings, threshold_type,
                                    spacing, deltas)
  dat <- expand.grid(id = seq_len(n_subjects), stimulus = c(0L, 1L))
  cbind(dat, as.data.frame(rsdt_rating(nrow(dat), n_trials, dat$stimulus,
                                       d, thr, sdratio = sdratio,
                                       dist = dist)))
}

test_that("sdt_rating check_data validates response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)

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

test_that("sdt_rating check_data refuses non-numeric count columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  dat <- data.frame(r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30),
                    r4 = c(5, 55), stimulus = c(0L, 1L))

  for (conv in list(as.character, as.factor)) {
    bad <- dat
    bad$r1 <- conv(bad$r1)
    expect_error(check_data(model, bad, formula),
                 "Response column 'r1' must be numeric")
  }
})

test_that("sdt_rating check_data refuses an all-NA column as NA, not a type", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  dat <- data.frame(r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30),
                    r4 = NA, stimulus = c(0L, 1L))
  expect_error(check_data(model, dat, bmf(d ~ 1, criterion ~ 1, spacing ~ 1)),
               "must not contain NA counts")
})

test_that("sdt_rating check_data still only warns on non-integer counts", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  dat <- data.frame(r1 = c(10.5, 5), r2 = c(20, 10), r3 = c(15, 30),
                    r4 = c(5, 55), stimulus = c(0L, 1L))
  expect_warning(check_data(model, dat, bmf(d ~ 1, criterion ~ 1, spacing ~ 1)),
                 "Response column 'r1' should contain integer counts")
})

test_that("sdt_rating check_data rejects missing response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10),
    stimulus = c(0L, 1L)
  )
  expect_error(check_data(model, invalid_data, formula), "missing in the data")
})

test_that("sdt_rating check_data rejects negative counts", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(-1, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 1L)
  )
  expect_error(check_data(model, invalid_data, formula), "non-negative")
})

test_that("sdt_rating check_data refuses an NA or constant stimulus", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  counts <- data.frame(r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55))

  expect_error(check_data(model, cbind(counts, stimulus = c(0L, NA)), formula),
               "1 of 2 values are NA")
  expect_error(check_data(model, cbind(counts, stimulus = c(1L, 1L)), formula),
               "stimulus")
})

test_that("sdt_rating check_data validates stimulus coding", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)

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
  counts <- rsdt_rating(10, 50, rep(c(0L, 1L), 5), d = 1.5,
                        thresholds = c(-0.5, 0, 0.5))
  expect_true(is.matrix(counts))
  expect_equal(dim(counts), c(10L, 4L))
  expect_equal(colnames(counts), paste0("r", 1:4))
  expect_true(all(rowSums(counts) == 50))
})

test_that("rsdt_rating recycles vectorized parameters per observation", {
  counts <- rsdt_rating(6, c(50, 100), c(0L, 1L), d = rnorm(6, 1.5, 0.2),
                        thresholds = c(-0.5, 0, 0.5), sdratio = c(1, 1.3))
  expect_equal(dim(counts), c(6L, 4L))
  expect_true(all(rowSums(counts) == c(50, 100)))

  thr <- rbind(c(-0.5, 0, 0.5), c(-1, 0, 1))
  counts2 <- rsdt_rating(2, 100, c(0L, 1L), d = 1.5, thresholds = thr)
  expect_equal(dim(counts2), c(2L, 4L))
})

test_that("rsdt_rating validates inputs", {
  expect_error(rsdt_rating(c(2, 3), 50, 1L, d = 1,
                           thresholds = c(-0.5, 0, 0.5)),
               "single positive integer")
  expect_error(rsdt_rating(2, 50, c(0L, 2L), d = 1,
                           thresholds = c(-0.5, 0, 0.5)),
               "0 \\(noise\\) or 1 \\(signal\\)")
  expect_error(rsdt_rating(2, 50, c(0L, 1L), d = 1,
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
    dat <- do.call(sim_rating, c(list(3, 50, d = 1.5, criterion = 0), cfg))
    rcols <- paste0("r", seq_len(cfg$n_ratings))
    expect_equal(nrow(dat), 6, info = cfg$threshold_type)
    expect_true(all(rowSums(dat[, rcols]) == 50), info = cfg$threshold_type)
  }
})

test_that("dsdt_rating returns valid density", {
  dens <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                      d = 1.5, thresholds = c(-0.5, 0, 0.5))
  expect_true(dens > 0)
  expect_true(dens <= 1)
})

test_that("dsdt_rating log matches log of density", {
  dens <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                      d = 1.5, thresholds = c(-0.5, 0, 0.5))
  ld <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                    d = 1.5, thresholds = c(-0.5, 0, 0.5), log = TRUE)
  expect_equal(log(dens), ld, tolerance = 1e-10)
})

test_that("dsdt_rating validates inputs", {
  expect_error(
    dsdt_rating(counts = c(-1, 20, 30, 40), stimulus = 1,
                d = 1.5, thresholds = c(-0.5, 0, 0.5)),
    "non-negative"
  )
  expect_error(
    dsdt_rating(counts = c(10, 20, 30), stimulus = 1,
                d = 1.5, thresholds = c(-0.5, 0, 0.5)),
    "K - 1"
  )
  expect_error(
    dsdt_rating(counts = c(NA, 20, 30, 40), stimulus = 1,
                d = 1.5, thresholds = c(-0.5, 0, 0.5)),
    "must not contain NA"
  )
})

test_that("dsdt_rating works for all distributions", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    dens <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                        d = 1.5, thresholds = c(-0.5, 0, 0.5), dist = di)
    expect_true(dens > 0, info = paste("dist:", di))
  }
})

test_that("dsdt_rating with sdratio != 1 changes density", {
  d_ev <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                      d = 1.5, thresholds = c(-0.5, 0, 0.5), sdratio = 1)
  d_uv <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                      d = 1.5, thresholds = c(-0.5, 0, 0.5), sdratio = 1.3)
  expect_false(d_ev == d_uv)
})

test_that("dsdt_rating matches dmultinom and is vectorized over rows", {
  thr <- c(-0.5, 0, 0.5)
  probs <- bmm:::.sdt_category_probs(thr, 1.5, 1, 1, "normal")
  expect_equal(dsdt_rating(c(10, 20, 30, 40), 1, 1.5, thr),
               dmultinom(c(10, 20, 30, 40), prob = probs), tolerance = 1e-12)

  counts <- rbind(c(10, 20, 30, 40), c(40, 30, 20, 10))
  dens <- dsdt_rating(counts, stimulus = c(1L, 0L), d = 1.5,
                      thresholds = thr, log = TRUE)
  expect_length(dens, 2)
  expect_true(all(is.finite(dens)))
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
  d <- matrix(c(1, 1.5, 0.5, 2), nrow = 2)
  out <- sdt_rating_logmu(2L, 4L, 1L, 1L, d, criterion = 0.1,
                          spacing = 0.3, sdratio = 0.2, stimulus = c(0, 1))
  expect_equal(dim(out), dim(d))

  thr <- bmm:::.sdt_make_thresholds(0.1, 4L, "parsimonious", 0.3)
  stim <- rep_len(c(0, 1), 4)
  ref <- vapply(seq_along(d), function(j) {
    log(bmm:::.sdt_category_probs(thr, as.vector(d)[j], exp(0.2),
                                  stim[j], "normal")[2])
  }, numeric(1))
  expect_equal(as.vector(out), ref, tolerance = 1e-12)
})

test_that("log_ratio requires at least 4 rating categories", {
  expect_error(sdt_rating(c("r1", "r2", "r3"), "stimulus",
                          threshold_type = "log_ratio"),
               "at least 4 rating categories")
})

test_that("K = 3 softmax reduces to a single interval centred on criterion", {
  thr <- bmm:::.sdt_make_thresholds(0.1, 3L, "softmax", spacing = 0.2)
  expect_equal(thr, 0.1 + c(-0.5, 0.5) * exp(0.2), tolerance = 1e-12)

  dat <- sim_rating(3, 50, d = 1.2, criterion = 0, n_ratings = 3,
                    spacing = 0.2, threshold_type = "softmax")
  model <- sdt_rating(paste0("r", 1:3), "stimulus", threshold_type = "softmax")
  code <- stancode(bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
                   data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("bmf2bf maps every rating category to its own logit", {
  # names that do not sort in category order, so a mapping by name or by sorted
  # position would show up as a shifted category index
  cats <- c("surenew", "new", "old", "sureold")
  bf <- bmf2bf(sdt_rating(cats, "stimulus"), bmf(d ~ 1))
  expect_match(deparse(bf$formula, width.cutoff = 500),
               "sdt_rating_logmu(1, 4,", fixed = TRUE)
  for (k in 2:4) {
    pform <- deparse(bf$pforms[[paste0("mu", cats[k])]], width.cutoff = 500)
    expect_match(pform, paste0("sdt_rating_logmu(", k, ", 4,"), fixed = TRUE)
  }
})

test_that("a link on d or criterion is inverted inside the rating formula", {
  # the multinomial family has no link of its own for d and criterion, so a
  # link the user sets reaches the kernel only through the formula
  bf <- bmf2bf(sdt_rating(paste0("r", 1:4), "stimulus",
                          links = list(d = "log", criterion = "softplus")),
               bmf(d ~ 1))
  for (f in c(list(bf$formula), bf$pforms[paste0("mur", 2:4)])) {
    txt <- paste(deparse(f, width.cutoff = 500), collapse = "")
    expect_match(txt, "exp(d)", fixed = TRUE)
    expect_match(txt, "log1p_exp(criterion)", fixed = TRUE)
  }
  bf <- bmf2bf(sdt_rating(paste0("r", 1:4), "stimulus"), bmf(d ~ 1))
  expect_false(grepl("exp(d)", deparse(bf$formula, width.cutoff = 500), fixed = TRUE))
})

test_that("sdt_rating produces valid stancode with parsimonious thresholds", {
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
  expect_true(grepl("sdt_thresholds_parsimonious_rating", code, fixed = TRUE))
})

test_that("sdt_rating stancode loads the shared distribution dispatch", {
  # the noise-distribution dispatch lives in sdt_dist_funs.stan; the rating
  # likelihood relies on the log-scale dispatchers from that shared chunk.
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(grepl("sdt_log_cumprob", code, fixed = TRUE))
  expect_true(grepl("sdt_log_one_minus_cumprob", code, fixed = TRUE))
  expect_true(grepl("log_diff_exp", code, fixed = TRUE))
})

test_that("sdt_rating produces valid stancode with log_distance thresholds", {
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0.5, 0.5), threshold_type = "log_distance")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  formula <- bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
  expect_true(grepl("sdt_thresholds_log_distance_rating", code, fixed = TRUE))
})

test_that("sdt_rating produces valid stancode with log_ratio thresholds", {
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0, 0), threshold_type = "log_ratio")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_ratio")
  formula <- bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_thresholds_log_ratio_rating", code, fixed = TRUE))
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
})

test_that("sdt_rating produces valid stancode with softmax thresholds", {
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 6,
                    spacing = 0.3, deltas = c(0, 0, 0),
                    threshold_type = "softmax")
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "softmax")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1,
                 delta1 ~ 1, delta2 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_thresholds_softmax_rating", code, fixed = TRUE))
})

test_that("sdt_rating handles K=6 parsimonious", {
  dat <- sim_rating(2, 50, d = 1.5, criterion = 0, n_ratings = 6,
                    spacing = 0.3)
  model <- sdt_rating(paste0("r", 1:6), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_rating handles predictors on d", {
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ condition, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


# The Stan kernels run through a fixed_param generated-quantities program, as
# the ranking parity tests do. sig_figs = 17 round-trips a double. Every row
# builds K thresholds with its own builder and returns one category's log
# probability, so both the builders and the kernel meet their R counterparts.
rating_stan_logmu <- function(grid, deltas, K) {
  sc <- system.file("stan_chunks", package = "bmm")
  funs <- paste(read_lines2(file.path(sc, "sdt_dist_funs.stan")),
                read_lines2(file.path(sc, "sdt_rating_funs.stan")), sep = "\n")
  program <- paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n  int N; int K; array[N] int cat; array[N] int thresh;\n",
    "  vector[N] criterion; vector[N] spacing; array[N, K - 2] real deltas;\n",
    "  vector[N] d; vector[N] sdratio; vector[N] stimulus; array[N] int dist;\n}\n",
    "generated quantities {\n  vector[N] lp;\n  for (i in 1:N) {\n",
    "    vector[K - 1] thr = sdt_make_thresholds_rating(criterion[i], spacing[i],\n",
    "      deltas[i], K, thresh[i]);\n",
    "    lp[i] = sdt_rating_logmu_cat(cat[i], thr, d[i], sdratio[i], stimulus[i], dist[i]);\n",
    "  }\n}\n"
  )
  data <- c(list(N = nrow(grid), K = K, deltas = deltas),
            as.list(grid[c("cat", "thresh", "criterion", "spacing", "d",
                           "sdratio", "stimulus", "dist")]))
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1,
    iter_warmup = 0, refresh = 0, show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#",
                         check.names = FALSE)
  as.numeric(csv[1, paste0("lp.", seq_len(nrow(grid)))])
}

test_that("the Stan rating kernel matches its R counterpart in both tails", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  # criterion = 12 puts every threshold far into the upper tail, where the log
  # cdf of both bounds of an interval rounds to 0; -12 does the same below.
  # K = 5 and K = 6 exercise the odd and even placement of the criterion; K = 7
  # is the smallest odd K with log_ratio intervals beyond the first on each side.
  for (K in c(5L, 6L, 7L)) {
    grid <- expand.grid(cat = seq_len(K), thresh = seq_along(bmm:::.sdt_threshold_types),
                        criterion = c(-12, -0.5, 0.3, 12), spacing = c(-0.5, 0.1),
                        d = c(0, 2), sdratio = c(-0.6, 0, 0.9),
                        stimulus = c(0, 1), dist = seq_along(bmm:::.sdt_dists))
    deltas <- matrix(seq(-0.6, 0.6, length.out = K - 2L), nrow(grid), K - 2L, byrow = TRUE)
    stan <- rating_stan_logmu(grid, deltas, K)

    r <- vapply(seq_len(nrow(grid)), function(i) {
      g <- grid[i, ]
      type <- bmm:::.sdt_threshold_types[g$thresh]
      n_deltas <- if (type == "softmax") K - 3L else K - 2L
      thr <- bmm:::.sdt_make_thresholds(
        g$criterion, K, type, g$spacing,
        if (type %in% c("log_distance", "log_ratio", "softmax")) deltas[1, seq_len(n_deltas)]
      )
      bmm:::.sdt_category_log_probs(thr, g$d, exp(g$sdratio), g$stimulus,
                                    names(bmm:::.sdt_dists)[g$dist])[g$cat]
    }, numeric(1))

    expect_true(all(is.finite(stan)), info = paste("K =", K))
    expect_equal(stan, r, tolerance = 1e-10, info = paste("K =", K))
  }
})

test_that("an upper-tail rating category keeps its mass", {
  # under gumbel_min the log cdf rounds to 0 from eta ~ 6.6, so an interval
  # taken as the difference of two log cdfs came out -Inf. Its survival
  # function exp(-exp(x)) gives the interval's mass in closed form.
  lo <- c(4, 7.5, 12)
  hi <- lo + 0.5
  expected <- -exp(lo) + log1m_exp(exp(lo) - exp(hi))
  got <- vapply(seq_along(lo), function(i) {
    bmm:::.sdt_category_log_probs(c(lo[i] - 1, lo[i], hi[i]), 0, 1, 0,
                                  "gumbel_min")[3]
  }, numeric(1))
  expect_equal(got, expected, tolerance = 1e-12)
})


############################################################################# !
# UV-SDT TESTS (sdratio overridable fixed parameter)                     ####
############################################################################# !

test_that("sdt_rating with sdratio ~ 1 produces valid stancode", {
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5, sdratio = 1.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
  expect_true(grepl("sdratio", code))
})

test_that("sdt_rating UV-SDT with log_distance thresholds produces valid stancode", {
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0.5, 0.5), threshold_type = "log_distance",
                    sdratio = 1.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  formula <- bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1,
                 sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdratio", code))
})


############################################################################# !
# PIPELINE INTEGRATION TESTS (mock backend)                              ####
############################################################################# !

test_that("sdt_rating integrates with the bmm pipeline via mock backend", {
  dat <- sim_rating(6, 80, d = 1.2, criterion = 0, n_ratings = 4,
                    spacing = 0.4)
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", dist = di)
    expect_silent(
      bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
          dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
    )
  }
})

test_that("sdt_rating UV-SDT integrates with the bmm pipeline via mock backend", {
  dat <- sim_rating(6, 80, d = 1.2, criterion = 0, n_ratings = 4,
                    spacing = 0.4, sdratio = 1.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_silent(
    bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
        dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_rating log_distance integrates with the bmm pipeline via mock backend", {
  # log_distance has variable-arity delta parameters, so this exercises the
  # generated multinomial formula and Stan function for that case.
  dat <- sim_rating(5, 60, d = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0.5, 0.5), threshold_type = "log_distance")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  expect_silent(
    bmm(bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1),
        dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_rating posterior_predict draws joint multinomial counts (real fit)", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)))

  dat <- sim_rating(6, 120, d = 1.4, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  fit <- bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1), dat, model,
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

  dat <- sim_rating(6, 150, d = 1.4, criterion = 0, n_ratings = 6,
                    deltas = c(-0.4, -0.4, -0.4, -0.4),
                    threshold_type = "log_distance")
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "log_distance")
  fit <- bmm(bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1,
                 delta3 ~ 1, delta4 ~ 1),
             dat, model, backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  expect_true(all(is.finite(brms::log_lik(fit))))
})
