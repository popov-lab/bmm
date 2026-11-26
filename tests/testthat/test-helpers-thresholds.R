# The ordered-threshold helpers shared by sdt_rating(), sdt_cdp() and
# likert_cumulative(): the threshold builders, the ordinal kernel and the noise
# registry. The .sdt_* names are aliases until no branch reads them.

test_that("the .sdt_* names are aliases of the shared helpers", {
  expect_identical(bmm:::.sdt_threshold_types, bmm:::.threshold_types)
  expect_identical(bmm:::.sdt_threshold_parameter_parts, bmm:::.threshold_parameter_parts)
  expect_identical(bmm:::.sdt_log_ratio_delta_labels, bmm:::.log_ratio_delta_labels)
  expect_identical(bmm:::.sdt_threshold_delta_names, bmm:::.threshold_delta_names)
  expect_identical(bmm:::.sdt_make_thresholds, bmm:::.make_thresholds)
  expect_identical(bmm:::.sdt_assemble_thresholds, bmm:::.assemble_thresholds)
  expect_identical(bmm:::.sdt_log_ratio_widths, bmm:::.log_ratio_widths)
  expect_identical(bmm:::.sdt_dists, bmm:::.noise_dists)
  expect_identical(bmm:::.sdt_dist_id, bmm:::.noise_dist_id)
  expect_identical(bmm:::.sdt_dist_names, bmm:::.noise_dist_names)
})

test_that("the noise registry positions are the Stan dist_type codes", {
  expect_identical(bmm:::.noise_dist_names, c("normal", "gumbel_min", "gumbel_max", "logistic"))
  expect_identical(bmm:::.noise_dist_id(c("logistic", "normal")), c(4L, 1L))
})

test_that(".ordinal_category_log_probs differences the cdf at the standardized thresholds", {
  # against the probability-scale closed form, which is exact away from the tails
  thr <- c(-1.2, -0.3, 0.4, 1.5)
  for (dist in bmm:::.noise_dist_names) {
    cdf <- bmm:::.noise_dists[[dist]]$cdf
    for (location in c(-0.7, 0, 1.1)) {
      for (scale in c(1, 1.6)) {
        eta <- (thr - location) / scale
        expected <- log(diff(c(0, cdf(eta), 1)))
        got <- bmm:::.ordinal_category_log_probs(thr, location, scale, dist)
        expect_equal(got, expected, tolerance = 1e-13, info = paste(dist, location, scale))
        expect_length(got, 5L)
      }
    }
  }
})

test_that(".ordinal_category_log_probs is vectorized over rows and recycles its inputs", {
  thr <- rbind(c(-1, 0, 1), c(-0.5, 0.2, 0.9), c(-2, -1, 3))
  location <- c(0.3, -0.4, 1)
  scale <- c(1, 1.3, 0.8)
  out <- bmm:::.ordinal_category_log_probs(thr, location, scale, "logistic")
  expect_equal(dim(out), c(3L, 4L))
  for (i in 1:3) {
    expect_equal(out[i, ], bmm:::.ordinal_category_log_probs(thr[i, ], location[i], scale[i], "logistic"))
  }
  # a threshold vector recycles across the rows the other arguments define
  rec <- bmm:::.ordinal_category_log_probs(c(-1, 0, 1), location, 1, "normal")
  expect_equal(dim(rec), c(3L, 4L))
  expect_equal(rec[2, ], bmm:::.ordinal_category_log_probs(c(-1, 0, 1), -0.4, 1, "normal"))
  # the categories of every row exhaust the probability
  expect_equal(log(rowSums(exp(out))), rep(0, 3), tolerance = 1e-12)
})

test_that(".ordinal_category_log_probs keeps the mass of categories deep in either tail", {
  # at location 12 every threshold sits in the lower tail of the latent
  # variable, at -12 in the upper tail; the sum over categories is still 1 and
  # every category keeps a finite log probability
  thr <- c(-1.5, -0.5, 0.5, 1.5)
  for (dist in bmm:::.noise_dist_names) {
    for (location in c(-12, 12)) {
      out <- bmm:::.ordinal_category_log_probs(thr, location, 1, dist)
      expect_true(all(is.finite(out)), info = paste(dist, location))
      expect_equal(log(sum(exp(out))), 0, tolerance = 1e-12, info = paste(dist, location))
    }
  }
})

test_that(".sdt_category_log_probs is the ordinal kernel at the d_a shift and the signal SD", {
  thr <- c(-0.8, 0.1, 0.9)
  d <- 1.4
  sdratio <- 1.3
  noise <- bmm:::.sdt_category_log_probs(thr, d, sdratio, 0, "normal")
  signal <- bmm:::.sdt_category_log_probs(thr, d, sdratio, 1, "normal")
  expect_equal(noise, bmm:::.ordinal_category_log_probs(thr, 0, 1, "normal"))
  expect_equal(signal, bmm:::.ordinal_category_log_probs(
    thr, d * sqrt((1 + sdratio^2) / 2), sdratio, "normal"
  ))
})

# The Stan kernel runs through a fixed_param generated-quantities program, as
# the rating parity test does. Every row builds its thresholds with
# make_thresholds() and returns one category's log probability, so the shared
# builders and the kernel meet their R twins. sig_figs = 17 round-trips a double.
ordinal_stan_logp <- function(grid, deltas, K) {
  sc <- system.file("stan_chunks", package = "bmm")
  chunks <- c("sdt_dist_functions.stan", "threshold_functions.stan", "ordinal_functions.stan")
  funs <- paste(vapply(file.path(sc, chunks), read_lines2, character(1)), collapse = "\n")
  program <- paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n  int N; int K; array[N] int cat; array[N] int thresh;\n",
    "  vector[N] criterion; vector[N] spacing; array[N, K - 2] real deltas;\n",
    "  vector[N] location; vector[N] scale; array[N] int dist;\n}\n",
    "generated quantities {\n  vector[N] lp;\n  for (i in 1:N) {\n",
    "    vector[K - 1] thr = make_thresholds(criterion[i], spacing[i], deltas[i], K, thresh[i]);\n",
    "    lp[i] = ordinal_log_prob_cat(cat[i], thr, location[i], scale[i], dist[i]);\n",
    "  }\n}\n"
  )
  data <- c(list(N = nrow(grid), K = K, deltas = deltas),
            as.list(grid[c("cat", "thresh", "criterion", "spacing", "location", "scale", "dist")]))
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1,
    iter_warmup = 0, refresh = 0, show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  as.numeric(csv[1, paste0("lp.", seq_len(nrow(grid)))])
}

test_that("the Stan ordinal kernel matches its R twin in both tails", {
  skip_unless_stan_tests()
  skip_on_cran()

  # location 8 puts every threshold into the lower tail of the latent variable,
  # -8 into the upper tail, where the log cdf of both bounds of an interval
  # rounds to 0. K = 5 and 6 exercise the odd and even placement of the
  # criterion; K = 7 is the smallest odd K with log_ratio intervals beyond the
  # first on each side.
  for (K in c(5L, 6L, 7L)) {
    grid <- expand.grid(cat = seq_len(K), thresh = seq_along(bmm:::.threshold_types),
                        criterion = c(-0.4, 0.6), spacing = c(-0.3, 0.2),
                        location = c(-8, -0.9, 0, 1.7, 8), scale = c(1, 1.8),
                        dist = seq_along(bmm:::.noise_dists))
    deltas <- matrix(seq(-0.5, 0.5, length.out = K - 2L), nrow(grid), K - 2L, byrow = TRUE)
    stan <- ordinal_stan_logp(grid, deltas, K)

    r <- vapply(seq_len(nrow(grid)), function(i) {
      g <- grid[i, ]
      type <- bmm:::.threshold_types[g$thresh]
      n_deltas <- if (type == "softmax") K - 3L else K - 2L
      thr <- bmm:::.make_thresholds(
        g$criterion, K, type, g$spacing,
        if (type %in% c("log_distance", "log_ratio", "softmax")) deltas[1, seq_len(n_deltas)]
      )
      bmm:::.ordinal_category_log_probs(thr, g$location, g$scale,
                                        bmm:::.noise_dist_names[g$dist])[g$cat]
    }, numeric(1))

    expect_true(all(is.finite(stan)), info = paste("K =", K))
    expect_equal(stan, r, tolerance = 1e-10, info = paste("K =", K))
  }
})

test_that("the log-likelihoods of the threshold models at fixed parameters are unchanged", {
  # pinned on 2026-11-05 before the helpers were lifted out of the SDT files;
  # the data are the help-page simulations under fixed seeds and the ranking
  # dataset. A change here is a change of the likelihood, not of the pin.
  dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
  dat <- cbind(dat, withr::with_seed(11, rsdt_rating(nrow(dat), 200, dat$stimulus, d = 1.5,
                                                     thresholds = c(0.25, 0.75, 1.25))))
  thr <- bmm:::.make_thresholds(0.3, 4L, "parsimonious", spacing = -0.2)
  expect_equal(sum(dsdt_rating(as.matrix(dat[, c("r1", "r2", "r3", "r4")]), dat$stimulus, d = 1.2,
                               thresholds = thr, log = TRUE)),
               -1647.6631528994869, tolerance = 1e-12)

  expect_equal(sum(dsdt_ranking(as.matrix(meyer_grant_jakob_2025[, paste0("rank", 1:5)]),
                                meyer_grant_jakob_2025$set_size, d = 1.1, log = TRUE)),
               -1965.9469193693587, tolerance = 1e-12)

  dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
  thr_cdp <- c(-1.1, -0.5, 0, 0.6, 1.3)
  dat <- cbind(dat, withr::with_seed(12, rsdt_cdp(nrow(dat), 200, dat$stimulus, dfam = 0.8, drec = 1.0,
                                                  thresholds = thr_cdp, rcrit = 0.7, n_new = 3)))
  counts <- grep("^id$|^stimulus$", names(dat), invert = TRUE, value = TRUE)
  expect_equal(sum(dsdt_cdp(as.matrix(dat[, counts]), dat$stimulus, dfam = 0.7, drec = 0.9,
                            thresholds = thr_cdp + 0.1, rcrit = 0.6, n_new = 3, log = TRUE)),
               -953.2166914451393, tolerance = 1e-12)

  dat <- expand.grid(contrast = c(0.01, 0.02, 0.04, 0.08, 0.16), id = 1:10)
  dat$n_correct <- withr::with_seed(13, rpsychometric_function(
    nrow(dat), dat$contrast, midpoint = 0.04, width = 2, guess = 0.5, lapse = 0.02,
    n_trials = 40L, sigmoid = "weibull"))
  expect_equal(sum(dpsychometric_function(dat$n_correct, dat$contrast, midpoint = 0.05, width = 1.8,
                                          guess = 0.5, lapse = 0.03, n_trials = 40, sigmoid = "weibull",
                                          log = TRUE)),
               -107.15958146403861, tolerance = 1e-12)
})
