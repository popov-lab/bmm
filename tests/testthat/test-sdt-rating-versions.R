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
  f <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  expect_equal(names(update_model_fixed_parameters(m, f)$fixed_parameters),
               c("sdratio", "Ro", "Rn"))

  f1 <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1)
  expect_equal(names(update_model_fixed_parameters(m, f1)$fixed_parameters),
               c("sdratio", "Rn"))

  f2 <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1)
  expect_equal(names(update_model_fixed_parameters(m, f2)$fixed_parameters),
               "sdratio")
})

test_that("the metad version refuses an odd number of rating categories", {
  # at odd K the criterion is the centre of the middle category, so there is no
  # old/new boundary for the type-1 normalisation to split at
  expect_error(sdt_rating(paste0("r", 1:5), "stimulus", version = "metad"),
               "even number of rating categories")
  thr5 <- c(-1, -0.3, 0.3, 1)
  expect_error(dsdt_metad(c(5, 15, 25, 30, 25), 1L, 1.5, thr5, metad = 1),
               "even number of rating categories")
  expect_error(rsdt_metad(2, 50, 1L, 1.5, thr5, metad = 1),
               "even number of rating categories")
  expect_s3_class(sdt_rating(paste0("r", 1:5), "stimulus", version = "dpsdt"),
                  "sdt_rating_dpsdt")
})

test_that("Ro, Rn and logmratio keep their identity links", {
  for (par in c("Ro", "Rn")) {
    expect_error(
      sdt_rating(paste0("r", 1:4), "stimulus", version = "dpsdt",
                 links = stats::setNames(list("log"), par)),
      paste0("link of '", par, "' cannot be changed")
    )
  }
  expect_error(
    sdt_rating(paste0("r", 1:4), "stimulus", version = "metad",
               links = list(logmratio = "log")),
    "link of 'logmratio' cannot be changed"
  )
})

test_that("a link on d reaches the dpsdt and metad kernels through the formula", {
  for (v in c("dpsdt", "metad")) {
    bf <- bmf2bf(sdt_rating(paste0("r", 1:4), "stimulus", version = v,
                            links = list(d = "log")), bmf(d ~ 1))
    for (f in c(list(bf$formula), bf$pforms[paste0("mur", 2:4)])) {
      txt <- paste(deparse(f, width.cutoff = 500), collapse = "")
      expect_match(txt, paste0("sdt_", v, "_logmu"), fixed = TRUE)
      expect_match(txt, "exp(d)", fixed = TRUE)
    }
  }
})

############################################################################# !
# CATEGORY PROBABILITIES: REDUCTION IDENTITIES & STRUCTURE               ####
############################################################################# !

test_that("dpsdt with no recollection reduces to the standard category probs", {
  thr <- bmm:::.sdt_make_thresholds(0, 5, "parsimonious", spacing = 0.3)
  for (stim in c(0L, 1L)) {
    base <- bmm:::.sdt_category_probs(thr, 1.4, 1.2, stim, "normal")
    off  <- bmm:::.sdt_dpsdt_category_probs(thr, 1.4, 1.2, stim, "normal",
                                            -100, -100)
    expect_equal(off, base, tolerance = 1e-12)
  }
})

test_that("metad equal to d reduces to the standard category probs", {
  thr <- bmm:::.sdt_make_thresholds(0, 6, "parsimonious", spacing = 0.3)
  for (stim in c(0L, 1L)) {
    base <- bmm:::.sdt_category_probs(thr, 1.4, 1.2, stim, "normal")
    md   <- bmm:::.sdt_metad_category_probs(thr, 1.4, 1.2, stim, "normal", 1.4)
    expect_equal(md, base, tolerance = 1e-12)
  }
})

test_that("metad holds each side of the criterion to the type-1 rate", {
  # the old/new boundary is criterion itself (threshold K/2), whatever the
  # threshold type, so the "new" side carries F((criterion - shift) / scale)
  sr <- 1.3
  d <- 1.5
  crit <- 0.4
  shift <- d * bmm:::.sdt_rms_scale(sr) / 2
  for (tt in c("parsimonious", "log_distance", "softmax")) {
    thr <- bmm:::.sdt_make_thresholds(crit, 6L, tt, spacing = 0.2,
                                      deltas = if (tt == "softmax") c(0.3, -0.2, 0.1)
                                               else c(-0.4, 0.2, 0.1, 0.5))
    ps <- bmm:::.sdt_metad_category_probs(thr, d, sr, 1L, "normal", 0.7)
    pn <- bmm:::.sdt_metad_category_probs(thr, d, sr, 0L, "normal", 0.7)
    expect_equal(sum(ps[1:3]), pnorm((crit - shift) / sr), tolerance = 1e-12,
                 info = tt)
    expect_equal(sum(pn[1:3]), pnorm(crit + shift), tolerance = 1e-12, info = tt)
  }
})

# Largest elementwise error relative to max(1, |y|): expect_equal()'s tolerance
# is relative to the mean over the whole vector, which a -1e6 log probability
# would dominate. Matching infinities count as equal.
max_rel_err <- function(x, y) {
  keep <- !(is.infinite(x) & x == y)
  max(0, abs(x - y)[keep] / pmax(1, abs(y[keep])))
}

test_that("dpsdt and metad keep their category mass in both tails", {
  # criterion = -12 or 12 puts every threshold far into one tail, where the
  # probability scale rounds the outer categories to 0 and the old kernels
  # clamped them to the double epsilon
  for (crit in c(-12, 12)) {
    for (dist in names(bmm:::.sdt_dists)) {
      for (stim in 0:1) {
        thr <- bmm:::.sdt_make_thresholds(crit, 6L, "equidistant", spacing = 0)
        base <- bmm:::.sdt_category_log_probs(thr, 1.5, 1.2, stim, dist)
        # recollection exactly off: the model's default of -100 floors the
        # loaded category at log p = -100, above a far-tail familiarity mass
        dp_off <- bmm:::.sdt_dpsdt_category_log_probs(thr, 1.5, 1.2, stim, dist,
                                                      -Inf, -Inf)
        md_ideal <- bmm:::.sdt_metad_category_log_probs(thr, 1.5, 1.2, stim,
                                                        dist, 1.5)
        info <- paste(crit, dist, stim)
        expect_lt(max_rel_err(dp_off, base), 1e-12, label = info)
        expect_lt(max_rel_err(md_ideal, base), 1e-12, label = info)

        dp <- bmm:::.sdt_dpsdt_category_log_probs(thr, 1.5, 1.2, stim, dist,
                                                  0.5, -0.5)
        md <- bmm:::.sdt_metad_category_log_probs(thr, 1.5, 1.2, stim, dist,
                                                  0.8)
        expect_equal(matrixStats::logSumExp(dp), 0, tolerance = 1e-10,
                     info = info)
        expect_equal(matrixStats::logSumExp(md[is.finite(md)]), 0,
                     tolerance = 1e-10, info = info)
      }
    }
  }
})

test_that("dpsdt category probs match an independent hand computation (K=4)", {
  thr <- c(-0.6, 0.1, 0.7)
  d <- 1.5
  Ro <- 0.35
  base_s <- diff(c(0, pnorm(thr - d / 2), 1))            # signal, sr = 1
  hand_s <- (1 - Ro) * base_s
  hand_s[4] <- hand_s[4] + Ro
  got_s <- bmm:::.sdt_dpsdt_category_probs(thr, d, 1, 1L, "normal",
                                           qlogis(Ro), -Inf)
  expect_equal(got_s, hand_s, tolerance = 1e-12)
})

test_that("recollection moves mass to the most-confident category and sums to 1", {
  thr <- bmm:::.sdt_make_thresholds(0, 5, "parsimonious", spacing = 0.3)
  base_s <- bmm:::.sdt_category_probs(thr, 1.4, 1, 1L, "normal")
  ps     <- bmm:::.sdt_dpsdt_category_probs(thr, 1.4, 1, 1L, "normal",
                                           qlogis(0.4), -Inf)
  expect_gt(ps[5], base_s[5])
  expect_equal(sum(ps), 1, tolerance = 1e-12)

  base_n <- bmm:::.sdt_category_probs(thr, 1.4, 1, 0L, "normal")
  pn     <- bmm:::.sdt_dpsdt_category_probs(thr, 1.4, 1, 0L, "normal",
                                           -Inf, qlogis(0.4))
  expect_gt(pn[1], base_n[1])
})

############################################################################# !
# DISTRIBUTION & SIMULATION FUNCTIONS                                    ####
############################################################################# !

test_that("rsdt_dpsdt / rsdt_metad return rating count matrices", {
  thr <- bmm:::.sdt_make_thresholds(0, 4L, "parsimonious", 0.3)
  counts_dp <- rsdt_dpsdt(8, 50, rep(c(0L, 1L), 4), d = 1.4, thresholds = thr,
                          Ro = 0.3, Rn = 0.1)
  expect_true(is.matrix(counts_dp))
  expect_equal(colnames(counts_dp), paste0("r", 1:4))
  expect_true(all(rowSums(counts_dp) == 50))

  m <- rsdt_metad(8, 50, rep(c(0L, 1L), 4), d = 1.4, thresholds = thr,
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

test_that("dsdt_dpsdt() and rsdt_dpsdt() default to Rn = 0, the one-sided model", {
  thr <- c(-0.5, 0, 0.5)
  for (stim in c(0L, 1L)) {
    counts <- if (stim == 1L) c(2, 8, 20, 70) else c(70, 20, 8, 2)
    expect_equal(dsdt_dpsdt(counts, stim, 1.5, thr, Ro = 0.3),
                 dsdt_dpsdt(counts, stim, 1.5, thr, Ro = 0.3, Rn = 0),
                 tolerance = 1e-12)
  }
  # noise rows: without recall-to-reject the bottom category keeps its
  # familiarity mass; 1e5 trials put the sampling error near 1e-3
  p_bottom <- bmm:::.sdt_dpsdt_category_probs(thr, 1.5, 1, 0L, "normal",
                                              qlogis(0.3), -Inf)[1]
  counts_n <- rsdt_dpsdt(1, 1e5, 0L, d = 1.5, thresholds = thr, Ro = 0.3)
  expect_lt(abs(counts_n[1, 1] / 1e5 - p_bottom), 0.01)
})

test_that("dsdt version densities match dmultinom and vectorize over rows", {
  thr <- c(-0.5, 0, 0.5)
  p_dp <- bmm:::.sdt_dpsdt_category_probs(thr, 1.5, 1, 1L, "normal",
                                          qlogis(0.3), -Inf)
  expect_equal(dsdt_dpsdt(c(2, 8, 20, 70), 1L, 1.5, thr, Ro = 0.3, Rn = 0),
               dmultinom(c(2, 8, 20, 70), prob = p_dp), tolerance = 1e-12)

  p_md <- bmm:::.sdt_metad_category_probs(thr, 1.5, 1, 1L, "normal", 1.0)
  expect_equal(dsdt_metad(c(5, 15, 25, 55), 1L, 1.5, thr, metad = 1.0),
               dmultinom(c(5, 15, 25, 55), prob = p_md), tolerance = 1e-12)

  counts <- rbind(c(2, 8, 20, 70), c(70, 20, 8, 2))
  dens <- dsdt_dpsdt(counts, c(1L, 0L), 1.5, thr, Ro = c(0.3, 0.2),
                     Rn = c(0, 0.1), log = TRUE)
  expect_length(dens, 2)
  expect_true(all(is.finite(dens)))
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
                                              qlogis(ro), qlogis(rn))
    ref_dp <- t(vapply(1:3, function(i) {
      bmm:::.sdt_dpsdt_category_probs(thr[i, ], dp[i], sr[i], stim, "normal",
                                      qlogis(ro[i]), qlogis(rn[i]))
    }, numeric(4)))
    expect_equal(vec_dp, ref_dp, tolerance = 1e-12, info = paste("dpsdt", stim))

    vec_md <- bmm:::.sdt_metad_category_probs(thr, dp, sr, stim, "normal", md)
    ref_md <- t(vapply(1:3, function(i) {
      bmm:::.sdt_metad_category_probs(thr[i, ], dp[i], sr[i], stim, "normal",
                                      md[i])
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
                                             d = 1.4, thresholds = thr,
                                             Ro = 0.3, Rn = 0.1)))
  md <- sdt_rating(paste0("r", 1:4), "stimulus", version = "dpsdt")
  mm <- sdt_rating(paste0("r", 1:4), "stimulus", version = "metad")
  expect_silent(bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1),
                    dat, md, backend = "mock", mock_fit = 1, rename = FALSE))
  expect_silent(bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1, logmratio ~ 1),
                    dat, mm, backend = "mock", mock_fit = 1, rename = FALSE))
})

test_that("the generated Stan wrapper carries the version-specific call", {
  md <- sdt_rating(paste0("r", 1:4), "stimulus", version = "dpsdt")
  code <- bmm:::.sdt_rating_logmu_stan(md)
  expect_match(code, "real sdt_dpsdt_logmu\\(")
  expect_match(code, "real Ro, real Rn, real stimulus")
  expect_match(code, "sdt_dpsdt_logmu_cat\\(cat, thr, d, sdratio, stimulus, dist_type, Ro, Rn\\)")

  mm <- sdt_rating(paste0("r", 1:4), "stimulus", version = "metad")
  code_md <- bmm:::.sdt_rating_logmu_stan(mm)
  # the estimated parameter is the log M-ratio; meta-d' is derived in the call
  expect_match(code_md, "real logmratio, real stimulus")
  expect_match(code_md, "sdt_metad_logmu_cat\\(cat, thr, d, exp\\(logmratio\\) \\* d")
})

test_that("metad R companion with logmratio = 0 reduces to the standard model", {
  args <- list(K = 6L, dist = 1L, thresh = 1L, d = matrix(1.4),
               criterion = matrix(0), spacing = matrix(0.3),
               sdratio = matrix(0), stimulus = matrix(1))
  for (k in 1:6) {
    md_k <- do.call(sdt_metad_logmu, c(list(cat = k), args, list(logmratio = matrix(0))))
    st_k <- do.call(sdt_rating_logmu, c(list(cat = k), args))
    expect_equal(as.numeric(md_k), as.numeric(st_k), tolerance = 1e-12)
  }
})


# The Stan kernels run through a fixed_param generated-quantities program, as
# the rating parity test in test-model_sdt_rating.R does. sig_figs = 17
# round-trips a double.
versions_stan_logmu <- function(grid, deltas) {
  sc <- system.file("stan_chunks", package = "bmm")
  chunks <- c("sdt_dist_funs.stan", "sdt_rating_funs.stan",
              "sdt_dpsdt_funs.stan", "sdt_metad_funs.stan")
  funs <- paste(vapply(file.path(sc, chunks), read_lines2, character(1)),
                collapse = "\n")
  program <- paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n  int N; int K; array[N] int cat; array[N, K - 2] real deltas;\n",
    "  vector[N] criterion; vector[N] d; vector[N] sdratio; vector[N] stimulus;\n",
    "  array[N] int dist; vector[N] Ro; vector[N] Rn; vector[N] metad;\n}\n",
    "generated quantities {\n  vector[N] dp; vector[N] md;\n  for (i in 1:N) {\n",
    "    vector[K - 1] thr = sdt_make_thresholds_rating(criterion[i], 0,\n",
    "      deltas[i], K, 3);\n",
    "    dp[i] = sdt_dpsdt_logmu_cat(cat[i], thr, d[i], sdratio[i], stimulus[i],\n",
    "      dist[i], Ro[i], Rn[i]);\n",
    "    md[i] = sdt_metad_logmu_cat(cat[i], thr, d[i], metad[i], sdratio[i],\n",
    "      stimulus[i], dist[i]);\n",
    "  }\n}\n"
  )
  data <- c(list(N = nrow(grid), K = 6L, deltas = deltas),
            as.list(grid[c("cat", "criterion", "d", "sdratio", "stimulus",
                           "dist", "Ro", "Rn", "metad")]))
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1,
    iter_warmup = 0, refresh = 0, show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#",
                         check.names = FALSE)
  idx <- seq_len(nrow(grid))
  list(dp = as.numeric(csv[1, paste0("dp.", idx)]),
       md = as.numeric(csv[1, paste0("md.", idx)]))
}

test_that("the Stan dpsdt and metad kernels match their R counterparts in both tails", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  grid <- expand.grid(cat = 1:6, criterion = c(-12, -0.5, 0.3, 12),
                      sdratio = c(-0.6, 0, 0.9), stimulus = c(0, 1),
                      dist = seq_along(bmm:::.sdt_dists))
  n <- nrow(grid)
  grid$d <- rep_len(c(0.4, 1.3, 2.5), n)
  grid$Ro <- rep_len(c(-100, -1.2, 0.4, 3), n)
  grid$Rn <- rep_len(c(-100, 0.8, -2.5), n)
  grid$metad <- rep_len(c(0.2, 1.1, 2.9, 0.7, 1.8), n)
  deltas <- matrix(c(-0.6, 0.2, -0.1, 0.5), n, 4, byrow = TRUE)
  stan <- versions_stan_logmu(grid, deltas)

  r <- t(vapply(seq_len(n), function(i) {
    g <- grid[i, ]
    thr <- bmm:::.sdt_make_thresholds(g$criterion, 6L, "log_distance",
                                      deltas = deltas[1, ])
    dist <- names(bmm:::.sdt_dists)[g$dist]
    c(bmm:::.sdt_dpsdt_category_log_probs(thr, g$d, exp(g$sdratio), g$stimulus,
                                          dist, g$Ro, g$Rn)[g$cat],
      bmm:::.sdt_metad_category_log_probs(thr, g$d, exp(g$sdratio), g$stimulus,
                                          dist, g$metad)[g$cat])
  }, numeric(2)))

  expect_true(all(is.finite(stan$dp)))
  expect_true(all(is.finite(stan$md)))
  expect_lt(max_rel_err(stan$dp, r[, 1]), 1e-10)
  expect_lt(max_rel_err(stan$md, r[, 2]), 1e-10)
})

############################################################################# !
# VERSION-AWARE POST-PROCESSING (mock fits)                              ####
############################################################################# !

test_that("roc_sdt reflects dual-process recollection (higher AUC, lifted curve)", {
  fit_std <- fake_rating_fit(n_ratings = 6L)
  fit_dp  <- fake_rating_fit(n_ratings = 6L, version = "dpsdt")
  draws_std <- list(d = 1.2, criterion = 0, spacing = 0)
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

test_that("the dpsdt summary curve runs from (0, Ro) to (1 - Rn, 1), as auc_sdt() integrates", {
  fit <- fake_rating_fit(n_ratings = 6L, version = "dpsdt")
  ro <- 0.4
  rn <- 0.2
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(
      list(d = 1.2, criterion = 0, spacing = 0, Ro = qlogis(ro), Rn = qlogis(rn))),
    ranef = function(...) list(),
    .package = "brms")
  summ <- attr(roc_sdt(fit, n_points = 1000), "summary")
  n <- nrow(summ)
  # the mock draws are all equal, so mean(plogis(Ro)) over draws is ro
  expect_lt(abs(summ$Hit_mean[1] - ro), 1e-10)
  expect_equal(summ$FA[1], 0)
  expect_lt(abs(summ$FA[n - 1] - (1 - rn)), 1e-10)
  expect_equal(summ$Hit_mean[n - 1], 1)
  expect_equal(c(summ$FA[n], summ$Hit_mean[n]), c(1, 1))
  # plot, summary and AUC are one curve: the trapezoid over the summary curve
  # is the integrated area
  trapezoid <- sum(diff(summ$FA) * (summ$Hit_mean[-1] + summ$Hit_mean[-n])) / 2
  expect_lt(abs(trapezoid - mean(auc_sdt(fit)$AUC)), 1e-3)
})

test_that("auc_sdt() on a dpsdt fit integrates the recollection-lifted curve", {
  fit <- fake_rating_fit(n_ratings = 6L, version = "dpsdt")
  d_true <- 1.2
  ro <- 0.4
  rn <- 0.2
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(
      list(d = d_true, sdratio = 0, criterion = 0, spacing = 0,
           Ro = qlogis(ro), Rn = qlogis(rn))),
    ranef = function(...) list(),
    variables = function(...) c("b_d_Intercept", "b_sdratio_Intercept"),
    .package = "brms"
  )
  auc <- auc_sdt(fit)
  # FA' = (1 - Rn) FA and Hit' = Ro + (1 - Ro) Hit up to FA' = 1 - Rn, then
  # Hit' = 1, so the area follows from the familiarity AUC Phi(d / sqrt(2))
  oracle <- (1 - rn) * (ro + (1 - ro) * stats::pnorm(d_true / sqrt(2))) + rn
  expect_lt(abs(mean(auc$AUC) - oracle), 1e-3)
})

test_that("default dpsdt roc_sdt (recollection off) matches the standard roc", {
  fit_std <- fake_rating_fit(n_ratings = 6L)
  fit_dp  <- fake_rating_fit(n_ratings = 6L, version = "dpsdt")
  draws_std <- list(d = 1.2, criterion = 0, spacing = 0)
  draws_off <- c(draws_std, list(Ro = -100, Rn = -100))

  local_mocked_bindings(ranef = function(...) list(), .package = "brms")
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(draws_std), .package = "brms")
  roc_std <- roc_sdt(fit_std, n_points = 40)

  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(draws_off), .package = "brms")
  roc_dp <- roc_sdt(fit_dp, n_points = 40)

  summ_dp  <- attr(roc_dp, "summary")
  summ_std <- attr(roc_std, "summary")
  # the dpsdt curve ends at (1 - Rn, 1) before (1, 1); with recollection off
  # that extra node sits at (1, 1) and every other node is the standard curve
  n <- nrow(summ_dp)
  expect_equal(n, nrow(summ_std) + 1L)
  expect_equal(summ_dp$Hit_mean[-(n - 1L)], summ_std$Hit_mean, tolerance = 1e-8)
  expect_lt(abs(summ_dp$FA[n - 1L] - 1), 1e-40)
  expect_lt(abs(summ_dp$Hit_mean[1L]), 1e-40)
})

test_that("latent_sdt reports the response-process parameters as an attribute", {
  fit_dp <- fake_rating_fit(n_ratings = 6L, version = "dpsdt")
  local_mocked_bindings(ranef = function(...) list(), .package = "brms")
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(
      list(d = 1.2, criterion = 0, spacing = 0,
           Ro = qlogis(0.4), Rn = qlogis(0.2))),
    .package = "brms")
  lat <- latent_sdt(fit_dp)
  extra <- attr(lat, "extra")
  expect_false(is.null(extra))
  expect_setequal(extra$parameter, c("Ro", "Rn"))
  expect_equal(extra$mean[extra$parameter == "Ro"], 0.4, tolerance = 1e-6)
})
