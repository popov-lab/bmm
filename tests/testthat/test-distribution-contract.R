# The recycling contract of the r*() and d*() functions (#444, ?bmm_distributions):
# r*(n, ...) returns n draws, draw i from the i-th element (row) of every
# recycled argument; d*(x, ...) returns one value per element of the longest
# argument, value i from the i-th elements alone. Every r*() pairing test
# interleaves parameter groups and compares each group with a reference
# outside the function's own code (a closed form, quadrature, or the single-
# trial density of the matching d*()). Mutants of the samplers (vector
# reversed, recycling removed, one group shifted 10 %) were run against these
# tests when they were written; a Monte Carlo test is not expected to see a
# 10 % shift, so each sampler also gets a deterministic check with n not a
# multiple of the parameter length. Non-multiple lengths recycle silently, as
# rnorm() does; the cswald tests use multiples only until #406 lands. The
# circular mixtures are circmix's (#434).

pick <- function(x, i) {
  if (is.matrix(x)) x[(i - 1) %% nrow(x) + 1, ] else x[(i - 1) %% length(x) + 1]
}
obs_length <- function(x) if (is.matrix(x)) nrow(x) else length(x)

# the vectorised call equals the calls on the i-th elements alone, without
# recycling warnings
expect_elementwise <- function(f, args, fixed = list()) {
  n <- max(vapply(args, obs_length, integer(1)))
  vec <- NULL
  expect_no_warning(vec <- do.call(f, c(args, fixed)))
  each <- vapply(seq_len(n), function(i) {
    do.call(f, c(lapply(args, pick, i = i), fixed))
  }, numeric(1))
  expect_length(vec, n)
  expect_equal(vec, each)
}

# draw i belongs to group (i - 1) %% G + 1
group_of <- function(n, G) rep_len(seq_len(G), n)
expect_mean_near <- function(x, reference, se) expect_lt(abs(mean(x) - reference), 5 * se)
se_prop <- function(p, n) sqrt(p * (1 - p) / n)

# Wiener process with drift v and diffusion s between the barriers 0 and a,
# started at z (Cox & Miller, 1965): the probability of the upper barrier,
# the mean decision time, and the mean decision time given each barrier
wiener_p_upper <- function(v, a, z, s = 1) (1 - exp(-2 * v * z / s^2)) / (1 - exp(-2 * v * a / s^2))
wiener_mean_dt <- function(v, a, z, s = 1) (a * wiener_p_upper(v, a, z, s) - z) / v
wiener_mean_dt_upper <- function(v, a, z, s = 1) (a / tanh(a * v / s^2) - z / tanh(z * v / s^2)) / v
wiener_mean_dt_lower <- function(v, a, z, s = 1) wiener_mean_dt_upper(-v, a, a - z, s)

# Gaussian SDT with d_a sensitivity d, thresholds from the noise mean, and the
# signal standard deviation sdratio: category k is the interval between
# thresholds k - 1 and k of the noise axis
sdt_category_probs <- function(d, thresholds, stimulus, sdratio = 1) {
  shift <- d * sqrt((1 + sdratio^2) / 2) * stimulus
  scale <- if (stimulus == 1) sdratio else 1
  diff(pnorm((c(-Inf, thresholds, Inf) - shift) / scale))
}
unit_counts <- function(k, K) replace(numeric(K), k, 1)

# --- sdm ----------------------------------------------------------------------

test_that("dsdm() recycles x and the parameters to the longest argument", {
  expect_elementwise(dsdm, list(x = c(0.1, -0.5), mu = c(0, 1, 2), c = c(2, 3, 4, 5), kappa = 3))
  expect_elementwise(dsdm, list(x = 0.3, mu = 0, c = c(2, 3), kappa = c(1, 5, 20)),
                     fixed = list(log = TRUE))
  expect_elementwise(dsdm, list(x = c(0.1, -0.5), mu = c(0, 1, 2), c = 3, kappa = 3),
                     fixed = list(parametrization = "bessel"))
})

test_that("rsdm() draws each observation from its own location, strength and precision", {
  withr::local_seed(4441)
  mu <- c(-2, 0, 2)
  c_ <- c(2, 4, 3)
  kappa <- c(2, 10, 5)
  n <- 3000
  g <- group_of(n, 3)
  x <- rsdm(n, mu = mu, c = c_, kappa = kappa)
  expect_length(x, n)
  for (k in 1:3) {
    # E[cos(x - mu)] under the group's density, by quadrature of dsdm()
    ref <- integrate(function(t) cos(t - mu[k]) * dsdm(t, mu[k], c_[k], kappa[k]), -pi, pi)$value
    obs <- cos(x[g == k] - mu[k])
    expect_mean_near(obs, ref, sd(obs) / sqrt(length(obs)))
  }
  # n not a multiple of the parameter length: draw i takes element i, and a
  # sharp density puts every draw on its location's side of zero
  x <- NULL
  expect_no_warning(x <- rsdm(5, mu = c(-2, 2), c = 10, kappa = 500))
  expect_equal(sign(x), rep_len(c(-1, 1), 5))
})

test_that("rejection_sampling() is deprecated and still samples", {
  withr::local_seed(4441)
  x <- NULL
  expect_warning(
    x <- rejection_sampling(10, stats::dunif, 1, stats::runif),
    "deprecated since bmm 1.5.0"
  )
  expect_length(x, 10)
})

# --- ddm ----------------------------------------------------------------------

test_that("dddm() recycles the observations and the parameters to the longest argument", {
  expect_elementwise(dddm, list(rt = c(0.6, 0.7), response = c(1, 0), drift = 1:3,
                                bound = c(1, 1.5), ndt = 0.2, zr = c(0.4, 0.5, 0.6)))
  expect_elementwise(dddm, list(rt = 0.6, response = 1, drift = c(1, 2, 3), bound = 1, ndt = 0.2),
                     fixed = list(log = FALSE))
  expect_elementwise(dddm, list(rt = c(0.6, 0.7, 0.8), response = c("upper", "lower", "upper"),
                                drift = 1, bound = c(1, 1.5), ndt = 0.2))
})

test_that("rddm() returns n draws, draw i from the i-th parameter elements", {
  withr::local_seed(4442)
  drift <- c(0.5, 2, 1)
  bound <- c(1, 1.5, 2)
  ndt <- c(0.2, 0.3, 0.4)
  zr <- c(0.3, 0.5, 0.7)
  n <- 6000
  g <- group_of(n, 3)
  out <- rddm(n, drift, bound, ndt, zr)
  expect_equal(nrow(out), n)
  for (k in 1:3) {
    z <- zr[k] * bound[k]
    p_up <- wiener_p_upper(drift[k], bound[k], z)
    rows <- out[g == k, ]
    expect_mean_near(rows$response, p_up, se_prop(p_up, nrow(rows)))
    expect_mean_near(rows$rt, ndt[k] + wiener_mean_dt(drift[k], bound[k], z), sd(rows$rt) / sqrt(nrow(rows)))
  }
  expect_equal(nrow(rddm(2, drift = c(1, 2, 3), bound = 1, ndt = 0.2)), 2)
  # n not a multiple of the parameter length: draw i takes element i, and a
  # strong drift decides the boundary
  expect_equal(rddm(7, drift = c(8, -8, 8), bound = 1, ndt = 0.2)$response, rep_len(c(1, 0, 1), 7))
})

test_that("rddm(1, <vectors>) keeps one draw per element until 1.6.0, with a warning", {
  withr::local_seed(4442)
  out <- NULL
  expect_warning(out <- rddm(1, drift = c(1, 2, 3), bound = 1, ndt = 0.2), "deprecated since bmm 1.5.0")
  expect_equal(nrow(out), 3)
  expect_no_warning(expect_equal(nrow(rddm(1, drift = 1, bound = 1, ndt = 0.2)), 1))
})

# --- cswald (tests only while #406 is open; multiples of each other) ----------

test_that("dcswald() recycles the observations and the parameters to the longest argument", {
  expect_elementwise(dcswald, list(rt = 0.6, response = 1, drift = c(1, 2, 3), bound = c(1.5, 1, 2), ndt = 0.2))
  expect_elementwise(dcswald, list(rt = c(0.6, 0.7, 0.8), response = c(1, 0, 1), drift = c(1, 2, 3),
                                   bound = 1.5, ndt = c(0.2, 0.3, 0.1), zr = c(0.4, 0.5, 0.6), s = c(1, 1.2, 0.9)),
                     fixed = list(version = "crisk", log = FALSE))
})

test_that("rcswald() returns n draws, draw i from the i-th parameter elements", {
  withr::local_seed(4443)
  drift <- c(0.5, 2, 1)
  bound <- c(1, 1.5, 2)
  ndt <- c(0.2, 0.3, 0.4)
  zr <- c(0.3, 0.5, 0.7)
  s <- c(1, 1.5, 0.8)
  n <- 6000
  g <- group_of(n, 3)
  out <- rcswald(n, drift, bound, ndt, zr, s)
  expect_equal(nrow(out), n)
  for (k in 1:3) {
    z <- zr[k] * bound[k]
    p_up <- wiener_p_upper(drift[k], bound[k], z, s[k])
    rows <- out[g == k, ]
    expect_mean_near(rows$response, p_up, se_prop(p_up, nrow(rows)))
    expect_mean_near(rows$rt, ndt[k] + wiener_mean_dt(drift[k], bound[k], z, s[k]), sd(rows$rt) / sqrt(nrow(rows)))
  }
})

# --- ezdm ---------------------------------------------------------------------

test_that("dezdm() recycles the summary statistics and the parameters to the longest argument", {
  expect_elementwise(dezdm, list(mean_rt = c(0.5, 0.55), var_rt = c(0.02, 0.025), n_upper = c(80, 75),
                                 n_trials = 100, drift = c(2, 1.5, 1), bound = 1.5, ndt = c(0.3, 0.25)))
  expect_elementwise(dezdm, list(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100, drift = 2,
                                 bound = c(1, 1.5, 2), ndt = 0.3, s = c(1, 1.2)),
                     fixed = list(log = FALSE))
  mean_rt <- rbind(c(0.45, 0.55), c(0.5, 0.6))
  var_rt <- rbind(c(0.018, 0.025), c(0.02, 0.03))
  expect_elementwise(dezdm, list(mean_rt = mean_rt, var_rt = var_rt, n_upper = c(80, 60, 70), n_trials = 100,
                                 drift = c(2, 1.5, 1), bound = 1.5, ndt = 0.3, zr = c(0.55, 0.45)),
                     fixed = list(version = "4par"))
})

test_that("rezdm() returns n draws, draw i from the i-th parameter elements", {
  withr::local_seed(4444)
  drift <- c(1, 2, 0.5)
  bound <- c(1, 1.5, 2)
  ndt <- c(0.2, 0.3, 0.4)
  n_trials <- c(50, 100, 80)
  s <- c(1, 1.2, 0.9)
  n <- 1500
  g <- group_of(n, 3)
  out <- rezdm(n, n_trials, drift, bound, ndt, s = s)
  expect_equal(nrow(out), n)
  for (k in 1:3) {
    rows <- out[g == k, ]
    expect_equal(unique(rows$n_trials), n_trials[k])
    # a symmetric start point: P(upper) = plogis(v a / s^2), mean decision time a / (2 v) tanh(a v / (2 s^2))
    p_up <- plogis(drift[k] * bound[k] / s[k]^2)
    expect_mean_near(rows$n_upper / n_trials[k], p_up, se_prop(p_up, n_trials[k] * nrow(rows)))
    m_rt <- ndt[k] + bound[k] / (2 * drift[k]) * tanh(bound[k] * drift[k] / (2 * s[k]^2))
    expect_mean_near(rows$mean_rt, m_rt, sd(rows$mean_rt) / sqrt(nrow(rows)))
  }
})

test_that("rezdm(version = '4par') pairs the start point and both boundary means per draw", {
  withr::local_seed(4444)
  drift <- c(0.5, 1, 0.8)
  bound <- c(1, 1.2, 1.5)
  ndt <- c(0.2, 0.3, 0.4)
  zr <- c(0.4, 0.5, 0.6)
  n_trials <- c(100, 150, 120)
  n <- 1500
  g <- group_of(n, 3)
  out <- rezdm(n, n_trials, drift, bound, ndt, zr = zr, version = "4par")
  expect_equal(nrow(out), n)
  for (k in 1:3) {
    rows <- out[g == k, ]
    z <- zr[k] * bound[k]
    p_up <- wiener_p_upper(drift[k], bound[k], z)
    expect_mean_near(rows$n_upper / n_trials[k], p_up, se_prop(p_up, n_trials[k] * nrow(rows)))
    upper <- rows$mean_rt_upper
    lower <- rows$mean_rt_lower
    expect_mean_near(upper[!is.na(upper)], ndt[k] + wiener_mean_dt_upper(drift[k], bound[k], z),
                     sd(upper, na.rm = TRUE) / sqrt(sum(!is.na(upper))))
    expect_mean_near(lower[!is.na(lower)], ndt[k] + wiener_mean_dt_lower(drift[k], bound[k], z),
                     sd(lower, na.rm = TRUE) / sqrt(sum(!is.na(lower))))
  }
})

# --- SDT ----------------------------------------------------------------------

thr3 <- rbind(c(0.2, 0.7, 1.2), c(0.3, 0.8, 1.3), c(0.4, 0.9, 1.4))
counts2 <- rbind(c(2, 8, 20, 70), c(30, 40, 20, 10))

test_that("the SDT d*() functions recycle counts, thresholds and parameters to the longest argument", {
  expect_elementwise(dsdt_yn, list(n_old = c(50, 60), n_trials = 100, stimulus = c(1, 0, 1), d = 1:3,
                                   criterion = c(0.5, 0.7)))
  expect_elementwise(dsdt_yn, list(n_old = 50, n_trials = c(100, 80), stimulus = 1, d = c(1, 2),
                                   bias = c(0.1, -0.1, 0.2), sdratio = c(1, 1.5)),
                     fixed = list(dist = "logistic", log = TRUE))
  expect_elementwise(dsdt_mafc, list(n_correct = c(50, 60), n_trials = c(100, 80, 90), m = c(2, 4), d = 1:3))
  expect_elementwise(dsdt_rating, list(counts = counts2, stimulus = c(1, 0, 1), d = c(1, 1.5, 2),
                                       thresholds = thr3, sdratio = c(1, 1.2)))
  expect_elementwise(dsdt_rating, list(counts = rbind(c(2, 8, 20, 70)), stimulus = 1, d = 1.5, thresholds = thr3),
                     fixed = list(log = TRUE))
  # a vector of counts is one observation, scored under each row of thresholds
  expect_equal(
    dsdt_rating(c(2, 8, 20, 70), 1, d = 1.5, thresholds = thr3),
    vapply(1:3, function(i) dsdt_rating(c(2, 8, 20, 70), 1, d = 1.5, thresholds = thr3[i, ]), numeric(1))
  )
  expect_elementwise(dsdt_dpsdt, list(counts = counts2, stimulus = c(1, 0, 1), d = c(1, 1.5, 2),
                                      thresholds = thr3, Ro = c(0.2, 0.3, 0.4), Rn = c(0, 0.1)))
  expect_elementwise(dsdt_metad, list(counts = counts2, stimulus = c(1, 0, 1), d = c(1, 1.5, 2),
                                      thresholds = thr3, metad = c(0.8, 1, 1.2), sdratio = c(1, 1.2)))
  expect_elementwise(dsdt_ranking, list(counts = rbind(c(50, 30, 15, 5), c(20, 30, 30, 20)), m = 4,
                                        d = c(0.5, 1, 2), sdratio = c(1, 1.5)),
                     fixed = list(dist = "normal"))
  expect_elementwise(dsdt_ranking, list(counts = rbind(c(50, 30, 15, 5), c(60, 40, 0, 0)), m = c(4, 2),
                                        d = c(0.5, 1, 2)))
  thr5 <- rbind(c(-1.1, -0.5, 0, 0.6, 1.3), c(-1, -0.4, 0.1, 0.7, 1.4), c(-1.2, -0.6, -0.1, 0.5, 1.2))
  counts_cdp <- rbind(c(10, 20, 30, 10, 15, 5, 5, 3, 2), c(30, 20, 10, 5, 5, 10, 2, 3, 15))
  expect_elementwise(dsdt_cdp, list(counts = counts_cdp, stimulus = c(1, 0, 1), dfam = c(0.8, 1, 0.6),
                                    drec = c(1, 0.8), thresholds = thr5, rcrit = c(0.5, 0.4, 0.6),
                                    sigmar = c(0, 0.1), rho = c(0, 0.2, -0.1)),
                     fixed = list(n_new = 3))
  counts_cdp_k <- rbind(c(10, 20, 30, 10, 15, 5, 5, 3, 2, 1, 2, 1), c(30, 20, 10, 5, 5, 10, 2, 3, 15, 2, 1, 1))
  expect_elementwise(dsdt_cdp, list(counts = counts_cdp_k, stimulus = 1, dfam = 0.8, drec = 1,
                                    thresholds = rbind(c(-1.1, -0.5, 0, 0.6, 1.3)), rcrit = 0.5,
                                    kcrit = c(0.3, 0.2, 0.4)),
                     fixed = list(n_new = 3))
})

test_that("rsdt_yn() pairs sensitivity, criterion, variance ratio and stimulus per draw", {
  withr::local_seed(4445)
  d <- c(0.5, 1.5, 2.5)
  criterion <- c(0.2, 0.8, 0.5)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  sdratio <- c(1, 1.5, 1)
  n_trials <- c(50, 100, 80)
  n <- 1500
  g <- group_of(n, 3)
  n_old <- rsdt_yn(n, n_trials, stimulus, d, criterion, sdratio = sdratio)
  expect_length(n_old, n)
  for (k in 1:3) {
    p_old <- pnorm((d[k] * sqrt((1 + sdratio[k]^2) / 2) * stimulus[k] - criterion[k]) / sdratio[k]^stimulus[k])
    expect_mean_near(n_old[g == k] / n_trials[k], p_old, se_prop(p_old, n_trials[k] * sum(g == k)))
  }
  # n not a multiple of the parameter length: draw i takes element i, and a
  # criterion far from the signal mean makes every response certain
  expect_no_warning(n_old <- rsdt_yn(7, 100, 1, d = 10, criterion = c(-5, 15, -5)))
  expect_equal(n_old, rep_len(c(100L, 0L, 100L), 7))
})

test_that("rsdt_mafc() pairs sensitivity and the number of alternatives per draw", {
  withr::local_seed(4445)
  d <- c(0.5, 1.5, 1)
  m <- c(2, 4, 3)
  n_trials <- c(50, 100, 80)
  n <- 1500
  g <- group_of(n, 3)
  n_correct <- rsdt_mafc(n, n_trials, m, d)
  expect_length(n_correct, n)
  for (k in 1:3) {
    p_correct <- integrate(function(x) dnorm(x - d[k]) * pnorm(x)^(m[k] - 1), -Inf, Inf)$value
    expect_mean_near(n_correct[g == k] / n_trials[k], p_correct, se_prop(p_correct, n_trials[k] * sum(g == k)))
  }
})

# category proportions of a group against the reference probabilities
expect_category_props <- function(counts, n_trials, probs) {
  for (j in seq_along(probs)) {
    expect_mean_near(counts[, j] / n_trials, probs[j], se_prop(probs[j], n_trials * nrow(counts)))
  }
}

test_that("rsdt_rating() pairs sensitivity, stimulus, thresholds and variance ratio per draw", {
  withr::local_seed(4445)
  d <- c(0.5, 1.5, 1)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  sdratio <- c(1, 1.5, 1.2)
  n_trials <- 100
  n <- 1500
  g <- group_of(n, 3)
  counts <- rsdt_rating(n, n_trials, stimulus, d, thresholds = thr3, sdratio = sdratio)
  expect_equal(dim(counts), c(n, 4))
  for (k in 1:3) {
    expect_category_props(counts[g == k, ], n_trials, sdt_category_probs(d[k], thr3[k, ], stimulus[k], sdratio[k]))
  }
  # n not a multiple of the threshold rows: row i takes threshold row i, and
  # thresholds far from the mean put every response in one category
  counts <- rsdt_rating(7, 100, 1, d = 0, thresholds = rbind(c(-9, -8, -7), c(7, 8, 9), c(-9, 8, 9)))
  expect_equal(max.col(counts), rep_len(c(4L, 1L, 2L), 7))
})

test_that("rsdt_dpsdt() pairs recollection with the rating parameters per draw", {
  withr::local_seed(4445)
  d <- c(0.5, 1.5, 1)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  Ro <- c(0.2, 0.5, 0.3)
  Rn <- c(0, 0.2, 0.1)
  n_trials <- 100
  n <- 1500
  g <- group_of(n, 3)
  counts <- rsdt_dpsdt(n, n_trials, stimulus, d, thresholds = thr3, Ro = Ro, Rn = Rn)
  for (k in 1:3) {
    p <- sdt_category_probs(d[k], thr3[k, ], stimulus[k])
    # recollection loads the most confident category of the stimulus
    probs <- if (stimulus[k] == 1) (1 - Ro[k]) * p + Ro[k] * unit_counts(4, 4) else (1 - Rn[k]) * p + Rn[k] * unit_counts(1, 4)
    expect_category_props(counts[g == k, ], n_trials, probs)
  }
})

test_that("rsdt_metad() pairs metacognitive sensitivity with the rating parameters per draw", {
  withr::local_seed(4445)
  d <- c(1, 1, 1.5)
  metad <- c(0.5, 1.5, 1)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  n_trials <- 100
  n <- 1500
  g <- group_of(n, 3)
  counts <- rsdt_metad(n, n_trials, stimulus, d, thresholds = thr3, metad = metad)
  for (k in 1:3) {
    probs <- vapply(1:4, function(j) {
      dsdt_metad(unit_counts(j, 4), stimulus[k], d[k], thr3[k, ], metad = metad[k])
    }, numeric(1))
    expect_category_props(counts[g == k, ], n_trials, probs)
  }
})

test_that("rsdt_ranking() pairs sensitivity and set size per draw", {
  withr::local_seed(4445)
  d <- c(0.8, 1.5, 1)
  m <- c(2, 4, 3)
  n_trials <- 100
  n <- 1500
  g <- group_of(n, 3)
  counts <- rsdt_ranking(n, n_trials, m, d)
  expect_equal(dim(counts), c(n, 4))
  expect_true(all(counts[g == 1, 3:4] == 0))
  # Gumbel-min with two items: the target is ranked first with probability plogis(d)
  expect_mean_near(counts[g == 1, 1] / n_trials, plogis(d[1]), se_prop(plogis(d[1]), n_trials * sum(g == 1)))
  for (k in 2:3) {
    probs <- vapply(seq_len(m[k]), function(j) {
      dsdt_ranking(unit_counts(j, m[k]), m[k], d[k])
    }, numeric(1))
    expect_category_props(counts[g == k, seq_len(m[k])], n_trials, probs)
  }
})

test_that("rsdt_cdp() pairs familiarity, recollection, criterion and thresholds per draw", {
  withr::local_seed(4445)
  dfam <- c(0.8, 1.2, 0.5)
  drec <- c(1, 0.5, 1.5)
  rcrit <- c(0.5, 0.3, 0.7)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  thr5 <- rbind(c(-1.1, -0.5, 0, 0.6, 1.3), c(-1, -0.4, 0.1, 0.7, 1.4), c(-1.2, -0.6, -0.1, 0.5, 1.2))
  n_trials <- 100
  n <- 900
  g <- group_of(n, 3)
  counts <- rsdt_cdp(n, n_trials, stimulus, dfam, drec, thresholds = thr5, rcrit = rcrit, n_new = 3)
  expect_equal(dim(counts), c(n, 9))
  for (k in 1:3) {
    probs <- vapply(1:9, function(j) {
      dsdt_cdp(unit_counts(j, 9), stimulus[k], dfam[k], drec[k], thr5[k, ], rcrit = rcrit[k], n_new = 3)
    }, numeric(1))
    expect_category_props(counts[g == k, ], n_trials, probs)
  }
})
