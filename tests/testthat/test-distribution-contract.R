# The recycling contract of the r*() and d*() functions (#444, ?bmm_distributions):
# r*(n, ...) returns n draws, draw i from the i-th element (row) of every
# recycled argument; d*(x, ...) returns one value per element of the longest
# argument, value i from the i-th elements alone. Every r*() pairing test
# interleaves parameter groups and compares each group with a reference
# outside the function's own code (a closed form, quadrature, or the single-
# trial density of the matching d*()). Every per-observation argument of the
# sdm, ddm, ezdm and SDT d*() and r*() functions is varied in at least one
# test with a vector (matrix) whose recycled form reads differently reversed:
# a mutant that reverses any one of them after recycling (the rezdm() ones
# once per version) fails a test in this file or in test-distributions.R.
# Mutants of the samplers (vector reversed, recycling removed, one group
# shifted 10 %) were run against these tests when they were written; a Monte
# Carlo test is not expected to see a 10 % shift, so each sampler also gets a
# deterministic check: certain outcomes, n = 7 and parameter lengths 2, 3
# and 4. Recycled to 7, a vector of length 3
# or 4 is chosen so that reversing it changes the expected draws; one of
# length 2 reads the same reversed (a b a b a b a) and is caught only when
# shifted. Non-multiple lengths recycle silently, as rnorm() does; the cswald
# tests use multiples only until #406 lands. The circular mixtures are
# circmix's (#434).

pick <- function(x, i) {
  if (is.matrix(x)) x[(i - 1) %% nrow(x) + 1, ] else x[(i - 1) %% length(x) + 1]
}
obs_length <- function(x) if (is.matrix(x)) nrow(x) else length(x)

# the vectorised call equals the calls on the i-th elements alone, without
# recycling warnings. A case can fail only if every varied argument changes
# when reversed (a b a at n = 3 does not) and every reference value is finite
# and non-zero (a density that underflows to 0 needs log = TRUE).
expect_elementwise <- function(f, args, fixed = list()) {
  n <- max(vapply(args, obs_length, integer(1)))
  palindromes <- Filter(function(a) {
    rows <- lapply(seq_len(n), pick, x = args[[a]])
    length(unique(rows)) > 1 && identical(rev(rows), rows)
  }, names(args))
  expect(length(palindromes) == 0,
         paste("recycled to", n, "these read the same reversed:", toString(palindromes)))
  vec <- NULL
  expect_no_warning(vec <- do.call(f, c(args, fixed)))
  each <- vapply(seq_len(n), function(i) {
    do.call(f, c(lapply(args, pick, i = i), fixed))
  }, numeric(1))
  expect_length(vec, n)
  expect(all(is.finite(each) & each != 0), "a reference value is zero or not finite")
  # as ratios, so that each element is compared on its own scale:
  # expect_equal(vec, each) scales its tolerance to the mean over the vector,
  # and one large density hides an error of many orders in a small one
  expect_equal(vec / each, rep(1, n))
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

# the deprecated `parametrization` (#558) warns on every call
dsdm_bessel <- function(...) {
  withCallingHandlers(dsdm(..., parametrization = "bessel"), warning = function(w) {
    if (grepl("`parametrization = ", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}

test_that("dsdm() recycles x and the parameters to the longest argument", {
  expect_elementwise(dsdm, list(x = c(0.1, -0.5), mu = c(0, 1, 2), c = c(2, 3, 4, 5), kappa = 3))
  expect_elementwise(dsdm, list(x = 0.3, mu = 0, c = c(2, 3), kappa = c(1, 5, 20, 10)),
                     fixed = list(log = TRUE))
  expect_elementwise(dsdm_bessel, list(x = c(0.1, -0.5), mu = c(0, 1, 2, -1), c = 3, kappa = 3))
  # parameter sets that agree in their first four digits get their own integral;
  # x longer than c, so that sets repeat and dsdm() looks them up by key
  expect_elementwise(dsdm, list(x = c(0.1, 0.2, 0.3, 0.4), mu = 0, c = c(3, 3.0001), kappa = 2))
})

test_that("dsdm() recycles alike for every activation and scaling", {
  cells <- expand.grid(activation = c("vonmises", "laplace"), scaling = c("density", "peak", "range"),
                       stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cells))) {
    cell <- as.list(cells[i, ])
    expect_elementwise(dsdm, list(x = c(0.1, -0.5), mu = c(0, 1, 2), c = c(2, 3, 4, 5), kappa = c(1, 5, 20)),
                       fixed = c(cell, log = TRUE))
    # sets repeat, and two of them differ in kappa only: under the peak scaling
    # their c* is the same, so the key must hold kappa as well
    expect_elementwise(dsdm, list(x = seq(0.1, 0.8, by = 0.1), mu = 0, c = c(3, 3.0001), kappa = c(1, 5, 2, 4)),
                       fixed = c(cell, log = TRUE))
  }
})

test_that("rsdm() draws each observation from its own location, strength and precision", {
  withr::local_seed(4441)
  mu <- c(-2, 0, 2)
  c_ <- c(2, 4, 3)
  kappa <- c(2, 10, 5)
  n <- 3000
  g <- group_of(n, 3)
  # a density sharp enough below to keep every draw on its location's side, and
  # wide enough for the sampler: c* is about 89 (von Mises) and 100 (Laplace)
  sharp <- list(vonmises = list(c = 10, kappa = 500), laplace = list(c = 10, kappa = 20))
  for (activation in c("vonmises", "laplace")) {
    cell <- list(activation = activation, scaling = "density")
    x <- do.call(rsdm, c(list(n, mu = mu, c = c_, kappa = kappa), cell))
    expect_length(x, n)
    for (k in 1:3) {
      # E[cos(x - mu)] under the group's density, by quadrature of dsdm()
      dens <- function(t) do.call(dsdm, c(list(t, mu[k], c_[k], kappa[k]), cell))
      ref <- integrate(function(t) cos(t - mu[k]) * dens(t), -pi, pi)$value
      obs <- cos(x[g == k] - mu[k])
      expect_mean_near(obs, ref, sd(obs) / sqrt(length(obs)))
    }
    # c longer than kappa: c* pairs them draw by draw
    x <- do.call(rsdm, c(list(n, mu = 0, c = c(0.5, 10), kappa = 5), cell))
    for (k in 1:2) {
      dens <- function(t) do.call(dsdm, c(list(t, 0, c(0.5, 10)[k], 5), cell))
      ref <- integrate(function(t) cos(t) * dens(t), -pi, pi)$value
      obs <- cos(x[group_of(n, 2) == k])
      expect_mean_near(obs, ref, sd(obs) / sqrt(length(obs)))
    }
    # n not a multiple of the parameter length: draw i takes element i, and a
    # sharp density puts every draw on its location's side of zero
    x <- NULL
    expect_no_warning(x <- do.call(rsdm, c(list(5, mu = c(-2, 2)), sharp[[activation]], cell)))
    expect_equal(sign(x), rep_len(c(-1, 1), 5))
    # n = 2, as in posterior_predict(ndraws = 2), pairs as well
    x <- replicate(10, do.call(rsdm, c(list(2, mu = c(-2, 2)), sharp[[activation]], cell)))
    expect_equal(sign(x), matrix(c(-1, 1), 2, 10))
  }
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
  expect_elementwise(dddm, list(rt = c(0.6, 0.7), response = c(1, 0), drift = 1:4,
                                bound = c(1, 1.5), ndt = c(0.2, 0.3, 0.25), zr = c(0.4, 0.5, 0.6)))
  expect_elementwise(dddm, list(rt = 0.6, response = 1, drift = c(1, 2, 3), bound = 1, ndt = 0.2),
                     fixed = list(log = FALSE))
  expect_elementwise(dddm, list(rt = c(0.6, 0.7, 0.8, 0.65), response = c("upper", "upper", "lower"),
                                drift = 1, bound = c(1, 1.5), ndt = 0.2))
  expect_elementwise(dddm, list(rt = 0.6, response = c(1, 0), drift = 1, bound = 1, ndt = 0.2))
  expect_elementwise(dddm, list(rt = c(0.6, 0.7, 0.8, 0.9), response = c("upper", "lower"),
                                drift = c(1, 2), bound = 1, ndt = 0.2))
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
  expect_warning(out <- rddm(1, drift = c(1, 2, 3), bound = 1, ndt = 0.2),
                 "`rddm(1, drift = c(1, 2, 3), ...)` is deprecated since bmm 1.5.0", fixed = TRUE)
  expect_equal(nrow(out), 3)
  expect_warning(out <- rddm(1, drift = 1, bound = 1, ndt = 0.2, zr = c(0.4, 0.6)),
                 "Use `rddm(2, zr = c(0.4, 0.6), ...)` instead", fixed = TRUE)
  expect_equal(nrow(out), 2)
  expect_no_warning(expect_equal(nrow(rddm(1, drift = 1, bound = 1, ndt = 0.2)), 1))
})

test_that("the rddm() deprecation warning shows a long literal argument by its length", {
  # do.call() passes values, not the code that made them; deparsed whole, 2000
  # drifts outran warning.length and cut off the suggested call (100 here keep
  # the test fast)
  withr::local_seed(4443)
  expect_warning(do.call(rddm, list(n = 1, drift = seq(0.1, 3, length.out = 100), bound = 1, ndt = 0.2)),
                 "Use `rddm(100, drift = <100 values>, ...)` instead", fixed = TRUE)
  expect_warning(do.call(rddm, list(n = seq(0.1, 1, length.out = 50), drift = 1, bound = 1, ndt = 0.2)),
                 "`rddm(<50 values>, ...)` is deprecated since bmm 1.5.0", fixed = TRUE)
  expect_warning(rddm(1, drift = sapply(1:2, function(i) {
    i
  }), bound = 1, ndt = 0.2), "rddm(1, drift = sapply(1:2, function(i) { i }), ...)", fixed = TRUE)
  # a literal is shown whole up to 60 characters of its deparse (this one has
  # 60) and by its length above (61)
  expect_warning(do.call(rddm, list(n = 1, drift = c(1, 11:24), bound = 1, ndt = 0.2)),
                 "Use `rddm(15, drift = c(1, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24), ...)`",
                 fixed = TRUE)
  expect_warning(do.call(rddm, list(n = 1, drift = c(10, 11:24), bound = 1, ndt = 0.2)),
                 "Use `rddm(15, drift = <15 values>, ...)`", fixed = TRUE)
})

test_that("the rddm() deprecation warning shows long typed code as typed", {
  withr::local_seed(4443)
  expect_warning(
    rddm(1, drift = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.5), bound = 1, ndt = 0.2),
    "Use `rddm(13, drift = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.5), ...)` instead",
    fixed = TRUE
  )
  expect_warning(rddm(c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.5), drift = 1, bound = 1,
                      ndt = 0.2),
                 "`rddm(c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.5), ...)` is deprecated",
                 fixed = TRUE)
})

# --- cswald (tests only while #406 is open; multiples of each other) ----------

test_that("dcswald() recycles the observations and the parameters to the longest argument", {
  expect_elementwise(dcswald, list(rt = 0.6, response = 1, drift = c(1, 2, 3), bound = c(1.5, 1, 2), ndt = 0.2))
  expect_elementwise(dcswald, list(rt = c(0.6, 0.7, 0.8), response = c(0, 1, 1), drift = c(1, 2, 3),
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
                                 n_trials = 100, drift = c(2, 1.5, 1, 1.8), bound = 1.5, ndt = c(0.3, 0.25)))
  expect_elementwise(dezdm, list(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100, drift = 2,
                                 bound = c(1, 1.5, 2), ndt = 0.3, s = c(1, 1.2, 0.9)),
                     fixed = list(log = FALSE))
  mean_rt <- rbind(c(0.45, 0.55), c(0.5, 0.6))
  var_rt <- rbind(c(0.018, 0.025), c(0.02, 0.03))
  expect_elementwise(dezdm, list(mean_rt = mean_rt, var_rt = var_rt, n_upper = c(80, 60, 70, 75), n_trials = 100,
                                 drift = c(2, 1.5, 1), bound = 1.5, ndt = 0.3, zr = c(0.55, 0.45)),
                     fixed = list(version = "4par"))
})

test_that("dezdm() recycles parameters of non-multiple lengths to a longer set of observations", {
  # the counts are checked after recycling, so non-multiple lengths raise no warning
  expect_elementwise(dezdm, list(mean_rt = 0.6, var_rt = 0.04, n_upper = c(20, 22, 18, 21), n_trials = c(40, 42, 45),
                                 drift = c(0.8, 1, 1.2), bound = 1.4, ndt = c(0.25, 0.3)))
  # zr is not a parameter of version "3par" and does not count towards the length
  expect_elementwise(dezdm, list(mean_rt = c(0.6, 0.65), var_rt = 0.04, n_upper = 20, n_trials = 40, drift = 0.8,
                                 bound = 1.4, ndt = 0.25),
                     fixed = list(zr = numeric(0)))
  expect_elementwise(dezdm, list(mean_rt = rbind(c(0.55, 0.62)), var_rt = rbind(c(0.03, 0.05)),
                                 n_upper = c(20, 22, 18, 21), n_trials = c(40, 42, 45), drift = 0.8,
                                 bound = c(1.2, 1.5, 1.8), ndt = 0.25, zr = c(0.4, 0.6, 0.5)),
                     fixed = list(version = "4par"))
  expect_elementwise(dezdm, list(mean_rt = rbind(c(0.55, 0.62)), var_rt = rbind(c(0.03, 0.05)),
                                 n_upper = c(20, 22, 18, 21), n_trials = 40, drift = 0.8, bound = 1.4,
                                 ndt = 0.25, zr = 0.45, s = c(1, 1.3, 0.8)),
                     fixed = list(version = "4par"))
  expect_elementwise(dezdm, list(mean_rt = rbind(c(0.55, 0.62)), var_rt = rbind(c(0.03, 0.05), c(0.06, 0.02)),
                                 n_upper = 20, n_trials = 40, drift = 0.8, bound = 1.4, ndt = 0.25),
                     fixed = list(version = "4par"))
})

test_that("dezdm() reads a length-2 vector beside a matrix as one 4par observation", {
  pair_mean <- c(0.55, 0.62)
  pair_var <- c(0.03, 0.05)
  var_rt <- rbind(c(0.03, 0.05), c(0.06, 0.02))
  expect_elementwise(dezdm, list(var_rt = var_rt, n_upper = 20, n_trials = 40,
                                 drift = 0.8, bound = 1.4, ndt = 0.25),
                     fixed = list(mean_rt = pair_mean, version = "4par"))
  expect_elementwise(dezdm, list(mean_rt = rbind(c(0.55, 0.62), c(0.5, 0.6)), n_upper = 20, n_trials = 40,
                                 drift = 0.8, bound = 1.4, ndt = 0.25),
                     fixed = list(var_rt = pair_var, version = "4par"))
  expect_equal(
    dezdm(pair_mean, var_rt, 20, 40, drift = 0.8, bound = 1.4, ndt = 0.25, version = "4par"),
    dezdm(rbind(c(0.55, 0.62), c(0.55, 0.62)), var_rt, 20, 40, drift = 0.8, bound = 1.4, ndt = 0.25, version = "4par")
  )
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
  # n not a multiple of the parameter lengths: draw i takes element i, a
  # strong drift decides every response and ndt sets the whole part of mean_rt
  out <- NULL
  expect_no_warning(out <- rezdm(7, c(30, 40, 50, 60), drift = c(20, 20, -20), bound = 2, ndt = c(0.1, 5)))
  expect_equal(out$n_trials, c(30, 40, 50, 60, 30, 40, 50))
  expect_equal(out$n_upper, c(30, 40, 0, 60, 30, 0, 50))
  expect_equal(floor(out$mean_rt), c(0, 5, 0, 5, 0, 5, 0))
})

test_that("rezdm(version = '4par') pairs the start point, s and both boundary means per draw", {
  withr::local_seed(4444)
  drift <- c(0.5, 1, 0.8)
  bound <- c(1, 1.2, 1.5)
  ndt <- c(0.2, 0.3, 0.4)
  zr <- c(0.4, 0.5, 0.6)
  s <- c(1, 1.3, 0.8)
  n_trials <- c(100, 150, 120)
  n <- 1500
  g <- group_of(n, 3)
  out <- rezdm(n, n_trials, drift, bound, ndt, zr = zr, s = s, version = "4par")
  expect_equal(nrow(out), n)
  for (k in 1:3) {
    rows <- out[g == k, ]
    z <- zr[k] * bound[k]
    p_up <- wiener_p_upper(drift[k], bound[k], z, s[k])
    expect_mean_near(rows$n_upper / n_trials[k], p_up, se_prop(p_up, n_trials[k] * nrow(rows)))
    upper <- rows$mean_rt_upper
    lower <- rows$mean_rt_lower
    expect_mean_near(upper[!is.na(upper)], ndt[k] + wiener_mean_dt_upper(drift[k], bound[k], z, s[k]),
                     sd(upper, na.rm = TRUE) / sqrt(sum(!is.na(upper))))
    expect_mean_near(lower[!is.na(lower)], ndt[k] + wiener_mean_dt_lower(drift[k], bound[k], z, s[k]),
                     sd(lower, na.rm = TRUE) / sqrt(sum(!is.na(lower))))
  }
  # n not a multiple of the parameter lengths: draw i takes element i, a
  # strong drift sends every response to one boundary and ndt sets the whole
  # part of that boundary's mean_rt
  out <- NULL
  expect_no_warning(out <- rezdm(7, c(30, 40, 50, 60), drift = c(20, 20, -20), bound = 2, ndt = c(0.1, 5),
                                 version = "4par"))
  expect_equal(out$n_upper, c(30, 40, 0, 60, 30, 0, 50))
  upper <- out$n_upper > 0
  expect_equal(floor(ifelse(upper, out$mean_rt_upper, out$mean_rt_lower)), c(0, 5, 0, 5, 0, 5, 0))
})

# --- SDT ----------------------------------------------------------------------

thr3 <- rbind(c(0.2, 0.7, 1.2), c(0.3, 0.8, 1.3), c(0.4, 0.9, 1.4))
counts2 <- rbind(c(2, 8, 20, 70), c(30, 40, 20, 10))

test_that("the binomial and rating SDT d*() functions recycle counts, thresholds and parameters", {
  expect_elementwise(dsdt_yn, list(n_old = c(50, 60), n_trials = 100, stimulus = c(1, 1, 0), d = c(1, 2, 3, 1.5),
                                   criterion = c(0.4, 0.7, 0.5)))
  expect_elementwise(dsdt_yn, list(n_old = 50, n_trials = c(100, 80, 90, 70), stimulus = 1, d = c(1, 2),
                                   bias = c(0.1, -0.1, 0.2), sdratio = c(1, 1.5)),
                     fixed = list(dist = "logistic", log = TRUE))
  expect_elementwise(dsdt_mafc, list(n_correct = c(50, 60), n_trials = c(100, 80, 90), m = c(2, 4),
                                     d = c(1, 2, 3, 1.5)))
  expect_elementwise(dsdt_rating, list(counts = counts2, stimulus = c(1, 1, 0), d = c(1, 1.5, 2, 0.8),
                                       thresholds = thr3, sdratio = c(1, 1.2)))
  expect_elementwise(dsdt_rating, list(counts = rbind(c(2, 8, 20, 70)), stimulus = 1, d = 1.5, thresholds = thr3),
                     fixed = list(log = TRUE))
  # a vector of counts is one observation, scored under each row of thresholds
  expect_equal(
    dsdt_rating(c(2, 8, 20, 70), 1, d = 1.5, thresholds = thr3),
    vapply(1:3, function(i) dsdt_rating(c(2, 8, 20, 70), 1, d = 1.5, thresholds = thr3[i, ]), numeric(1))
  )
  expect_elementwise(dsdt_dpsdt, list(counts = counts2, stimulus = c(1, 1, 0), d = c(1, 1.5, 2, 0.8),
                                      thresholds = thr3, Ro = c(0.2, 0.3, 0.4), Rn = c(0, 0.1),
                                      sdratio = c(1, 1.3, 0.8)))
  # Rn acts on noise trials only, so three of them take three different values
  expect_elementwise(dsdt_dpsdt, list(counts = counts2, stimulus = c(0, 0, 1, 0), d = c(1, 1.5, 2, 0.8),
                                      thresholds = thr3, Ro = c(0.2, 0.3, 0.4), Rn = c(0.05, 0.3, 0.1, 0.2)))
  expect_elementwise(dsdt_metad, list(counts = counts2, stimulus = c(1, 1, 0), d = c(1, 1.5, 2, 0.8),
                                      thresholds = thr3, metad = c(0.8, 1, 1.2), sdratio = c(1, 1.2)))
})

test_that("dsdt_ranking() and dsdt_cdp() recycle counts, thresholds and parameters", {
  expect_elementwise(dsdt_ranking, list(counts = rbind(c(50, 30, 15, 5), c(20, 30, 30, 20)), m = 4,
                                        d = c(0.5, 1, 2, 1.5), sdratio = c(1, 1.5)),
                     fixed = list(dist = "normal"))
  expect_elementwise(dsdt_ranking, list(counts = rbind(c(50, 30, 15, 5), c(60, 40, 0, 0)), m = c(4, 2),
                                        d = c(0.5, 1, 2, 1.5)))
  thr5 <- rbind(c(-1.1, -0.5, 0, 0.6, 1.3), c(-1, -0.4, 0.1, 0.7, 1.4), c(-1.2, -0.6, -0.1, 0.5, 1.2))
  counts_cdp <- rbind(c(10, 20, 30, 10, 15, 5, 5, 3, 2), c(30, 20, 10, 5, 5, 10, 2, 3, 15))
  expect_elementwise(dsdt_cdp, list(counts = counts_cdp, stimulus = c(1, 1, 0), dfam = c(0.8, 1, 0.6, 1.2),
                                    drec = c(1, 0.8), thresholds = thr5, rcrit = c(0.5, 0.4, 0.6),
                                    sigmar = c(0, 0.1), rho = c(0, 0.2, -0.1)),
                     fixed = list(n_new = 3))
  counts_cdp_k <- rbind(c(10, 20, 30, 10, 15, 5, 5, 3, 2, 1, 2, 1), c(30, 20, 10, 5, 5, 10, 2, 3, 15, 2, 1, 1))
  expect_elementwise(dsdt_cdp, list(counts = counts_cdp_k, stimulus = 1, dfam = 0.8, drec = 1,
                                    thresholds = rbind(c(-1.1, -0.5, 0, 0.6, 1.3)), rcrit = 0.5,
                                    kcrit = c(0.3, 0.2, 0.4, 0.1)),
                     fixed = list(n_new = 3, log = TRUE))
})

test_that("each per-observation argument alone sets the length of a d*() result", {
  # every argument of length 3 in turn, the others scalar (one-row matrices):
  # an argument left out of the length gives fewer than 3 values
  thr <- rbind(c(0.2, 0.7, 1.2))
  thr5 <- rbind(c(-1.1, -0.5, 0, 0.6, 1.3))
  c4 <- rbind(c(2, 8, 20, 70), c(30, 40, 20, 10), c(10, 20, 30, 40))
  c9 <- rbind(c(10, 20, 30, 10, 15, 5, 5, 3, 2), c(30, 20, 10, 5, 5, 10, 2, 3, 15),
              c(5, 5, 5, 20, 20, 20, 10, 10, 5))
  rating <- list(counts = c4, stimulus = c(0, 1, 1), d = c(0.5, 1, 1.5),
                 thresholds = rbind(thr, thr + 0.1, thr - 0.1), sdratio = c(1, 1.3, 0.8))
  case <- function(f, base, long, fixed = list()) list(f = f, base = base, long = long, fixed = fixed)
  cases <- list(
    dsdm = case(dsdm, list(x = 0.3, mu = 0, c = 3, kappa = 3),
                list(x = c(0.1, 0.3, -0.5), mu = c(0, 0.5, 1), c = c(2, 3, 4), kappa = c(1, 3, 8))),
    dddm = case(dddm, list(rt = 0.6, response = 1, drift = 1, bound = 1, ndt = 0.2, zr = 0.5),
                list(rt = c(0.5, 0.6, 0.8), response = c(0, 1, 1), drift = c(0.5, 1, 2), bound = c(1, 1.2, 1.5),
                     ndt = c(0.2, 0.25, 0.3), zr = c(0.4, 0.5, 0.6))),
    dezdm3 = case(dezdm, list(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100, drift = 2, bound = 1.5,
                              ndt = 0.3, s = 1),
                  list(mean_rt = c(0.45, 0.5, 0.55), var_rt = c(0.02, 0.025, 0.03), n_upper = c(70, 80, 75),
                       n_trials = c(100, 110, 120), drift = c(1.5, 2, 2.5), bound = c(1.2, 1.5, 1.8),
                       ndt = c(0.25, 0.3, 0.35), s = c(1, 1.1, 0.9))),
    dezdm4 = case(dezdm, list(mean_rt = rbind(c(0.45, 0.55)), var_rt = rbind(c(0.02, 0.03)), n_upper = 80,
                              n_trials = 100, drift = 1, bound = 1.5, ndt = 0.3, zr = 0.5, s = 1),
                  list(mean_rt = rbind(c(0.45, 0.55), c(0.5, 0.6), c(0.48, 0.52)),
                       var_rt = rbind(c(0.02, 0.03), c(0.025, 0.035), c(0.03, 0.02)), n_upper = c(70, 80, 75),
                       n_trials = c(100, 110, 120), drift = c(0.5, 1, 1.5), bound = c(1.2, 1.5, 1.8),
                       ndt = c(0.25, 0.3, 0.35), zr = c(0.4, 0.5, 0.6), s = c(1, 1.1, 0.9)),
                  list(version = "4par")),
    dsdt_yn = case(dsdt_yn, list(n_old = 50, n_trials = 100, stimulus = 1, d = 1, criterion = 0.5, sdratio = 1),
                   list(n_old = c(40, 50, 60), n_trials = c(100, 90, 80), stimulus = c(0, 1, 1), d = c(0.5, 1, 1.5),
                        criterion = c(0.2, 0.5, 0.8), sdratio = c(1, 1.3, 0.8))),
    dsdt_yn_bias = case(dsdt_yn, list(n_old = 50, n_trials = 100, stimulus = 1, d = 1, bias = 0.1),
                        list(bias = c(-0.2, 0.1, 0.3))),
    dsdt_mafc = case(dsdt_mafc, list(n_correct = 50, n_trials = 100, m = 2, d = 1),
                     list(n_correct = c(40, 50, 60), n_trials = c(100, 90, 80), m = c(2, 3, 4), d = c(0.5, 1, 1.5))),
    dsdt_rating = case(dsdt_rating, list(counts = c4[1, , drop = FALSE], stimulus = 1, d = 1, thresholds = thr,
                                         sdratio = 1),
                       rating),
    dsdt_dpsdt = case(dsdt_dpsdt, list(counts = c4[1, , drop = FALSE], stimulus = 0, d = 1, thresholds = thr,
                                       Ro = 0.2, Rn = 0.1, sdratio = 1),
                      c(rating, list(Ro = c(0.1, 0.2, 0.3), Rn = c(0.05, 0.1, 0.2)))),
    dsdt_metad = case(dsdt_metad, list(counts = c4[1, , drop = FALSE], stimulus = 1, d = 1, thresholds = thr,
                                       metad = 1, sdratio = 1),
                      c(rating, list(metad = c(0.8, 1, 1.2)))),
    dsdt_ranking = case(dsdt_ranking, list(counts = c4[1, , drop = FALSE], m = 4, d = 1, sdratio = 1),
                        list(counts = c4, m = c(4, 4, 4), d = c(0.5, 1, 1.5), sdratio = c(1, 1.3, 0.8)),
                        list(dist = "normal")),
    dsdt_cdp = case(dsdt_cdp, list(counts = c9[1, , drop = FALSE], stimulus = 1, dfam = 0.8, drec = 1,
                                   thresholds = thr5, rcrit = 0.5, sigmar = 0, rho = 0),
                    list(counts = c9, stimulus = c(0, 1, 1), dfam = c(0.6, 0.8, 1), drec = c(0.8, 1, 1.2),
                         thresholds = rbind(thr5, thr5 + 0.1, thr5 - 0.1), rcrit = c(0.4, 0.5, 0.6),
                         sigmar = c(0, 0.1, -0.1), rho = c(0, 0.2, -0.2)),
                    list(n_new = 3, log = TRUE)),
    dsdt_cdp_kcrit = case(dsdt_cdp, list(counts = rbind(c(10, 20, 30, 10, 15, 5, 5, 3, 2, 1, 2, 1)), stimulus = 1,
                                         dfam = 0.8, drec = 1, thresholds = thr5, rcrit = 0.5, kcrit = 0.3),
                          list(kcrit = c(0.2, 0.3, 0.4)), list(n_new = 3, log = TRUE))
  )
  # two expectations over all 69 cases rather than expect_elementwise() for
  # each, whose per-expectation overhead tripled the time of this test
  vec <- each <- list()
  expect_no_warning(for (nm in names(cases)) {
    cs <- cases[[nm]]
    for (a in names(cs$long)) {
      args <- replace(cs$base, a, cs$long[a])
      vec[[paste(nm, a)]] <- do.call(cs$f, c(args, cs$fixed))
      each[[paste(nm, a)]] <- vapply(1:3, function(i) do.call(cs$f, c(lapply(args, pick, i = i), cs$fixed)), 0)
    }
  })
  expect_identical(lengths(vec), setNames(rep(3L, 69), names(vec)))
  expect_equal(unlist(vec) / unlist(each), rep(1, 3 * 69), ignore_attr = TRUE)
})

test_that("a matrix in a one-number-per-observation argument gives one value per element", {
  # only counts, thresholds and the dezdm() 4par pairs are read by rows; every
  # other argument of every d*() function, in turn a 2 x 2 matrix, holds one
  # observation per element, in column-major order, as in dnorm() and as
  # before bmm 1.5.0. Read by rows, it would give 2 values, not 4.
  thr <- c(0.2, 0.7, 1.2)
  thr5 <- c(-1.1, -0.5, 0, 0.6, 1.3)
  rating_c <- c(2, 8, 20, 70)
  case <- function(f, base, mats, fixed = list()) list(f = f, base = base, mats = mats, fixed = fixed)
  stim <- c(1, 0, 1, 1)
  sdr <- c(1, 1.3, 0.8, 1.1)
  dee <- c(0.5, 1, 1.5, 0.8)
  cases <- list(
    dsdm = case(dsdm, list(x = 0.3, mu = 0, c = 3, kappa = 3),
                list(x = c(0.1, 0.3, -0.5, 0.8), mu = c(0, 0.5, -0.4, 1), c = c(2, 3, 4, 2.5),
                     kappa = c(1, 3, 8, 2))),
    dddm = case(dddm, list(rt = 0.6, response = 1, drift = 1, bound = 1, ndt = 0.2, zr = 0.5),
                list(rt = c(0.5, 0.6, 0.8, 0.7), response = stim, drift = dee, bound = c(1, 1.2, 1.5, 0.9),
                     ndt = c(0.2, 0.25, 0.3, 0.15), zr = c(0.4, 0.5, 0.6, 0.45))),
    dezdm3 = case(dezdm, list(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100, drift = 2,
                              bound = 1.5, ndt = 0.3, s = 1),
                  list(mean_rt = c(0.45, 0.5, 0.55, 0.48), var_rt = c(0.02, 0.025, 0.03, 0.022),
                       n_upper = c(70, 80, 75, 85), n_trials = c(100, 110, 120, 95), drift = c(1.5, 2, 2.5, 1.8),
                       bound = c(1.2, 1.5, 1.8, 1.4), ndt = c(0.25, 0.3, 0.35, 0.28), s = c(1, 1.1, 0.9, 1.05))),
    dezdm4 = case(dezdm, list(mean_rt = c(0.45, 0.55), var_rt = c(0.02, 0.03), n_upper = 80, n_trials = 100,
                              drift = 1, bound = 1.5, ndt = 0.3, zr = 0.5, s = 1),
                  list(n_upper = c(70, 80, 75, 85), n_trials = c(100, 110, 120, 95), drift = dee,
                       bound = c(1.2, 1.5, 1.8, 1.4), ndt = c(0.25, 0.3, 0.35, 0.28), zr = c(0.4, 0.5, 0.6, 0.45),
                       s = c(1, 1.1, 0.9, 1.05)),
                  list(version = "4par")),
    dsdt_yn = case(dsdt_yn, list(n_old = 50, n_trials = 100, stimulus = 1, d = 1, criterion = 0.5, sdratio = 1),
                   list(n_old = c(40, 50, 60, 45), n_trials = c(100, 90, 80, 95), stimulus = stim, d = dee,
                        criterion = c(0.2, 0.5, 0.8, 0.3), sdratio = sdr)),
    dsdt_yn_bias = case(dsdt_yn, list(n_old = 50, n_trials = 100, stimulus = 1, d = 1, bias = 0.1),
                        list(bias = c(-0.2, 0.1, 0.3, 0.2))),
    dsdt_mafc = case(dsdt_mafc, list(n_correct = 50, n_trials = 100, m = 3, d = 1),
                     list(n_correct = c(40, 50, 60, 45), n_trials = c(100, 90, 80, 95), m = c(2, 3, 4, 5), d = dee)),
    dsdt_rating = case(dsdt_rating, list(counts = rating_c, stimulus = 1, d = 1, thresholds = thr, sdratio = 1),
                       list(stimulus = stim, d = dee, sdratio = sdr), list(log = TRUE)),
    dsdt_dpsdt = case(dsdt_dpsdt, list(counts = rating_c, stimulus = 1, d = 1, thresholds = thr, Ro = 0.2,
                                       Rn = 0.1, sdratio = 1),
                      list(stimulus = stim, d = dee, sdratio = sdr, Ro = c(0.1, 0.2, 0.3, 0.15),
                           Rn = c(0.05, 0.1, 0.2, 0.15)),
                      list(log = TRUE)),
    dsdt_metad = case(dsdt_metad, list(counts = rating_c, stimulus = 1, d = 1, thresholds = thr, metad = 1,
                                       sdratio = 1),
                      list(stimulus = stim, d = dee, metad = c(0.8, 1, 1.2, 0.9), sdratio = sdr), list(log = TRUE)),
    dsdt_ranking = case(dsdt_ranking, list(counts = c(30, 40, 0, 0), m = 4, d = 1, sdratio = 1),
                        list(m = c(2, 3, 4, 3), d = dee, sdratio = sdr), list(dist = "normal", log = TRUE)),
    dsdt_cdp = case(dsdt_cdp, list(counts = c(10, 20, 30, 10, 15, 5, 5, 3, 2), stimulus = 1, dfam = 0.8,
                                   drec = 1, thresholds = thr5, rcrit = 0.5, sigmar = 0, rho = 0),
                    list(stimulus = stim, dfam = c(0.6, 0.8, 1, 0.7), drec = c(0.8, 1, 1.2, 0.9),
                         rcrit = c(0.4, 0.5, 0.6, 0.45), sigmar = c(0, 0.1, -0.1, 0.05), rho = c(0, 0.2, -0.2, 0.1)),
                    list(n_new = 3, log = TRUE)),
    dsdt_cdp_kcrit = case(dsdt_cdp, list(counts = c(10, 20, 30, 10, 15, 5, 5, 3, 2, 1, 2, 1), stimulus = 1,
                                         dfam = 0.8, drec = 1, thresholds = thr5, rcrit = 0.5, kcrit = 0.3),
                          list(kcrit = c(0.2, 0.3, 0.4, 0.25)), list(n_new = 3, log = TRUE))
  )
  vec <- each <- list()
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    for (a in names(cs$mats)) {
      value <- matrix(cs$mats[[a]], 2)
      vec[[paste(nm, a)]] <- as.vector(do.call(cs$f, c(replace(cs$base, a, list(value)), cs$fixed)))
      each[[paste(nm, a)]] <- vapply(value, function(v) do.call(cs$f, c(replace(cs$base, a, v), cs$fixed)), 0)
    }
  }
  expect_identical(lengths(vec), setNames(rep(4L, 58), names(vec)))
  expect(all(is.finite(unlist(each)) & unlist(each) != 0), "a reference value is zero or not finite")
  expect_equal(unlist(vec) / unlist(each), rep(1, 4 * 58), ignore_attr = TRUE)
})

test_that("dsdm() and dsdt_yn() keep the dim, names and class they returned before bmm 1.5.0", {
  # dsdm() from any argument as long as the result, x first, as dnorm();
  # dsdt_yn() from n_trials, or from n_old when it is longer, as lchoose() did
  x <- matrix(c(0.1, 0.3, -0.5, 0.8), 2, dimnames = list(c("a", "b"), NULL))
  out <- dsdm(x, mu = 0, c = 3, kappa = 3)
  expect_identical(dimnames(out), dimnames(x))
  expect_equal(as.vector(out), dsdm(c(0.1, 0.3, -0.5, 0.8), mu = 0, c = 3, kappa = 3))
  expect_named(dsdm(c(a = 0.1, b = 0.3), mu = 0, c = 3, kappa = c(1, 3)), c("a", "b"))
  expect_named(dsdm(c(x = 0.1, y = 0.2), mu = c(a = 0, b = 1)), c("x", "y"))
  expect_named(dsdm(c(0.1, 0.2), mu = c(a = 0, b = 1), c = 2, kappa = 3), c("a", "b"))
  expect_null(names(dsdm(c(0.1, 0.2), mu = c(a = 0), kappa = c(1, 2))))
  expect_named(dsdt_yn(c(a = 30, b = 10), 50, c(1, 0), 1, 0.5), c("a", "b"))
  expect_identical(dim(dsdt_yn(matrix(c(30, 10, 20, 5), 2), 50, 1, 1, 0.5)), c(2L, 2L))
  expect_identical(dim(dsdm(0.3, kappa = matrix(c(1, 3, 8, 2), 2))), c(2L, 2L))
  expect_null(dim(dsdt_yn(30, 50, 1, matrix(c(1, 2, 3, 4), 2), 0.5)))
  expect_named(dsdt_yn(c(30, 10), c(a = 50, b = 60), 1, 1, 0.5), c("a", "b"))
  expect_named(dsdt_yn(c(x = 30, y = 10), c(a = 50, b = 60), 1, 1, 0.5), c("a", "b"))
  expect_null(names(dsdt_yn(c(x = 30, y = 10), c(50, 60), 1, 1, 0.5)))
  expect_null(dim(dsdt_yn(matrix(c(7, 3), 2, 1), c(10, 11), c(1, 0), 1, 0.2)))
  expect_null(attributes(dsdt_yn(ts(c(7, 3)), c(10, 11), c(1, 0), 1, 0.2)))
  x_ts <- ts(c(0.1, 0.3, -0.5, 0.8), start = 2000)
  expect_identical(tsp(dsdm(x_ts, mu = 0, c = 3, kappa = c(1, 3))), c(2000, 2003, 1))
  expect_s3_class(dsdm(structure(c(0.1, 0.3), class = "foo"), mu = 0, c = 3, kappa = 3), "foo")
  expect_identical(tsp(dsdt_yn(ts(c(30, 10, 20, 5), start = 2000), 50, 1, 1, 0.5)), c(2000, 2003, 1))
})

test_that("dsdm() returns a data frame of densities for a data-frame x, as before bmm 1.5.0", {
  # reference values from bmm 1.4 (develop c13dd2c5), where the arithmetic on
  # the data frame paired its elements, column by column, with the parameters
  expect_equal(dsdm(data.frame(v = c(0.1, 0.2))), data.frame(v = c(0.6255801, 0.5592532)), tolerance = 1e-6)
  expect_equal(dsdm(data.frame(v = c(0.1, 0.2)), log = TRUE), data.frame(v = c(-0.4690759, -0.5811530)),
               tolerance = 1e-6)
  x <- data.frame(v = c(0.1, 0.2), w = c(0.3, -0.4), row.names = c("r1", "r2"))
  expect_equal(dsdm(x, kappa = c(2, 5, 7)),
               data.frame(v = c(0.4205239, 0.6845533), w = c(0.5156991, 0.3339982), row.names = c("r1", "r2")),
               tolerance = 1e-6)
  # the parameter longer than the cells is cut without a warning
  expect_no_warning(out <- dsdm(x, kappa = c(2, 5, 7, 1.5, 9)))
  expect_equal(out,
               data.frame(v = c(0.4205239, 0.6845533), w = c(0.5156991, 0.3009992), row.names = c("r1", "r2")),
               tolerance = 1e-6)
  expect_identical(dsdm(data.frame(v = numeric(0))), data.frame(v = numeric(0)))
  # an empty parameter empties the result, as an empty x does; bmm 1.4 stopped
  expect_identical(dsdm(data.frame(v = c(0.1, 0.2)), mu = numeric(0)), data.frame(v = numeric(0)))
  expect_identical(dsdm(data.frame(v = c(0.1, 0.2), w = 3:4), mu = NULL), data.frame(v = numeric(0), w = numeric(0)))
})

test_that("dsdm() gives NA with a warning for a factor x, as before bmm 1.5.0", {
  expect_warning(out <- dsdm(factor(1)), "not meaningful for factors")
  expect_identical(out, NA_real_)
  expect_warning(out <- dsdm(factor(c(1, 2)), kappa = c(a = 2, b = 5)), "not meaningful for factors")
  expect_identical(out, c(a = NA_real_, b = NA_real_))
})

test_that("dsdm() keeps the outcomes of bmm 1.4 for data-frame columns that are not numbers", {
  # reference values from bmm 1.4 (develop c13dd2c5): a factor column beside a
  # numeric one gives NA with a warning, not densities at its level codes
  x <- data.frame(a = c(0.1, 0.2), b = factor(c(0.3, 0.4)))
  expect_warning(out <- dsdm(x), "not meaningful for factors")
  expect_equal(out, data.frame(a = c(0.6255800913, 0.5592531962), b = c(NA_real_, NA_real_)), tolerance = 1e-9)
  # bmm 1.4 stopped on a date column, which must not run as day numbers
  expect_error(dsdm(data.frame(a = c(0.1, 0.2), b = as.Date(c("2020-01-01", "2020-01-03")))), "must be numeric")
  # a matrix column spreads into one column per cell, and a complex column runs
  x <- data.frame(id = 1:2)
  x$m <- matrix(c(0.1, 0.2, 0.3, 0.4), 2)
  expect_equal(dsdm(x),
               data.frame(id = c(0.1084691396508, 0.0703997091249),
                          m.1 = c(0.6255800912855, 0.5592531962204),
                          m.2 = c(0.4703425690427, 0.3787859868264)),
               tolerance = 1e-9)
  # a matrix column kept by I() keeps its shape, and the column beside it its own cells
  x <- data.frame(a = I(matrix(c(0.1, 0.2, 0.3, 0.4), 2)), b = c(-1, 1))
  expect_equal(dsdm(x),
               data.frame(a = I(matrix(c(0.6255800912855, 0.5592531962204, 0.4703425690427, 0.3787859868264), 2)),
                          b = c(0.1084691396508, 0.1084691396508)),
               tolerance = 1e-9)
  expect_equal(dsdm(data.frame(a = complex(real = c(0.1, 0.2), imaginary = 0))),
               data.frame(a = complex(real = c(0.6255800912855, 0.5592531962204), imaginary = 0)),
               tolerance = 1e-9)
})

test_that("dsdm() normalises each observation by its own parameter set when sets repeat", {
  # x longer than every parameter, so that sets repeat and are looked up by
  # key; references by quadrature of the sqrtexp numerator, not by dsdm().
  # The (mu, kappa) sets of the first period are 1 2 3 2 1 4: not periodic
  # in the order of their first appearance
  expect_equal(dsdm(seq(-0.7, 0.7, length.out = 8), mu = c(0, 1), kappa = c(1, 1, 2)),
               c(0.2191874724561, 0.1365186194367, 0.3700701487885, 0.1701868481098,
                 0.2799925229630, 0.2266224266039, 0.2453847761430, 0.2673206318095),
               tolerance = 1e-8)
  # kappa values that agree in their first four digits get their own integral;
  # normalised by kappa = 2, the second value would be 0.4003006905
  expect_equal(dsdm(c(0.1, 0.2, 0.3), mu = 0, kappa = c(2, 2.0001)),
               c(0.4205238558571, 0.4002974625841, 0.3700701487885), tolerance = 1e-8)
  # sets 1 2 3 2 1 4 again, but mu differs between positions 4 and 6, so an
  # integral bounded by the mu at a position, not at the set's first
  # appearance, gives a wrong value
  expect_equal(dsdm(seq(0.1, 0.6, 0.1), mu = c(0, 0, 1), kappa = c(1, 2)),
               c(0.2799925229630, 0.4002857123791, 0.2191874724561, 0.3339982008881,
                 0.2453847761430, 0.3339982008881),
               tolerance = 1e-8)
})

test_that("dsdt_yn() refuses n_old above n_trials, also after recycling", {
  expect_error(dsdt_yn(c(5, 15), c(20, 10), 1, d = 1, criterion = 0), "n_old must not exceed n_trials")
  expect_error(dsdt_yn(15, c(20, 10, 30), 1, d = 1, criterion = 0), "n_old must not exceed n_trials")
})

test_that("non-multiple lengths recycle without a warning in dsdm(), rsdt_yn() and rddm()", {
  withr::local_seed(316)
  expect_no_warning(dsdm(c(0.1, -0.5), mu = c(0, 1, 2), c = 3, kappa = 3))
  expect_no_warning(rsdt_yn(7, c(10, 20, 30), stimulus = c(1, 0), d = 1, criterion = 0.2))
  expect_no_warning(rddm(7, drift = 1, bound = 1, ndt = 0.2, zr = c(0.4, 0.6)))
})

test_that("counts given as a data frame read as the matrix of its rows", {
  counts <- rbind(c(50, 30, 15, 5), c(20, 30, 30, 20), c(40, 40, 10, 10))
  thr <- c(0.2, 0.7, 1.2)
  expect_equal(dsdt_ranking(as.data.frame(counts), m = 4, d = 1), dsdt_ranking(counts, m = 4, d = 1))
  expect_equal(dsdt_ranking(as.data.frame(counts)[1, ], m = 4, d = c(0.5, 1)),
               dsdt_ranking(counts[1, , drop = FALSE], m = 4, d = c(0.5, 1)))
  expect_equal(dsdt_rating(as.data.frame(counts), 1, d = 1, thresholds = thr),
               dsdt_rating(counts, 1, d = 1, thresholds = thr))
  # a count column read as text is refused, not converted
  text_counts <- data.frame(a = c(2, 30), b = c(8, 40), c = c(20, 20), d = c("70", "10"))
  expect_error(dsdt_rating(text_counts, 1, d = 1, thresholds = thr),
               "counts must be numbers; the data frame has non-numeric columns: 'd'", fixed = TRUE)
  text_counts$d <- factor(text_counts$d)
  expect_error(dsdt_ranking(text_counts, m = 4, d = 1), "counts must be numbers")
})

test_that("rsdt_yn() pairs sensitivity, criterion, variance ratio and stimulus per draw", {
  withr::local_seed(4445)
  d <- c(0.5, 1.5, 2.5)
  criterion <- c(0.2, 0.8, 0.5)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  sdratio <- c(1, 1.5, 0.6)
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
  # the same for bias, measured from the midpoint, which is d / 2 = 5 above the noise mean
  expect_no_warning(n_old <- rsdt_yn(7, 100, 1, d = 10, bias = c(-10, 20, -10)))
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
  # n not a multiple of the parameter lengths: draw i takes element i, and a
  # sensitivity of +-15 makes every trial correct or every trial wrong, so m
  # cannot change the draw and is left to the test above
  n_correct <- NULL
  expect_no_warning(n_correct <- rsdt_mafc(7, c(10, 20, 30, 40), m = c(2, 5), d = c(15, 15, -15)))
  expect_equal(n_correct, c(10, 20, 0, 40, 10, 0, 30))
})

# category proportions of a group against the reference probabilities, each
# within 5 standard errors; one expectation for all categories keeps the
# samplers' tests within the CRAN time per test
expect_category_props <- function(counts, n_trials, probs) {
  se <- se_prop(probs, n_trials * nrow(counts))
  expect_lt(max(abs(colMeans(counts) / n_trials - probs) - 5 * se), 0)
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
  # n not a multiple of the threshold rows or n_trials: row i takes threshold
  # row i and n_trials element i, and thresholds far from the mean put every
  # response in one category
  counts <- rsdt_rating(7, c(10, 20, 30, 40), 1, d = 0, thresholds = rbind(c(-9, -8, -7), c(7, 8, 9), c(-9, 8, 9)))
  expect_equal(max.col(counts), rep_len(c(4L, 1L, 2L), 7))
  expect_equal(counts[cbind(1:7, max.col(counts))], c(10, 20, 30, 40, 10, 20, 30))
})

test_that("rsdt_dpsdt() pairs recollection with the rating parameters per draw", {
  withr::local_seed(4445)
  d <- c(0.5, 1.5, 1)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  Ro <- c(0.2, 0.5, 0.3)
  Rn <- c(0, 0.2, 0.1)
  sdratio <- c(0.6, 1.5, 2)
  n_trials <- 100
  n <- 1500
  g <- group_of(n, 3)
  counts <- rsdt_dpsdt(n, n_trials, stimulus, d, thresholds = thr3, Ro = Ro, Rn = Rn, sdratio = sdratio)
  for (k in 1:3) {
    p <- sdt_category_probs(d[k], thr3[k, ], stimulus[k], sdratio[k])
    # recollection loads the most confident category of the stimulus
    probs <- if (stimulus[k] == 1) (1 - Ro[k]) * p + Ro[k] * unit_counts(4, 4) else (1 - Rn[k]) * p + Rn[k] * unit_counts(1, 4)
    expect_category_props(counts[g == k, ], n_trials, probs)
  }
  # n not a multiple of the parameter lengths: draw i takes element i;
  # thresholds far above the familiarity put a response in category 1, far
  # below in category 4, and recollection with probability 1 moves it to 4
  counts <- NULL
  expect_no_warning(counts <- rsdt_dpsdt(7, c(10, 20, 30, 40), 1, d = 1,
                                         thresholds = rbind(c(20, 21, 22), c(-22, -21, -20)), Ro = c(0, 1, 0)))
  expect_equal(max.col(counts), c(1L, 4L, 1L, 4L, 4L, 4L, 1L))
  expect_equal(counts[cbind(1:7, max.col(counts))], c(10, 20, 30, 40, 10, 20, 30))
  # the same for Rn, which acts on noise trials only: familiarity puts every
  # response in category 4, recollection of a noise item moves it to 1
  counts <- NULL
  expect_no_warning(counts <- rsdt_dpsdt(7, 10, stimulus = c(0, 0, 1, 0), d = 1,
                                         thresholds = c(-22, -21, -20), Ro = 0, Rn = c(1, 0, 0)))
  expect_equal(max.col(counts), c(1L, 4L, 4L, 1L, 4L, 4L, 4L))
})

test_that("rsdt_metad() pairs metacognitive sensitivity with the rating parameters per draw", {
  withr::local_seed(4445)
  d <- c(1, 1, 1.5)
  metad <- c(0.5, 1.5, 1)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  sdratio <- c(1.4, 1, 0.7)
  n_trials <- 100
  n <- 1500
  g <- group_of(n, 3)
  counts <- rsdt_metad(n, n_trials, stimulus, d, thresholds = thr3, metad = metad, sdratio = sdratio)
  for (k in 1:3) {
    probs <- vapply(1:4, function(j) {
      dsdt_metad(unit_counts(j, 4), stimulus[k], d[k], thr3[k, ], metad = metad[k], sdratio = sdratio[k])
    }, numeric(1))
    expect_category_props(counts[g == k, ], n_trials, probs)
  }
  # n not a multiple of the parameter lengths: draw i takes element i. With
  # d = 20 and thresholds (-10, 10, 60), metad = 20 leaves the thresholds in
  # place (signal in category 3, noise in 2); metad = 150 moves them up by
  # (150 - 20) / 2 = 65 in the meta-d' space (signal in 4, noise in 1)
  counts <- NULL
  expect_no_warning(counts <- rsdt_metad(7, c(30, 50, 40, 60), stimulus = c(1, 1, 0, 1), d = 20,
                                         thresholds = c(-10, 10, 60), metad = c(20, 150, 20)))
  expect_equal(max.col(counts), c(3L, 4L, 2L, 3L, 4L, 3L, 2L))
  expect_equal(counts[cbind(1:7, max.col(counts))], c(30, 50, 40, 60, 30, 50, 40))
})

test_that("rsdt_ranking() pairs sensitivity, set size and variance ratio per draw", {
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
  # the variance ratio acts under normal noise only; with two items the target
  # is ranked first with probability pnorm(d / sqrt(2)) whatever the ratio
  sdratio <- c(0.6, 1.5, 2)
  counts <- rsdt_ranking(n, n_trials, m, d, sdratio = sdratio, dist = "normal")
  p_first <- pnorm(d[1] / sqrt(2))
  expect_mean_near(counts[g == 1, 1] / n_trials, p_first, se_prop(p_first, n_trials * sum(g == 1)))
  for (k in 2:3) {
    probs <- vapply(seq_len(m[k]), function(j) {
      dsdt_ranking(unit_counts(j, m[k]), m[k], d[k], sdratio = sdratio[k], dist = "normal")
    }, numeric(1))
    expect_category_props(counts[g == k, seq_len(m[k])], n_trials, probs)
  }
  # n not a multiple of the parameter lengths: draw i takes element i, and a
  # sensitivity of +-20 ranks the target first or last (rank m)
  for (dist in c("gumbel_min", "normal")) {
    counts <- NULL
    expect_no_warning(counts <- rsdt_ranking(7, c(10, 20, 30, 40), m = c(2, 4), d = c(20, 20, -20), dist = dist))
    expect_equal(max.col(counts), c(1L, 1L, 2L, 1L, 1L, 4L, 1L))
    expect_equal(counts[cbind(1:7, max.col(counts))], c(10, 20, 30, 40, 10, 20, 30))
  }
})

test_that("rsdt_cdp() pairs familiarity, recollection, criteria, correlation and thresholds per draw", {
  withr::local_seed(4445)
  dfam <- c(0.8, 1.2, 0.5)
  drec <- c(1, 0.5, 1.5)
  rcrit <- c(0.5, 0.3, 0.7)
  kcrit <- c(0.3, -0.4, 0.8)
  sigmar <- c(0.4, 0, -0.3)
  rho <- c(0.5, 0, -0.4)
  # not a palindrome, so a reversed pairing is visible
  stimulus <- c(0, 1, 1)
  thr5 <- rbind(c(-1.1, -0.5, 0, 0.6, 1.3), c(-1, -0.4, 0.1, 0.7, 1.4), c(-1.2, -0.6, -0.1, 0.5, 1.2))
  n_trials <- 100
  n <- 600
  g <- group_of(n, 3)
  # remember/know, then remember/know/guess
  counts <- rsdt_cdp(n, n_trials, stimulus, dfam, drec, thresholds = thr5, rcrit = rcrit, sigmar = sigmar,
                     rho = rho, n_new = 3)
  expect_equal(dim(counts), c(n, 9))
  for (k in 1:3) {
    # one call scores the nine single-trial outcomes (rows of diag(9)) under group k
    probs <- dsdt_cdp(diag(9), stimulus[k], dfam[k], drec[k], thr5[k, ], rcrit = rcrit[k], sigmar = sigmar[k],
                      rho = rho[k], n_new = 3)
    expect_category_props(counts[g == k, ], n_trials, probs)
  }
  counts <- rsdt_cdp(n, n_trials, stimulus, dfam, drec, thresholds = thr5, rcrit = rcrit, kcrit = kcrit,
                     sigmar = sigmar, rho = rho, n_new = 3)
  expect_equal(dim(counts), c(n, 12))
  for (k in 1:3) {
    probs <- dsdt_cdp(diag(12), stimulus[k], dfam[k], drec[k], thr5[k, ], rcrit = rcrit[k], kcrit = kcrit[k],
                      sigmar = sigmar[k], rho = rho[k], n_new = 3)
    expect_category_props(counts[g == k, ], n_trials, probs)
  }
  # n not a multiple of the parameter lengths: draw i takes element i. Far low
  # thresholds under a strong familiarity make every response the most
  # confident "old": remember (column 12) at rcrit -20, else know (9) at
  # kcrit 20 and guess (6) at kcrit 40
  counts <- NULL
  expect_no_warning(counts <- rsdt_cdp(7, c(30, 50, 40), 1, dfam = 30, drec = 0,
                                       thresholds = c(-28, -27, -26, -25, -24), rcrit = c(50, -20, 50),
                                       kcrit = c(20, 40, 40, 20), n_new = 3))
  expect_equal(max.col(counts), c(9L, 12L, 6L, 9L, 12L, 6L, 6L))
  expect_equal(counts[cbind(1:7, max.col(counts))], c(30, 50, 40, 30, 50, 40, 30))
})

# --- empty observations ---------------------------------------------------------

test_that("the d*() functions return an empty result for no observations", {
  thr <- c(0.2, 0.7, 1.2)
  expect_identical(dsdm(numeric(0), mu = c(0, 1)), numeric(0))
  expect_identical(dddm(numeric(0), numeric(0), drift = 1:2, bound = 1, ndt = 0.2), numeric(0))
  expect_identical(dezdm(numeric(0), numeric(0), numeric(0), 100, drift = 1:2, bound = 1, ndt = 0.3), numeric(0))
  expect_identical(dezdm(matrix(0, 0, 2), matrix(0, 0, 2), numeric(0), 100, drift = 1, bound = 1, ndt = 0.3,
                         zr = 0.5, version = "4par"), numeric(0))
  expect_identical(dsdt_yn(integer(0), 10, 1, d = 1:2, criterion = 0), numeric(0))
  expect_identical(dsdt_mafc(integer(0), 10, 4, d = 1:2), numeric(0))
  expect_identical(dsdt_rating(matrix(0, 0, 4), 1, d = 1, thresholds = thr), numeric(0))
  expect_identical(dsdt_dpsdt(matrix(0, 0, 4), 1, d = 1, thresholds = thr, Ro = 0.3), numeric(0))
  expect_identical(dsdt_metad(matrix(0, 0, 4), 1, d = 1, thresholds = thr, metad = 1), numeric(0))
  expect_identical(dsdt_ranking(matrix(0, 0, 4), m = 4, d = 1:2), numeric(0))
  expect_identical(dsdt_cdp(matrix(0, 0, 9), 1, dfam = 0.8, drec = 1, thresholds = c(-1.1, -0.5, 0, 0.6, 1.3),
                            rcrit = 0.5, n_new = 3), numeric(0))
})

test_that("a NULL per-observation argument is empty, as a misspelt data column", {
  # dsdm() returned numeric(0) for it before bmm 1.5.0; counts, thresholds and
  # the 4par pairs are one observation as a vector and stop when empty
  thr <- c(0.2, 0.7, 1.2)
  case <- function(f, base, args) list(f = f, base = base, args = args)
  cases <- list(
    dsdm = case(dsdm, list(x = 0.3, mu = 0, c = 3, kappa = 3), c("x", "mu", "c", "kappa")),
    dddm = case(dddm, list(rt = 0.6, response = 1, drift = 1, bound = 1, ndt = 0.2, zr = 0.5),
                c("rt", "response", "drift", "bound", "ndt", "zr")),
    dezdm3 = case(dezdm, list(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100, drift = 1, bound = 1,
                              ndt = 0.3, s = 1),
                  c("mean_rt", "var_rt", "n_upper", "n_trials", "drift", "bound", "ndt", "s")),
    dezdm4 = case(dezdm, list(mean_rt = c(0.45, 0.55), var_rt = c(0.02, 0.03), n_upper = 80, n_trials = 100,
                              drift = 1, bound = 1, ndt = 0.3, zr = 0.5, s = 1, version = "4par"),
                  c("n_upper", "n_trials", "drift", "bound", "ndt", "zr", "s")),
    dsdt_yn = case(dsdt_yn, list(n_old = 5, n_trials = 10, stimulus = 1, d = 1, criterion = 0, sdratio = 1),
                   c("n_old", "n_trials", "stimulus", "d", "sdratio")),
    dsdt_mafc = case(dsdt_mafc, list(n_correct = 5, n_trials = 10, m = 2, d = 1),
                     c("n_correct", "n_trials", "m", "d")),
    dsdt_rating = case(dsdt_rating, list(counts = c(2, 8, 20, 70), stimulus = 1, d = 1, thresholds = thr,
                                         sdratio = 1),
                       c("stimulus", "d", "sdratio")),
    dsdt_dpsdt = case(dsdt_dpsdt, list(counts = c(2, 8, 20, 70), stimulus = 1, d = 1, thresholds = thr,
                                       Ro = 0.2, Rn = 0.1, sdratio = 1),
                      c("stimulus", "d", "Ro", "Rn", "sdratio")),
    dsdt_metad = case(dsdt_metad, list(counts = c(2, 8, 20, 70), stimulus = 1, d = 1, thresholds = thr,
                                       metad = 1, sdratio = 1),
                      c("stimulus", "d", "metad", "sdratio")),
    dsdt_ranking = case(dsdt_ranking, list(counts = c(30, 40, 20, 10), m = 4, d = 1, sdratio = 1),
                        c("m", "d", "sdratio")),
    dsdt_cdp = case(dsdt_cdp, list(counts = c(10, 20, 30, 10, 15, 5, 5, 3, 2), stimulus = 1, dfam = 0.8,
                                   drec = 1, thresholds = c(-1.1, -0.5, 0, 0.6, 1.3), rcrit = 0.5, sigmar = 0,
                                   rho = 0, n_new = 3),
                    c("stimulus", "dfam", "drec", "rcrit", "sigmar", "rho"))
  )
  out <- list()
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    for (a in cs$args) out[[paste(nm, a)]] <- do.call(cs$f, replace(cs$base, a, list(NULL)))
  }
  expect_identical(out, setNames(rep(list(numeric(0)), 55), names(out)))
})

test_that("NULL in an optional argument means it is not given", {
  thr5 <- c(-1.1, -0.5, 0, 0.6, 1.3)
  expect_identical(dsdt_yn(c(5, 6), 10, 1, d = 1, criterion = NULL, bias = 0.2),
                   dsdt_yn(c(5, 6), 10, 1, d = 1, bias = 0.2))
  expect_length(dsdt_yn(c(5, 6), 10, 1, d = 1, criterion = 0.2, bias = NULL), 2)
  counts <- rbind(c(10, 20, 30, 10, 15, 5, 5, 3, 2), c(5, 10, 20, 30, 10, 5, 5, 3, 2))
  expect_length(dsdt_cdp(counts, 1, dfam = 0.8, drec = 1, thresholds = thr5, rcrit = 0.5, kcrit = NULL,
                         n_new = 3), 2)
  expect_identical(dezdm(0.5, 0.02, c(80, 70), 100, drift = 1, bound = 1, ndt = 0.3, zr = NULL),
                   dezdm(0.5, 0.02, c(80, 70), 100, drift = 1, bound = 1, ndt = 0.3))
  # rsdt_yn() gave NA draws for a NULL d before bmm 1.5.0, as rnorm() does for an empty mean
  expect_warning(draws <- rsdt_yn(3, 10, 1, d = NULL, criterion = 0), "NAs produced")
  expect_identical(draws, rep(NA_integer_, 3))
  expect_warning(draws <- rsdt_yn(3, 10, 1, d = 1, criterion = 0, sdratio = NULL), "NAs produced")
  expect_identical(draws, rep(NA_integer_, 3))
})

test_that("rsdt_yn() gives NA draws with a warning for a factor d, as before bmm 1.5.0", {
  # not draws from the level codes 1, 2
  warnings <- character()
  draws <- withCallingHandlers(
    rsdt_yn(6, 100, 1, d = factor(c(0.5, 3)), criterion = 0),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(draws, rep(NA_integer_, 6))
  expect_match(warnings, "not meaningful for factors", all = FALSE)
  expect_error(rsdt_yn(3, 100, 1, d = c("1", "2"), criterion = 0), "non-numeric argument")
})

test_that("the rsdt_*() samplers return integer counts, as before bmm 1.5.0", {
  withr::local_seed(4447)
  thr <- c(-0.5, 0.3, 1)
  draws <- list(
    yn = rsdt_yn(2, 100, 1, d = 1, criterion = 0),
    mafc = rsdt_mafc(2, 100, m = 3, d = 1),
    rating = rsdt_rating(2, 100, 1, d = 1, thresholds = thr),
    dpsdt = rsdt_dpsdt(2, 100, 1, d = 1, thresholds = thr, Ro = 0.2),
    metad = rsdt_metad(2, 100, 1, d = 1, thresholds = thr, metad = 0.8),
    ranking = rsdt_ranking(2, 100, m = 3, d = 1),
    cdp = rsdt_cdp(2, 100, 1, 0.8, 1, thr, 0.7, n_new = 1),
    cdp_guess = rsdt_cdp(2, 100, 1, 0.8, 1, thr, 0.7, kcrit = 0.2, n_new = 1)
  )
  expect_identical(vapply(draws, typeof, character(1)), setNames(rep("integer", 8), names(draws)))
})

test_that("dezdm(version = '4par') refuses an empty mean_rt or var_rt vector", {
  # a 4par vector is one observation, as a vector of counts is; a log_lik()
  # that loses the boundary indicators passes empty vectors and must not
  # return an empty log-likelihood
  for (pair in list(list(numeric(0), numeric(0)), list(c(0.5, 0.6), numeric(0)),
                    list(numeric(0), rbind(c(0.02, 0.03))))) {
    expect_error(dezdm(pair[[1]], pair[[2]], 80, 100, drift = 1, bound = 1, ndt = 0.3, version = "4par"),
                 "empty vector.*zero rows and 2 columns")
  }
  expect_error(dezdm(c(0.5, 0.6, 0.7), c(0.02, 0.03), 80, 100, drift = 1, bound = 1, ndt = 0.3, version = "4par"),
               "must be length 2")
})

test_that("dezdm(version = '4par') reads the pair shapes that bmm 1.4 read", {
  # reference values from bmm 1.4 (develop c13dd2c5)
  expect_equal(dezdm(c(0.5, 0.6), matrix(c(0.02, 0.03), 2, 1), 60, 100, 1, 1, 0.2, version = "4par"),
               -41.8388675903824, tolerance = 1e-10)
  expect_equal(dezdm(matrix(c(0.5, 0.55, 0.6, 0.62), 2), data.frame(u = c(0.02, 0.025), l = c(0.03, 0.035)),
                     c(60, 55), 100, 1, 1, 0.2, version = "4par"),
               c(-41.8388675903824, -69.1887379676426), tolerance = 1e-10)
  # a data-frame mean_rt is read as the same matrix; bmm 1.4 stopped
  expect_equal(dezdm(data.frame(u = c(0.5, 0.55), l = c(0.6, 0.62)), matrix(c(0.02, 0.025, 0.03, 0.035), 2),
                     c(60, 55), 100, 1, 1, 0.2, version = "4par"),
               c(-41.8388675903824, -69.1887379676426), tolerance = 1e-10)
  # s != 1 scales the drift in the logit of the accuracy pair by s^2
  expect_equal(dezdm(c(0.5, 0.6), c(0.02, 0.03), 60, 100, 1, 1, 0.2, zr = 0.4, s = 1.3, version = "4par"),
               -244.004534245268, tolerance = 1e-10)
})

test_that("dezdm(version = '4par') refuses a mean_rt or var_rt that is no pair and no 2-column matrix", {
  # bmm 1.4 refused these too; read as pairs they would be reshaped into rows
  refused <- list(
    list(matrix(c(0.5, 0.6, 0.55, 0.62, 0.52, 0.61), 2, 3), rbind(c(0.03, 0.05), c(0.04, 0.06))),
    list(c(0.5, 0.6), matrix(c(0.02, 0.03, 0.025, 0.035), 1, 4)),
    list(rbind(c(0.5, 0.6), c(0.55, 0.62)), matrix(c(0.02, 0.03, 0.025, 0.035), 4, 1)),
    list(c(0.5, 0.6), data.frame(v = c(0.02, 0.03, 0.025)))
  )
  for (pair in refused) {
    expect_error(dezdm(pair[[1]], pair[[2]], 20, 40, drift = 0.8, bound = 1.4, ndt = 0.25, version = "4par"),
                 "must be length 2 (one observation) or matrices with 2 columns", fixed = TRUE)
  }
})
