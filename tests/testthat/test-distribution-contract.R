# The recycling contract of the r*() and d*() functions (#444, ?bmm_distributions):
# r*(n, ...) returns n draws, draw i from the i-th element (row) of every
# recycled argument; d*(x, ...) returns one value per element of the longest
# argument, value i from the i-th elements alone. Every r*() pairing test
# interleaves parameter groups and compares each group with a reference
# outside the function's own code (a closed form, quadrature, or the single-
# trial density of the matching d*()). The mutants in local/reviews/444a/
# (vector reversed, recycling removed, one group shifted) record which of them
# each test catches; a Monte Carlo test is not expected to see a 10 % shift.
# Non-multiple lengths recycle silently, as rnorm() does; the cswald tests use
# multiples only until #406 lands. The circular mixtures are circmix's (#434).

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
  skip_on_cran()
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
  expect_no_warning(expect_length(rsdm(5, mu = c(0, 1), c = 3, kappa = 3), 5))
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
  skip_on_cran()
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
  expect_equal(nrow(rddm(5, drift = c(1, 2), bound = 1, ndt = 0.2)), 5)
  expect_equal(nrow(rddm(2, drift = c(1, 2, 3), bound = 1, ndt = 0.2)), 2)
})

test_that("rddm(1, <vectors>) keeps one draw per element until 1.6.0, with a warning", {
  withr::local_seed(4442)
  out <- NULL
  expect_warning(out <- rddm(1, drift = c(1, 2, 3), bound = 1, ndt = 0.2), "deprecated since bmm 1.5.0")
  expect_equal(nrow(out), 3)
  expect_no_warning(expect_equal(nrow(rddm(1, drift = 1, bound = 1, ndt = 0.2)), 1))
})
