nt_model <- function(...) {
  mixture3p("dev_rad", nt_features = "col_nt", set_size = "set_size", regex = TRUE, ...)
}

test_that("mixture3p() exposes the parameters each version implies", {
  expect_equal(names(nt_model()$parameters), c("mu", "kappa", "thetat", "thetant"))
  expect_equal(
    names(nt_model(version = "slot")$parameters), c("mu", "kappa", "K", "pnt")
  )
  expect_equal(
    names(nt_model(variable_precision = TRUE)$parameters),
    c("mu", "kappa", "tau", "thetat", "thetant")
  )
  expect_equal(
    class(nt_model(version = "slot_averaging")),
    c("bmmodel", "circular", "non_targets", "mixture3p", "mixture3p_slot_averaging")
  )
  expect_error(nt_model(version = "slots"), "should be one of")
  expect_error(nt_model(vp_nodes = 21), "at least 41")
})

# The previous implementation was a brms mixture whose weights were a softmax
# over the target, the active non-targets and guessing, with a -100 sentinel
# switching off the components a trial's set size did not have. That is what
# .mixture3p_softmax_weights() has to reproduce.
test_that("the simple version reproduces the softmax weights it replaces", {
  previous <- function(thetat, thetant, set_size, max_set_size) {
    theta <- c(thetat, rep(-100, max_set_size - 1), 0)
    active <- seq_len(set_size - 1)
    theta[1 + active] <- thetant + log(1 / (set_size - 1))
    exp(theta - matrixStats::logSumExp(theta))
  }
  for (set_size in c(1, 2, 4, 8)) {
    weights <- .mixture3p_softmax_weights(0.7, -0.4, set_size)
    expect_equal(
      c(exp(weights$logw), exp(weights$logw_guess)),
      previous(0.7, -0.4, set_size, 8)[c(seq_len(set_size), 9)],
      tolerance = 1e-12
    )
  }
})

test_that("dmixture3p() reproduces the mixture it replaces", {
  x <- seq(-pi, pi, length.out = 40)
  mu <- c(0, 2, -1.5)
  previous <- function(x, mu, kappa, p_mem, p_nt) {
    components <- vapply(mu, function(m) {
      brms::dvon_mises(x, mu = m, kappa = kappa, log = TRUE)
    }, numeric(length(x)))
    probs <- c(p_mem, rep(p_nt / (length(mu) - 1), length(mu) - 1))
    matrixStats::rowLogSumExps(cbind(
      sweep(components, 2, log(probs), `+`),
      log(1 - p_mem - p_nt) + stats::dunif(x, -pi, pi, log = TRUE)
    ))
  }
  for (kappa in c(1, 6, 30)) {
    expect_equal(
      dmixture3p(x, mu = mu, kappa = kappa, p_mem = 0.6, p_nt = 0.2, log = TRUE),
      previous(x, mu, kappa, 0.6, 0.2),
      tolerance = 1e-12
    )
  }
})

test_that("dmixture3p() is a normalised density", {
  total <- stats::integrate(
    function(x) dmixture3p(x, mu = c(0, 2, -1.5), kappa = 6, p_mem = 0.6, p_nt = 0.2),
    -pi, pi
  )
  expect_equal(total$value, 1, tolerance = 1e-6)
  varying <- stats::integrate(
    function(x) {
      dmixture3p(x, mu = c(0, 2, -1.5), kappa = 6, p_mem = 0.6, p_nt = 0.2, tau = 1.5)
    },
    -pi, pi
  )
  expect_equal(varying$value, 1, tolerance = 1e-6)
})

test_that("rmixture3p() draws from dmixture3p()", {
  compare <- function(n = 2e5, n_bins = 20, ...) {
    draws <- rmixture3p(n, ...)
    breaks <- seq(-pi, pi, length.out = n_bins + 1)
    observed <- as.numeric(table(cut(draws, breaks))) / n
    expected <- vapply(seq_len(n_bins), function(i) {
      stats::integrate(function(x) dmixture3p(x, ...), breaks[i], breaks[i + 1])$value
    }, numeric(1))
    expect_lt(max(abs(observed - expected)), 0.005)
  }
  compare(mu = c(0, 2, -1.5), kappa = 6, p_mem = 0.6, p_nt = 0.2)
  compare(mu = c(0, 2, -1.5), kappa = 6, p_mem = 0.6, p_nt = 0.2, tau = 1.5)
})

test_that("the capacity versions separate storage from swapping", {
  x <- seq(-pi, pi, length.out = 30)
  nt <- c(2, -1.5, 1, 0.5, -2, 2.5, -0.7)
  # below capacity nothing is guessed, so the density is the target and the
  # non-targets only, in the proportions pnt implies
  expect_equal(
    .dmixture3p_slot(x, 0, 6, K = 6, p_nt = 0.25, set_size = 3, nt = nt, tau = 0, nodes = 41L),
    .dmixture3p_slot(x, 0, 6, K = 99, p_nt = 0.25, set_size = 3, nt = nt, tau = 0, nodes = 41L),
    tolerance = 1e-12
  )
  # with pnt = 0 and full capacity it collapses to a plain von Mises
  expect_equal(
    .dmixture3p_slot(x, 0, 6, K = 6, p_nt = 0, set_size = 3, nt = nt, tau = 0, nodes = 41L),
    brms::dvon_mises(x, mu = 0, kappa = 6, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("a set size of one leaves no room for a swap", {
  x <- seq(-pi, pi, length.out = 30)
  nt <- c(2, -1.5)
  # whatever the swap parameter says, a single item cannot be confused with
  # another, so the density matches the two-parameter model
  for (p_nt in c(0, 0.3, 0.9)) {
    expect_equal(
      .dmixture3p_slot(x, 0, 6, K = 0.5, p_nt = p_nt, set_size = 1, nt = nt, tau = 0, nodes = 41L),
      dmixture2p(x, mu = 0, kappa = 6, p_mem = 0.5, log = TRUE),
      tolerance = 1e-12
    )
  }
})

test_that("configure_model() drops the set-size sentinel machinery", {
  code <- stancode(
    bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1), oberauer_lin_2017, nt_model()
  )
  expect_match(code, "mixture3p_simple_lpdf")
  expect_false(grepl("LureIdx", code, fixed = TRUE))
  expect_false(grepl("-100", code, fixed = TRUE))
  expect_false(grepl("inv_ss", code, fixed = TRUE))
  expect_match(code, "array[N] int vint1", fixed = TRUE)
  expect_match(code, "array[N] real vreal7", fixed = TRUE)
})

test_that("an intercept is still refused when set_size is a predictor", {
  expect_error(
    bmm(bmf(kappa ~ 1 + set_size, thetat ~ 1, thetant ~ 1), oberauer_lin_2017,
      nt_model(),
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "contains \\s*an intercept"
  )
})

test_that("the swap parameter is flagged as prior-only at set size 1", {
  expect_warning(
    bmm(bmf(kappa ~ 1, thetat ~ 1, thetant ~ 0 + set_size), oberauer_lin_2017,
      nt_model(),
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "set-size-1 level"
  )
})

test_that("every version runs through the bmm() pipeline", {
  dat <- oberauer_lin_2017
  expect_silent(bmm(
    bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1), dat, nt_model(),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  expect_silent(bmm(
    bmf(kappa ~ 1, tau ~ 1, thetat ~ 1, thetant ~ 1), dat,
    nt_model(variable_precision = TRUE),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  expect_silent(bmm(
    bmf(kappa ~ 1, K ~ 1, pnt ~ 1), dat, nt_model(version = "slot_averaging"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
})

test_that("the Stan and R likelihoods agree for every version", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the mixture3p likelihood"
  )

  tab <- .circmix_kappa_table()
  code <- paste0(
    "functions {\n",
    read_lines2(system.file("stan_chunks", "circmix_funs.stan", package = "bmm")),
    "\n",
    read_lines2(system.file("stan_chunks", "mixture3p_funs.stan", package = "bmm")),
    "\n}\n",
    "data {
      int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ;
      int N; int NT; vector[N] y; vector[N] kappa; vector[N] tau;
      vector[N] thetat; vector[N] thetant; vector[N] K; vector[N] pnt;
      array[N] int ss; matrix[N, NT] nt;
    }
    generated quantities {
      vector[N] out_simple; vector[N] out_slot; vector[N] out_sa;
      for (i in 1:N) {
        out_simple[i] = mixture3p_simple_core(y[i], 0.0, kappa[i], tau[i],
                          thetat[i], thetant[i], ss[i], to_vector(nt[i]), 41,
                          logk, dlogk, logJ_min, dlogJ);
        out_slot[i] = mixture3p_slot_core(y[i], 0.0, kappa[i], tau[i], K[i],
                        pnt[i], ss[i], to_vector(nt[i]), 41, logk, dlogk,
                        logJ_min, dlogJ);
        out_sa[i] = mixture3p_slot_averaging_core(y[i], 0.0, kappa[i], tau[i],
                      K[i], pnt[i], ss[i], to_vector(nt[i]), 41, logk, dlogk,
                      logJ_min, dlogJ);
      }
    }"
  )

  y <- c(-2.5, -0.8, 0, 0.4, 1.9, 3.0)
  kappa <- c(0.5, 2, 6, 6, 20, 80)
  tau <- c(0, 0, 0, 1, 0.5, 2)
  thetat <- c(-1, 0, 0.8, 1.5, 2, 0.3)
  thetant <- c(-2, -1, -0.5, 0, 0.5, -1.5)
  K <- c(1.2, 2.5, 3, 3.5, 6, 8)
  pnt <- c(0.05, 0.15, 0.25, 0.4, 0.1, 0.5)
  ss <- c(1L, 2L, 3L, 4L, 6L, 8L)
  nt <- matrix(seq(-3, 3, length.out = 6 * 7), nrow = 6)

  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = list(
      NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
      logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, N = length(y), NT = ncol(nt),
      y = y, kappa = kappa, tau = tau, thetat = thetat, thetant = thetant,
      K = K, pnt = pnt, ss = ss, nt = nt
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  stan <- function(variable) as.numeric(fit$draws(variable))

  # the R side takes one observation at a time, because the non-target
  # locations differ per row
  expect_equal(
    stan("out_simple"),
    vapply(seq_along(y), function(i) {
      .dmixture3p_simple(y[i], 0, kappa[i], thetat[i], thetant[i], ss[i],
        nt[i, ], tau[i], 41L)
    }, numeric(1)),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_slot"),
    vapply(seq_along(y), function(i) {
      .dmixture3p_slot(y[i], 0, kappa[i], K[i], pnt[i], ss[i], nt[i, ], tau[i], 41L)
    }, numeric(1)),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_sa"),
    vapply(seq_along(y), function(i) {
      .dmixture3p_slot_averaging(y[i], 0, kappa[i], K[i], pnt[i], ss[i],
        nt[i, ], tau[i], 41L)
    }, numeric(1)),
    tolerance = 1e-12
  )
})
