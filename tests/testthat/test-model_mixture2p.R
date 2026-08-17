test_that("mixture2p() validates its arguments", {
  expect_error(mixture2p(), "arguments are missing in mixture2p\\(\\): resp_error")
  expect_error(mixture2p("y", version = "slots"), "should be one of")
  expect_error(mixture2p("y", version = "slot"), "set_size argument is required")
  expect_error(
    mixture2p("y", variable_precision = NA),
    "must be either TRUE or FALSE"
  )
  expect_error(mixture2p("y", vp_nodes = 21), "odd number of at least 41")
  expect_error(mixture2p("y", vp_nodes = 42), "odd number of at least 41")
})

test_that("mixture2p() exposes the parameters each version implies", {
  expect_equal(names(mixture2p("y")$parameters), c("mu", "kappa", "thetat"))
  expect_equal(
    names(mixture2p("y", set_size = "ss", version = "slot")$parameters),
    c("mu", "kappa", "K")
  )
  expect_equal(
    names(mixture2p("y", variable_precision = TRUE)$parameters),
    c("mu", "kappa", "tau", "thetat")
  )
  expect_equal(mixture2p("y", variable_precision = TRUE)$links$tau, "log")
  expect_equal(
    class(mixture2p("y", set_size = "ss", version = "slot_averaging")),
    c("bmmodel", "circular", "mixture2p", "mixture2p_slot_averaging")
  )
})

# The likelihood before this model moved to a custom family was a two-component
# brms mixture, which is exactly p_mem * von Mises + (1 - p_mem) * uniform. This
# is the check that the default version still fits the same model.
test_that("dmixture2p() reproduces the two-component mixture it replaces", {
  x <- seq(-pi, pi, length.out = 50)
  previous <- function(x, mu, kappa, p_mem) {
    matrixStats::rowLogSumExps(cbind(
      log(p_mem) + brms::dvon_mises(x, mu = mu, kappa = kappa, log = TRUE),
      log(1 - p_mem) + stats::dunif(x, -pi, pi, log = TRUE)
    ))
  }
  for (kappa in c(0.5, 4, 30)) {
    for (p_mem in c(0.2, 0.8, 1)) {
      expect_equal(
        dmixture2p(x, mu = 0, kappa = kappa, p_mem = p_mem, log = TRUE),
        previous(x, 0, kappa, p_mem),
        tolerance = 1e-12
      )
    }
  }
})

test_that("dmixture2p() is a normalised density for every version", {
  integrates_to_one <- function(...) {
    total <- stats::integrate(function(x) dmixture2p(x, ...), -pi, pi)
    expect_equal(total$value, 1, tolerance = 1e-6)
  }
  integrates_to_one(kappa = 6, p_mem = 0.7)
  integrates_to_one(kappa = 6, p_mem = 0.7, tau = 1.5)
  integrates_to_one(kappa = 6, K = 3, set_size = 5, version = "slot")
  integrates_to_one(kappa = 6, K = 3, set_size = 2, version = "slot")
  integrates_to_one(kappa = 6, K = 3, set_size = 5, version = "slot_averaging")
  integrates_to_one(kappa = 6, K = 7, set_size = 2, version = "slot_averaging")
  integrates_to_one(kappa = 6, K = 7, set_size = 2, version = "slot_averaging", tau = 1)
})

test_that("the slot version remembers everything below capacity", {
  x <- seq(-pi, pi, length.out = 30)
  # p_mem is exactly 1 at set sizes below K, the boundary a softmax over mixture
  # weights cannot represent and the reason this version needs its own likelihood
  expect_equal(
    dmixture2p(x, kappa = 6, K = 4, set_size = 2, version = "slot", log = TRUE),
    brms::dvon_mises(x, mu = 0, kappa = 6, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    dmixture2p(x, kappa = 6, K = 4, set_size = 8, version = "slot", log = TRUE),
    dmixture2p(x, kappa = 6, p_mem = 0.5, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("slot averaging puts every item on the same number of slots when K divides evenly", {
  x <- seq(-pi, pi, length.out = 30)
  two_slots <- .circmix_kappa(2 * .circmix_J(6))
  expect_equal(
    dmixture2p(x, kappa = 6, K = 8, set_size = 4, version = "slot_averaging", log = TRUE),
    brms::dvon_mises(x, mu = 0, kappa = two_slots, log = TRUE),
    tolerance = 1e-8
  )
})

test_that("slot averaging guesses whenever fewer slots than items are available", {
  x <- seq(-pi, pi, length.out = 30)
  # K = 2, set_size = 4: half the items get one slot, half get none
  expect_equal(
    dmixture2p(x, kappa = 6, K = 2, set_size = 4, version = "slot_averaging", log = TRUE),
    dmixture2p(x, kappa = 6, p_mem = 0.5, log = TRUE),
    tolerance = 1e-8
  )
})

test_that("rmixture2p() draws from dmixture2p() for every version", {
  compare <- function(n = 2e5, n_bins = 20, ...) {
    draws <- rmixture2p(n, ...)
    breaks <- seq(-pi, pi, length.out = n_bins + 1)
    observed <- as.numeric(table(cut(draws, breaks))) / n
    expected <- vapply(seq_len(n_bins), function(i) {
      stats::integrate(function(x) dmixture2p(x, ...), breaks[i], breaks[i + 1])$value
    }, numeric(1))
    expect_lt(max(abs(observed - expected)), 0.005)
  }
  compare(kappa = 6, p_mem = 0.7)
  compare(kappa = 6, p_mem = 0.7, tau = 1.5)
  compare(kappa = 6, K = 3, set_size = 5, version = "slot")
  compare(kappa = 6, K = 7, set_size = 2, version = "slot_averaging")
  compare(kappa = 6, K = 7, set_size = 2, version = "slot_averaging", tau = 1)
})

test_that("configure_model() builds a custom family without the mixture machinery", {
  code <- stancode(bmf(kappa ~ 1, thetat ~ 1), oberauer_lin_2017, mixture2p("dev_rad"))
  expect_match(code, "mixture2p_simple_lpdf")
  expect_false(grepl("LureIdx", code, fixed = TRUE))
  expect_false(grepl("-100", code, fixed = TRUE))
  expect_false(grepl("log_sum_exp_theta", code, fixed = TRUE))
})

test_that("the set-size versions pass the set size as an integer covariate", {
  model <- mixture2p("dev_rad", set_size = "set_size", version = "slot")
  code <- stancode(bmf(kappa ~ 1, K ~ 1), oberauer_lin_2017, model)
  expect_match(code, "array[N] int vint1", fixed = TRUE)
  expect_match(code, "mixture2p_slot_lpdf(Y[n] | mu[n], kappa[n], K[n], vint1[n]", fixed = TRUE)
})

test_that("check_data() warns when the set size cannot identify the capacity", {
  model <- mixture2p("dev_rad", set_size = "set_size", version = "slot")
  constant_set_size <- subset(oberauer_lin_2017, set_size == 4)
  expect_warning(
    check_data(model, constant_set_size, bmf(kappa ~ 1, K ~ 1)),
    "not identified by the data"
  )
})

test_that("mu1 keeps working as a deprecated name for mu", {
  expect_warning(
    check_formula(
      mixture2p("dev_rad"), oberauer_lin_2017, bmf(mu1 = 0, kappa ~ 1, thetat ~ 1)
    ),
    "mu1.*renamed to.*mu"
  )
})

test_that("every version runs through the bmm() pipeline", {
  dat <- oberauer_lin_2017
  expect_silent(bmm(
    bmf(kappa ~ 1, thetat ~ 1), dat, mixture2p("dev_rad"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  expect_silent(bmm(
    bmf(kappa ~ 1, tau ~ 1, thetat ~ 1), dat,
    mixture2p("dev_rad", variable_precision = TRUE),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  expect_silent(bmm(
    bmf(kappa ~ 1, K ~ 1), dat,
    mixture2p("dev_rad", set_size = "set_size", version = "slot_averaging"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
})

test_that("default priors reach the parameters of each version", {
  priors <- default_prior(
    bmf(kappa ~ 1, thetat ~ 1), oberauer_lin_2017, mixture2p("dev_rad")
  )
  intercepts <- priors[priors$class == "Intercept", ]
  expect_equal(intercepts$prior[intercepts$dpar == "kappa"], "normal(2, 1)")
  expect_equal(intercepts$prior[intercepts$dpar == "thetat"], "logistic(0, 1)")
  # brms writes the response parameter's own priors with an empty dpar
  expect_equal(intercepts$prior[intercepts$dpar == ""], "constant(0)")
  # the parameters are distributional now, so their intercept priors carry a
  # real Intercept class instead of being stored as b/coef = "Intercept", which
  # is what let a user prior on class "b" silently fail to replace them
  expect_false(any(priors$class == "b" & priors$coef == "Intercept"))
})

test_that("the Stan and R likelihoods agree for every version", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the mixture2p likelihood"
  )

  tab <- .circmix_kappa_table()
  code <- paste0(
    "functions {\n",
    read_lines2(system.file("stan_chunks", "circmix_funs.stan", package = "bmm")),
    "\n",
    read_lines2(system.file("stan_chunks", "mixture2p_funs.stan", package = "bmm")),
    "\n}\n",
    "data {
      int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ;
      int N; vector[N] y; vector[N] kappa; vector[N] tau; vector[N] thetat;
      vector[N] K; array[N] int ss;
    }
    generated quantities {
      vector[N] out_simple; vector[N] out_slot; vector[N] out_sa;
      for (i in 1:N) {
        out_simple[i] = mixture2p_simple_core(y[i], 0.0, kappa[i], tau[i],
                          thetat[i], 41, logk, dlogk, logJ_min, dlogJ);
        out_slot[i] = mixture2p_slot_core(y[i], 0.0, kappa[i], tau[i], K[i],
                        ss[i], 41, logk, dlogk, logJ_min, dlogJ);
        out_sa[i] = mixture2p_slot_averaging_core(y[i], 0.0, kappa[i], tau[i],
                      K[i], ss[i], 41, logk, dlogk, logJ_min, dlogJ);
      }
    }"
  )

  y <- c(-2.5, -0.8, 0, 0.4, 1.9, 3.0)
  kappa <- c(0.5, 2, 6, 6, 20, 80)
  tau <- c(0, 0, 0, 1, 0.5, 2)
  thetat <- c(0.1, 0.5, 0.75, 0.9, 0.99, 0.6)
  K <- c(1.2, 2.5, 3, 3.5, 6, 8)
  ss <- c(1L, 2L, 3L, 4L, 6L, 8L)

  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = list(
      NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
      logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, N = length(y), y = y,
      kappa = kappa, tau = tau, thetat = thetat, K = K, ss = ss
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  stan <- function(variable) as.numeric(fit$draws(variable))

  expect_equal(
    stan("out_simple"),
    .dmixture2p_simple(y, 0, kappa, thetat, tau, 41L),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_slot"),
    .dmixture2p_slot(y, 0, kappa, K, ss, tau, 41L),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_sa"),
    .dmixture2p_slot_averaging(y, 0, kappa, K, ss, tau, 41L),
    tolerance = 1e-12
  )
})
