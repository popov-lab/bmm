n_rows <- function(y) length(y)

ld_mixture <- function(y, mu_nt = numeric(0), kappa = 5, w_mem = 1, w_guess = 0) {
  n <- n_rows(y)
  cosd <- cbind(cos(y), vapply(mu_nt, function(m) cos(y - m), numeric(n)))
  logw <- matrix(log(w_mem), nrow = n, ncol = length(w_mem), byrow = TRUE)
  .circmix_ld(cosd, logw, rep(log(w_guess), n), rep(kappa, n))
}

test_that(".circmix_kappa() inverts .circmix_J() over the tabulated range", {
  kappa <- exp(seq(log(0.02), log(4e4), length.out = 400))
  expect_equal(.circmix_kappa(.circmix_J(kappa)), kappa, tolerance = 1e-9)
})

test_that(".circmix_kappa() uses the exact asymptotics outside the table", {
  expect_equal(.circmix_kappa(1e-9), sqrt(2e-9), tolerance = 1e-12)
  expect_equal(.circmix_kappa(1e7), 1e7 + 0.5, tolerance = 1e-12)
})

test_that(".circmix_J() reproduces the limits of kappa * I1/I0", {
  expect_equal(.circmix_J(1e-4), 1e-8 / 2, tolerance = 1e-6)
  expect_equal(.circmix_J(1e4), 1e4 - 0.5, tolerance = 1e-6)
  expect_equal(.circmix_J(0), 0)
})

test_that(".circmix_ld() is a normalised density on the circle", {
  total <- stats::integrate(
    function(y) exp(ld_mixture(y, mu_nt = 1.2, kappa = 8, w_mem = c(0.5, 0.3), w_guess = 0.2)),
    -pi, pi
  )
  expect_equal(total$value, 1, tolerance = 1e-6)
})

test_that(".circmix_ld() reduces to a von Mises with one component and no guessing", {
  y <- seq(-pi, pi, length.out = 25)
  expect_equal(
    ld_mixture(y, kappa = 12, w_mem = 1, w_guess = 0),
    brms::dvon_mises(y, mu = 0, kappa = 12, log = TRUE),
    tolerance = 1e-10
  )
})

test_that(".circmix_het_ld() agrees with .circmix_ld() when the components share kappa", {
  y <- seq(-pi, pi, length.out = 25)
  cosd <- cbind(cos(y), cos(y - 1.2))
  logw <- matrix(log(c(0.5, 0.3)), nrow = length(y), ncol = 2, byrow = TRUE)
  guess <- rep(log(0.2), length(y))
  expect_equal(
    .circmix_het_ld(cosd, logw, guess, matrix(9, nrow = length(y), ncol = 2)),
    .circmix_ld(cosd, logw, guess, rep(9, length(y)))
  )
})

test_that(".circmix_vp_ld() returns the constant-precision density when tau is zero", {
  y <- seq(-pi, pi, length.out = 25)
  cosd <- matrix(cos(y))
  logw <- matrix(log(0.8), nrow = length(y))
  guess <- rep(log(0.2), length(y))
  expect_equal(
    .circmix_vp_ld(cosd, logw, guess, rep(10, length(y)), rep(0, length(y))),
    .circmix_ld(cosd, logw, guess, rep(10, length(y)))
  )
})

test_that(".circmix_vp_ld() is a normalised density on the circle", {
  total <- stats::integrate(
    function(y) {
      n <- length(y)
      exp(.circmix_vp_ld(
        matrix(cos(y)), matrix(log(0.8), nrow = n), rep(log(0.2), n),
        rep(10, n), rep(2, n)
      ))
    },
    -pi, pi
  )
  expect_equal(total$value, 1, tolerance = 1e-6)
})

test_that(".circmix_vp_ld() matches adaptive integration above the shape threshold", {
  vp_reference <- function(y, kappa, tau, p_mem) {
    shape <- .circmix_J(kappa) / tau
    integrand <- function(t) {
      J <- exp(t)
      k <- .circmix_kappa(J)
      density <- (p_mem * exp(k * cos(y) - .circmix_log_besselI0(k)) + (1 - p_mem)) / (2 * pi)
      density * stats::dgamma(J, shape = shape, scale = tau) * J
    }
    centre <- digamma(shape) + log(tau)
    width <- 30 * sqrt(trigamma(shape))
    log(stats::integrate(integrand, centre - width, centre + width,
      rel.tol = 1e-12, subdivisions = 8000L
    )$value)
  }

  for (shape in c(4, 2, 1)) {
    tau <- .circmix_J(10) / shape
    expect_equal(
      .circmix_vp_ld(matrix(cos(0.3)), matrix(log(0.8)), log(0.2), 10, tau),
      vp_reference(0.3, 10, tau, 0.8),
      tolerance = 1e-4
    )
  }
})

test_that(".circmix_vp_ld() refuses shapes the node count cannot integrate", {
  tau <- .circmix_J(10) / 0.5
  expect_error(
    .circmix_vp_ld(matrix(cos(0.3)), matrix(log(0.8)), log(0.2), 10, tau),
    "gamma shape"
  )
  expect_silent(
    .circmix_vp_ld(matrix(cos(0.3)), matrix(log(0.8)), log(0.2), 10, tau, nodes = 81L)
  )
})

test_that(".circmix_slots() gives a slot distribution that is continuous in K", {
  slot_weights <- function(K, set_size, max_slots) {
    allocation <- .circmix_slots(K, set_size)
    weights <- numeric(max_slots + 1)
    weights[allocation$slots + 1] <- 1 - allocation$extra
    weights[allocation$slots + 2] <- weights[allocation$slots + 2] + allocation$extra
    weights
  }
  expect_equal(slot_weights(6 - 1e-9, 3, 4), slot_weights(6, 3, 4), tolerance = 1e-8)
  expect_equal(slot_weights(3 - 1e-9, 3, 4), slot_weights(3, 3, 4), tolerance = 1e-8)
  expect_equal(.circmix_slots(3.5, 1), list(slots = 3, extra = 0.5))
})

# Compares the sampler against the density it is meant to draw from. n is large
# enough that a bin proportion has a standard error near 5e-4, so the tolerance
# below is about ten standard errors and does not need a seed to be stable.
expect_samples_match_density <- function(draws, density, n_bins = 20) {
  breaks <- seq(-pi, pi, length.out = n_bins + 1)
  observed <- as.numeric(table(cut(draws, breaks))) / length(draws)
  expected <- vapply(seq_len(n_bins), function(i) {
    stats::integrate(density, breaks[i], breaks[i + 1])$value
  }, numeric(1))
  expect_lt(max(abs(observed - expected)), 0.005)
}

test_that(".rcircmix() draws from the constant-precision mixture density", {
  n <- 2e5
  mu_nt <- 1.2
  draws <- .rcircmix(
    mu = matrix(c(0, mu_nt), nrow = n, ncol = 2, byrow = TRUE),
    logw = matrix(log(c(0.5, 0.3)), nrow = n, ncol = 2, byrow = TRUE),
    logw_guess = rep(log(0.2), n), kappa = rep(8, n), tau = rep(0, n)
  )
  expect_samples_match_density(draws, function(y) {
    exp(ld_mixture(y, mu_nt = mu_nt, kappa = 8, w_mem = c(0.5, 0.3), w_guess = 0.2))
  })
})

test_that(".rcircmix() draws from the variable-precision mixture density", {
  n <- 2e5
  draws <- .rcircmix(
    mu = matrix(0, nrow = n, ncol = 1),
    logw = matrix(log(0.8), nrow = n, ncol = 1),
    logw_guess = rep(log(0.2), n), kappa = rep(10, n), tau = rep(2, n)
  )
  expect_samples_match_density(draws, function(y) {
    exp(.circmix_vp_ld(
      matrix(cos(y)), matrix(log(0.8), nrow = length(y)),
      rep(log(0.2), length(y)), rep(10, length(y)), rep(2, length(y))
    ))
  })
})

test_that(".circmix_add_variable_precision() inserts tau next to kappa", {
  spec <- list(
    parameters = list(mu = "location", kappa = "precision", thetat = "weight"),
    links = list(mu = "tan_half", kappa = "log", thetat = "logit"),
    priors = list(kappa = list(main = "normal(2, 1)")),
    init_ranges = list(kappa = c(3, 8))
  )
  out <- .circmix_add_variable_precision(spec)
  expect_equal(names(out$parameters), c("mu", "kappa", "tau", "thetat"))
  expect_equal(names(out$links), c("mu", "kappa", "tau", "thetat"))
  expect_equal(out$links$tau, "log")
  expect_true(!is.null(out$priors$tau) && !is.null(out$init_ranges$tau))
})

test_that(".circmix_recycle() extends every argument to the longest", {
  out <- .circmix_recycle(x = 1:3, mu = 0, kappa = c(5, 6, 7))
  expect_equal(out$x, 1:3)
  expect_equal(out$mu, c(0, 0, 0))
  expect_equal(out$kappa, c(5, 6, 7))
})

test_that(".circmix_bounds() derives natural-scale bounds from the links", {
  bounds <- .circmix_bounds(list(mu = "tan_half", kappa = "log", thetat = "logit"))
  expect_equal(bounds$lb, c(NA, 0, 0))
  expect_equal(bounds$ub, c(NA, NA, 1))
})

test_that(".circmix_prep_tau() is zero when the model has no tau parameter", {
  expect_equal(.circmix_prep_tau(list(dpars = list(kappa = 1), ndraws = 4), 1), rep(0, 4))
})

test_that(".circmix_prep_nodes() falls back to the default node count", {
  expect_equal(.circmix_prep_nodes(list(family = list())), 41L)
  expect_equal(.circmix_prep_nodes(list(family = list(vp_nodes = 81L))), 81L)
})

test_that(".circmix_family_vars() emits only bare indices", {
  vars <- .circmix_family_vars(vint = TRUE, n_vreal = 3)
  indexed <- grep("\\[", vars, value = TRUE)
  expect_equal(indexed, c("vint1[n]", "vreal1[n]", "vreal2[n]", "vreal3[n]"))
  expect_true(all(grepl("^(vint|vreal)[0-9]+\\[n\\]$", indexed)))
  expect_true(all(.circmix_table_vars() %in% vars))
})

test_that(".circmix_aterm() builds the response formula for each covariate combination", {
  expect_equal(as.character(.circmix_aterm("dev_rad")), "dev_rad ~ 1")
  expect_equal(
    as.character(.circmix_aterm("dev_rad", vint = "ss_numeric", vreal = c("nt1", "nt2"))),
    "dev_rad | vint(ss_numeric) + vreal(nt1, nt2) ~ 1"
  )
})

test_that(".circmix_stan_wrapper() packs the vreal arguments and passes the table through", {
  wrapper <- .circmix_stan_wrapper(
    "imm_full",
    dpars = c("mu", "kappa", "tau", "c", "a", "s"),
    vint = "ss", vreal = list(nt = 2, dist = 2), nodes = 41L
  )
  expect_match(wrapper, "real imm_full_lpdf\\(real y, real mu, real kappa, real tau")
  expect_match(wrapper, "int ss, real nt1, real nt2, real dist1, real dist2")
  expect_match(wrapper, "to_vector\\(\\{nt1, nt2\\}\\), to_vector\\(\\{dist1, dist2\\}\\)", fixed = FALSE)
  expect_match(wrapper, "imm_full_core\\(")
  expect_match(wrapper, "41, circmix_logk, circmix_dlogk, circmix_logJ_min, circmix_dlogJ")
})

test_that(".circmix_stan_wrapper() can pass the core a value its family omits", {
  wrapper <- .circmix_stan_wrapper(
    "mixture2p_simple",
    dpars = c("mu", "kappa", "thetat"),
    core_dpars = c("mu", "kappa", "0.0", "thetat")
  )
  expect_match(wrapper, "real mixture2p_simple_lpdf\\(real y, real mu, real kappa, real thetat")
  expect_false(grepl("real 0.0", wrapper, fixed = TRUE))
  expect_match(wrapper, "mixture2p_simple_core\\(y, mu, kappa, 0.0, thetat, 41")
})

# The functions are evaluated through a fixed_param run rather than through
# cmdstanr::expose_functions(), which links against RcppParallel's libtbb and so
# fails on toolchains where that disagrees with the one cmdstan was built with.
# sig_figs is raised because cmdstan writes only 6 significant digits by default,
# which alone accounts for differences of ~1e-7.
circmix_stan_values <- function(data) {
  code <- paste0(
    "functions {\n",
    read_lines2(system.file("stan_chunks", "circmix_funs.stan", package = "bmm")),
    "\n}\n",
    "data {
      int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ;
      int NC; vector[NC] cosd; vector[NC] logw; real logw_guess; vector[NC] kappa_het;
      int NKAP; vector[NKAP] kappa;
      int NJ; vector[NJ] Jtest;
      int NTAU; vector[NTAU] tau; real kappa_vp; int nodes;
      real K_slots; int ss_slots;
    }
    generated quantities {
      vector[NKAP] out_J; vector[NKAP] out_ld; vector[NJ] out_kappa;
      real out_het; vector[NTAU] out_vp; vector[2] out_slots;
      for (i in 1:NKAP) {
        out_J[i] = circmix_J(kappa[i]);
        out_ld[i] = circmix_ld(cosd, logw, logw_guess, kappa[i]);
      }
      for (i in 1:NJ) {
        out_kappa[i] = circmix_kappa(Jtest[i], logk, dlogk, logJ_min, dlogJ);
      }
      out_het = circmix_het_ld(cosd, logw, logw_guess, kappa_het);
      for (i in 1:NTAU) {
        out_vp[i] = circmix_vp_ld(cosd, logw, logw_guess, kappa_vp, tau[i], nodes,
                                  logk, dlogk, logJ_min, dlogJ);
      }
      out_slots = circmix_slots(K_slots, ss_slots);
    }"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = data, fixed_param = TRUE, iter_sampling = 1, chains = 1,
    refresh = 0, show_messages = FALSE, sig_figs = 18
  )
  function(variable) as.numeric(fit$draws(variable))
}

test_that("the Stan and R implementations of the shared core agree", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the shared Stan functions"
  )

  tab <- .circmix_kappa_table()
  y <- 0.4
  cosd <- cos(y - c(0, 1.2, -2.0))
  logw <- log(c(0.4, 0.2, 0.15))
  logw_guess <- log(0.25)
  kappa <- c(0.05, 0.5, 2, 10, 100, 1000, 6e4)
  Jtest <- c(1e-7, 1e-3, 0.5, 9.486, 500, 5e4, 1e7)
  kappa_het <- c(3, 12, 40)
  tau <- c(0, 0.5, 2)

  stan <- circmix_stan_values(list(
    NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
    logJ_min = tab$logJ_min, dlogJ = tab$dlogJ,
    NC = length(cosd), cosd = cosd, logw = logw, logw_guess = logw_guess,
    kappa_het = kappa_het, NKAP = length(kappa), kappa = kappa,
    NJ = length(Jtest), Jtest = Jtest, NTAU = length(tau), tau = tau,
    kappa_vp = 10, nodes = 41L, K_slots = 3.5, ss_slots = 2L
  ))

  n <- length(kappa)
  as_rows <- function(x) matrix(x, nrow = n, ncol = length(x), byrow = TRUE)

  expect_equal(stan("out_J"), .circmix_J(kappa), tolerance = 1e-10)
  expect_equal(stan("out_kappa"), .circmix_kappa(Jtest), tolerance = 1e-10)
  expect_equal(
    stan("out_ld"),
    .circmix_ld(as_rows(cosd), as_rows(logw), rep(logw_guess, n), kappa),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_het"),
    .circmix_het_ld(
      matrix(cosd, nrow = 1), matrix(logw, nrow = 1), logw_guess,
      matrix(kappa_het, nrow = 1)
    ),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_vp"),
    vapply(tau, function(c_tau) {
      .circmix_vp_ld(matrix(cosd, nrow = 1), matrix(logw, nrow = 1), logw_guess, 10, c_tau)
    }, numeric(1)),
    tolerance = 1e-12
  )
  expect_equal(stan("out_slots"), unname(unlist(.circmix_slots(3.5, 2))))
})
