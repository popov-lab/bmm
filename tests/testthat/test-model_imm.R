imm_model <- function(version = "full", ...) {
  imm("dev_rad",
    nt_features = "col_nt",
    nt_distances = if (version == "abc") NULL else "dist_nt",
    set_size = "set_size", regex = TRUE, version = version, ...
  )
}

test_that("imm() exposes the parameters each version implies", {
  expect_equal(names(imm_model()$parameters), c("mu", "kappa", "a", "c", "s", "b"))
  expect_equal(names(imm_model("bsc")$parameters), c("mu", "kappa", "c", "s", "b"))
  expect_equal(names(imm_model("abc")$parameters), c("mu", "kappa", "a", "c", "b"))
  expect_equal(
    names(imm_model(variable_precision = TRUE)$parameters),
    c("mu", "kappa", "tau", "a", "c", "s", "b")
  )
  expect_error(imm_model(vp_nodes = 21), "at least 41")
})

test_that("the background activation is a visible fixed parameter", {
  model <- imm_model()
  expect_equal(model$fixed_parameters, list(mu = 0, b = 0))
  # a log link makes the fixed value of 0 an activation of 1, the reference the
  # remaining activations are weighed against
  expect_equal(model$links$b, "log")
})

# The previous implementation was a brms mixture whose weights were a softmax
# over log activations, with a -100 sentinel switching off components the
# trial's set size did not have. dimm() computes the same activations, so it is
# the reference the version-specific weights have to reproduce.
test_that("the version weights reproduce the activations dimm() computes", {
  x <- seq(-pi, pi, length.out = 30)
  nt <- c(2, -1.5, 1)
  dist <- c(0.5, 1.5, 2)
  expect_equal(
    .dimm_version(x, 0, 6, c = 5, a = 2, s = 2, b = 1, set_size = 4, nt = nt,
      dist = dist, tau = 0, nodes = 41L, version = "full"),
    dimm(x, mu = c(0, nt), dist = c(0, dist), c = 5, a = 2, b = 1, s = 2,
      kappa = 6, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("the reduced versions drop the term they are named for", {
  x <- seq(-pi, pi, length.out = 30)
  nt <- c(2, -1.5, 1)
  dist <- c(0.5, 1.5, 2)
  # bsc has no cue-independent activation, so it matches full with a -> 0
  expect_equal(
    .dimm_version(x, 0, 6, c = 5, a = NULL, s = 2, b = 1, set_size = 4, nt = nt,
      dist = dist, tau = 0, nodes = 41L, version = "bsc"),
    .dimm_version(x, 0, 6, c = 5, a = 1e-12, s = 2, b = 1, set_size = 4, nt = nt,
      dist = dist, tau = 0, nodes = 41L, version = "full"),
    tolerance = 1e-9
  )
  # abc gives the target c + a and every non-target a, which dimm() cannot be
  # coaxed into producing, since it applies the same gradient to all items
  activation <- c(5 + 2, rep(2, 3), 1)
  probs <- activation / sum(activation)
  locations <- c(0, nt)
  reference <- matrixStats::rowLogSumExps(cbind(
    vapply(seq_along(locations), function(i) {
      log(probs[i]) + brms::dvon_mises(x, mu = locations[i], kappa = 6, log = TRUE)
    }, numeric(length(x))),
    log(probs[5]) + stats::dunif(x, -pi, pi, log = TRUE)
  ))
  expect_equal(
    .dimm_version(x, 0, 6, c = 5, a = 2, s = NULL, b = 1, set_size = 4, nt = nt,
      dist = dist, tau = 0, nodes = 41L, version = "abc"),
    reference,
    tolerance = 1e-12
  )
})

test_that("a set size of one leaves only the target and the background", {
  x <- seq(-pi, pi, length.out = 30)
  expect_equal(
    .dimm_version(x, 0, 6, c = 5, a = 2, s = 2, b = 1, set_size = 1,
      nt = c(2, -1.5), dist = c(0.5, 1.5), tau = 0, nodes = 41L, version = "full"),
    dmixture2p(x, mu = 0, kappa = 6, p_mem = 7 / 8, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("dimm() is a normalised density with and without variable precision", {
  for (tau in c(0, 1.5)) {
    total <- stats::integrate(
      function(x) {
        dimm(x, mu = c(0, 2, -1.5), dist = c(0, 0.5, 2), c = 5, a = 2, b = 1,
          s = 2, kappa = 8, tau = tau)
      },
      -pi, pi
    )
    expect_equal(total$value, 1, tolerance = 1e-6)
  }
})

test_that("rimm() draws from dimm()", {
  compare <- function(n = 2e5, n_bins = 20, ...) {
    draws <- rimm(n, ...)
    breaks <- seq(-pi, pi, length.out = n_bins + 1)
    observed <- as.numeric(table(cut(draws, breaks))) / n
    expected <- vapply(seq_len(n_bins), function(i) {
      stats::integrate(function(x) dimm(x, ...), breaks[i], breaks[i + 1])$value
    }, numeric(1))
    expect_lt(max(abs(observed - expected)), 0.005)
  }
  compare(mu = c(0, 2, -1.5), dist = c(0, 0.5, 2), c = 5, a = 2, b = 1, s = 2, kappa = 8)
  compare(mu = c(0, 2, -1.5), dist = c(0, 0.5, 2), c = 5, a = 2, b = 1, s = 2, kappa = 8, tau = 1)
})

test_that("configure_model() drops the set-size sentinel machinery", {
  code <- stancode(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1), oberauer_lin_2017, imm_model())
  expect_match(code, "imm_full_lpdf")
  expect_false(grepl("LureIdx", code, fixed = TRUE))
  expect_false(grepl("-100", code, fixed = TRUE))
  expect_false(grepl("expS", code, fixed = TRUE))
  # the non-target features come first in the vreal block, the distances after
  expect_match(code, "array[N] real vreal14", fixed = TRUE)
})

test_that("the abc version asks for no distances", {
  code <- stancode(bmf(kappa ~ 1, c ~ 1, a ~ 1), oberauer_lin_2017, imm_model("abc"))
  expect_match(code, "array[N] real vreal7", fixed = TRUE)
  expect_false(grepl("vreal8", code, fixed = TRUE))
})

test_that("every version runs through the bmm() pipeline", {
  dat <- oberauer_lin_2017
  expect_silent(bmm(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1), dat, imm_model(),
    backend = "mock", mock_fit = 1, rename = FALSE))
  expect_silent(bmm(bmf(kappa ~ 1, c ~ 1, s ~ 1), dat, imm_model("bsc"),
    backend = "mock", mock_fit = 1, rename = FALSE))
  expect_silent(bmm(bmf(kappa ~ 1, c ~ 1, a ~ 1), dat, imm_model("abc"),
    backend = "mock", mock_fit = 1, rename = FALSE))
  expect_silent(bmm(bmf(kappa ~ 1, tau ~ 1, c ~ 1, a ~ 1, s ~ 1), dat,
    imm_model(variable_precision = TRUE),
    backend = "mock", mock_fit = 1, rename = FALSE))
})

test_that("check_data() still validates the non-target distances", {
  dat <- oberauer_lin_2017
  dat$dist_nt1 <- -1
  formula <- bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1)
  # check_model() expands the regex into the actual column names first
  model <- check_model(imm_model(), dat, formula)
  expect_error(
    check_data(model, dat, formula),
    "distances to the target need to be postive"
  )
})

test_that("the Stan and R likelihoods agree for every version", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the imm likelihood"
  )

  tab <- .circmix_kappa_table()
  code <- paste0(
    "functions {\n",
    read_lines2(system.file("stan_chunks", "circmix_funs.stan", package = "bmm")),
    "\n",
    read_lines2(system.file("stan_chunks", "imm_funs.stan", package = "bmm")),
    "\n}\n",
    "data {
      int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ;
      int N; int NT; vector[N] y; vector[N] kappa; vector[N] tau;
      vector[N] c; vector[N] a; vector[N] s; vector[N] b;
      array[N] int ss; matrix[N, NT] nt; matrix[N, NT] dist;
    }
    generated quantities {
      vector[N] out_full; vector[N] out_bsc; vector[N] out_abc;
      for (i in 1:N) {
        out_full[i] = imm_full_core(y[i], 0.0, kappa[i], tau[i], c[i], a[i],
                        s[i], b[i], ss[i], to_vector(nt[i]), to_vector(dist[i]),
                        41, logk, dlogk, logJ_min, dlogJ);
        out_bsc[i] = imm_bsc_core(y[i], 0.0, kappa[i], tau[i], c[i], s[i], b[i],
                       ss[i], to_vector(nt[i]), to_vector(dist[i]), 41, logk,
                       dlogk, logJ_min, dlogJ);
        out_abc[i] = imm_abc_core(y[i], 0.0, kappa[i], tau[i], c[i], a[i], b[i],
                       ss[i], to_vector(nt[i]), 41, logk, dlogk, logJ_min, dlogJ);
      }
    }"
  )

  y <- c(-2.5, -0.8, 0, 0.4, 1.9, 3.0)
  kappa <- c(0.5, 2, 6, 6, 20, 80)
  tau <- c(0, 0, 0, 1, 0.5, 2)
  cc <- c(0.5, 1, 3, 5, 8, 12)
  a <- c(0.1, 0.3, 1, 2, 0.5, 4)
  ss_grad <- c(0.2, 1, 2, 3, 0.5, 5)
  b <- c(1, 1, 1, 2, 0.5, 1)
  ss <- c(1L, 2L, 3L, 4L, 6L, 8L)
  nt <- matrix(seq(-3, 3, length.out = 6 * 7), nrow = 6)
  dist <- matrix(seq(0.1, 3, length.out = 6 * 7), nrow = 6)

  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = list(
      NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
      logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, N = length(y), NT = ncol(nt),
      y = y, kappa = kappa, tau = tau, c = cc, a = a, s = ss_grad, b = b,
      ss = ss, nt = nt, dist = dist
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  stan <- function(variable) as.numeric(fit$draws(variable))

  r_side <- function(version) {
    vapply(seq_along(y), function(i) {
      .dimm_version(y[i], 0, kappa[i],
        c = cc[i], a = if (version != "bsc") a[i], s = if (version != "abc") ss_grad[i],
        b = b[i], set_size = ss[i], nt = nt[i, ], dist = dist[i, ],
        tau = tau[i], nodes = 41L, version = version
      )
    }, numeric(1))
  }

  expect_equal(stan("out_full"), r_side("full"), tolerance = 1e-12)
  expect_equal(stan("out_bsc"), r_side("bsc"), tolerance = 1e-12)
  expect_equal(stan("out_abc"), r_side("abc"), tolerance = 1e-12)
})
