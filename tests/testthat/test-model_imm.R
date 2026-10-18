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

test_that("imm() and dimm() default to 81 quadrature nodes", {
  # the core refuses gamma shapes below 60 / vp_nodes, which at 41 nodes is 1.46
  # and reaches into real-data posteriors of mixture2p and mixture3p (#398)
  expect_equal(imm_model(variable_precision = TRUE)$vp_nodes, 81L)
  expect_equal(formals(imm)$vp_nodes, 81L)
  expect_equal(formals(dimm)$vp_nodes, 81L)
})

test_that("dimm() refuses node counts the model refuses", {
  # 21 nodes miss the density by about 2e-3, and a fractional count built a grid
  # that recycled with a warning
  expect_error(dimm(0.1, tau = 1, vp_nodes = 21), "odd number of at least 41")
  expect_error(dimm(0.1, tau = 1, vp_nodes = 42), "odd number of at least 41")
  expect_error(dimm(0.1, tau = 1, vp_nodes = 80.5), "odd number of at least 41")
  expect_error(dimm(0.1, tau = 1, vp_nodes = c(81, 121)), "odd number of at least 41")
  expect_equal(dimm(0.1, tau = 1, vp_nodes = 81), dimm(0.1, tau = 1, vp_nodes = 81L))
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

test_that("rimm() keeps its spread at high precision", {
  # a sampler on a grid of locations loses spread at large kappa (#434); b is
  # tiny so that no draw is a guess
  withr::local_seed(434)
  draws <- rimm(1e5, mu = 0, dist = 0, c = 1, a = 1, b = 1e-12, s = 1, kappa = 1e4)
  # 1 / sqrt(J(1e4)), the sd of a von Mises at this precision
  expect_equal(sd(draws), 0.0100003, tolerance = 0.01)
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
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
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

test_that("a formula naming mu1 still fits, with a warning", {
  expect_warning(
    fit <- bmm(bmf(mu1 ~ 1, kappa ~ 1, c ~ 1, a ~ 1), oberauer_lin_2017,
      imm_model("abc"),
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "Use `bmf\\(mu ~ 1\\)` instead"
  )
  expect_true("mu" %in% names(fit$formula$pforms))
  # mu is fixed to 0 unless a formula frees it, so a default intercept prior
  # on the location shows that the mu1 formula reached it
  mu_rows <- fit$prior[fit$prior$class == "Intercept" & fit$prior$dpar == "", ]
  expect_equal(mu_rows$prior, "normal(0, 0.5)")
  expect_false("mu" %in% names(fit$bmm$model$fixed_parameters))
})

test_that("imm cites van den Berg et al. only with variable precision", {
  expect_length(model_citation(imm_model()), 1)
  vp <- model_citation(imm_model(variable_precision = TRUE))
  expect_length(vp, 2)
  expect_identical(vp[2], model_citation(mixture2p("y", variable_precision = TRUE))[2])
})

test_that("an intercept next to set_size reaches imm without a set-size-1 pin", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  # set size 1 is left to the priors; check_formula.imm() says which coefficients
  fit <- suppressWarnings(
    bmm(bmf(kappa ~ 1 + set_size, c ~ 1, a ~ 1 + set_size, s ~ 0 + set_size),
      oberauer_lin_2017, imm_model(),
      backend = "mock", mock_fit = 1, rename = FALSE
    )
  )
  expect_false(any(grepl("constant", fit$prior$prior[fit$prior$dpar %in% c("a", "s")])))
})

test_that("each imm family refuses the expected response under its own name", {
  formulas <- list(
    full = bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1),
    bsc = bmf(kappa ~ 1, c ~ 1, s ~ 1),
    abc = bmf(kappa ~ 1, c ~ 1, a ~ 1)
  )
  for (version in names(formulas)) {
    model <- check_model(imm_model(version), oberauer_lin_2017, formulas[[version]])
    data <- check_data(model, oberauer_lin_2017, formulas[[version]])
    family <- configure_model(model, data, formulas[[version]])$formula$family
    # add_posterior_epred() finds the function for old fits by the family name
    expect_identical(
      family$posterior_epred,
      get(paste0("posterior_epred_", family$name), envir = asNamespace("bmm"))
    )
    expect_error(family$posterior_epred(NULL), "not defined for the imm model")
  }
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
  skip_unless_stan_tests()
  skip_on_cran()

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

# names of the coefficients the set-size-1 warnings name, from the warnings
# check_formula() gives for the formula
set_size_one_names <- function(formula, data = oberauer_lin_2017, version = "full") {
  model <- check_model(imm_model(version), data, formula)
  data <- check_data(model, data, formula)
  messages <- character(0)
  withCallingHandlers(
    check_formula(model, data, formula),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  messages <- grep("set size of 1", messages, value = TRUE)
  as.character(sort(unlist(regmatches(messages, gregexpr("(b|sd)_[A-Za-z0-9_.]+(?=')", messages, perl = TRUE)))))
}

test_that("set size 1 leaves a and s to their priors instead of pinning them", {
  formula <- bmf(kappa ~ 1, c ~ 1, a ~ 0 + set_size + (0 + set_size | ID),
    s ~ 0 + set_size + (0 + set_size | ID))
  pr <- suppressWarnings(default_prior(formula, oberauer_lin_2017, imm_model()))
  expect_false(any(grepl("constant", pr$prior[pr$dpar %in% c("a", "s")], fixed = TRUE)))
  expect_equal(unique(pr$prior[pr$class == "sd" & pr$dpar %in% c("a", "s") & pr$group == ""]),
    "exponential(1)")
})

test_that("an intercept next to set_size is accepted, with the set-size-1 warning", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  expect_warning(
    bmm(bmf(kappa ~ 1 + set_size, c ~ 1, a ~ 1 + set_size, s ~ 1 + set_size),
      oberauer_lin_2017, imm_model(),
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "set size of 1"
  )
})

test_that("a missing set size leaves the set-size-1 check to the complete rows", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  # the check read the rows complete in the formula's variables only, so one
  # NA set size reached svd() and stopped every fit with set-size-1 trials
  dat <- oberauer_lin_2017
  dat$set_size[5] <- NA
  fit_mock <- function(formula, version = "full") {
    withCallingHandlers(
      bmm(formula, dat, imm_model(version),
        backend = "mock", mock_fit = 1, rename = FALSE
      ),
      # brms dropping the row is the expected outcome; any other warning shows
      warning = function(w) {
        if (grepl("Rows containing NAs were excluded", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  expect_s3_class(fit_mock(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1)), "bmmfit")
  expect_s3_class(fit_mock(bmf(kappa ~ 1, c ~ 1, s ~ 1), "bsc"), "bmmfit")
})

test_that("the set-size-1 warning names exactly the coefficients no other trial informs", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  set_sizes <- paste0("set_size", 2:8)
  # s is absent at set size 1
  expect_equal(set_size_one_names(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1)), character(0))
  expect_equal(set_size_one_names(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 0 + set_size)),
    "b_s_set_size1")
  expect_equal(set_size_one_names(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1 + set_size)),
    sort(c("b_s_Intercept", paste0("b_s_", set_sizes))))
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 1, a ~ 1,
      s ~ 0 + set_size + (0 + set_size | ID))),
    c("b_s_set_size1", "sd_ID__s_set_size1")
  )
  # a enters set size 1 through c + a, so c from the other trials fixes it
  expect_equal(set_size_one_names(bmf(kappa ~ 1, c ~ 1, a ~ 0 + set_size, s ~ 1)),
    character(0))
  expect_equal(set_size_one_names(bmf(kappa ~ 1, c ~ 1, a ~ 1 + set_size, s ~ 1)),
    character(0))
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 1, a ~ 0 + set_size + (0 + set_size | ID), s ~ 1)),
    character(0)
  )
  expect_equal(set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size, a ~ 1, s ~ 1)),
    character(0))
  # with c free per set size too, set size 1 fixes only their sum
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size, a ~ 0 + set_size, s ~ 1)),
    c("b_a_set_size1", "b_c_set_size1")
  )
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size, a ~ 1 + set_size, s ~ 1)),
    sort(c("b_a_Intercept", paste0("b_a_", set_sizes), "b_c_set_size1"))
  )
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size,
      a ~ 0 + set_size + (0 + set_size | ID), s ~ 1)),
    c("b_a_set_size1", "b_c_set_size1")
  )
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size + (0 + set_size | ID),
      a ~ 0 + set_size + (0 + set_size | ID), s ~ 1)),
    c("b_a_set_size1", "b_c_set_size1", "sd_ID__a_set_size1", "sd_ID__c_set_size1")
  )
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size, a ~ 0 + set_size),
      version = "abc"),
    c("b_a_set_size1", "b_c_set_size1")
  )
  # bsc has no a, and c alone is the target's activation at set size 1
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size, s ~ 0 + set_size),
      version = "bsc"),
    "b_s_set_size1"
  )
  # a numeric set size draws its line through the set sizes above 1
  numeric_ss <- oberauer_lin_2017
  numeric_ss$set_size <- as.numeric(as.character(numeric_ss$set_size))
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 1 + set_size, a ~ 1 + set_size,
      s ~ 1 + set_size), numeric_ss),
    character(0)
  )
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size, a ~ 0 + set_size,
      s ~ 0 + set_size), numeric_ss),
    character(0)
  )
  # without set size 1 nothing is left to the prior
  no_single <- oberauer_lin_2017[oberauer_lin_2017$set_size != 1, ]
  expect_equal(
    set_size_one_names(bmf(kappa ~ 1, c ~ 0 + set_size, a ~ 0 + set_size,
      s ~ 0 + set_size), no_single),
    character(0)
  )
})

test_that("the generated likelihood passes each parameter to its place in the core", {
  skip_unless_stan_tests()
  skip_on_cran()
  # the cores are checked against R above; this checks the generated wrappers
  # that brms calls, whose argument order is built from the version tables
  cases <- list(
    full = list(version = "full", formula = bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1)),
    bsc = list(version = "bsc", formula = bmf(kappa ~ 1, c ~ 1, s ~ 1)),
    abc = list(version = "abc", formula = bmf(kappa ~ 1, c ~ 1, a ~ 1)),
    full_vp = list(version = "full", formula = bmf(kappa ~ 1, tau ~ 1, c ~ 1, a ~ 1, s ~ 1),
      variable_precision = TRUE)
  )
  values <- c(
    y = 0.7, mu = 0.1, kappa = 6, tau = 0.8, c = 3, a = 0.5, s = 1.5, b = 1.2,
    ss = 5, setNames(c(2, -1.5, 1, -0.4, 2.8, 0.3, -2), paste0("nt", 1:7)),
    setNames(c(0.3, 1.5, 2, 0.8, 2.5, 1.1, 0.6), paste0("dist", 1:7))
  )
  wrappers <- character(0)
  calls <- character(0)
  for (name in names(cases)) {
    case <- cases[[name]]
    model <- check_model(
      imm_model(case$version, variable_precision = isTRUE(case$variable_precision)),
      oberauer_lin_2017, case$formula
    )
    data <- check_data(model, oberauer_lin_2017, case$formula)
    family <- paste0("imm_", case$version)
    code <- configure_model(model, data, case$formula)$stanvars[[paste0(family, "_lpdf")]]$scode
    code <- gsub(paste0(family, "_lpdf"), paste0("w_", name, "_lpdf"), code, fixed = TRUE)
    wrappers <- c(wrappers, code)
    signature <- sub(".*_lpdf\\((.*?)\\).*", "\\1", strsplit(code, "\n")[[1]][1])
    arg_names <- sub(".* ", "", strsplit(signature, ", ")[[1]])
    arg_values <- ifelse(arg_names %in% names(values),
      format(values[arg_names], digits = 17), arg_names
    )
    arg_values[arg_names == "ss"] <- as.character(as.integer(values[["ss"]]))
    calls <- c(calls, glue::glue(
      "out[{length(calls) + 1}] = w_{name}_lpdf({arg_values[1]} | ",
      "{paste(arg_values[-1], collapse = ', ')});"
    ))
  }
  program <- paste0(
    "functions {\n",
    read_lines2(system.file("stan_chunks", "circmix_funs.stan", package = "bmm")), "\n",
    read_lines2(system.file("stan_chunks", "imm_funs.stan", package = "bmm")), "\n",
    paste(wrappers, collapse = "\n"), "}\n",
    "data { int NK; vector[NK] circmix_logk; vector[NK] circmix_dlogk;",
    " real circmix_logJ_min; real circmix_dlogJ; }\n",
    "generated quantities { vector[", length(calls), "] out;\n",
    paste(calls, collapse = "\n"), "\n}\n"
  )
  tab <- .circmix_kappa_table()
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = list(
      NK = length(tab$logkappa), circmix_logk = tab$logkappa,
      circmix_dlogk = tab$dlogkappa, circmix_logJ_min = tab$logJ_min,
      circmix_dlogJ = tab$dlogJ
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  r_side <- vapply(cases, function(case) {
    v <- case$version
    .dimm_version(values[["y"]], values[["mu"]], values[["kappa"]],
      c = values[["c"]], a = if (v != "bsc") values[["a"]],
      s = if (v != "abc") values[["s"]], b = values[["b"]], set_size = values[["ss"]],
      nt = values[paste0("nt", 1:7)], dist = values[paste0("dist", 1:7)],
      tau = if (isTRUE(case$variable_precision)) values[["tau"]] else 0,
      nodes = 81L, version = v
    )
  }, numeric(1))
  expect_equal(as.numeric(fit$draws("out")), unname(r_side), tolerance = 1e-10)
})

test_that("the link of the fixed background activation cannot be changed", {
  # b = 0 on the log link is the background of 1 the activations are weighed
  # against; on another link the same fixed value is another model
  expect_error(
    imm_model(links = list(b = "identity")),
    "The link of 'b' cannot be changed in imm"
  )
  expect_error(imm_model("abc", links = list(b = "softplus")), "cannot be changed")
  expect_silent(imm_model(links = list(c = "softplus")))
  expect_silent(imm_model(links = list(b = "log")))
})

test_that("the set-size-1 warnings of s and of a with c share one wording", {
  formula <- bmf(kappa ~ 1, c ~ 0 + set_size, a ~ 0 + set_size, s ~ 0 + set_size)
  model <- check_model(imm_model(), oberauer_lin_2017, formula)
  data <- check_data(model, oberauer_lin_2017, formula)
  messages <- character(0)
  withCallingHandlers(
    check_formula(model, data, formula),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(sort(gsub("\\s+", " ", messages)), c(
    paste(
      "Trials with a set size of 1 have no non-target, so they inform 'a' and 'c'",
      "only through their sum. The other trials do not determine 'b_a_set_size1',",
      "'b_c_set_size1' in the formulas for 'a' and 'c', so the posteriors of these",
      "coefficients depend on their priors. Give set size 1 a level of its own in a",
      "factor, so that only the set-size-1 coefficients depend on their priors."
    ),
    paste(
      "Trials with a set size of 1 have no non-target, so they say nothing about 's'.",
      "The other trials do not determine 'b_s_set_size1' in the formula for 's', so",
      "the posteriors of these coefficients depend on their priors. Give set size 1 a",
      "level of its own in a factor, so that only the set-size-1 coefficients depend",
      "on their priors."
    )
  ))
})

test_that("the deprecated IMM constructors warn and build the model imm() builds", {
  nt <- paste0("col_nt", 1:7)
  dist <- paste0("dist_nt", 1:7)
  cases <- list(
    full = list(old = function(...) IMMfull("dev_rad", nt, dist, ...), name = "IMMfull"),
    bsc = list(old = function(...) IMMbsc("dev_rad", nt, dist, ...), name = "IMMbsc"),
    abc = list(old = function(...) IMMabc("dev_rad", nt, ...), name = "IMMabc")
  )
  for (version in names(cases)) {
    case <- cases[[version]]
    current <- if (version == "abc") {
      imm("dev_rad", nt_features = nt, set_size = "set_size", version = "abc")
    } else {
      imm("dev_rad", nt_features = nt, nt_distances = dist, set_size = "set_size",
          version = version)
    }
    expect_warning(old <- case$old(set_size = "set_size"),
                   paste0("`", case$name, "()` is deprecated"), fixed = TRUE)
    expect_identical(class(old), class(current))
    expect_identical(old$resp_vars, current$resp_vars)
    expect_identical(old$other_vars, current$other_vars)
    # the argument name of bmm < 0.4 still reaches set_size, with a warning
    warnings <- character()
    withCallingHandlers(old <- case$old(setsize = "set_size"), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_true(any(grepl("'setsize' is deprecated", warnings, fixed = TRUE)), label = version)
    expect_identical(old$other_vars, current$other_vars)
  }
})
