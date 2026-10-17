nt_model <- function(...) {
  mixture3p("dev_rad", nt_features = "col_nt", set_size = "set_size", regex = TRUE, ...)
}

test_that("mixture3p() and dmixture3p() default to 81 quadrature nodes", {
  # 41 nodes put the reject threshold inside a real-data posterior: slot
  # averaging + VP on three subjects of oberauer_lin_2017 diverged on 165 of
  # 200 transitions at 41 nodes (#398)
  expect_equal(nt_model(variable_precision = TRUE)$vp_nodes, 81L)
  expect_equal(formals(mixture3p)$vp_nodes, 81L)
  expect_equal(formals(dmixture3p)$vp_nodes, 81L)
})

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

test_that("the capacity versions put effects priors on K and pnt", {
  for (version in c("slot", "slot_averaging")) {
    pr <- default_prior(
      bmf(kappa ~ 1, K ~ 1 + session, pnt ~ 1 + session), small_data(),
      mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size",
                version = version)
    )
    expect_equal(pr[pr$class == "b" & pr$coef == "" & pr$dpar == "K", ]$prior, "normal(0, 0.3)")
    expect_equal(pr[pr$class == "b" & pr$coef == "" & pr$dpar == "pnt", ]$prior, "normal(0, 0.5)")
  }
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

test_that("dmixture3p() refuses node counts the model refuses", {
  # 21 nodes miss the log density by about 1e-2 at a gamma shape near 3
  expect_error(
    dmixture3p(0.1, kappa = 4, tau = 1.5, vp_nodes = 21),
    "odd number of at least 41"
  )
  expect_error(dmixture3p(0.1, tau = 1, vp_nodes = 42), "odd number of at least 41")
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

test_that("rmixture3p() is exact at high precision", {
  # the grid-based sampler this replaces bounded the density on a grid and
  # missed the peak between grid points: an SD of 0.01067 at this mu (#434)
  # 1e-4 is about six Monte Carlo standard errors of the SD at 2e5 draws
  withr::local_seed(434)
  draws <- rmixture3p(2e5, mu = 0.0087, kappa = 1e4, p_mem = 1, p_nt = 0)
  expect_lt(abs(sd(draws) - 0.01), 1e-4)
  draws <- rmixture3p(2e5, mu = c(0.0087, 2, -1.5), kappa = 1e4, p_mem = 1, p_nt = 0)
  expect_lt(abs(sd(draws) - 0.01), 1e-4)
})

test_that("posterior_predict() samples the density log_lik() evaluates, for every version", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(398)
  n <- 2e5
  nt <- c(2, -1.5)
  # K / set_size = 1.5: the slot version stores every item, slot averaging
  # gives the reported item one or two slots
  pars <- list(mu = 0, kappa = 6, thetat = 0.6, thetant = -0.4, K = 4.5, pnt = 0.3)
  density <- function(x, version, tau) {
    exp(switch(version,
      simple = .dmixture3p_simple(x, 0, 6, 0.6, -0.4, 3, nt, tau, 41L),
      slot = .dmixture3p_slot(x, 0, 6, 4.5, 0.3, 3, nt, tau, 41L),
      slot_averaging = .dmixture3p_slot_averaging(x, 0, 6, 4.5, 0.3, 3, nt, tau, 41L)
    ))
  }
  breaks <- seq(-pi, pi, length.out = 21)
  for (tau in c(0, 1)) {
    dpars <- lapply(c(pars, if (tau > 0) list(tau = tau)), function(x) matrix(x, n, 1))
    prep <- structure(
      list(
        dpars = dpars, data = list(Y = 0, vint1 = 3L, vreal1 = nt[1], vreal2 = nt[2]),
        ndraws = n, family = list(vp_nodes = 41L)
      ),
      class = "brmsprep"
    )
    for (version in c("simple", "slot", "slot_averaging")) {
      draws <- .mixture3p_posterior_predict(version)(1, prep)
      observed <- as.numeric(table(cut(draws, breaks))) / n
      expected <- vapply(seq_len(20), function(i) {
        stats::integrate(density, breaks[i], breaks[i + 1], version = version, tau = tau)$value
      }, numeric(1))
      expect_lt(max(abs(observed - expected)), 0.005)
    }
  }
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
    "do not determine 'b_thetant_set_size1' in the formula for 'thetant'"
  )
})

test_that("a missing set size leaves the set-size-1 check to the complete rows", {
  # the check read the rows complete in the formula's variables only, so one
  # NA set size reached svd() and stopped every fit with set-size-1 trials
  dat <- oberauer_lin_2017
  dat$set_size[5] <- NA
  fit_mock <- function(formula, version = "simple") {
    withCallingHandlers(
      bmm(formula, dat, nt_model(version = version),
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
  expect_s3_class(fit_mock(bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)), "bmmfit")
  expect_s3_class(fit_mock(bmf(kappa ~ 1, K ~ 1, pnt ~ 1), "slot"), "bmmfit")
})

test_that("the set-size-1 warning names exactly the coefficients no other trial informs", {
  dat <- oberauer_lin_2017
  dat$ss_numeric <- as.numeric(as.character(dat$set_size))
  dat$ss_num <- dat$ss_numeric
  dat$single <- factor(ifelse(dat$ss_numeric == 1, "single", "multi"))
  named <- function(rhs, data = dat, version = "slot") {
    par <- if (version == "simple") "thetant" else "pnt"
    ff <- bmf(kappa ~ 1) + bmf(stats::as.formula(paste(par, rhs)))
    msg <- NULL
    withCallingHandlers(
      .mixture3p_warn_set_size_one(nt_model(version = version), data, ff),
      warning = function(w) {
        msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
    if (is.null(msg)) {
      return(character(0))
    }
    gsub("'", "", unlist(regmatches(msg, gregexpr("'(b|sd)_[^']+'", msg))))
  }
  levels_2_to_8 <- paste0("set_size", 2:8)

  # every coefficient is informed by the trials with non-targets
  expect_equal(named("~ 1"), character(0))
  expect_equal(named("~ 1 + (1 | ID)"), character(0))
  expect_equal(named("~ 0 + set_size", droplevels(dat[dat$ss_numeric > 1, ])), character(0))
  # a numeric set size is a line, fixed by set sizes 2 to 8
  expect_equal(named("~ 1 + ss_num"), character(0))
  expect_equal(named("~ 0 + ss_num"), character(0))

  expect_equal(named("~ 0 + set_size"), "b_pnt_set_size1")
  expect_equal(named("~ 0 + set_size", version = "simple"), "b_thetant_set_size1")
  # with set size 1 as the reference, every coefficient is a contrast with it
  expect_equal(named("~ 1 + set_size"), paste0("b_pnt_", c("Intercept", levels_2_to_8)))
  expect_equal(
    named("~ set_size + (set_size | ID)"),
    c(
      paste0("b_pnt_", c("Intercept", levels_2_to_8)),
      paste0("sd_ID__pnt_", c("Intercept", levels_2_to_8))
    )
  )
  expect_equal(
    named("~ 0 + set_size + (0 + set_size || ID)"),
    c("b_pnt_set_size1", "sd_ID__pnt_set_size1")
  )
  # a predictor that only varies at set size 1 is just as uninformed
  expect_equal(named("~ 0 + single"), "b_pnt_singlesingle")
  # two set sizes, one of them 1: the line through them has one point left
  expect_equal(
    named("~ 1 + ss_num", droplevels(dat[dat$ss_numeric %in% c(1, 4), ])),
    c("b_pnt_Intercept", "b_pnt_ss_num")
  )
  # a variable missing from the data is left to brms, which names it
  expect_equal(named("~ 0 + typo"), character(0))
  expect_error(
    bmm(bmf(kappa ~ 1, K ~ 1, pnt ~ 0 + typo), oberauer_lin_2017, nt_model(version = "slot"),
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "can neither be found in 'data'"
  )
})

test_that("every version runs through the bmm() pipeline", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
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

test_that("every version's family refuses the expected response by name", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  dat <- oberauer_lin_2017
  formulas <- list(
    simple = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1),
    slot = bmf(kappa ~ 1, K ~ 1, pnt ~ 1),
    slot_averaging = bmf(kappa ~ 1, K ~ 1, pnt ~ 1)
  )
  for (version in names(formulas)) {
    fit <- bmm(formulas[[version]], dat, nt_model(version = version),
      backend = "mock", mock_fit = 1, rename = FALSE
    )
    family <- fit$formula$family
    expect_identical(
      family$posterior_epred,
      get(paste0("posterior_epred_mixture3p_", version), envir = asNamespace("bmm"))
    )
    expect_error(brms::posterior_epred(fit), "not defined for the mixture3p model")
    # a fit saved without the function gets it back by family name
    fit$formula$family$posterior_epred <- NULL
    expect_identical(add_posterior_epred(fit)$formula$family$posterior_epred, family$posterior_epred)
  }
})

test_that("the capacity versions also cite the slot model", {
  expect_length(model_citation(nt_model()), 1)
  expect_match(model_citation(nt_model(version = "slot"))[2], "Zhang, W., & Luck")
})

test_that("a formula for mu1 still reaches mu, with a warning", {
  expect_warning(
    fit <- bmm(bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1, mu1 ~ 1), oberauer_lin_2017,
      nt_model(), backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "renamed to 'mu'"
  )
  # mu is fixed to 0 unless a formula frees it, so a default intercept prior
  # on the location shows that the mu1 formula reached it
  mu_rows <- fit$prior[fit$prior$class == "Intercept" & fit$prior$dpar == "", ]
  expect_equal(mu_rows$prior, "normal(0, 0.5)")
})

test_that("the Stan and R likelihoods agree for every version", {
  skip_unless_stan_tests()
  skip_on_cran()

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

test_that("the generated likelihood passes each parameter to its place in the core", {
  skip_unless_stan_tests()
  skip_on_cran()
  # the cores are checked against R above; this checks the generated wrappers
  # that brms calls, whose argument order is built from the version tables
  dat <- oberauer_lin_2017
  cases <- list(
    simple = list(
      version = "simple", formula = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
    ),
    simple_vp = list(
      version = "simple", variable_precision = TRUE,
      formula = bmf(kappa ~ 1, tau ~ 1, thetat ~ 1, thetant ~ 1)
    ),
    slot = list(version = "slot", formula = bmf(kappa ~ 1, K ~ 1, pnt ~ 1)),
    slot_averaging = list(
      version = "slot_averaging", formula = bmf(kappa ~ 1, K ~ 1, pnt ~ 1)
    )
  )
  values <- c(
    y = 0.7, mu = 0.1, kappa = 6, tau = 0.8, thetat = 0.9, thetant = -0.4,
    K = 4.5, pnt = 0.3, ss = 5,
    setNames(c(2, -1.5, 1, -0.4, 2.8, 0.3, -2), paste0("nt", 1:7))
  )
  wrappers <- character(0)
  calls <- character(0)
  for (name in names(cases)) {
    case <- cases[[name]]
    model <- check_model(
      nt_model(version = case$version, variable_precision = isTRUE(case$variable_precision)),
      dat, case$formula
    )
    data <- check_data(model, dat, case$formula)
    family <- paste0("mixture3p_", case$version)
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
    read_lines2(system.file("stan_chunks", "mixture3p_funs.stan", package = "bmm")), "\n",
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
  nt <- values[paste0("nt", seq_len(values[["ss"]] - 1))]
  r_side <- c(
    simple = .dmixture3p_simple(
      values[["y"]], values[["mu"]], values[["kappa"]], values[["thetat"]],
      values[["thetant"]], values[["ss"]], nt, 0, 81L
    ),
    simple_vp = .dmixture3p_simple(
      values[["y"]], values[["mu"]], values[["kappa"]], values[["thetat"]],
      values[["thetant"]], values[["ss"]], nt, values[["tau"]], 81L
    ),
    slot = .dmixture3p_slot(
      values[["y"]], values[["mu"]], values[["kappa"]], values[["K"]],
      values[["pnt"]], values[["ss"]], nt, 0, 81L
    ),
    slot_averaging = .dmixture3p_slot_averaging(
      values[["y"]], values[["mu"]], values[["kappa"]], values[["K"]],
      values[["pnt"]], values[["ss"]], nt, 0, 81L
    )
  )
  expect_equal(as.numeric(fit$draws("out")), unname(r_side), tolerance = 1e-10)
})
