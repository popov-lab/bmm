vm_density <- function(x, kappa) {
  exp(kappa * cos(x) - .circmix_log_besselI0(kappa)) / (2 * pi)
}

# one target component at mu = 0 with retrieval weight p_mem, plus guessing
psame_simple <- function(probe, kappa, p_mem, criterion = 0, p_s = p_mem) {
  .circmix_cd_psame(
    matrix(probe), matrix(log(p_mem), nrow = length(probe)),
    log1p(-p_mem), kappa, criterion, p_s
  )
}

test_that(".cd_gl_rule() integrates polynomials of degree up to 63 exactly", {
  gl <- .cd_gl_rule()
  expect_length(gl$x, 32)
  for (k in c(0, 1, 2, 7, 20, 40, 62)) {
    expect_equal(sum(gl$w * gl$x^k), (1 + (-1)^k) / (k + 1), tolerance = 1e-12)
  }
})

test_that(".circmix_cd_arc_mass() gives the full von Mises mass on the full circle", {
  d <- matrix(c(0, 1.3, -2.9), nrow = 1)
  for (kappa in c(0.5, 8, 60)) {
    expect_equal(
      as.vector(.circmix_cd_arc_mass(d, pi, kappa)), rep(1, 3),
      tolerance = 1e-10
    )
  }
})

test_that(".circmix_cd_crit_angle() puts the boundary where the likelihood ratio equals the criterion", {
  llr <- function(x, kappa, p_s) {
    -log(2 * pi) - log(p_s * vm_density(x, kappa) + (1 - p_s) / (2 * pi))
  }
  for (case in list(c(3, 0.6, 0), c(8, 0.8, 0.4), c(20, 0.3, -0.2), c(8, 0.5, 0.1))) {
    hw <- .circmix_cd_crit_angle(case[1], .circmix_cd_offset(case[3], case[2]))
    expect_gt(hw, 0)
    expect_lt(hw, pi)
    expect_equal(llr(hw, case[1], case[2]), case[3], tolerance = 1e-10)
  }
})

test_that("p_s drops out of the boundary at criterion 0 (Lin & Oberauer, Appendix B)", {
  hw <- vapply(c(0.05, 0.2, 0.55, 0.9, 1), function(p_s) {
    .circmix_cd_crit_angle(8, .circmix_cd_offset(0, p_s))
  }, numeric(1))
  expect_equal(hw, rep(acos(.circmix_log_besselI0(8) / 8), 5))
  expect_false(isTRUE(all.equal(
    .circmix_cd_crit_angle(8, .circmix_cd_offset(0.4, 0.2)),
    .circmix_cd_crit_angle(8, .circmix_cd_offset(0.4, 0.9))
  )))
})

test_that(".circmix_cd_offset() resolves the degenerate decisions", {
  # nothing is stored: the likelihood ratio is 1 and the criterion decides
  expect_equal(.circmix_cd_offset(c(-0.3, 0.3), 0), c(Inf, -Inf))
  # a criterion this liberal accepts every retrieved value as "same"
  expect_equal(.circmix_cd_offset(-log(0.5), 0.4), -Inf)
  expect_equal(.circmix_cd_crit_angle(8, c(Inf, -Inf)), c(0, pi))
  expect_equal(psame_simple(0.3, 8, 0.7, criterion = -0.5, p_s = 0), 0)
  expect_equal(psame_simple(0.3, 8, 0.7, criterion = 0.5, p_s = 0), 1)
})

# range of the log-likelihood ratio over a fine grid of retrieved values, from
# the decision density itself
llr_range <- function(kappa, w, p_s) {
  x <- seq(-pi, pi, length.out = 20001)
  density <- rowSums(vapply(seq_along(kappa), function(s) {
    w[s] * vm_density(x, kappa[s])
  }, numeric(length(x))))
  range(-log(1 - p_s + 2 * pi * p_s * density))
}

test_that(".circmix_cd_bounded_criterion() keeps 0 and stays inside the evidence range", {
  one <- function(criterion, p_s = 0.6, kappa = 8) {
    .circmix_cd_bounded_criterion(criterion, p_s, kappa, 0, FALSE)
  }
  expect_identical(one(0), 0)
  expect_equal(one(1e-7) / 1e-7, 1, tolerance = 1e-6)
  expect_equal(one(-1e-7) / -1e-7, 1, tolerance = 1e-6)
  expect_true(all(diff(one(seq(-5, 5, by = 0.1))) > 0))
  # nothing stored: the ratio is 1 everywhere, and the criterion decides alone
  expect_identical(one(c(-0.4, 0.4), p_s = 0), c(-0.4, 0.4))

  for (case in list(
    list(kappa = 8, w = 1, p_s = 0.6), list(kappa = 1, w = 1, p_s = 0.2),
    list(kappa = 40, w = 1, p_s = 1), list(kappa = c(6, 14), w = c(0.65, 0.35), p_s = 1)
  )) {
    marginal <- llr_range(case$kappa, case$w, case$p_s)
    # up to three times each end of the range; tanh rounds to 1 past about 18
    criteria <- c(c(3, 1, 0.2) * marginal[1], c(0.2, 1, 3) * marginal[2])
    n <- length(criteria)
    kappa <- matrix(case$kappa, n, length(case$kappa), byrow = TRUE)
    lw <- matrix(log(case$w), n, length(case$w), byrow = TRUE)
    bounded <- .circmix_cd_bounded_criterion(criteria, case$p_s, kappa, lw, FALSE)
    expect_true(all(bounded > marginal[1] & bounded < marginal[2]))
    expect_equal(bounded[6], marginal[2] * tanh(3), tolerance = 1e-6)
    # the rich observer has an arc per state, so the narrowest single-state range binds
    narrowest <- c(
      max(vapply(case$kappa, function(k) llr_range(k, 1, case$p_s)[1], numeric(1))),
      min(vapply(case$kappa, function(k) llr_range(k, 1, case$p_s)[2], numeric(1)))
    )
    bounded <- .circmix_cd_bounded_criterion(criteria, case$p_s, kappa, lw, TRUE)
    expect_true(all(bounded > narrowest[1] & bounded < narrowest[2]))
  }
})

test_that("a criterion beyond a trial's evidence leaves both responses possible", {
  # unbounded, criterion 0.5 lies past -log(1 - p_s) and the arc covers the circle
  p <- psame_simple(c(0.3, pi), 8, 0.7, criterion = 0.5, p_s = 0.1)
  expect_true(all(p > 1e-3 & p < 1 - 1e-3))
  p <- psame_simple(c(0.3, pi), 8, 0.7, criterion = -3, p_s = 0.6)
  expect_true(all(p > 1e-3 & p < 1 - 1e-3))
  # a held slot: p_s = 1, where -log(1 - p_s) is infinite but the range is not
  branch <- .circmix_slot_averaging_branches(c(4, 3.3), 2, c(6, 6))
  for (rich in c(TRUE, FALSE)) {
    p <- .circmix_cd_slot_averaging_psame(
      matrix(c(pi, 0), 2), matrix(0, 2), branch, c(10, -10), rich
    )
    expect_true(all(p > .Machine$double.eps & p < 1 - .Machine$double.eps))
  }
})

test_that(".circmix_cd_psame() matches a fine trapezoid rule over the arc", {
  for (case in list(c(0.6, 8, 0.7, 0.3), c(-2, 2, 0.4, 0), c(0.1, 40, 0.9, -0.4))) {
    probe <- case[1]
    kappa <- case[2]
    p_mem <- case[3]
    hw <- .circmix_cd_crit_angle(kappa, .circmix_cd_level(case[4], p_mem, kappa, 0, FALSE))
    x <- seq(probe - hw, probe + hw, length.out = 2e5 + 1)
    f <- p_mem * vm_density(x, kappa) + (1 - p_mem) / (2 * pi)
    trapezoid <- sum((f[-1] + f[-length(f)]) / 2) * (x[2] - x[1])
    expect_equal(psame_simple(probe, kappa, p_mem, case[4]), trapezoid, tolerance = 1e-8)
  }
})

test_that("an identical probe is called 'same' well above chance", {
  # the corrected decision rule: the same probe is called "same" well above 0.5
  expect_gt(psame_simple(0, 8, 0.8), 0.65)
  # and "change" grows with the distance of the probe from the target
  p <- psame_simple(seq(0, pi, length.out = 20), 8, 0.8)
  expect_true(all(diff(p) < 0))
})

test_that("non-targets raise P('same') for probes near them (intrusion cost)", {
  probe <- 1.2
  d <- matrix(c(probe, probe - 1.2), nrow = 1)
  with_lure <- .circmix_cd_psame(d, matrix(log(c(0.6, 0.3)), 1), log(0.1), 8, 0, 0.6)
  without_lure <- .circmix_cd_psame(
    d[, 1, drop = FALSE], matrix(log(0.6), 1), log(0.4), 8, 0, 0.6
  )
  expect_gt(with_lure, without_lure + 0.2)
})

# two precision states, as slot averaging gives a held item
two_states <- function(kappa_lo = 6, kappa_hi = 14, extra = 0.35) {
  list(kappa = matrix(c(kappa_lo, kappa_hi), 1), lw = matrix(log(c(1 - extra, extra)), 1))
}

test_that(".circmix_cd_mix_psame() matches adaptive integration over the arc", {
  d <- matrix(c(0.5, 1.4), nrow = 1)
  logw <- matrix(log(c(0.6, 0.3)), 1)
  states <- two_states()
  w <- exp(states$lw)
  arc <- function(kappa, hw) {
    x_mass <- function(loc) {
      stats::integrate(function(x) vm_density(x - loc, kappa), -hw, hw,
        rel.tol = 1e-12
      )$value
    }
    0.1 * hw / pi + sum(c(0.6, 0.3) * vapply(-d, x_mass, numeric(1)))
  }
  for (criterion in c(0, 0.4, -0.3)) {
    offset <- .circmix_cd_level(criterion, 0.6, states$kappa, states$lw, TRUE)
    rich <- sum(w * vapply(1:2, function(s) {
      arc(states$kappa[s], .circmix_cd_crit_angle(states$kappa[s], offset))
    }, numeric(1)))
    offset <- .circmix_cd_level(criterion, 0.6, states$kappa, states$lw, FALSE)
    marginal <- function(h) {
      sum(w * exp(states$kappa * cos(h) - .circmix_log_besselI0(states$kappa)))
    }
    hw <- stats::uniroot(function(h) log(marginal(h)) - offset, c(1e-9, pi - 1e-9),
      tol = 1e-14
    )$root
    limited <- sum(w * vapply(states$kappa, function(k) arc(k, hw), numeric(1)))
    for (rich_obs in c(TRUE, FALSE)) {
      expect_equal(
        .circmix_cd_mix_psame(d, logw, log(0.1), states$kappa, states$lw,
          criterion, 0.6, rich_obs
        ),
        if (rich_obs) rich else limited,
        tolerance = 1e-9
      )
    }
  }
})

test_that("both observers coincide when there is a single precision state", {
  d <- matrix(c(0.5, 1.4), nrow = 1)
  logw <- matrix(log(c(0.6, 0.3)), 1)
  for (criterion in c(0, 0.4, -0.3)) {
    simple <- .circmix_cd_psame(d, logw, log(0.1), 10, criterion, 0.6)
    for (rich in c(TRUE, FALSE)) {
      expect_equal(
        .circmix_cd_mix_psame(d, logw, log(0.1), matrix(10), matrix(0), criterion, 0.6, rich),
        simple
      )
    }
  }
  # K a multiple of the set size: every item holds exactly two slots
  branch <- .circmix_slot_averaging_branches(4, 2, 6)
  expect_equal(branch$extra, 0)
  for (rich in c(TRUE, FALSE)) {
    expect_equal(
      .circmix_cd_slot_averaging_psame(d, logw, branch, 0.3, rich),
      .circmix_cd_psame(d, logw, -Inf, branch$kappa_lo, 0.3, 1),
      tolerance = 1e-10
    )
  }
})

test_that("the knowledge-limited boundary lies between the per-state boundaries", {
  for (extra in c(0.1, 0.5, 0.9)) {
    states <- two_states(extra = extra)
    for (criterion in c(0, 0.4, -0.3)) {
      offset <- .circmix_cd_offset(criterion, 0.6)
      per_state <- .circmix_cd_crit_angle(states$kappa, offset)
      start <- .circmix_cd_crit_angle(
        exp(matrixStats::rowLogSumExps(states$lw + log(states$kappa))), offset
      )
      shared <- .circmix_cd_limited_hw(states$kappa, states$lw, offset, start)
      expect_gte(shared, min(per_state))
      expect_lte(shared, max(per_state))
    }
  }
})

test_that("P('same') falls with the distance of the probe under both observers", {
  probes <- seq(0, pi, length.out = 15)
  states <- two_states()
  for (rich in c(TRUE, FALSE)) {
    p <- .circmix_cd_mix_psame(
      matrix(probes), matrix(log(0.8), nrow = 15), rep(log(0.2), 15),
      states$kappa[rep(1, 15), ], states$lw[rep(1, 15), ], 0.2, 0.8, rich
    )
    expect_true(all(diff(p) < 0))
  }
})

test_that("p_s drops out at criterion 0 under both observers", {
  d <- matrix(c(0.5, 1.4), nrow = 1)
  logw <- matrix(log(c(0.6, 0.3)), 1)
  states <- two_states()
  for (rich in c(TRUE, FALSE)) {
    expect_equal(
      .circmix_cd_mix_psame(d, logw, log(0.1), states$kappa, states$lw, 0, 0.2, rich),
      .circmix_cd_mix_psame(d, logw, log(0.1), states$kappa, states$lw, 0, 0.9, rich)
    )
  }
})

test_that(".circmix_cd_slot_averaging_psame() judges an item without a slot at kappa_hi", {
  d <- matrix(0.4, nrow = 1)
  logw <- matrix(0, nrow = 1)
  branch <- .circmix_slot_averaging_branches(1.5, 4, 6)
  expect_false(branch$held)
  extra <- branch$extra
  hw <- .circmix_cd_crit_angle(branch$kappa_hi, .circmix_cd_offset(0, extra))
  expected <- (1 - extra) * hw / pi + extra *
    .circmix_cd_arc_mass(d, hw, branch$kappa_hi)[1, 1]
  for (rich in c(TRUE, FALSE)) {
    expect_equal(
      .circmix_cd_slot_averaging_psame(d, logw, branch, 0, rich),
      expected,
      tolerance = 1e-12
    )
  }
})

test_that(".cd_bernoulli_ld() scores 1 as 'change' and clamps the probability", {
  expect_equal(.cd_bernoulli_ld(c(1, 0), c(0.3, 0.3)), log(c(0.7, 0.3)))
  expect_true(all(is.finite(.cd_bernoulli_ld(c(1, 0), c(1, 0)))))
  expect_equal(.cd_bernoulli_ld(1, c(0.2, 0.4, 0.6)), log(c(0.8, 0.6, 0.4)))
})

test_that(".cd_add_criterion() adds a fixed, identity-linked criterion", {
  spec <- .cd_add_criterion(list(
    parameters = list(mu = "m", kappa = "k"), links = list(mu = "tan_half", kappa = "log"),
    fixed_parameters = list(mu = 0), priors = list(), init_ranges = list()
  ))
  expect_equal(names(spec$parameters), c("mu", "kappa", "criterion"))
  expect_equal(spec$links$criterion, "identity")
  expect_equal(spec$fixed_parameters, list(mu = 0, criterion = 0))
  expect_named(spec$priors$criterion, c("main", "effects", "sd"))
  expect_equal(spec$init_ranges$criterion, c(-0.2, 0.2))
})

test_that(".cd_check_binary_response() accepts 0/1 and logical responses only", {
  dat <- data.frame(r = c(0, 1, 1))
  expect_identical(.cd_check_binary_response(dat, "r"), dat)
  expect_warning(
    out <- .cd_check_binary_response(data.frame(r = c(TRUE, FALSE)), "r"),
    "logical"
  )
  expect_equal(out$r, c(1L, 0L))
  expect_error(.cd_check_binary_response(data.frame(r = c(0, 2)), "r"), "0 \\('same'\\)")
  expect_error(.cd_check_binary_response(data.frame(r = c("a", "b")), "r"), "0 \\('same'\\)")
  expect_error(.cd_check_binary_response(data.frame(r = c(0, NA)), "r"), "missing")
})

cd_stub_model <- function() {
  structure(
    list(resp_vars = list(response = "r", probe = "p", target = "t"), other_vars = list()),
    class = c("bmmodel", "change_detection", "cd_stub")
  )
}

test_that("check_data.change_detection() centres the probe on the target", {
  dat <- data.frame(r = c(0, 1, 1), p = c(0.5, 3, -3), t = c(0.2, -3, 3))
  out <- check_data(cd_stub_model(), dat, bmf(kappa ~ 1))
  expect_equal(out$probe_centered, wrap(dat$p - dat$t))
  expect_error(check_data(cd_stub_model(), dat[c("r", "p")], bmf(kappa ~ 1)), "'t'")
  expect_warning(
    check_data(cd_stub_model(), transform(dat, p = p * 100), bmf(kappa ~ 1)),
    "degrees"
  )
})

test_that("check_data.change_detection() refuses a missing probe or target", {
  dat <- data.frame(r = c(0, 1, 1), p = c(0.5, 3, -3), t = c(0.2, -3, 3))
  for (column in c("p", "t")) {
    holed <- dat
    holed[[column]][2] <- NA
    expect_error(
      check_data(cd_stub_model(), holed, bmf(kappa ~ 1)),
      glue::glue("'{column}' contains missing values")
    )
  }
  expect_error(
    check_data(cd_stub_model(), transform(dat, r = c(0, NA, 1)), bmf(kappa ~ 1)),
    "missing"
  )
})

test_that("check_data.change_detection() declares the column it builds", {
  model <- cd_stub_model()
  dat <- data.frame(r = c(0, 1, 1), p = c(0.5, 3, -3), t = c(0.2, -3, 3))
  built <- setdiff(names(check_data(model, dat, bmf(kappa ~ 1))), names(dat))
  expect_identical(built_data_columns(model), built)
})

test_that("revert_check_data.change_detection() rebuilds what check_data() needs", {
  model <- cd_stub_model()
  dat <- data.frame(r = c(0, 1, 1), p = c(0.5, 3, -3), t = c(0.2, -3, 3))
  checked <- check_data(model, dat, bmf(kappa ~ 1))
  for (keep in list(c("r", "probe_centered"), c("r", "probe_centered", "t"),
                    c("r", "probe_centered", "p"))) {
    stored <- revert_check_data(model, as.data.frame(checked)[keep])
    expect_setequal(attr(stored, "rebuilt"), setdiff(c("p", "t"), keep))
    expect_equal(
      check_data(model, stored, bmf(kappa ~ 1))$probe_centered,
      checked$probe_centered,
      tolerance = 1e-14
    )
  }
})

test_that("change-detection models print their response units", {
  expect_named(response_annotations(cd_stub_model()), c("response", "probe", "target"))
  expect_equal(
    model_group("Visual working memory (change detection)"),
    "Detection, recognition and confidence judgments"
  )
})

test_that(".cd_custom_family() appends criterion and declares an integer response", {
  model <- list(
    links = list(mu = "tan_half", kappa = "log", thetat = "logit", criterion = "identity"),
    variable_precision = FALSE, knowledge = "rich"
  )
  stub_lik <- function(i, prep) NULL
  stub_pred <- function(i, prep, ...) NULL
  family <- .cd_custom_family(
    model, "demo_cd", "thetat",
    log_lik = stub_lik, posterior_predict = stub_pred,
    posterior_epred = posterior_epred_undefined("demo")
  )
  expect_equal(family$dpars, c("mu", "kappa", "thetat", "criterion"))
  expect_equal(family$type, "int")
  expect_equal(family$knowledge, "rich")
  expect_equal(
    family$vars,
    c("vreal1[n]", .circmix_table_vars(), "cd_gl_x", "cd_gl_w")
  )
})

test_that(".cd_prep_nt() drops the probe column", {
  prep <- list(data = list(vreal3 = c(3, 0), vreal1 = c(1, 0), vreal2 = c(2, 0)))
  expect_equal(unname(.cd_prep_nt(prep, 1)), c(2, 3))
  expect_equal(.cd_prep_probe(prep, 1), 1)
})

test_that(".cd_data() passes the kappa(J) table in the order of its vars", {
  expect_identical(head(names(.cd_data()), 4), .circmix_table_vars())
})

test_that(".cd_stan_wrapper() emits an lpmf with a scalar probe and the observer", {
  wrapper <- .cd_stan_wrapper(
    "demo_cd", c("mu", "kappa", "thetat", "criterion"),
    vint = "ss", vreal = list(probe = 1, nt = 2), rich = TRUE
  )
  expect_match(
    wrapper,
    paste0(
      "real demo_cd_lpmf\\(int y, real mu, real kappa, real thetat, real criterion, ",
      "int ss, real probe, real nt1, real nt2, data vector circmix_logk"
    )
  )
  expect_match(wrapper, "data vector cd_gl_x, data vector cd_gl_w\\) \\{")
  expect_match(
    wrapper,
    paste0(
      "demo_cd_core\\(y, mu, kappa, thetat, criterion, ss, probe, ",
      "to_vector\\(\\{nt1, nt2\\}\\), 1, circmix_logk"
    )
  )
  # set size 1 only: no non-target columns, an empty vector for the core
  wrapper <- .cd_stan_wrapper("demo_cd", c("mu", "kappa"), vreal = list(probe = 1, nt = 0))
  expect_match(wrapper, "real probe, data vector circmix_logk")
  expect_match(wrapper, "probe, rep_vector\\(0, 0\\), 0, circmix_logk")
})

cd_chunks <- function() {
  paste(vapply(
    c("circmix_funs.stan", "cd_funs.stan", "circmix_cd_funs.stan"),
    function(f) read_lines2(system.file("stan_chunks", f, package = "bmm")),
    character(1)
  ), collapse = "\n")
}

test_that("a change-detection family built from the helpers reaches brms and parses", {
  skip_unless_stan_tests()
  skip_on_cran()
  model <- list(
    links = list(mu = "tan_half", kappa = "log", thetat = "logit", criterion = "identity"),
    variable_precision = FALSE, knowledge = "limited"
  )
  family <- .cd_custom_family(
    model, "demo_cd", "thetat",
    vint = TRUE, n_vreal = 2,
    log_lik = function(i, prep) NULL,
    posterior_predict = function(i, prep, ...) NULL,
    posterior_epred = posterior_epred_undefined("demo")
  )
  core <- brms::stanvar(scode = "
  real demo_cd_core(int y, real mu, real kappa, real thetat, real criterion,
                    int ss, real probe, vector nt, int rich, data vector logk,
                    data vector dlogk, data real logJ_min, data real dlogJ,
                    data vector gl_x, data vector gl_w) {
    vector[2] w = [log(thetat), log1m(thetat)]';
    real p = circmix_cd_psame([probe - mu]', head(w, 1), w[2], kappa, criterion,
                              thetat, gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }", block = "functions", name = "demo_cd_core")
  stanvars <- .cd_model_stanvars(
    model, family, character(0),
    vint = "ss", vreal = list(probe = 1, nt = 1)
  ) + core
  dat <- data.frame(r = c(0L, 1L, 1L), probe_centered = c(0.1, 2, -1), nt1 = 1, ss = 2L)
  formula <- brms::bf(r | vint(ss) + vreal(probe_centered, nt1) ~ 1, kappa ~ 1,
    thetat ~ 1, criterion ~ 1,
    family = family
  )
  code <- brms::stancode(formula, dat,
    stanvars = stanvars,
    prior = brms::prior(constant(0), class = "Intercept", dpar = "criterion")
  )
  expect_match(code, "demo_cd_lpmf\\(Y\\[n\\] \\| mu\\[n\\], kappa\\[n\\], thetat\\[n\\], criterion\\[n\\], vint1\\[n\\], vreal1\\[n\\], vreal2\\[n\\]")
  expect_match(code, "Intercept_criterion = 0")
  standata <- brms::standata(formula, dat, stanvars = stanvars)
  expect_length(standata$cd_gl_x, 32)
  expect_equal(as.numeric(standata$vreal1), dat$probe_centered)
  stan_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
  expect_silent(stan_model$check_syntax(quiet = TRUE))
})

# Evaluates the change-detection functions through a fixed_param run, as
# circmix_stan_values() in test-helpers-circmix.R does
cd_stan_values <- function(data) {
  code <- paste0(
    "functions {\n", cd_chunks(), "\n}\n",
    "data {
      int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ;
      vector[32] gl_x; vector[32] gl_w;
      int NC; array[NC] real crit; array[NC] real ps;
      int ND; vector[ND] d; vector[ND] logw; real logw_guess;
      int NS; vector[NS] kappa_s; vector[NS] lw_s;
      int NB; array[NB] real K; array[NB] int ss;
    }
    generated quantities {
      vector[NC] out_offset; vector[NC] out_crit; vector[ND] out_mass;
      vector[NC] out_psame; vector[NC] out_limited; array[2] vector[NC] out_mix;
      array[2, NB] vector[NC] out_sa; array[2] vector[NC] out_bound;
      array[2] real out_bern;
      for (c in 1:NC) {
        out_offset[c] = circmix_cd_offset(crit[c], ps[c]);
        out_crit[c] = circmix_cd_crit_angle(7, out_offset[c]);
        out_psame[c] = circmix_cd_psame(d, logw, logw_guess, 7, crit[c], ps[c], gl_x, gl_w);
        out_limited[c] = circmix_cd_limited_hw(kappa_s, lw_s, out_offset[c], 1.0);
        for (r in 1:2) {
          out_bound[r, c] = circmix_cd_bounded_criterion(crit[c], ps[c], kappa_s, lw_s, r - 1);
          out_mix[r, c] = circmix_cd_mix_psame(d, logw, logw_guess, kappa_s, lw_s,
                                               crit[c], ps[c], r - 1, gl_x, gl_w);
          for (b in 1:NB) {
            vector[4] branch = circmix_slot_averaging_branches(K[b], ss[b], 7, logk,
                                                               dlogk, logJ_min, dlogJ);
            out_sa[r, b, c] = circmix_cd_slot_averaging_psame(
              d, logw - log_sum_exp(logw), branch, crit[c], r - 1, gl_x, gl_w);
          }
        }
      }
      out_mass = circmix_cd_arc_mass(d, 0.8, 7, gl_x, gl_w);
      out_bern[1] = cd_bernoulli_lpmf(1 | 0.3);
      out_bern[2] = cd_bernoulli_lpmf(0 | 0.3);
    }"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = data, fixed_param = TRUE, iter_sampling = 1, chains = 1,
    refresh = 0, show_messages = FALSE, sig_figs = 18
  )
  function(variable) as.numeric(fit$draws(variable))
}

test_that("the Stan and R change-detection functions agree", {
  skip_unless_stan_tests()
  skip_on_cran()
  tab <- .circmix_kappa_table()
  gl <- .cd_gl_rule()
  # past -log(1 - p_s), nothing stored, a held slot, the liberal side
  crit <- c(0, 0.4, -0.3, -0.5, 1.5, 3, -2)
  ps <- c(0.6, 0.6, 0.6, 0, 0.3, 1, 0.9)
  d <- c(0.4, 0.4 - 1.2, 0.4 + 2)
  logw <- log(c(0.5, 0.2, 0.15))
  logw_guess <- log(0.15)
  kappa_s <- c(3, 9, 20)
  lw_s <- log(c(0.2, 0.5, 0.3))
  # held with a fractional slot, no slot held, a whole K / set size held and not
  K <- c(3.5, 1.5, 4, 0)
  ss <- c(2L, 4L, 2L, 3L)

  stan <- cd_stan_values(list(
    NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
    logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, gl_x = gl$x, gl_w = gl$w,
    NC = length(crit), crit = crit, ps = ps, ND = length(d), d = d, logw = logw,
    logw_guess = logw_guess, NS = length(kappa_s), kappa_s = kappa_s, lw_s = lw_s,
    NB = length(K), K = K, ss = ss
  ))

  n <- length(crit)
  rows <- function(x) matrix(x, nrow = n, ncol = length(x), byrow = TRUE)
  offset <- .circmix_cd_offset(crit, ps)
  expect_equal(stan("out_offset"), offset, tolerance = 1e-12)
  bound <- vapply(c(FALSE, TRUE), function(rich) {
    .circmix_cd_bounded_criterion(crit, ps, rows(kappa_s), rows(lw_s), rich)
  }, numeric(n))
  expect_equal(stan("out_bound"), as.vector(t(bound)), tolerance = 1e-12)
  expect_equal(stan("out_crit"), .circmix_cd_crit_angle(7, offset), tolerance = 1e-12)
  expect_equal(
    stan("out_mass"),
    as.vector(.circmix_cd_arc_mass(matrix(d, 1), 0.8, 7)),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_psame"),
    .circmix_cd_psame(rows(d), rows(logw), rep(logw_guess, n), 7, crit, ps),
    tolerance = 1e-12
  )
  limited <- .circmix_cd_limited_hw(rows(kappa_s), rows(lw_s), offset, rep(1, n))
  expect_equal(stan("out_limited"), limited, tolerance = 1e-12)

  # cmdstanr flattens arrays with the first index fastest, as R does
  mix <- vapply(c(FALSE, TRUE), function(rich) {
    .circmix_cd_mix_psame(rows(d), rows(logw), rep(logw_guess, n), rows(kappa_s),
      rows(lw_s), crit, ps, rich
    )
  }, numeric(n))
  expect_equal(stan("out_mix"), as.vector(t(mix)), tolerance = 1e-12)

  sa <- array(NA_real_, c(2, length(K), n))
  for (r in 1:2) {
    for (b in seq_along(K)) {
      branch <- .circmix_slot_averaging_branches(rep(K[b], n), rep(ss[b], n), rep(7, n))
      sa[r, b, ] <- .circmix_cd_slot_averaging_psame(
        rows(d), rows(logw - matrixStats::logSumExp(logw)), branch, crit, r == 2
      )
    }
  }
  expect_equal(stan("out_sa"), as.vector(sa), tolerance = 1e-12)
  expect_equal(stan("out_bern"), log(c(0.7, 0.3)), tolerance = 1e-12)
})

test_that("an infinite level gives kappa a finite gradient in circmix_cd_crit_angle()", {
  skip_unless_stan_tests()
  skip_on_cran()
  # the division by kappa must not run for an infinite level: its partial in
  # kappa is infinite, and reverse mode multiplies it by the unused node's zero
  # adjoint, which is NaN although the function returns 0 or pi
  code <- paste0(
    "functions {\n", cd_chunks(), "\n}\n",
    "data { real level; }
    parameters { real log_kappa; }
    model { target += circmix_cd_crit_angle(exp(log_kappa), level); }"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  for (level in c(-Inf, Inf)) {
    gradient <- tryCatch(
      model$diagnose(
        data = list(level = level), init = list(list(log_kappa = log(8))),
        error = 1e300
      )$gradients()$model,
      error = function(e) NA_real_
    )
    expect_identical(gradient, 0, label = paste("gradient at level", level))
  }
})

test_that("a criterion past the cliff keeps a finite gradient that points back", {
  skip_unless_stan_tests()
  skip_on_cran()
  # unbounded, criterion 2 at p_s 0.3 lies past -log(1 - p_s) = 0.36: the arc
  # covers the circle, a "change" response sits on the probability floor
  # (log 2.2e-16 = -36) and the criterion has no gradient to leave it
  code <- paste0(
    "functions {\n", cd_chunks(), "\n}\n",
    "data { int y; real probe; real p_s; vector[32] gl_x; vector[32] gl_w; }
    parameters { real log_kappa; real criterion; }
    model {
      target += cd_bernoulli_lpmf(y | circmix_cd_psame([probe]', [log(p_s)]', log1m(p_s),
                                  exp(log_kappa), criterion, p_s, gl_x, gl_w));
    }"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  gl <- .cd_gl_rule()
  diagnosed <- model$diagnose(
    data = list(y = 1L, probe = pi, p_s = 0.3, gl_x = gl$x, gl_w = gl$w),
    init = list(list(log_kappa = log(8), criterion = 2)), error = 1e300
  )
  expect_gt(diagnosed$lp(), log(1e-4))
  gradient <- diagnosed$gradients()
  expect_true(all(is.finite(gradient$model)))
  # parameters in declaration order: log_kappa, criterion
  expect_lt(gradient$model[2], 0)
})

test_that("an estimated criterion is announced as bounded", {
  model <- structure(
    list(fixed_parameters = list(mu = 0, criterion = 0)),
    class = c("change_detection", "cd_stub")
  )
  expect_silent(check_formula(model, NULL, bmf(kappa ~ 1)))
  model$fixed_parameters$criterion <- NULL
  expect_message(check_formula(model, NULL, bmf(kappa ~ 1, criterion ~ 1)), "bounds it")
})
