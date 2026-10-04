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
  expect_equal(psame_simple(0.3, 8, 0.7, criterion = 2, p_s = 0.1), 1)
  expect_equal(psame_simple(0.3, 8, 0.7, criterion = -0.5, p_s = 0), 0)
})

test_that(".circmix_cd_psame() matches a fine trapezoid rule over the arc", {
  for (case in list(c(0.6, 8, 0.7, 0.3), c(-2, 2, 0.4, 0), c(0.1, 40, 0.9, -0.4))) {
    probe <- case[1]
    kappa <- case[2]
    p_mem <- case[3]
    hw <- .circmix_cd_crit_angle(kappa, .circmix_cd_offset(case[4], p_mem))
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

vp_states_reference <- function(kappa, tau) {
  shape <- .circmix_J(kappa) / tau
  list(
    density = function(J) stats::dgamma(J, shape = shape, scale = tau),
    lower = stats::qgamma(1e-12, shape = shape, scale = tau),
    upper = stats::qgamma(1 - 1e-12, shape = shape, scale = tau)
  )
}

test_that(".circmix_cd_vp_psame() matches adaptive integration over J", {
  d <- matrix(c(0.5, 1.4), nrow = 1)
  logw <- matrix(log(c(0.6, 0.3)), 1)
  criterion <- 0.2
  p_s <- 0.6
  offset <- .circmix_cd_offset(criterion, p_s)
  arc <- function(kappa, hw) {
    .circmix_cd_psame_at(d, logw, log(0.1), kappa, hw)
  }
  for (shape in c(4, 2, 1)) {
    tau <- .circmix_J(10) / shape
    ref <- vp_states_reference(10, tau)
    integrate_J <- function(f) {
      stats::integrate(Vectorize(function(J) f(J) * ref$density(J)),
        ref$lower, ref$upper,
        rel.tol = 1e-10, subdivisions = 2000L
      )$value
    }
    rich <- integrate_J(function(J) {
      kappa <- .circmix_kappa(J)
      arc(kappa, .circmix_cd_crit_angle(kappa, offset))
    })
    marginal <- function(h) {
      integrate_J(function(J) {
        kappa <- .circmix_kappa(J)
        exp(kappa * cos(h) - .circmix_log_besselI0(kappa))
      })
    }
    hw <- stats::uniroot(function(h) log(marginal(h)) - offset, c(1e-6, pi - 1e-6),
      tol = 1e-12
    )$root
    limited <- integrate_J(function(J) arc(.circmix_kappa(J), hw))

    expect_equal(.circmix_cd_vp_psame(d, logw, log(0.1), 10, tau, criterion, p_s, TRUE),
      rich,
      tolerance = 1e-4
    )
    expect_equal(.circmix_cd_vp_psame(d, logw, log(0.1), 10, tau, criterion, p_s, FALSE),
      limited,
      tolerance = 1e-4
    )
  }
})

test_that("both observers coincide when there is a single precision state", {
  d <- matrix(c(0.5, 1.4), nrow = 1)
  logw <- matrix(log(c(0.6, 0.3)), 1)
  for (criterion in c(0, 0.4, -0.3)) {
    simple <- .circmix_cd_psame(d, logw, log(0.1), 10, criterion, 0.6)
    expect_equal(.circmix_cd_vp_psame(d, logw, log(0.1), 10, 0, criterion, 0.6, TRUE), simple)
    expect_equal(.circmix_cd_vp_psame(d, logw, log(0.1), 10, 0, criterion, 0.6, FALSE), simple)
  }
  # K a multiple of the set size: every item holds exactly two slots
  branch <- .circmix_slot_averaging_branches(4, 2, 6)
  for (rich in c(TRUE, FALSE)) {
    expect_equal(
      .circmix_cd_slot_averaging_psame(d, logw, branch, 0, 0.3, rich),
      .circmix_cd_psame(d, logw, -Inf, branch$kappa_lo, 0.3, 1),
      tolerance = 1e-10
    )
  }
})

test_that("the knowledge-limited boundary lies between the per-state boundaries", {
  for (tau in c(0.3, 1, 3)) {
    states <- .circmix_cd_vp_states(10, tau, 41L)
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
  for (rich in c(TRUE, FALSE)) {
    p <- .circmix_cd_vp_psame(
      matrix(probes), matrix(log(0.8), nrow = 15), rep(log(0.2), 15),
      rep(8, 15), rep(1, 15), 0.2, 0.8, rich
    )
    expect_true(all(diff(p) < 0))
  }
})

test_that("p_s drops out at criterion 0 under both observers", {
  d <- matrix(c(0.5, 1.4), nrow = 1)
  logw <- matrix(log(c(0.6, 0.3)), 1)
  for (rich in c(TRUE, FALSE)) {
    expect_equal(
      .circmix_cd_vp_psame(d, logw, log(0.1), 10, 1, 0, 0.2, rich),
      .circmix_cd_vp_psame(d, logw, log(0.1), 10, 1, 0, 0.9, rich)
    )
  }
})

test_that(".circmix_cd_slot_averaging_psame() judges an item without a slot at kappa_hi", {
  d <- matrix(0.4, nrow = 1)
  logw <- matrix(0, nrow = 1)
  branch <- .circmix_slot_averaging_branches(1.5, 4, 6)
  expect_false(branch$held)
  extra <- exp(branch$log_w_hi)
  hw <- .circmix_cd_crit_angle(branch$kappa_hi, .circmix_cd_offset(0.3, extra))
  expected <- (1 - extra) * hw / pi + extra *
    .circmix_cd_arc_mass(d, hw, branch$kappa_hi)[1, 1]
  for (rich in c(TRUE, FALSE)) {
    expect_equal(
      .circmix_cd_slot_averaging_psame(d, logw, branch, 0, 0.3, rich),
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
    links = list(mu = "tan_half", kappa = "log", tau = "log", thetat = "logit", criterion = "identity"),
    vp_nodes = 41L, variable_precision = TRUE, knowledge = "rich"
  )
  stub_lik <- function(i, prep) NULL
  stub_pred <- function(i, prep, ...) NULL
  family <- .cd_custom_family(
    model, "demo_cd", "thetat",
    log_lik = stub_lik, posterior_predict = stub_pred,
    posterior_epred = posterior_epred_undefined("demo")
  )
  expect_equal(family$dpars, c("mu", "kappa", "tau", "thetat", "criterion"))
  expect_equal(family$core_dpars, family$dpars)
  expect_equal(family$type, "int")
  expect_equal(family$knowledge, "rich")
  expect_equal(
    family$vars,
    c("vreal1[n]", .circmix_table_vars(), "cd_gl_x", "cd_gl_w")
  )
})

test_that(".circmix_prep_nt() can skip the probe column", {
  prep <- list(data = list(vreal3 = c(3, 0), vreal1 = c(1, 0), vreal2 = c(2, 0)))
  expect_equal(unname(.circmix_prep_nt(prep, 1, skip = 1)), c(2, 3))
  expect_equal(.cd_prep_probe(prep, 1), 1)
})

test_that(".circmix_stan_wrapper() emits an lpmf with a scalar probe and literals", {
  wrapper <- .circmix_stan_wrapper(
    "demo_cd", c("mu", "kappa", "thetat", "criterion"),
    core_dpars = c("mu", "kappa", "0.0", "thetat", "criterion"),
    vint = "ss", vreal = list(probe = 1, nt = 2), type = "int",
    vreal_scalars = "probe", literals = "0",
    data = c(.circmix_table_data(), .cd_gl_data())
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
      "demo_cd_core\\(y, mu, kappa, 0.0, thetat, criterion, ss, probe, ",
      "to_vector\\(\\{nt1, nt2\\}\\), 41, 0, circmix_logk"
    )
  )
})

cd_chunks <- function() {
  paste(vapply(
    c("circmix_funs.stan", "cd_funs.stan", "circmix_cd_funs.stan"),
    function(f) read_lines2(system.file("stan_chunks", f, package = "bmm")),
    character(1)
  ), collapse = "\n")
}

test_that("a change-detection family built from the helpers reaches brms and parses", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to check the generated Stan"
  )
  model <- list(
    links = list(mu = "tan_half", kappa = "log", thetat = "logit", criterion = "identity"),
    vp_nodes = 41L, variable_precision = FALSE, knowledge = "limited"
  )
  family <- .cd_custom_family(
    model, "demo_cd", "thetat",
    vint = TRUE, n_vreal = 2,
    log_lik = function(i, prep) NULL,
    posterior_predict = function(i, prep, ...) NULL,
    posterior_epred = posterior_epred_undefined("demo")
  )
  core <- brms::stanvar(scode = "
  real demo_cd_core(int y, real mu, real kappa, real tau, real thetat,
                    real criterion, int ss, real probe, vector nt, int nodes,
                    int rich, data vector logk, data vector dlogk,
                    data real logJ_min, data real dlogJ, data vector gl_x,
                    data vector gl_w) {
    vector[2] w = [log(thetat), log1m(thetat)]';
    real p = circmix_cd_vp_psame([probe - mu]', head(w, 1), w[2], kappa, tau,
                                 criterion, thetat, rich, nodes, logk, dlogk,
                                 logJ_min, dlogJ, gl_x, gl_w);
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
      int NT; array[NT] real tau;
      int NB; array[NB] real K; array[NB] int ss;
    }
    generated quantities {
      vector[NC] out_offset; vector[NC] out_crit; vector[ND] out_mass;
      vector[NC] out_psame; vector[NC] out_limited; array[2] vector[NC] out_mix;
      array[2, NT] vector[NC] out_vp; array[2, NT, NB] vector[NC] out_sa;
      array[2] real out_bern;
      for (c in 1:NC) {
        out_offset[c] = circmix_cd_offset(crit[c], ps[c]);
        out_crit[c] = circmix_cd_crit_angle(7, out_offset[c]);
        out_psame[c] = circmix_cd_psame(d, logw, logw_guess, 7, crit[c], ps[c], gl_x, gl_w);
        out_limited[c] = circmix_cd_limited_hw(kappa_s, lw_s, out_offset[c], 1.0);
        for (r in 1:2) {
          out_mix[r, c] = circmix_cd_mix_psame(d, logw, logw_guess, kappa_s, lw_s,
                                               crit[c], ps[c], r - 1, gl_x, gl_w);
          for (t in 1:NT) {
            out_vp[r, t, c] = circmix_cd_vp_psame(d, logw, logw_guess, 7, tau[t],
                                                  crit[c], ps[c], r - 1, 41, logk,
                                                  dlogk, logJ_min, dlogJ, gl_x, gl_w);
            for (b in 1:NB) {
              vector[5] branch = circmix_slot_averaging_branches(K[b], ss[b], 7, logk,
                                                                 dlogk, logJ_min, dlogJ);
              out_sa[r, t, b, c] = circmix_cd_slot_averaging_psame(
                d, logw - log_sum_exp(logw), branch, tau[t], crit[c], r - 1, 41,
                logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
            }
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
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the change-detection functions"
  )
  tab <- .circmix_kappa_table()
  gl <- .cd_gl_rule()
  crit <- c(0, 0.4, -0.3, -0.5, 1.5)
  ps <- c(0.6, 0.6, 0.6, 0, 0.3)
  d <- c(0.4, 0.4 - 1.2, 0.4 + 2)
  logw <- log(c(0.5, 0.2, 0.15))
  logw_guess <- log(0.15)
  kappa_s <- c(3, 9, 20)
  lw_s <- log(c(0.2, 0.5, 0.3))
  tau <- c(0, 0.5, 2)
  K <- c(3.5, 1.5, 4)
  ss <- c(2L, 4L, 2L)

  stan <- cd_stan_values(list(
    NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
    logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, gl_x = gl$x, gl_w = gl$w,
    NC = length(crit), crit = crit, ps = ps, ND = length(d), d = d, logw = logw,
    logw_guess = logw_guess, NS = length(kappa_s), kappa_s = kappa_s, lw_s = lw_s,
    NT = length(tau), tau = tau, NB = length(K), K = K, ss = ss
  ))

  n <- length(crit)
  rows <- function(x) matrix(x, nrow = n, ncol = length(x), byrow = TRUE)
  offset <- .circmix_cd_offset(crit, ps)
  expect_equal(stan("out_offset"), offset, tolerance = 1e-12)
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

  vp <- array(NA_real_, c(2, length(tau), n))
  sa <- array(NA_real_, c(2, length(tau), length(K), n))
  for (r in 1:2) {
    for (t in seq_along(tau)) {
      vp[r, t, ] <- .circmix_cd_vp_psame(rows(d), rows(logw), rep(logw_guess, n),
        rep(7, n), rep(tau[t], n), crit, ps, r == 2
      )
      for (b in seq_along(K)) {
        branch <- .circmix_slot_averaging_branches(rep(K[b], n), rep(ss[b], n), rep(7, n))
        sa[r, t, b, ] <- .circmix_cd_slot_averaging_psame(
          rows(d), rows(logw - matrixStats::logSumExp(logw)), branch,
          rep(tau[t], n), crit, r == 2
        )
      }
    }
  }
  expect_equal(stan("out_vp"), as.vector(vp), tolerance = 1e-12)
  expect_equal(stan("out_sa"), as.vector(sa), tolerance = 1e-12)
  expect_equal(stan("out_bern"), log(c(0.7, 0.3)), tolerance = 1e-12)
})
