# simulate_data() for every model. The round trip runs each `simulate` case of
# the model-case registry (helper-model-cases.R) through check_data() and a
# mock fit. The statistic tests compare one statistic of n = 1e5 simulated
# observations of one cell or display with a value computed here from literal
# native values, within 5 Monte Carlo SE, so a wrong link conversion, a missing
# `2 * bound` or a dropped softmax reference fails.
#
# simulate_data() calls the r*() functions, so a bug that r*() shares with
# d*() cancels out wherever the reference is a d*() function. For m3 and mpt,
# r*() and d*() even share the probability vector; their references are
# therefore the closed forms of the model, written out here.

expect_within_mc <- function(sim, expected, label) {
  se <- sqrt(stats::var(sim) / length(sim))
  expect_lt(abs(mean(sim) - expected) / se, 5, label = label)
}

circular_mean_cos <- function(density) {
  stats::integrate(function(x) cos(x) * density(x), -pi, pi)$value
}

# P(upper boundary) of a Wiener process with drift v, boundary separation a,
# relative start point z and diffusion constant 1
p_upper <- function(v, a, z) (1 - exp(-2 * v * a * z)) / (1 - exp(-2 * v * a))

n_sim <- 1e5

test_that("simulate_data() output passes check_data() and a mock fit for every simulate case", {
  # about 30 cases at 0.3 s each: a mock fit per case is what proves the layout
  skip_on_cran()
  withr::local_seed(20261102)
  for (name in names(cases <- model_test_cases("simulate"))) {
    case <- cases[[name]]
    data <- simulate_data(case$model, case$pars, case$design, formula = case$formula)
    response <- case$response %||% unlist(case$model$resp_vars, use.names = FALSE)
    expect_in(c(names(case$design), response), names(data))
    expect_identical(nrow(data), nrow(case$design), label = name)
    expect_no_warning(suppressMessages(check_data(case$model, data, case$formula)),
                      message = name)
    # the custom m3 version warns that its default priors follow its links
    fit <- suppressWarnings(suppressMessages(bmm(case$formula, data, case$model,
                                                 backend = "mock", mock_fit = 1, rename = FALSE)))
    expect_identical(parameter_info(fit)$parameter,
                     parameter_info(case$model, case$formula)$parameter, label = name)
    expect_mapequal(lapply(attr(data, "simulate_pars"), unique), case$simulate_pars)
  }
})

test_that("scale = 'native' gives the data of the link values and of the generator", {
  model <- mixture2p("y")
  link <- withr::with_seed(5, simulate_data(model, list(kappa = 2, thetat = 0), n_trials = 50))
  native <- withr::with_seed(5, simulate_data(model, list(kappa = exp(2), thetat = 0.5),
                                              scale = "native", n_trials = 50))
  direct <- withr::with_seed(5, rmixture2p(50, mu = 0, kappa = exp(2), p_mem = 0.5))
  expect_identical(link$y, direct)
  expect_identical(native$y, direct)
})

test_that("the simulate_pars attribute holds the link values whatever the scale", {
  model <- mixture3p("y", nt_features = c("nt1", "nt2"), set_size = 3)
  design <- data.frame(nt1 = c(1, -1), nt2 = c(2, -2), cond = c("a", "b"))
  data <- simulate_data(model, list(kappa = exp(1), thetat = c(0.6, 0.5), thetant = 0.2),
                        design, scale = "native")
  expect_identical(data$cond, design$cond)
  expect_equal(attr(data, "simulate_pars"), data.frame(
    mu1 = c(0, 0), kappa = c(1, 1), thetat = log(c(0.6, 0.5) / c(0.2, 0.3)),
    thetant = log(0.2 / c(0.2, 0.3))
  ))
  expect_identical(attr(data, "simulate_call"), list(model = class(model), scale = "native"))
})

test_that("simulate_data() refuses what it cannot simulate", {
  model <- mixture2p("y")
  expect_error(simulate_data(list(name = "x"), list()), "has no simulate_data\\(\\) method")
  expect_error(simulate_data(model, list(kappa = 2, thetat = 0, p = 1), n_trials = 5),
               "Unknown parameter.*'p'")
  expect_error(simulate_data(model, list(kappa = 2), n_trials = 5), "No value.*'thetat'")
  expect_error(simulate_data(model, list(kappa = 1:2, thetat = 0), n_trials = 5),
               "'kappa'.*length 2")
  expect_error(simulate_data(model, list(kappa = TRUE, thetat = 0), n_trials = 5),
               "'kappa' must be numeric\\.$")
  expect_error(simulate_data(model, list(kappa = 2, thetat = 0), NULL, "link", 100),
               "must be named")
  expect_error(simulate_data(model, list(kappa = 2, thetat = 0)), "n_trials")
  expect_error(simulate_data(model, list(kappa = c(2, NA), thetat = NA), data.frame(x = 1:2)),
               "'kappa', 'thetat' contain NA")
  expect_error(simulate_data(model, list(kappa = 2, thetat = 0), data.frame(x = 1:4), n_trials = 10),
               "number of rows: drop `n_trials =`")
  expect_error(simulate_data(m3(resp_cats = c("corr", "npl"), num_options = c(1, 4), version = "ss"),
                             list(c = 1), data.frame(n_trials = 5), n_trials = 5),
               "read the trials of each row\\s+from a column")
  expect_identical(nrow(simulate_data(model, list(kappa = 2, thetat = 0), data.frame(x = 1:4),
                                      n_trials = NULL)), 4L)
  expect_error(simulate_data(sdt_yn("k", "stim", "n"), list(d = 1, criterion = 0)),
               "'stim', 'n'")
  expect_error(simulate_data(sdt_yn("k", "stim", "n"), list(d = 1, criterion = 0, bias = 0),
                             data.frame(stim = 0:1, n = 10)),
               "one of 'criterion', 'bias'")
  custom <- m3(resp_cats = c("corr", "npl"), num_options = c(1, 4), version = "custom",
               links = list(a = "identity"))
  expect_error(simulate_data(custom, list(a = 1), data.frame(n_trials = 5)), "formula")
  custom <- m3(resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5), version = "custom",
               links = list(a = "identity", c = "identity", c0 = "identity", c1 = "identity"))
  formula <- bmf(corr ~ b + a + c * cond, other ~ b + a, npl ~ b, a ~ 1, c ~ 1)
  for (cond in list(factor(c("low", "high")), c("low", "high"))) {
    expect_error(simulate_data(custom, list(a = 1, c = 1), data.frame(n_trials = 5, cond = cond),
                               formula = formula),
                 "design column\\(s\\) 'cond', so they must be numeric.*0/1 indicators")
  }
  expect_identical(nrow(simulate_data(custom, list(a = 1, c = 1), data.frame(n_trials = 5, cond = c(TRUE, FALSE)),
                                      formula = formula)), 2L)
  expect_error(simulate_data(custom, list(a = 1), data.frame(n_trials = 5, cond = 1),
                             formula = bmf(corr ~ b + a + c * cond, other ~ b + a, npl ~ b, a ~ 1)),
               "'c' in your activation formula\\(s\\) is neither a design column nor a model parameter")
  expect_error(simulate_data(custom, list(a = 1, c0 = 0, c1 = 0.5), data.frame(n_trials = 5, ss = 1, w = 9),
                             formula = bmf(corr ~ b + a + w, other ~ b + a, npl ~ b, a ~ 1,
                                           w ~ c0 + c1 * ss, c0 ~ 1, c1 ~ 1)),
               "'w' in the activation formula\\(s\\) is predicted by a non-linear formula of its own")
  expect_error(simulate_data(sdt_yn("k", "stim", "n"), list(d = 1, criterion = 0),
                             data.frame(stim = 0:1, n = c(10.7, 20))),
               "'n' hold numbers of trials: they must be whole numbers of 1 or more")
  expect_error(simulate_data(sdt_yn("k", "stim", "n"), list(d = 1, criterion = 0),
                             data.frame(stim = 0:1, n = c("10", "20"))),
               "'n' hold numbers of trials")
  expect_error(simulate_data(mpt(mpt_2htm_trees(), tree_id = "item"), list(D = 0, g = 0),
                             data.frame(item = c("old", "new"), n_trials = 10.7)),
               "'n_trials' hold numbers of trials")
  simplex <- mpt(mpt_tree("t", list(A = "gA", B = "gB", C = "gC")), simplex = c("gA", "gB", "gC"))
  expect_error(simulate_data(simplex, list(gA = 0.2, gB = 0.3, gC = 0.4), data.frame(n_trials = 5)),
               "sum to 1")
})

test_that("simulate_data() names the design column that holds NA", {
  expect_error(simulate_data(sdt_yn("k", "stim", "n"), list(d = 1, criterion = 0),
                             data.frame(stim = 0:1, n = c(10, NA))),
               "column\\(s\\) 'n' contain NA")
  model <- mixture3p("y", nt_features = c("nt1", "nt2"), set_size = "ss")
  pars <- list(kappa = 2, thetat = 0, thetant = 0)
  expect_error(simulate_data(model, pars, data.frame(ss = c(2, NA), nt1 = 1, nt2 = NA)),
               "column\\(s\\) 'ss' contain NA")
  expect_error(simulate_data(model, pars, data.frame(ss = c(3, 3, 2), nt1 = c(1, 1, 1), nt2 = c(1, NA, NA))),
               "'nt2' contain NA within the set size, in row\\(s\\) 2:")
  expect_error(simulate_data(model, pars, data.frame(ss = 3, nt1 = 1, nt2 = c(1, rep(NA, 7)))),
               "in row\\(s\\) 2, 3, 4, 5, 6, \\.\\.\\. \\(7 in all\\):")
  # past a row's set size NA is the layout of displays of different sizes
  data <- simulate_data(model, pars, data.frame(ss = c(1, 2, 3), nt1 = c(NA, 1, 1), nt2 = c(NA, NA, 1)))
  expect_false(anyNA(data$y))
  imm_model <- imm("y", nt_features = "nt1", nt_distances = "d1", set_size = 2, version = "full")
  expect_error(simulate_data(imm_model, list(c = 1, a = 0, s = 0, kappa = 2), data.frame(nt1 = 1, d1 = NA)),
               "'d1' contain NA within the set size")
})

test_that("the simulate_pars of a result simulate the same values again", {
  model <- sdt_yn("k", "stim", "n")
  design <- data.frame(stim = rep(0:1, 3), n = 20)
  data <- simulate_data(model, list(d = c(0.5, 1, 1.5, 0.5, 1, 1.5), criterion = 0.2), design)
  again <- simulate_data(model, attr(data, "simulate_pars"), design)
  expect_identical(attr(again, "simulate_pars"), attr(data, "simulate_pars"))
  expect_error(simulate_data(model, attr(data, "simulate_pars")[1:2, ], design),
               "one row or one row per row of the design \\(6\\)")
})

test_that("count models give each row its own number of trials and options", {
  model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"),
              choice_rule = "softmax", version = "ss")
  design <- data.frame(n_trials = c(1, 1000, 40, 7), n_corr = 1, n_other = c(0, 4, 0, 4), n_npl = 5)
  data <- simulate_data(model, list(c = 3, a = 3), design)
  expect_equal(unname(rowSums(data[c("corr", "other", "npl")])), c(1, 1000, 40, 7))
  expect_identical(data$other[design$n_other == 0], c(0L, 0L))
  mpt_model <- mpt(mpt_2htm_trees(), tree_id = "item")
  mpt_data <- simulate_data(mpt_model, list(D = 0, g = 0),
                            data.frame(item = c("old", "new"), n_trials = c(1, 1000)))
  expect_equal(unname(rowSums(mpt_data[c("old", "new")])), c(1, 1000))
})

test_that("the one-display generators get each trial's own values", {
  model <- mixture3p("y", nt_features = "nt1", set_size = 2)
  withr::local_seed(11)
  data <- simulate_data(model, list(kappa = 3, thetat = rep_len(c(5, -5), 2000), thetant = -5),
                        data.frame(nt1 = rep(2, 2000)))
  near_target <- abs(data$y) < 0.5
  expect_gt(mean(near_target[c(TRUE, FALSE)]), 0.9)
  expect_lt(mean(near_target[c(FALSE, TRUE)]), 0.2)
})

test_that("ezdm drops the summaries bmm() rejects and keeps those below ndt", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials")
  withr::local_seed(3)
  # with three trials and ndt = exp(-1) = 0.37, this seed draws mean RTs at or
  # below 0 and mean RTs between 0 and ndt
  expect_warning(
    data <- simulate_data(model, list(drift = 0, bound = 1.5, ndt = -1),
                          data.frame(n_trials = rep(3, 2000))),
    "Dropped [0-9]+ of 2000 rows.*not positive"
  )
  dropped <- attr(data, "dropped_rows")
  expect_gt(length(dropped), 0)
  expect_identical(nrow(data) + length(dropped), 2000L)
  expect_identical(nrow(attr(data, "simulate_pars")), nrow(data))
  expect_true(all(data$mean_rt > 0))
  expect_true(any(data$mean_rt < exp(-1)))
})

test_that("simulated statistics of the circular models match their densities", {
  skip_on_cran()
  withr::local_seed(1)
  y <- simulate_data(mixture2p("y"), list(kappa = 2, thetat = 0.4), n_trials = n_sim)$y
  expect_within_mc(cos(y), circular_mean_cos(function(x) dmixture2p(x, 0, 7.389056, 0.5986877)),
                   "mixture2p")

  y <- simulate_data(sdm("y"), list(c = 2, kappa = 1.75), n_trials = n_sim)$y
  expect_within_mc(cos(y), circular_mean_cos(function(x) dsdm(x, 0, 7.389056, 5.754603)), "sdm")

  display <- data.frame(set_size = rep(3, n_sim), nt1 = 1.5, nt2 = -2, d1 = 1, d2 = 2)
  model <- mixture3p("y", nt_features = c("nt1", "nt2"), set_size = "set_size")
  y <- simulate_data(model, list(kappa = 2, thetat = 0.5, thetant = -0.5), display)$y
  # softmax(c(0.5, -0.5, 0)): the guessing weight is the reference
  expect_within_mc(cos(y), circular_mean_cos(function(x) {
    dmixture3p(x, c(0, 1.5, -2), 7.389056, p_mem = 0.5064804, p_nt = 0.1863237)
  }), "mixture3p")

  model <- imm("y", nt_features = c("nt1", "nt2"), nt_distances = c("d1", "d2"),
               set_size = "set_size", version = "full")
  y <- simulate_data(model, list(kappa = 2, a = 0, c = 0.5, s = 0), display)$y
  expect_within_mc(cos(y), circular_mean_cos(function(x) {
    dimm(x, c(0, 1.5, -2), c(0, 1, 2), c = 1.648721, a = 1, b = 1, s = 1, kappa = 7.389056)
  }), "imm")

  model <- imm("y", nt_features = c("nt1", "nt2"), nt_distances = c("d1", "d2"),
               set_size = "set_size", version = "bsc")
  y <- simulate_data(model, list(kappa = 2, c = 0.5, s = 0), display)$y
  expect_within_mc(cos(y), circular_mean_cos(function(x) {
    dimm(x, c(0, 1.5, -2), c(0, 1, 2), c = 1.648721, a = 0, b = 1, s = 1, kappa = 7.389056)
  }), "imm bsc")

  # abc: weights c + a for the target, a for each non-target, 1 for guessing;
  # E[cos] of a von Mises at mu is A(kappa) cos(mu), A(7.389056) = 0.9296366
  model <- imm("y", nt_features = c("nt1", "nt2"), set_size = "set_size", version = "abc")
  y <- simulate_data(model, list(kappa = 2, c = 0.5, a = -0.5), display)$y
  expect_within_mc(cos(y), 0.9296366 * (2.255252 + 0.6065307 * (cos(1.5) + cos(-2))) / 4.468313,
                   "imm abc")

  # tcc: P(option k) = integral of dnorm(z - d f_k) prod_j pnorm(z - d f_j) dz,
  # by quadrature with d = e^0.5 on a wheel of 12 with f = exp(-(0:6) / 2)
  model <- tcc("y", similarity = exp(-(0:6) / 2), n_options = 12)
  y <- simulate_data(model, list(d = 0.5), n_trials = n_sim)$y
  expect_true(all(abs(y * 6 / pi - round(y * 6 / pi)) < 1e-12))
  expect_within_mc(y == 0, 0.3484979, "tcc P(target)")
  expect_within_mc(cos(y), 0.5483012, "tcc")
  # Gumbel noise: the softmax of d f over the options
  model <- tcc("y", similarity = exp(-(0:6) / 2), n_options = 12, noise = "gumbel")
  y <- simulate_data(model, list(d = 0.5), n_trials = n_sim)$y
  expect_within_mc(y == 0, 0.2254180, "tcc gumbel P(target)")
  expect_within_mc(cos(y), 0.3218072, "tcc gumbel")
  # option sets on the same wheel: {0, pi/2} with Gaussian noise has
  # P(target) = pnorm(d (1 - f(3)) / sqrt(2)); four alternatives with Gumbel
  # noise the softmax of d f over the offsets 0, 3, 6 and -3
  model <- tcc("y", similarity = exp(-(0:6) / 2), n_options = 12, options = c("o1", "o2"))
  y <- simulate_data(model, list(d = 0.5), data.frame(o1 = rep(0, n_sim), o2 = pi / 2))$y
  expect_true(all(y %in% c(0, pi / 2)))
  expect_within_mc(y == 0, 0.8174505, "tcc options {0, pi/2} P(target)")
  model <- tcc("y", similarity = exp(-(0:6) / 2), n_options = 12, n_alternatives = 4, noise = "gumbel")
  y <- simulate_data(model, list(d = 0.5), n_trials = n_sim)$y
  expect_within_mc(y == 0, 0.5667801, "tcc 4 alternatives gumbel P(target)")
  # correlated noise, corr = exp(-k / 1.5) at k steps: P(target) is the
  # orthant probability of the differences X_0 - X_j (mvtnorm::pmvnorm);
  # motor noise of one step: the independent P convolved with the rounded
  # normal; the independent P(target) is 0.3484979
  model <- tcc("y", similarity = exp(-(0:6) / 2), n_options = 12, correlation = exp(-(0:6) / 1.5))
  y <- simulate_data(model, list(d = 0.5), n_trials = n_sim)$y
  expect_within_mc(y == 0, 0.3896782, "tcc correlated P(target)")
  model <- tcc("y", similarity = exp(-(0:6) / 2), n_options = 12, motor_sd = pi / 6)
  y <- simulate_data(model, list(d = 0.5), n_trials = n_sim)$y
  expect_within_mc(y == 0, 0.2086131, "tcc motor P(target)")
  # with units = "bins" motor_sd counts option steps: the same one-step cell
  model <- tcc("y", similarity = exp(-(0:6) / 2), n_options = 12, units = "bins", motor_sd = 1)
  y <- simulate_data(model, list(d = 0.5), n_trials = n_sim)$y
  expect_within_mc(y == 0, 0.2086131, "tcc bins motor P(target)")
  bins <- simulate_data(tcc("y", similarity = exp(-(0:6) / 2), n_options = 12, units = "bins"),
                        list(d = 0.5), n_trials = 50)$y
  expect_true(all(bins %in% -5:6))
})

test_that("simulated decisions of the RT models match the closed forms", {
  skip_on_cran()
  withr::local_seed(2)
  data <- simulate_data(ddm("rt", "response"), list(drift = 1, bound = 0, ndt = -1.2, zr = 0.4),
                        n_trials = n_sim)
  expect_within_mc(data$response, p_upper(v = 1, a = 1, z = 0.5986877), "ddm")

  # the simple version's bound is half the boundary separation
  data <- simulate_data(cswald("rt", "response", version = "simple"),
                        list(drift = 0, bound = -0.2876821, ndt = -1.2), n_trials = n_sim)
  expect_within_mc(data$response, p_upper(v = 1, a = 1.5, z = 0.5), "cswald simple")

  data <- simulate_data(cswald("rt", "response", version = "crisk"),
                        list(drift = 0.8, bound = 0, ndt = -1.2, zr = 0.3), n_trials = n_sim)
  expect_within_mc(data$response, p_upper(v = 0.8, a = 1, z = 0.5744425), "cswald crisk")

  design <- data.frame(n_trials = rep(100, n_sim))
  data <- simulate_data(ezdm("mean_rt", "var_rt", "n_upper", "n_trials"),
                        list(drift = 1, bound = 0.1823216, ndt = -1.2), design)
  expect_within_mc(data$n_upper / 100, stats::plogis(1 * 1.2), "ezdm 3par")

  data <- simulate_data(ezdm(c("m_up", "m_lo"), c("v_up", "v_lo"), "n_upper", "n_trials",
                             version = "4par"),
                        list(drift = 1, bound = 0.1823216, ndt = -1.2, zr = -0.3), design)
  expect_within_mc(data$n_upper / 100, p_upper(v = 1, a = 1.2, z = 0.4255575), "ezdm 4par")
})

test_that("simulated responses of the signal detection models match their densities", {
  skip_on_cran()
  withr::local_seed(3)
  signal <- data.frame(stimulus = rep(1, n_sim), n_trials = 1)

  k <- simulate_data(sdt_yn("k", "stimulus", "n_trials"),
                     list(d = 1, criterion = 0.5, sdratio = 0.2231436), signal)$k
  expect_within_mc(k, dsdt_yn(1, 1, 1, d = 1, criterion = 0.5, sdratio = 1.25), "sdt_yn")

  k <- simulate_data(sdt_mafc("k", "n_trials", m = 4), list(d = 1), signal)$k
  expect_within_mc(k, dsdt_mafc(1, 1, 4, d = 1), "sdt_mafc")

  # parsimonious thresholds: criterion + exp(spacing) * log(k / (K - k))
  thresholds <- c(-0.8862944, 0.0945349, 0.9054651, 1.8862944)
  rating <- function(version, ...) {
    sdt_rating(paste0("r", 1:5), "stimulus", version = version, ...)
  }
  base <- list(d = 1, criterion = 0.5, spacing = 0, sdratio = 0.2)
  k <- simulate_data(rating("standard"), base, signal)$r5
  expect_within_mc(k, dsdt_rating(c(0, 0, 0, 0, 1), 1, d = 1, thresholds = thresholds,
                                  sdratio = 1.221403), "sdt_rating standard")
  k <- simulate_data(rating("dpsdt"), c(base, Ro = 0.3), signal)$r5
  expect_within_mc(k, dsdt_dpsdt(c(0, 0, 0, 0, 1), 1, d = 1, thresholds = thresholds,
                                 Ro = 0.5744425, Rn = 0, sdratio = 1.221403), "sdt_rating dpsdt")
  metad <- sdt_rating(paste0("r", 1:4), "stimulus", version = "metad")
  k <- simulate_data(metad, c(base, logmratio = -0.2231436), signal)$r4
  expect_within_mc(k, dsdt_metad(c(0, 0, 0, 1), 1, d = 1, thresholds = c(-0.5986123, 0.5, 1.5986123),
                                 metad = 0.8, sdratio = 1.221403), "sdt_rating metad")

  k <- simulate_data(sdt_ranking(paste0("rank", 1:4), m = 4, dist = "normal"),
                     list(d = 1, sdratio = 0.2), signal)$rank1
  expect_within_mc(k, dsdt_ranking(c(1, 0, 0, 0), 4, d = 1, sdratio = 1.221403, dist = "normal"),
                   "sdt_ranking")

  k <- simulate_data(sdt_cdp("", "stimulus", n_new = 2, n_old = 2),
                     list(dfam = 1, drec = 1, criterion = 0, spacing = 0, rcrit = 0), signal)
  expect_within_mc(k$remember4, dsdt_cdp(c(0, 0, 0, 0, 0, 1), 1, dfam = 1, drec = 1,
                                         thresholds = c(-1.0986123, 0, 1.0986123), rcrit = 0,
                                         n_new = 2), "sdt_cdp")
})

test_that("simulated counts of the m3 and mpt models match their closed forms", {
  skip_on_cran()
  withr::local_seed(4)
  model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
              choice_rule = "softmax", version = "ss")
  k <- simulate_data(model, list(c = 2, a = 1), data.frame(n_trials = rep(1, n_sim)))$corr
  # softmax activations b + a + c, b + a, b with b = 0, times the number of options
  expect_within_mc(k, exp(3) / (exp(3) + 4 * exp(1) + 5), "m3")

  model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
              choice_rule = "simple", version = "ss")
  k <- simulate_data(model, list(c = 1, a = -1), data.frame(n_trials = rep(1, n_sim)))$corr
  # activations b + a + c, b + a, b on the native scale, with the fixed b = 0.1
  expect_within_mc(k, 3.18616 / (3.18616 + 4 * 0.4678794 + 5 * 0.1), "m3 simple")

  model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
              choice_rule = "softmax", version = "custom", links = list(a = "identity", c = "identity"))
  formula <- bmf(corr ~ b + a + c * ss, other ~ b + a, npl ~ b, a ~ 1, c ~ 1)
  data <- simulate_data(model, list(a = 1, c = 0.5), data.frame(n_trials = 1, ss = rep(c(1, 4), n_sim / 2)),
                        formula = formula)
  # corr activation b + a + c * ss with b = 0: 1.5 at ss = 1, 3 at ss = 4
  expect_within_mc(data$corr[data$ss == 1], exp(1.5) / (exp(1.5) + 4 * exp(1) + 5), "m3 custom ss 1")
  expect_within_mc(data$corr[data$ss == 4], exp(3) / (exp(3) + 4 * exp(1) + 5), "m3 custom ss 4")

  model <- mpt(mpt_2htm_trees(), tree_id = "item")
  k <- simulate_data(model, list(D = 0.5, g = -0.3),
                     data.frame(item = rep("old", n_sim), n_trials = 1))$old
  expect_within_mc(k, 0.6224593 + (1 - 0.6224593) * 0.4255575, "mpt")
})
