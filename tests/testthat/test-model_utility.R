rhs_string <- function(bmmform, par) {
  f <- bmmform[[par]]
  deparse1(f[[length(f)]])
}

welfare_model <- function(choice_rule = "simple") {
  utility(
    resp_cats = c("nkeep", "ningroup", "nuniversal"),
    payoffs = rbind(
      nkeep      = c(b = 1.0, wi = 0.0, wo = 0.0),
      ningroup   = c(b = 0.5, wi = 1.0, wo = 0.0),
      nuniversal = c(b = 0.3, wi = 0.6, wo = 0.9)
    ),
    choice_rule = choice_rule,
    numeraire = c(b = 1)
  )
}

welfare_data <- function(n = 30) {
  data.frame(
    nkeep = rpois(n, 8), ningroup = rpois(n, 10), nuniversal = rpois(n, 12)
  )
}

# ---- constructor structure -------------------------------------------------

test_that("utility() builds the payoff (welfare) model with the right structure", {
  m <- welfare_model()
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "utility")
  expect_s3_class(m, "utility_linear")
  expect_equal(m$version, "linear")
  expect_equal(names(m$parameters), c("b", "wi", "wo"))
  expect_equal(m$links, list(wi = "identity", wo = "identity"))
  expect_equal(m$fixed_parameters, list(b = 1))
  expect_true(inherits(m, supported_models(print_call = FALSE)))
})

test_that("utility() derives the version token from utility_fn x weighting", {
  base_args <- list(resp_cats = c("corr", "other"), value_cols = c(corr = "V"))
  expect_equal(do.call("utility", base_args)$version, "linear")
  expect_equal(do.call("utility", c(base_args, utility_fn = "power"))$version, "power")
  expect_equal(do.call("utility", c(base_args, weighting = "prelec"))$version, "prelec")
  expect_equal(
    do.call("utility", c(base_args, utility_fn = "power", weighting = "prelec"))$version,
    "powerprelec"
  )
})

test_that("value_cols parameters carry the correct links", {
  m <- utility(
    resp_cats = c("corr", "other"), value_cols = c(corr = "V"),
    utility_fn = "power", weighting = "prelec"
  )
  expect_equal(names(m$parameters), c("b", "gamma", "rho", "alpha"))
  expect_equal(m$links$gamma, "identity")
  expect_equal(m$links$rho, "log")
  expect_equal(m$links$alpha, "log")
})

# ---- default priors --------------------------------------------------------

test_that("default priors follow the round-5-correct schema", {
  mw <- welfare_model()
  expect_equal(mw$default_priors$wi, list(main = "normal(0, 1)", effects = "normal(0, 0.5)"))

  mv <- utility(
    resp_cats = c("corr", "other"), value_cols = c(corr = "V"),
    utility_fn = "power", weighting = "prelec"
  )
  expect_equal(mv$default_priors$gamma, list(main = "normal(0, 1)", effects = "normal(0, 0.5)"))
  # rho / alpha are normal on the log-link predictor, NOT lognormal
  expect_equal(mv$default_priors$rho, list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)"))
  expect_equal(mv$default_priors$alpha, list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)"))
})

test_that("numeraire fixes b to a constant prior", {
  dat <- welfare_data()
  prior <- suppressWarnings(default_prior(bmf(wi ~ 1, wo ~ 1), dat, welfare_model()))
  b_row <- prior[prior$nlpar == "b" & prior$coef == "Intercept", "prior"]
  expect_equal(b_row, "constant(1)")

  # explicit non-unit numeraire
  m2 <- utility(
    resp_cats = c("nkeep", "ningroup", "nuniversal"),
    payoffs = rbind(
      nkeep = c(b = 1, wi = 0, wo = 0), ningroup = c(b = 0.5, wi = 1, wo = 0),
      nuniversal = c(b = 0.3, wi = 0.6, wo = 0.9)
    ),
    numeraire = c(b = 2)
  )
  expect_equal(m2$fixed_parameters$b, 2)
})

# ---- activation-formula generation -----------------------------------------

test_that("payoff activations use per-category coefficients and drop zero terms", {
  af <- construct_utility_act_funs(welfare_model())
  expect_equal(rhs_string(af, "nkeep"), "1 * b")
  expect_equal(rhs_string(af, "ningroup"), "0.5 * b + 1 * wi")
  expect_equal(rhs_string(af, "nuniversal"), "0.3 * b + 0.6 * wi + 0.9 * wo")
})

test_that("payoff cells accept per-trial data columns", {
  m <- utility(
    resp_cats = c("a", "b", "c"),
    payoffs = data.frame(
      row.names = c("a", "b", "c"),
      b = c(1, 0.5, 0.3), wi = c(0, 1, "wi_trial"), wo = c(0, 0, "wo_trial")
    ),
    choice_rule = "simple", numeraire = c(b = 1)
  )
  af <- construct_utility_act_funs(m)
  expect_equal(rhs_string(af, "c"), "0.3 * b + wi_trial * wi + wo_trial * wo")
})

test_that("value_cols activations build linear, power, and prelec terms", {
  lin <- construct_utility_act_funs(
    utility(resp_cats = c("corr", "npl"), value_cols = c(corr = "V"))
  )
  expect_equal(rhs_string(lin, "corr"), "b + gamma * V")
  expect_equal(rhs_string(lin, "npl"), "b")

  pw <- construct_utility_act_funs(
    utility(resp_cats = c("corr", "npl"), value_cols = c(corr = "V"), utility_fn = "power")
  )
  expect_equal(rhs_string(pw, "corr"), "b + gamma * V^rho")

  pr <- construct_utility_act_funs(
    utility(
      resp_cats = c("corr", "npl"), num_options = c("nc", "nn"),
      value_cols = c(corr = "V"), weighting = "prelec"
    )
  )
  expect_match(rhs_string(pr, "corr"), "gamma * V", fixed = TRUE)
  expect_match(rhs_string(pr, "npl"), "-log(nn/n_opt_total)", fixed = TRUE)
})

test_that("log-linked parameters are exp-transformed in the checked formula", {
  m <- utility(resp_cats = c("corr", "npl"), value_cols = c(corr = "V"), utility_fn = "power")
  form <- check_formula(m, welfare_data(), bmf(gamma ~ 1, rho ~ 1))
  expect_match(rhs_string(form, "corr"), "V^exp(rho)", fixed = TRUE)
  expect_true(is_nl(form)[["corr"]])
})

# ---- choice kernel ---------------------------------------------------------

test_that("the choice kernel differs by rule and drops log(n) for prelec", {
  opt <- c(corr = "n_corr")
  expect_equal(utility_choice_kernel("softmax", "none", "corr", opt), "(corr + log(n_corr))")
  expect_equal(utility_choice_kernel("simple", "none", "corr", opt), "log(corr * n_corr)")
  expect_equal(utility_choice_kernel("softmax", "prelec", "corr", opt), "(corr)")
})

# ---- full formula assembly (mock, no Stan) ---------------------------------

test_that("configure_model builds a multinomial family with mu<cat> dpars", {
  m <- welfare_model()
  dat <- welfare_data()
  fit <- suppressWarnings(
    bmm(bmf(wi ~ 1, wo ~ 1), dat, m, backend = "mock", mock = 1, rename = FALSE)
  )
  expect_equal(fit$formula$family$family, "multinomial")
  expect_setequal(fit$formula$family$dpars, c("munkeep", "muningroup", "mununiversal"))
})

# ---- data-property validation ----------------------------------------------

test_that("missing per-trial payoff / value columns error with the column name", {
  mp <- utility(
    resp_cats = c("a", "b"),
    payoffs = data.frame(row.names = c("a", "b"), b = c(1, 0.5), wi = c(0, "wicol")),
    choice_rule = "simple", numeraire = c(b = 1)
  )
  dat <- data.frame(a = rpois(20, 5), b = rpois(20, 5))
  expect_error(check_data(mp, dat, bmf(wi ~ 1)), "wicol")

  mv <- utility(resp_cats = c("corr", "other"), value_cols = c(corr = "V"))
  expect_error(
    check_data(mv, data.frame(corr = rpois(20, 5), other = rpois(20, 5)), bmf(gamma ~ 1)),
    "V"
  )
})

# ---- identifiability guards ------------------------------------------------

test_that("G3 stops when a value column is constant, passes when it varies", {
  m <- utility(resp_cats = c("corr", "other"), value_cols = c(corr = "V"))
  n <- 30
  const <- data.frame(corr = rpois(n, 5), other = rpois(n, 5), V = 2)
  vary <- data.frame(corr = rpois(n, 5), other = rpois(n, 5), V = sample(c(1, 3), n, TRUE))
  expect_error(check_data(m, const, bmf(gamma ~ 1)), "gamma")
  expect_no_error(check_data(m, vary, bmf(gamma ~ 1)))
})

test_that("G2 stops when the choice-set proportions are constant", {
  m <- utility(
    resp_cats = c("corr", "other"), num_options = c("nc", "no"),
    value_cols = c(corr = "V"), weighting = "prelec"
  )
  n <- 30
  V <- sample(c(1, 3), n, TRUE)
  const <- data.frame(corr = rpois(n, 5), other = rpois(n, 5), V = V, nc = 1L, no = 1L)
  vary <- data.frame(corr = rpois(n, 5), other = rpois(n, 5), V = V, nc = 1L,
                     no = sample(c(2L, 6L), n, TRUE))
  expect_error(check_data(m, const, bmf(gamma ~ 1, alpha ~ 1)), "alpha")
  expect_no_error(check_data(m, vary, bmf(gamma ~ 1, alpha ~ 1)))
})

test_that("G1 warns on a narrow power-utility value range", {
  m <- utility(resp_cats = c("corr", "other"), value_cols = c(corr = "V"), utility_fn = "power")
  n <- 30
  narrow <- data.frame(corr = rpois(n, 5), other = rpois(n, 5), V = sample(c(1, 2), n, TRUE))
  wide <- data.frame(corr = rpois(n, 5), other = rpois(n, 5), V = sample(c(1, 5), n, TRUE))
  expect_warning(check_data(m, narrow, bmf(gamma ~ 1, rho ~ 1)), "rho")
  expect_no_warning(check_data(m, wide, bmf(gamma ~ 1, rho ~ 1)))
})

test_that("G4 warns for the simple rule with identity-linked parameters", {
  expect_warning(
    bmm(bmf(wi ~ 1, wo ~ 1), welfare_data(), welfare_model("simple"),
        backend = "mock", mock = 1, rename = FALSE),
    "simple"
  )
  expect_no_warning(
    bmm(bmf(wi ~ 1, wo ~ 1), welfare_data(), welfare_model("softmax"),
        backend = "mock", mock = 1, rename = FALSE)
  )
})

# ---- boundary validation ---------------------------------------------------

test_that("the constructor validates its arguments at the boundary", {
  pm <- rbind(a = c(b = 1, wi = 0), bb = c(b = 0.5, wi = 1))

  # exactly one of payoffs / value_cols
  expect_error(
    utility(resp_cats = c("a", "bb"), payoffs = pm, value_cols = c(a = "V", bb = "W")),
    "exactly one"
  )
  # G5: payoff path forbids non-default utility_fn / weighting
  expect_error(utility(resp_cats = c("a", "bb"), payoffs = pm, utility_fn = "power"), "defaults")
  # prelec requires softmax
  expect_error(
    utility(resp_cats = c("a", "bb"), value_cols = c(a = "V"),
            weighting = "prelec", choice_rule = "simple"),
    "softmax"
  )
  # numeraire column must be present
  expect_error(
    utility(resp_cats = c("a", "bb"), payoffs = rbind(a = c(wi = 1), bb = c(wi = 2))),
    "numeraire"
  )
  # payoffs needs dimnames
  expect_error(utility(resp_cats = c("a", "bb"), payoffs = matrix(1, 2, 2)), "row names")
})
