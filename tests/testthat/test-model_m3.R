test_that("construct_m3_act_funs works with simple m3", {
  model <- m3(
    resp_cats = c("correct", "other", "npl"),
    num_options = c(1, 4, 5),
    version = "ss"
  )
  expect_equal(
    construct_m3_act_funs(model, warnings = FALSE),
    bmf(correct ~ b + a + c, other ~ b + a, npl ~ b),
    ignore_formula_env = TRUE
  )
})

test_that("construct_m3_act_funs works with complex span m3", {
  model <- m3(
    resp_cats = c("correct", "dist_context", "other", "dist_other", "npl"),
    num_options = c(1, 4, 5, 4, 5),
    version = "cs"
  )
  expect_equal(
    construct_m3_act_funs(model, warnings = FALSE),
    bmf(
      correct ~ b + a + c,
      dist_context ~ b + f * a + f * c,
      other ~ b + a,
      dist_other ~ b + f * a,
      npl ~ b
    ),
    ignore_formula_env = TRUE
  )
})

test_that("construct_m3_act_funs gives error for other models", {
  model <- m3(
    resp_cats = c("correct", "dist_context", "other", "dist_other", "npl"),
    num_options = c(1, 4, 5, 4, 5),
    version = "custom"
  )
  expect_error(construct_m3_act_funs(model), "can only be generated for")

  model <- sdm("dev_rad")
  expect_error(construct_m3_act_funs(model), "can only be generated for")
})

test_that("m3 compiles for the simple_span / simple choice rule", {
  formula <- bmf(
    c ~ 1 + cond + (1 + cond || ID),
    a ~ 1 + cond + (1 + cond || ID)
  )

  my_model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "simple",
    version = "ss"
  )

  expect_silent(bmm(
    formula = formula,
    data = oberauer_lewandowsky_2019_e1,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  ))
})

test_that("m3 compiles for the simple_span / softmax choice rule", {
  formula <- bmf(
    c ~ 1 + cond + (1 + cond || ID),
    a ~ 1 + cond + (1 + cond || ID)
  )

  my_model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "softmax",
    version = "ss"
  )

  expect_silent(bmm(
    formula = formula,
    data = oberauer_lewandowsky_2019_e1,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  ))
})

test_that("m3 compiles for the complex_span / simple choice rule", {
  data <- oberauer_lewandowsky_2019_e1
  data$distother <- data$dist
  data$n_dist_other <- data$n_dist
  formula <- bmf(
    c ~ 1 + cond + (1 + cond || ID),
    a ~ 1 + cond + (1 + cond || ID),
    f ~ 1
  )

  my_model <- m3(
    resp_cats = c("corr", "dist", "other", "distother", "npl"),
    num_options = c("n_corr", "n_dist", "n_other", "n_dist_other", "n_npl"),
    choice_rule = "simple",
    version = "cs"
  )

  expect_silent(bmm(
    formula = formula,
    data = data,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  ))
})

test_that("m3 compiles for the complex_span / softmax choice rule", {
  data <- oberauer_lewandowsky_2019_e1
  data$distother <- data$dist
  data$n_dist_other <- data$n_dist
  formula <- bmf(
    c ~ 1 + cond + (1 + cond || ID),
    a ~ 1 + cond + (1 + cond || ID),
    f ~ 1
  )

  my_model <- m3(
    resp_cats = c("corr", "dist", "other", "distother", "npl"),
    num_options = c("n_corr", "n_dist", "n_other", "n_dist_other", "n_npl"),
    choice_rule = "softmax",
    version = "cs"
  )

  expect_silent(bmm(
    formula = formula,
    data = data,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  ))
})

test_that("m3 compiles for the custom model / simple choice rule", {
  formula <- bmf(
    corr ~ b + a + c,
    other ~ b + a,
    dist ~ b + d,
    npl ~ b,
    c ~ 1 + cond + (1 + cond || ID),
    a ~ 1 + cond + (1 + cond || ID),
    d ~ 1 + (1 || ID)
  )

  my_links <- list(c = "log", a = "log", d = "log")

  my_priors <- list(
    c = list(main = "normal(2, 0.5)", effects = "normal(0, 0.5)"),
    a = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)"),
    d = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)")
  )

  my_model <- m3(
    resp_cats = c("corr", "other", "dist", "npl"),
    num_options = c("n_corr", "n_other", "n_dist", "n_npl"),
    choice_rule = "simple",
    links = my_links,
    default_priors = my_priors
  )

  expect_silent(bmm(
    formula = formula,
    data = oberauer_lewandowsky_2019_e1,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  ))
})

test_that("m3 compiles for the custom model / softmax choice rule", {
  formula <- bmf(
    corr ~ b + a + c,
    other ~ b + a,
    dist ~ b + d,
    npl ~ b,
    c ~ 1 + cond + (1 + cond || ID),
    a ~ 1 + cond + (1 + cond || ID),
    d ~ 1 + (1 || ID)
  )

  my_links <- list(c = "log", a = "log", d = "log")

  my_priors <- list(
    c = list(main = "normal(2, 0.5)", effects = "normal(0, 0.5)"),
    a = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)"),
    d = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)")
  )

  my_model <- m3(
    resp_cats = c("corr", "other", "dist", "npl"),
    num_options = c("n_corr", "n_other", "n_dist", "n_npl"),
    choice_rule = "softmax",
    links = my_links,
    default_priors = my_priors
  )

  expect_silent(bmm(
    formula = formula,
    data = oberauer_lewandowsky_2019_e1,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  ))
})


test_that("m3 custom accepts softplus links and generates default priors", {
  formula <- bmf(
    corr ~ b + a + c,
    other ~ b + a,
    npl ~ b,
    c ~ 1,
    a ~ 1
  )

  softplus_main <- list(simple = "normal(2, 1)", softmax = "normal(1, 1)")

  for (rule in c("simple", "softmax")) {
    my_model <- m3(
      resp_cats = c("corr", "other", "npl"),
      num_options = c(1, 2, 5),
      choice_rule = rule,
      links = list(c = "softplus", a = "softplus")
    )

    fit <- suppressWarnings(bmm(
      formula = formula,
      data = oberauer_lewandowsky_2019_e1,
      model = my_model,
      backend = "mock",
      mock_fit = 1,
      rename = F
    ))

    for (par in c("c", "a")) {
      expect_equal(
        fit$bmm$model$default_priors[[par]],
        list(main = softplus_main[[rule]], effects = "normal(0, 0.5)", sd = "exponential(1)")
      )
    }
  }
})

test_that("m3 works with num_options as a numeric vector", {
  formula <- bmf(
    c ~ 1 + (1 | ID),
    a ~ 1 + (1 | ID)
  )

  my_model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 2, 5),
    choice_rule = "simple",
    version = "ss"
  )

  fit <- bmm(
    formula = formula,
    data = oberauer_lewandowsky_2019_e1,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  )

  nopts <- my_model$other_vars$num_options
  expect_named(nopts, paste0("n_opt_", my_model$resp_vars$resp_cats))
  expect_equal(unlist(unique(fit$data[names(nopts)])), nopts)
})

test_that("m3_custom version works with variables contained in data in the activation formulas", {
  my_data <- data.frame(
    corr = c(5, 6, 7, 8),
    other = c(1, 2, 3, 4),
    npl = c(1, 2, 3, 4),
    time = c(1, 2, 1, 2),
    id = c(1, 1, 2, 2)
  )

  formula <- bmf(
    corr ~ b + a + cstart + cslope * time,
    other ~ b + a,
    npl ~ b,
    a ~ 1,
    cstart ~ 1,
    cslope ~ 1
  )

  my_model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 2, 3),
    choice_rule = "softmax",
    version = "custom"
  )

  my_model$links <- list(
    a = "log",
    cstart = "log",
    cslope = "log"
  )

  my_model$default_priors <- list(
    a = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)"),
    cstart = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)"),
    cslope = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)")
  )

  expect_silent(bmm(
    formula = formula,
    data = my_data,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  ))
})


test_that("m3 with numerical vector as num_options containing 0 returns error", {
  formula <- bmf(
    c ~ 1 + (1 | ID),
    a ~ 1 + (1 | ID)
  )

  my_model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 0, 5),
    choice_rule = "simple",
    version = "ss"
  )

  expect_error(bmm(
    formula = formula,
    data = oberauer_lewandowsky_2019_e1,
    model = my_model,
    backend = "mock",
    mock_fit = 1,
    rename = F
  ), "not identified")
})

m3_num_options_fit <- function(num_options, choice_rule = "simple") {
  suppressWarnings(bmm(
    bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b, c ~ 1, a ~ 1, d ~ 1),
    oberauer_lewandowsky_2019_e1,
    m3(
      resp_cats = c("corr", "other", "dist", "npl"), num_options = num_options,
      choice_rule = choice_rule, links = list(c = "log", a = "log", d = "log")
    ),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
}

test_that("num_options named after the response categories are matched by name (#449)", {
  for (choice_rule in c("simple", "softmax")) {
    unnamed <- m3_num_options_fit(c(1, 4, 5, 5), choice_rule)
    by_category <- m3_num_options_fit(c(npl = 5, other = 4, corr = 1, dist = 5), choice_rule)
    expect_equal(by_category$formula, unnamed$formula)
    expect_equal(brms::standata(by_category), brms::standata(unnamed))
  }
})

test_that("character num_options named after the response categories are matched by name (#457)", {
  by_position <- c("n_corr", "n_other", "n_dist", "n_npl")
  by_category <- c(other = "n_other", npl = "n_npl", corr = "n_corr", dist = "n_dist")
  for (choice_rule in c("simple", "softmax")) {
    unnamed <- m3_num_options_fit(by_position, choice_rule)
    named <- m3_num_options_fit(by_category, choice_rule)
    expect_equal(named$formula, unnamed$formula)
    expect_equal(brms::standata(named), brms::standata(unnamed))
  }
})

test_that("the ss version matches character num_options named after the categories (#457)", {
  fit_ss <- function(num_options) {
    bmm(
      bmf(c ~ 1, a ~ 1),
      oberauer_lewandowsky_2019_e1,
      m3(
        resp_cats = c("corr", "other", "npl"), num_options = num_options,
        choice_rule = "simple", version = "ss"
      ),
      backend = "mock", mock_fit = 1, rename = FALSE
    )
  }
  unnamed <- fit_ss(c("n_corr", "n_other", "n_npl"))
  named <- fit_ss(c(npl = "n_npl", corr = "n_corr", other = "n_other"))
  expect_equal(named$formula, unnamed$formula)
  expect_equal(brms::standata(named), brms::standata(unnamed))
})

test_that("character num_options that are partly named, duplicated or incomplete give an error (#457)", {
  cats <- c("corr", "other", "npl")
  expect_error(m3(cats, num_options = c(corr = "n_corr", "n_other", "n_npl")), "all elements")
  expect_error(m3(cats, num_options = c(k = "n_corr", k = "n_other", j = "n_npl")), "only once")
  expect_error(
    m3(cats, num_options = c(corr = "n_corr", other = "n_other", dist = "n_npl")),
    "one element for each"
  )
})

test_that("m3 refuses NA among numeric num_options (#457)", {
  expect_error(m3(c("corr", "other", "npl"), num_options = c(1, NA, 3)), "num_options")
})

test_that("num_options names already taken by a column or parameter give an error", {
  expect_error(m3_num_options_fit(c(a = 1, b = 4, c = 5, d = 5)), "'a', 'b', 'c', 'd'")
  expect_error(m3_num_options_fit(c(nTrials = 1, k2 = 4, k3 = 5, k4 = 5)), "'nTrials'")
  expect_error(m3_num_options_fit(c(Idx_dist = 1, k2 = 4, k3 = 5, k4 = 5)), "'Idx_dist'")
  expect_error(m3_num_options_fit(c(ID = 1, k2 = 4, k3 = 5, k4 = 5)), "'ID'")
})

test_that("m3 rejects num_options it cannot map onto the response categories", {
  cats <- c("corr", "other", "npl")
  expect_error(m3(cats, num_options = c(corr = 1, 4, 5)), "all elements")
  expect_error(m3(cats, num_options = c(k = 1, k = 4, j = 5)), "only once")
  expect_error(m3(cats, num_options = c(corr = 1, other = 4, dist = 5)), "one element for each")
  expect_error(m3(cats, num_options = cats), "response category column")
})

test_that("softmax default priors give a and c equal main means (c - a centered at 0)", {
  for (v in c("ss", "cs")) {
    p <- m3(
      resp_cats = if (v == "ss") c("corr", "other", "npl") else
        c("corr", "dist_context", "other", "dist_other", "npl"),
      num_options = if (v == "ss") c(1, 2, 3) else c(1, 2, 2, 2, 3),
      choice_rule = "softmax", version = v
    )$default_priors
    expect_identical(p$a$main, p$c$main)
  }
})

test_that("no activation effect prior is wider than the shared normal(0,0.5)", {
  for (v in c("ss", "cs")) {
    for (cr in c("simple", "softmax")) {
      p <- m3(
        resp_cats = if (v == "ss") c("corr", "other", "npl") else
          c("corr", "dist_context", "other", "dist_other", "npl"),
        num_options = if (v == "ss") c(1, 2, 3) else c(1, 2, 2, 2, 3),
        choice_rule = cr, version = v
      )$default_priors
      expect_identical(p$a$effects, "normal(0,0.5)")
      expect_identical(p$c$effects, "normal(0,0.5)")
    }
  }
})
