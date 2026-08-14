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
        list(main = softplus_main[[rule]], effects = "normal(0, 0.5)")
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
