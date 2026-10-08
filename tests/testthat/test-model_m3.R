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

test_that("m3_custom refuses an activation symbol that is neither a column nor a parameter (#495)", {
  my_data <- data.frame(corr = c(5, 6), other = c(1, 2), npl = c(1, 2))
  my_model <- m3(c("corr", "other", "npl"), num_options = c(1, 2, 3))
  my_model$links <- list(a = "log", cstart = "log", cslope = "log", time = "log")
  formula <- bmf(
    corr ~ b + a + cstart + cslope * time,
    other ~ b + a,
    npl ~ b,
    a ~ 1,
    cstart ~ 1,
    cslope ~ 1
  )

  expect_error(
    bmm(formula, my_data, my_model, backend = "mock", mock_fit = 1, rename = FALSE),
    "'time' in your activation formula\\(s\\) is neither a data column nor a model parameter"
  )

  fit <- suppressWarnings(bmm(
    formula + bmf(time ~ 1), my_data, my_model,
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  expect_true("time" %in% names(fit$bmm$model$parameters))
})

test_that("m3_custom refuses an unknown symbol in an activation without formula parameters (#495)", {
  expect_error(
    bmm(
      bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + dd, npl ~ b, c ~ 1, a ~ 1),
      oberauer_lewandowsky_2019_e1,
      m3(c("corr", "other", "dist", "npl"), c(1, 4, 5, 5), links = list(c = "log", a = "log")),
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "'dd' in your activation formula\\(s\\) is neither a data column nor a model parameter"
  )

  err <- expect_error(suppressWarnings(bmm(
    bmf(corr ~ b + a + c, other ~ b + a, npl ~ b, c ~ 1 + cnd, a ~ 1),
    oberauer_lewandowsky_2019_e1,
    m3(c("corr", "other", "npl"), c(1, 4, 5), links = list(c = "log", a = "log")),
    backend = "mock", mock_fit = 1, rename = FALSE
  )), "'cnd'")
  expect_no_match(conditionMessage(err), "activation formula")
})

test_that("m3_custom activations can use numeric num_options by their names (#495)", {
  cats <- c("corr", "other", "dist", "npl")
  counts <- c(1, 4, 5, 5)
  dat <- as.data.frame(oberauer_lewandowsky_2019_e1)
  dat[c("n_corr", "n_other", "n_dist", "n_npl")] <- rep(counts, each = nrow(dat))
  fit_with <- function(count_name, num_options) {
    suppressWarnings(bmm(
      bmf(corr ~ b + a + c, other ~ b + a, npl ~ b, c ~ 1, a ~ 1, d ~ 1) +
        bmf(as.formula(paste("dist ~ b + d *", count_name))),
      dat,
      m3(cats, num_options, links = list(c = "log", a = "log", d = "log")),
      backend = "mock", mock_fit = 1, rename = FALSE
    ))
  }

  from_columns <- fit_with("n_dist", c("n_corr", "n_other", "n_dist", "n_npl"))
  unnamed <- fit_with("n_opt_dist", counts)
  user_named <- fit_with("k3", setNames(counts, c("k1", "k2", "k3", "k4")))
  category_named <- fit_with("n_opt_dist", setNames(counts, cats))
  expect_equal(brms::standata(unnamed), brms::standata(from_columns))
  expect_equal(brms::standata(user_named), brms::standata(from_columns))
  expect_equal(brms::standata(category_named), brms::standata(from_columns))
  expect_named(unnamed$bmm$model$parameters, c("b", "c", "a", "d"), ignore.order = TRUE)

  expect_error(
    fit_with("n_opt_dist", c("n_corr", "n_other", "n_dist", "n_npl")),
    "'n_opt_dist' in your activation formula\\(s\\) is neither a data column nor a model parameter"
  )
})

test_that("m3_custom linear activations can use nTrials and Idx_ columns (#495)", {
  fit_with <- function(activation) {
    suppressWarnings(bmm(
      bmf(corr ~ b + a + c, other ~ b + a, npl ~ b, c ~ 1, a ~ 1) + bmf(activation),
      oberauer_lewandowsky_2019_e1,
      m3(c("corr", "other", "dist", "npl"), c(1, 4, 5, 5), links = list(c = "log", a = "log")),
      backend = "mock", mock_fit = 1, rename = FALSE
    ))
  }

  for (activation in c(dist ~ b + nTrials, dist ~ b + Idx_dist)) {
    fit <- fit_with(activation)
    expect_named(fit$bmm$model$parameters, c("b", "c", "a"), ignore.order = TRUE)
    expect_true("C_dist_1" %in% names(brms::standata(fit)))
  }
  # Y is a matrix column, and as a predictor it breaks the Stan code
  expect_error(fit_with(dist ~ b + Y), "'Y' in your activation formula.*`Y` is reserved")
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
  expect_error(m3(c("corr", "other", "npl"), num_options = c(1, NA, 3)), "missing values")
})

test_that("num_options names already taken by a column or parameter give an error", {
  expect_error(m3_num_options_fit(c(a = 1, b = 4, c = 5, d = 5)), "'a', 'b', 'c', 'd'")
  expect_error(m3_num_options_fit(c(nTrials = 1, k2 = 4, k3 = 5, k4 = 5)), "'nTrials'")
  expect_error(m3_num_options_fit(c(Idx_dist = 1, k2 = 4, k3 = 5, k4 = 5)), "'Idx_dist'")
  expect_error(m3_num_options_fit(c(ID = 1, k2 = 4, k3 = 5, k4 = 5)), "'ID'")
})

m3_data_fit <- function(data) {
  bmm(
    bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b, c ~ 1, a ~ 1, d ~ 1),
    data,
    m3(
      resp_cats = c("corr", "other", "dist", "npl"),
      num_options = c("n_corr", "n_other", "n_dist", "n_npl"),
      choice_rule = "simple", links = list(c = "log", a = "log", d = "log"),
      default_priors = list(
        c = list(main = "normal(2, 0.5)", effects = "normal(0, 0.5)"),
        a = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)"),
        d = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)")
      )
    ),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
}

test_that("m3 warns about missing counts only where the category has options (#496)", {
  dat <- oberauer_lewandowsky_2019_e1
  expect_true(all(is.na(dat$dist) == (dat$n_dist == 0)))
  expect_silent(m3_data_fit(dat))

  dat$other[1:3] <- NA
  dat$dist[dat$n_dist > 0][1:2] <- NA
  expect_warning(m3_data_fit(dat), "contain 5 missing value\\(s\\)")
})

test_that("m3 refuses data columns it would overwrite (#496)", {
  for (col in c("nTrials", "Y", "Idx_corr")) {
    dat <- oberauer_lewandowsky_2019_e1
    dat[[col]] <- 1
    expect_error(m3_data_fit(dat), glue::glue("'{col}' would be overwritten"))
  }
})

test_that("a response category may be called Y or nTrials, also in the stored frame (#496)", {
  for (cat in c("Y", "nTrials")) {
    dat <- oberauer_lewandowsky_2019_e1
    names(dat)[names(dat) == "corr"] <- cat
    cats <- c(cat, "other", "npl")
    fit <- bmm(
      bmf(c ~ 1, a ~ 1), dat,
      m3(cats, num_options = c("n_corr", "n_other", "n_npl"), choice_rule = "simple", version = "ss"),
      backend = "mock", mock_fit = 1, rename = FALSE
    )
    stored <- check_stored_data(fit$bmm$model, fit$data, fit$bmm$user_formula)
    expect_equal(unname(stored$Y[, cats]), unname(as.matrix(dat[cats])))
  }
})

test_that("m3 refuses missing option counts (#496)", {
  dat <- oberauer_lewandowsky_2019_e1
  dat$n_other[1:3] <- NA
  dat$n_npl[1] <- NA
  expect_error(m3_data_fit(dat), "missing values: n_other \\(3 NA\\), n_npl \\(1 NA\\)")
})

test_that("m3 rejects num_options it cannot map onto the response categories", {
  cats <- c("corr", "other", "npl")
  expect_error(m3(cats, num_options = c(corr = 1, 4, 5)), "all elements")
  expect_error(m3(cats, num_options = c(k = 1, k = 4, j = 5)), "only once")
  expect_error(m3(cats, num_options = c(corr = 1, other = 4, dist = 5)), "one element for each")
  expect_error(m3(cats, num_options = cats), "response category column")
})

test_that("softmax and gaussian default priors give a and c equal main means (c - a centered at 0)", {
  for (v in c("ss", "cs")) {
    for (cr in c("softmax", "gaussian")) {
      p <- m3(
        resp_cats = if (v == "ss") c("corr", "other", "npl") else
          c("corr", "dist_context", "other", "dist_other", "npl"),
        num_options = if (v == "ss") c(1, 2, 3) else c(1, 2, 2, 2, 3),
        choice_rule = cr, version = v
      )$default_priors
      expect_identical(p$a$main, p$c$main)
    }
  }
})

test_that("no activation effect prior is wider than the shared normal(0,0.5)", {
  for (v in c("ss", "cs")) {
    for (cr in c("simple", "softmax", "gaussian")) {
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

test_that("m3 runs through the pipeline with the gaussian choice rule in all versions", {
  ss <- m3(c("corr", "other", "npl"), c("n_corr", "n_other", "n_npl"),
           choice_rule = "gaussian", version = "ss")
  expect_silent(bmm(bmf(c ~ 1 + cond + (1 | ID), a ~ 1 + cond + (1 | ID)),
                    oberauer_lewandowsky_2019_e1, ss, backend = "mock", mock_fit = 1, rename = FALSE))

  cs <- m3(c("corr", "dc", "other", "do", "npl"), c(1, 1, 4, 4, 10),
           choice_rule = "gaussian", version = "cs")
  dat <- data.frame(corr = c(50, 40), dc = c(5, 8), other = c(25, 30), do = c(10, 12), npl = c(10, 10))
  expect_silent(bmm(bmf(c ~ 1, a ~ 1, f ~ 1), dat, cs, backend = "mock", mock_fit = 1, rename = FALSE))

  custom <- m3(c("corr", "other", "dist", "npl"), c("n_corr", "n_other", "n_dist", "n_npl"),
               choice_rule = "gaussian")
  custom$links <- list(c = "identity", a = "identity", d = "identity")
  f <- bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b, c ~ 1 + cond, a ~ 1 + cond, d ~ 1)
  expect_warning(bmm(f, oberauer_lewandowsky_2019_e1, custom, backend = "mock", mock_fit = 1, rename = FALSE),
                 "Default priors for each parameter")
})

test_that("the gaussian choice rule fixes b at 0, accepts any case, and puts log P on each category", {
  model <- m3(c("corr", "other", "npl"), c(1, 4, 10), choice_rule = "Gaussian", version = "ss")
  expect_identical(model$other_vars$choice_rule, "gaussian")
  expect_identical(model$fixed_parameters$b, 0)
  expect_identical(model$default_priors$c$main, "normal(2,1)")

  code <- stancode(bmf(c ~ 1, a ~ 1), data.frame(corr = 50, other = 30, npl = 20), model)
  expect_match(code, "real m3_gauss_logp(int k, real A1, real A2, real A3, real n1, real n2, real n3)", fixed = TRUE)
  for (k in 1:3) expect_match(code, paste0("m3_gauss_logp(", k, " , nlp_corr[n] , nlp_other[n] , nlp_npl[n]"), fixed = TRUE)
  expect_false(grepl("multinomial_logit_lpmf(Y[n] | log", code, fixed = TRUE))
  expect_error(m3(c("corr", "other", "npl"), c(1, 4, 10), choice_rule = "probit"), "gaussian")
})

test_that("the gaussian rule refuses fractional option counts", {
  cats <- c("corr", "other", "npl")
  opts <- c("n_corr", "n_other", "n_npl")
  f <- bmf(c ~ 1, a ~ 1)
  dat <- data.frame(corr = c(50, 40), other = c(20, 30), npl = c(30, 30),
                    n_corr = 1, n_other = c(0.3, 4), n_npl = 10)
  gauss <- m3(cats, opts, choice_rule = "gaussian", version = "ss")
  expect_error(check_data(gauss, dat, f), "whole numbers")
  expect_s3_class(check_data(m3(cats, opts, choice_rule = "softmax", version = "ss"), dat, f), "data.frame")
  dat$n_other[1] <- 0
  expect_silent(check_data(gauss, dat, f))
})

test_that("custom m3 with the gaussian rule centres default identity-link priors on the Gaussian scale", {
  model <- m3(c("corr", "other", "npl"), c(1, 4, 10), choice_rule = "gaussian")
  model$links <- list(c = "identity", a = "identity")
  f <- bmf(corr ~ b + a + c, other ~ b + a, npl ~ b, c ~ 1, a ~ 1)
  model <- suppressWarnings(check_model(model, data.frame(corr = 1, other = 1, npl = 1), f))
  expect_identical(model$default_priors$c$main, "normal(2, 1)")
  expect_identical(model$default_priors$a$main, "normal(2, 1)")
})

# 40 quadrature nodes hold |delta log P| below 2e-5 where log P > -5, and below
# 1.2e-4 for categories with 30 options
test_that("Gaussian-rule probabilities match closed forms and integrate()", {
  # two categories with one option each: the difference of two N(0, 1) draws has variance 2
  A <- c(-1, 0, 0.7, 2.5)
  for (a in A) {
    expect_lt(abs(m3_gauss_logp(1, a, 0, 1, 1) - stats::pnorm(a / sqrt(2), log.p = TRUE)), 2e-5)
  }
  # equal activations: every option is equally likely to win, so P(k) = n_k / sum(n)
  n <- c(1, 4, 10, 3)
  lp <- vapply(1:4, function(k) do.call(m3_gauss_logp, as.list(c(k, rep(0.8, 4), n))), numeric(1))
  expect_lt(max(abs(lp - log(n / sum(n)))), 2e-5)

  ref_logp <- function(k, A, n) {
    others <- setdiff(seq_along(A), k)
    f <- function(z) {
      out <- stats::dnorm(z) * stats::pnorm(z)^(n[k] - 1)
      for (j in others) out <- out * stats::pnorm(z + A[k] - A[j])^n[j]
      out
    }
    log(n[k]) + log(stats::integrate(f, -Inf, Inf, rel.tol = 1e-12)$value)
  }
  cells <- list(list(A = c(3, 1.2, 0), n = c(1, 4, 10)), list(A = c(1, 1.5, 0), n = c(1, 8, 30)),
                list(A = c(2, 1.5, 0.8, 1, 0), n = c(1, 1, 4, 4, 10)))
  for (cl in cells) {
    for (k in seq_along(cl$A)) {
      lp <- do.call(m3_gauss_logp, as.list(c(k, cl$A, cl$n)))
      expect_lt(abs(lp - ref_logp(k, cl$A, cl$n)), if (max(cl$n) >= 30) 1.2e-4 else 2e-5)
    }
  }
})

test_that("m3_gauss_logp keeps the draws-by-observation shape and gives absent categories -100", {
  A <- matrix(c(0.5, 1, 2, 3, 1.5, 0), 2)
  one <- matrix(1, 2, 3)
  out <- m3_gauss_logp(1, A, A * 0, A * 0 - 1, one, one * 4, one * 1e-4)
  expect_identical(dim(out), c(2L, 3L))
  expect_true(all(out < 0))
  expect_identical(as.vector(m3_gauss_logp(3, A, A * 0, A * 0, one, one, one * 1e-4)), rep(-100, 6))
})

test_that("the Stan Gaussian-rule wrapper matches its R companion", {
  skip_on_cran() # compiles and runs a Stan program
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  model <- m3(c("corr", "other", "npl"), c(1, 4, 10), choice_rule = "gaussian", version = "ss")
  grid <- expand.grid(k = 1:3, A1 = c(-2, 0.5, 4, 40), A2 = c(0, 1.5), n2 = c(1e-4, 1, 8), n3 = c(1, 30))
  grid$A3 <- 0
  grid$n1 <- 1
  program <- paste0(
    "functions {\n", m3_gaussian_stanvars(model)[[1]]$scode, "\n}\n",
    "data {\n  int N; array[N] int k; vector[N] A1; vector[N] A2; vector[N] A3;\n",
    "  vector[N] n1; vector[N] n2; vector[N] n3;\n}\n",
    "generated quantities {\n  vector[N] lp;\n  for (i in 1:N)\n",
    "    lp[i] = m3_gauss_logp(k[i], A1[i], A2[i], A3[i], n1[i], n2[i], n3[i]);\n}\n"
  )
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(list(N = nrow(grid)), as.list(grid)), fixed_param = TRUE, chains = 1,
    iter_sampling = 1, iter_warmup = 0, refresh = 0, show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  stan <- as.numeric(csv[1, paste0("lp.", seq_len(nrow(grid)))])
  r <- vapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    m3_gauss_logp(g$k, g$A1, g$A2, g$A3, g$n1, g$n2, g$n3)
  }, numeric(1))
  expect_true(all(is.finite(stan)))
  expect_lt(max(abs(stan - r)), 1e-6)
})
