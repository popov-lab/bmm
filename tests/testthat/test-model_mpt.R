test_that("mpt_tree constructs and prints a tree", {
  tree <- mpt_tree("old", list(
    old = "D + (1 - D) * g",
    new = "(1 - D) * (1 - g)"
  ))
  expect_s3_class(tree, "mpt_tree")
  expect_equal(tree$name, "old")
  expect_equal(names(tree$branches), c("old", "new"))
  expect_output(print(tree), "P\\(old\\) = D")
})

test_that("mpt_tree validates its inputs", {
  expect_error(mpt_tree("t", list("D + g")), "named after its response category")
  expect_error(
    mpt_tree("t", list(a = "D + g", a = "1 - D - g")),
    "unique within a tree"
  )
  expect_error(mpt_tree("t", list(a = "D + * g", b = "x")), "Cannot parse")
  expect_error(mpt_tree("t", list(a = 0.5, b = "x")), "character strings")
})

test_that("mpt stores its derived state once", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  expect_equal(model$other_vars$link, "logit")
  expect_equal(model$other_vars$indicators$tree, c(withdist = "Idx_withdist", nodist = "Idx_nodist"))
  expect_equal(model$other_vars$indicators$possible, c(dist = "Poss_dist"))
  expect_null(model$other_vars$simplex_raw)
  expect_setequal(names(model$parameters), c("Pm", "Pb"))
  expect_equal(model$links$Pm, "logit")
  expect_equal(model$default_priors$Pm$main, "logistic(0, 1)")
  expect_equal(model$default_priors$Pm$effects, "logistic(0, 1)")

  single <- mpt(mpt_tree("t", list(A = "gA", B = "gB", C = "gC")),
    simplex = c("gA", "gB", "gC"), links = "probit"
  )
  expect_null(single$other_vars$indicators$tree)
  expect_null(single$other_vars$indicators$possible)
  expect_equal(single$other_vars$simplex_raw, c(gA = "gAraw", gB = "gBraw"))
  expect_equal(single$links$gA, "identity")
  expect_equal(single$default_priors$gAraw$main, "normal(0, 1)")
  expect_equal(single$default_priors$gAraw$effects, "normal(0, 1)")

  covariate_tree <- mpt_tree("main", list(
    correct = "Pb + (1 - Pb) * Pi * GcorrPi",
    other = "(1 - Pb) * Pi * (1 - GcorrPi) + (1 - Pb) * (1 - Pi)"
  ))
  with_covariate <- mpt(covariate_tree, covariates = "GcorrPi")
  expect_setequal(names(with_covariate$parameters), c("Pb", "Pi"))
  expect_equal(with_covariate$other_vars$covariates, "GcorrPi")
})

test_that("an empty formula fits every parameter with an intercept", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  fit <- suppressMessages(bmm(
    bmf(), mpt_2htm_data(), model, backend = "mock", mock_fit = 1, rename = FALSE
  ))
  expect_s3_class(fit, "bmmfit")
  expect_setequal(names(fit$bmm$model$parameters), c("D", "g"))
})

test_that("update() re-runs the mpt data preparation for new data", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat10 <- mpt_2htm_data(n_id = 10)
  dat20 <- mpt_2htm_data(n_id = 20)

  # a mock stanfit with just enough structure for brms::update.brmsfit to
  # reach its data-revalidation path without any compilation or sampling
  methods::setClass("bmm_mock_stanfit", representation(sim = "list"))
  mockfit <- methods::new("bmm_mock_stanfit", sim = list(
    warmup = 1000, iter = 2000, chains = 1, thin = 1,
    samples = list(structure(list(), args = list(control = list())))
  ))
  fit <- suppressMessages(bmm(
    bmf(D ~ 1, g ~ 1), dat10, model,
    backend = "mock", mock_fit = mockfit, rename = FALSE
  ))
  expect_equal(nrow(fit$data), nrow(dat10))
  expect_true(all(c("Idx_old", "Idx_new") %in% names(fit$data)))

  fit2 <- suppressMessages(
    update(fit, newdata = dat20, testmode = TRUE, recompile = FALSE)
  )
  expect_equal(nrow(fit2$data), nrow(dat20))
  expect_equal(sum(fit2$data$Idx_old), nrow(dat20) / 2)
})

test_that("mpt errors on invalid parameter and category names", {
  bad_par <- mpt_tree("t", list(a = "d_A + g", b = "1 - d_A - g"))
  expect_error(mpt(bad_par), "underscores or dots")
  bad_cat <- mpt_tree("t", list(cat_a = "D", cat_b = "1 - D"))
  expect_error(mpt(bad_cat), "underscores or dots")
})

test_that("mpt errors on inconsistent trees and a missing tree_id", {
  trees <- list(
    mpt_tree("t1", list(a = "p", b = "1 - p")),
    mpt_tree("t2", list(a = "p", c = "1 - p"))
  )
  expect_error(mpt(trees, tree_id = "cond"), "same response categories")
  expect_error(mpt(mpt_2htm_trees()), "require the tree_id argument")
  dup_trees <- list(
    mpt_tree("t1", list(a = "p", b = "1 - p")),
    mpt_tree("t1", list(a = "p", b = "1 - p"))
  )
  expect_error(mpt(dup_trees, tree_id = "cond"), "unique")
})

test_that("mpt refuses a zero branch and points to impossible categories", {
  zero_padded <- list(
    mpt_tree("a", list(x = "D", y = "1 - D", z = "0")),
    mpt_tree("b", list(x = "1 - 1", y = "g", z = "1 - g"))
  )
  expect_error(mpt(zero_padded, tree_id = "t"), "constant 0")
  expect_error(mpt(zero_padded, tree_id = "t"), "'z' in tree 'a', 'x' in tree 'b'")
  expect_error(mpt(zero_padded, tree_id = "t"), "mpt_tree(impossible = )", fixed = TRUE)
  expect_error(mpt(mpt_tree("c", list(x = "(0)", y = "D + (1 - D)"))), "constant 0")
  # not the literal 0 this guard reads, but 0 at every interior test point,
  # which the range check refuses
  zero_product <- mpt_tree("d", list(x = "0 * D", y = "1 - 0 * D"))
  expect_error(mpt(zero_product), "category 'x' in tree 'd' is 0 .*outside \\(0, 1\\]")

  declared <- list(
    mpt_tree("a", list(x = "D", y = "1 - D"), impossible = "z"),
    mpt_tree("b", list(y = "g", z = "1 - g"), impossible = "x")
  )
  model <- mpt(declared, tree_id = "t")
  expect_equal(model$other_vars$indicators$possible, c(z = "Poss_z", x = "Poss_x"))
})

test_that("mpt errors on name collisions and reserved names", {
  tree_collision <- mpt_tree("t", list(D = "D + g", other = "1 - D - g"))
  expect_error(mpt(tree_collision), "both a parameter and a response category")
  tree_reserved <- mpt_tree("t", list(Y = "p", other = "1 - p"))
  expect_error(mpt(tree_reserved), "reserved")
})

test_that("mpt validates simplex groups", {
  trees <- mpt_2htm_trees()
  expect_error(
    mpt(trees, tree_id = "item_type", simplex = c("g", "x")),
    "Unknown"
  )
  expect_error(
    mpt(trees, tree_id = "item_type", simplex = "g"),
    "at least two parameters"
  )
  three <- mpt_tree("t", list(A = "a", B = "b", C = "c"))
  expect_error(
    mpt(three, simplex = c("a", "b", "c"), restrictions = "a = b"),
    "'a' are also restricted, and the restriction removed them"
  )
  expect_error(mpt(three), "belong in the simplex argument")
})

test_that("the links of simplex parameters cannot be changed", {
  tree <- mpt_tree("src", list(
    a = "gA", b = "gB", n = "gN"
  ))
  model <- mpt(tree, simplex = c("gA", "gB", "gN"))
  dat <- data.frame(a = 10, b = 5, n = 5)
  model$links$gA <- "logit"
  expect_error(
    suppressMessages(check_model(model, dat, bmf(gA ~ 1))),
    "cannot be changed"
  )
})

test_that("check_formula generates linked category formulas", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  formula <- bmf(D ~ 1 + (1 | id), g ~ 1)
  model_checked <- check_model(model, dat, formula)
  dat_checked <- check_data(model_checked, dat, formula)
  formula_checked <- check_formula(model_checked, dat_checked, formula)

  expect_setequal(names(formula_checked), c("old", "new", "D", "g"))
  old_rhs <- paste(deparse(formula_checked$old[[3]]), collapse = " ")
  expect_true(grepl("Idx_old", old_rhs, fixed = TRUE))
  expect_true(grepl("inv_logit(D)", old_rhs, fixed = TRUE))
  expect_true(all(is_nl(formula_checked)[c("old", "new")]))
})

test_that("check_formula uses Phi for the probit link", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type", links = "probit")
  dat <- mpt_2htm_data()
  formula <- bmf(D ~ 1, g ~ 1)
  model_checked <- check_model(model, dat, formula)
  dat_checked <- check_data(model_checked, dat, formula)
  formula_checked <- check_formula(model_checked, dat_checked, formula)
  old_rhs <- paste(deparse(formula_checked$old[[3]]), collapse = " ")
  expect_true(grepl("Phi(D)", old_rhs, fixed = TRUE))
  expect_false(grepl("pnorm", old_rhs, fixed = TRUE))
})

test_that("formulas for response categories are rejected", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  formula <- bmf(old ~ 1, D ~ 1, g ~ 1)
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE),
    "cannot be predicted directly"
  )
})

test_that("mpt compiles for a multi-tree binary-category model", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  formula <- bmf(D ~ 1 + (1 | id), g ~ 1)
  expect_silent(bmm(
    formula = formula, data = dat, model = model,
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
})

test_that("mpt compiles for a single-tree multinomial model", {
  pcm_tree <- mpt_tree("study", list(
    C = "cp + (1 - cp) * rp * rp",
    E = "2 * (1 - cp) * rp * (1 - rp)",
    U = "(1 - cp) * (1 - rp) * (1 - rp)"
  ))
  model <- mpt(pcm_tree)
  dat <- data.frame(id = factor(1:10))
  counts <- t(rmultinom(10, 40, c(0.5, 0.3, 0.2)))
  colnames(counts) <- c("C", "E", "U")
  dat <- cbind(dat, counts)
  formula <- bmf(cp ~ 1 + (1 | id), rp ~ 1 + (1 | id))
  expect_silent(bmm(
    formula = formula, data = dat, model = model,
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
})

test_that("mpt compiles with design-fixed covariates", {
  tree <- mpt_tree("main", list(
    correct = "D + (1 - D) * Gcorr",
    incorrect = "(1 - D) * (1 - Gcorr)"
  ))
  model <- mpt(tree, covariates = "Gcorr")
  dat <- data.frame(
    id = factor(1:10), Gcorr = 0.25,
    correct = rbinom(10, 40, 0.7), incorrect = 0
  )
  dat$incorrect <- 40 - dat$correct
  expect_silent(bmm(
    bmf(D ~ 1 + (1 | id)), dat, model,
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
})

test_that("generated indicator and covariate columns do not trigger the clash warning", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  expect_no_warning(bmm(
    bmf(D ~ 1 + (1 | id), g ~ 1), dat, model,
    backend = "mock", mock_fit = 1, rename = FALSE
  ))

  cov_model <- mpt(
    mpt_tree("main", list(
      correct = "D + (1 - D) * Gcorr",
      incorrect = "(1 - D) * (1 - Gcorr)"
    )),
    covariates = "Gcorr"
  )
  cov_dat <- data.frame(
    id = factor(1:10), Gcorr = 0.25, correct = 28L, incorrect = 12L
  )
  expect_no_warning(bmm(
    bmf(D ~ 1 + (1 | id)), cov_dat, cov_model,
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
})

test_that("a fixed non-linear sub-parameter reaches the constant prior", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  dat$ptime <- rep(c(1, 2), length.out = nrow(dat))
  formula <- bmf(
    D ~ inv_logit(Dmax) * (1 - exp(-exp(rate) * ptime)),
    Dmax = 1, rate ~ 1, g ~ 1
  )
  prior <- suppressWarnings(suppressMessages(default_prior(formula, dat, model)))
  constant_row <- prior[grepl("constant", prior$prior), ]
  expect_equal(constant_row$nlpar, "Dmax")
  expect_equal(constant_row$prior, "constant(1)")
})

test_that("non-linear parameter formulas bypass the link transformation", {
  tree <- mpt_tree("main", list(
    correct = "D + (1 - D) * 0.25",
    incorrect = "(1 - D) * 0.75"
  ))
  model <- mpt(tree)
  dat <- data.frame(
    id = factor(1:10), ptime = rep(c(0.5, 1), 5),
    correct = rbinom(10, 40, 0.6)
  )
  dat$incorrect <- 40 - dat$correct
  formula <- bmf(
    D ~ inv_logit(Dmax) * (1 - exp(-exp(rate) * ptime)),
    Dmax ~ 1,
    rate ~ 1
  )
  model_checked <- suppressMessages(check_model(model, dat, formula))
  expect_equal(model_checked$links$D, "identity")
  expect_setequal(
    names(model_checked$parameters), c("D", "Dmax", "rate")
  )
  expect_equal(model_checked$links$Dmax, "identity")
  expect_equal(model_checked$default_priors$Dmax$main, "normal(0, 1)")
  expect_null(model_checked$default_priors[["D"]])

  dat_checked <- check_data(model_checked, dat, formula)
  formula_checked <- suppressMessages(
    check_formula(model_checked, dat_checked, formula)
  )
  correct_rhs <- paste(deparse(formula_checked$correct[[3]]), collapse = " ")
  expect_false(grepl("inv_logit(D)", correct_rhs, fixed = TRUE))
})

test_that("the message on non-linear sub-parameters names the sd prior they get", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  dat$x <- rep(c(0, 1), length.out = nrow(dat))
  formula <- bmf(D ~ inv_logit(a + b * x), a ~ 1 + (1 | id), b ~ 1, g ~ 1)
  priors <- suppressWarnings(suppressMessages(default_prior(formula, dat, model)))
  sd_prior <- priors$prior[priors$class == "sd" & priors$nlpar == "a" & priors$group == ""]
  expect_length(sd_prior, 1)
  expect_message(check_model(model, dat, formula), sd_prior, fixed = TRUE)
})

test_that("a symbol of a non-linear formula that is neither column nor parameter errors", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data(n_id = 4)
  formula <- bmf(D ~ inv_logit(a + b * x), a ~ 1, b ~ 1, g ~ 1)
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE),
    "'x' .* neither a data column nor a model parameter"
  )
  # the implicit intercept of a sub-parameter without a formula is gone too
  expect_error(
    check_model(model, dat, bmf(D ~ inv_logit(a + c), a ~ 1, g ~ 1)),
    "neither a data column"
  )
  dat$x <- rep(c(0, 1), length.out = nrow(dat))
  fit <- suppressMessages(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
  expect_setequal(names(fit$bmm$model$parameters), c("D", "g", "a", "b"))
})

test_that("links set after construction are checked like any other model's", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  formula <- bmf(D ~ 1, g ~ 1)

  model$links$D <- "cloglog"
  expect_error(check_model(model, dat, formula), "Unknown link function")
  model$links$D <- "log"
  expect_error(check_model(model, dat, formula), "Unknown link function")
})

test_that("model_docs() names no scale that contradicts a switched link", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  model$links$D <- "probit"
  expect_no_match(paste(model_docs(model), collapse = "\n"), "logit scale")
})

test_that("a per-parameter probit link gets the probit default prior", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  model$links$D <- "probit"
  dat <- mpt_2htm_data()
  prior <- suppressMessages(default_prior(bmf(D ~ 1, g ~ 1), dat, model))
  intercept <- prior[prior$coef == "Intercept", ]
  expect_equal(intercept$prior[intercept$nlpar == "D"], "normal(0, 1)")
  expect_equal(intercept$prior[intercept$nlpar == "g"], "logistic(0, 1)")
})

test_that("random-effects SDs get a default prior on the scale of the link", {
  sd_prior <- function(model, par) {
    pr <- suppressMessages(default_prior(
      bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id)), mpt_2htm_data(), model
    ))
    pr$prior[pr$class == "sd" & pr$coef == "" & pr$group == "" & pr$nlpar == par]
  }
  logit <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  probit <- mpt(mpt_2htm_trees(), tree_id = "item_type", links = "probit")
  for (par in c("D", "g")) {
    expect_equal(sd_prior(logit, par), "exponential(1)")
    expect_equal(sd_prior(probit, par), "exponential(2)")
  }
})

test_that("condition effects under reference coding get the intercept's prior scale", {
  dat <- mpt_2htm_data()
  dat$cond <- rep(c("a", "b"), length.out = nrow(dat))
  effect_prior <- function(model) {
    pr <- suppressMessages(default_prior(bmf(D ~ 1 + cond, g ~ 1), dat, model))
    pr$prior[pr$class == "b" & pr$coef == "" & pr$nlpar == "D"]
  }
  expect_equal(effect_prior(mpt(mpt_2htm_trees(), tree_id = "item_type")), "logistic(0, 1)")
  expect_equal(
    effect_prior(mpt(mpt_2htm_trees(), tree_id = "item_type", links = "probit")),
    "normal(0, 1)"
  )
})

test_that("the links a non-linear formula switches off survive a second check", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  dat$ptime <- rep(c(1, 2), length.out = nrow(dat))
  formula <- bmf(
    D ~ inv_logit(Dmax) * (1 - exp(-exp(rate) * ptime)),
    Dmax ~ 1, rate ~ 1, g ~ 1
  )
  checked <- suppressMessages(check_model(model, dat, formula))
  rechecked <- suppressMessages(check_model(checked, dat, formula))
  expect_equal(rechecked$links, checked$links)
})

test_that("a parameter made linear again in a re-check gets back its link and prior", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  dat$x <- rep(c(0, 1), length.out = nrow(dat))
  dat$cond <- rep(c("a", "b"), each = 2, length.out = nrow(dat))
  nl_formula <- bmf(D ~ inv_logit(a + b * x), a ~ 1, b ~ 1, g ~ 1)
  linear_formula <- bmf(D ~ 1 + cond, g ~ 1)

  checked <- suppressMessages(check_model(model, dat, nl_formula))
  rechecked <- check_model(checked, dat, linear_formula)
  expect_setequal(names(rechecked$parameters), c("D", "g"))
  expect_equal(rechecked$links$D, "logit")
  expect_identical(rechecked$default_priors$D, .mpt_latent_prior("logit"))

  model$links$D <- "probit"
  checked <- suppressMessages(check_model(model, dat, nl_formula))
  rechecked <- check_model(checked, dat, linear_formula)
  expect_equal(rechecked$links$D, "probit")
  expect_identical(rechecked$default_priors$D, .mpt_latent_prior("probit"))
})

test_that("a stick-breaking component made linear again in a re-check gets back its prior", {
  tree <- mpt_tree("t", list(A = "D * gA", B = "D * gB", C = "D * gC", N = "1 - D"))
  model <- mpt(tree, simplex = c("gA", "gB", "gC"))
  dat <- data.frame(id = factor(1:8), x = rep(0:1, 4), A = 10L, B = 10L, C = 10L, N = 10L)

  checked <- suppressMessages(check_model(model, dat, bmf(gAraw ~ a + b * x, a ~ 1, b ~ 1)))
  rechecked <- check_model(checked, dat, bmf(gAraw ~ 1 + x))
  expect_equal(rechecked$links$gAraw, "identity")
  expect_identical(rechecked$default_priors$gAraw, .mpt_latent_prior("logit"))
})

test_that("update() to a linear formula restores the link and prior of a non-linear parameter", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  dat$x <- rep(c(0, 1), length.out = nrow(dat))
  dat$cond <- rep(c("a", "b"), each = 2, length.out = nrow(dat))
  methods::setClass("bmm_mock_stanfit", representation(sim = "list"))
  mockfit <- methods::new("bmm_mock_stanfit", sim = list(
    warmup = 1000, iter = 2000, chains = 1, thin = 1,
    samples = list(structure(list(), args = list(control = list())))
  ))
  fit <- suppressWarnings(suppressMessages(bmm(
    bmf(D ~ inv_logit(a + b * x), a ~ 1, b ~ 1, g ~ 1), dat, model,
    backend = "mock", mock_fit = mockfit, rename = FALSE
  )))

  up <- suppressMessages(update(
    fit, formula. = bmf(D ~ 1 + cond, g ~ 1), newdata = dat,
    testmode = TRUE, recompile = FALSE
  ))
  expect_setequal(names(up$bmm$model$parameters), c("D", "g"))
  expect_equal(up$bmm$model$links$D, "logit")
  intercept <- up$prior$nlpar == "D" & up$prior$coef == "Intercept"
  expect_equal(up$prior$prior[intercept], .mpt_latent_prior("logit")$main)
})

test_that("printing an mpt model lists trees, restrictions and the identifiability bound", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  expect_output(print(model), "MPT tree 'old':")
  expect_output(print(model), "P\\(new\\) = \\(1 - D\\) \\* \\(1 - g\\)")
  expect_output(
    print(model),
    "2 free parameter\\(s\\), 2 degrees of freedom"
  )
  expect_no_match(capture.output(print(model)), "not identified")

  trees <- list(
    mpt_tree("old", list(yes = "Do + (1 - Do) * g", no = "(1 - Do) * (1 - g)")),
    mpt_tree("new", list(yes = "(1 - Dn) * g", no = "Dn + (1 - Dn) * (1 - g)"))
  )
  unidentified <- mpt(trees, tree_id = "item_type")
  expect_output(print(unidentified), "3 free parameter\\(s\\), 2 degrees of freedom")
  expect_output(print(unidentified), "not identified")

  dat <- mpt_2htm_data()
  dat$x <- rep(0:1, length.out = nrow(dat))
  nonlinear <- suppressMessages(check_model(
    model, dat, bmf(D ~ inv_logit(a + b * x), a ~ 1, b ~ 1, g ~ 1)
  ))
  expect_output(print(nonlinear), "2 free parameter\\(s\\)")
  expect_no_match(capture.output(print(nonlinear)), "not identified")

  fixed <- suppressMessages(check_model(model, dat, bmf(D ~ 1, g = 0.5)))
  expect_output(print(fixed), "1 free parameter\\(s\\)")

  restricted <- mpt(trees, tree_id = "item_type", restrictions = c("Dn = Do", "g = 0.5"))
  expect_output(print(restricted), "Restrictions: Dn = Do, g = 0.5")

  simplex <- mpt(
    mpt_tree("t", list(A = "gA", B = "gB", C = "gC")),
    simplex = c("gA", "gB", "gC")
  )
  expect_output(print(simplex), "2 free parameter\\(s\\), 2 degrees of freedom")
  expect_output(print(simplex), "'gAraw', 'gBraw' pass through inv_logit\\(\\)")
  expect_match(
    parameter_info(simplex)$description[4],
    "transformed by the inverse logit into the stick proportion"
  )
})

test_that("restrictions are substituted into the trees before parameters are identified", {
  trees <- list(
    mpt_tree("old", list(yes = "Do + (1 - Do) * g", no = "(1 - Do) * (1 - g)")),
    mpt_tree("new", list(yes = "(1 - Dn) * g", no = "Dn + (1 - Dn) * (1 - g)"))
  )
  equated <- mpt(trees, tree_id = "item_type", restrictions = "Dn = Do")
  expect_setequal(names(equated$parameters), c("Do", "g"))
  expect_equal(equated$other_vars$trees$new$branches$yes, quote((1 - Do) * g))
  expect_equal(equated$other_vars$restrictions, list(Dn = quote(Do)))

  fixed <- mpt(trees, tree_id = "item_type", restrictions = c("g = 0.5", "Dn = Do"))
  expect_equal(names(fixed$parameters), "Do")
  expect_equal(
    dmpt(c(yes = 6, no = 4), pars = c(Do = 0.2), mpt_model = fixed, tree = "old", log = FALSE),
    dmultinom(c(6, 4), prob = c(0.2 + 0.8 * 0.5, 0.8 * 0.5))
  )
  expect_equal(
    mpt(trees, tree_id = "item_type", restrictions = list(g = 0.5, Dn = "Do"))$other_vars$trees,
    fixed$other_vars$trees
  )
  for (restricted in list(equated, fixed)) {
    branch_errors <- .mpt_tree_branch_errors(
      restricted$other_vars$trees, names(restricted$parameters), list()
    )
    expect_true(all(is.na(branch_errors)))
  }
  one_symbol <- mpt_tree("t", list(A = "a", B = "(1 - a) * b", C = "(1 - a) * (1 - b)"))
  expect_error(mpt(one_symbol, restrictions = "a = 0"), "strictly between 0 and 1")
  for (r in c("g = 1.5", "g = -0.2", "g = 1")) {
    expect_error(mpt(trees, tree_id = "item_type", restrictions = r), "strictly between 0 and 1")
  }
  expect_error(
    mpt(trees, tree_id = "item_type", restrictions = list(g = 1)),
    "strictly between 0 and 1"
  )
  expect_equal(
    mpt(trees, tree_id = "item_type", restrictions = list("Dn = Do", "g = 0.5"))$other_vars$trees,
    mpt(trees, tree_id = "item_type", restrictions = c("Dn = Do", "g = 0.5"))$other_vars$trees
  )
  expect_error(
    mpt(mpt_tree("t", list(A = "a", B = "1 - a")), restrictions = "a = 0.3"),
    "no latent parameters"
  )
  expect_error(
    mpt(trees, tree_id = "item_type", restrictions = "Do = Dn = FE"),
    "person effects.*random effects"
  )

  expect_error(
    mpt(trees, tree_id = "item_type", restrictions = "Do > Dn"),
    "Order constraints"
  )
  expect_error(
    mpt(trees, tree_id = "item_type", restrictions = "Dn = Dold"),
    "do not appear"
  )
  expect_error(
    mpt(trees, tree_id = "item_type", restrictions = c("Dn = Do", "Do = Dn")),
    "circular"
  )
  expect_error(
    mpt(trees, tree_id = "item_type", restrictions = c("g = 0.5", "g = 0.4")),
    "more than once"
  )
  expect_error(
    mpt(trees, tree_id = "item_type", restrictions = "g = 0.00001"),
    "scientific"
  )

  cov_tree <- mpt_tree("t", list(a = "D + (1 - D) * gc", b = "(1 - D) * (1 - gc)"))
  expect_error(
    mpt(cov_tree, covariates = "gc", restrictions = "gc = 0.5"),
    "Covariates"
  )
})

mpt_printed <- function(model) {
  gsub("\\s+", " ", paste(capture.output(print(model)), collapse = " "))
}

test_that("print() names parameters that only enter as a product although the count passes", {
  # D and r enter only as D * r, so the four parameters fill four degrees of
  # freedom but leave one combination open
  model <- mpt(list(
    mpt_tree("a", list(x = "D * r", y = "(1 - D * r) * g", z = "(1 - D * r) * (1 - g)")),
    mpt_tree("b", list(x = "g", y = "(1 - g) * h", z = "(1 - g) * (1 - h)"))
  ), tree_id = "tree")
  printed <- mpt_printed(model)
  expect_match(printed, "4 free parameter\\(s\\), 4 degrees of freedom")
  expect_match(printed, "not identified: .* has rank 3 for 4 free parameters")
  expect_match(printed, "1 combination\\(s\\) of 'D', 'r' cannot be estimated")
})

test_that("print() reports full Jacobian rank for identified models", {
  expect_match(
    mpt_printed(mpt(mpt_2htm_trees(), tree_id = "item_type")),
    "Jacobian rank 2 of 2 at interior test values: locally identified"
  )

  # Do and Dn are identified once the guessing rate varies across two bias trees
  bias_trees <- lapply(1:2, function(i) {
    list(
      mpt_tree(paste0("old", i), list(
        yes = glue("Do + (1 - Do) * g{i}"), no = glue("(1 - Do) * (1 - g{i})")
      )),
      mpt_tree(paste0("new", i), list(
        yes = glue("(1 - Dn) * g{i}"), no = glue("Dn + (1 - Dn) * (1 - g{i})")
      ))
    )
  })
  bias <- mpt(unlist(bias_trees, recursive = FALSE), tree_id = "tree")
  expect_match(mpt_printed(bias), "Jacobian rank 4 of 4 at interior test values")

  # pair-clustering model, pairs tree
  pair_clustering <- mpt(mpt_tree("pairs", list(
    adjacent = "c * r",
    apart = "(1 - c) * u * u",
    one = "2 * (1 - c) * u * (1 - u)",
    none = "c * (1 - r) + (1 - c) * (1 - u) * (1 - u)"
  )))
  printed <- mpt_printed(pair_clustering)
  expect_match(printed, "Jacobian rank 3 of 3 at interior test values")
  expect_no_match(printed, "not identified")
})

test_that("the identifiability check counts a parameter fixed in the formula as known", {
  trees <- list(
    mpt_tree("old", list(yes = "Do + (1 - Do) * g", no = "(1 - Do) * (1 - g)")),
    mpt_tree("new", list(yes = "(1 - Dn) * g", no = "Dn + (1 - Dn) * (1 - g)"))
  )
  model <- mpt(trees, tree_id = "item_type")
  printed <- mpt_printed(model)
  expect_match(printed, "More free parameters than degrees of freedom")
  expect_match(printed, "1 combination\\(s\\) of all free parameters cannot")

  dat <- mpt_2htm_data()
  expect_warning(
    check_model(model, dat, bmf(Do ~ 1 + (1 | id), Dn ~ 1, g ~ 1)),
    "rank 2 for 3 free parameters.*the posterior follows the prior"
  )
  fixed <- expect_no_warning(check_model(model, dat, bmf(Do ~ 1, Dn ~ 1, g = 0.5)))
  expect_match(mpt_printed(fixed), "Jacobian rank 2 of 2 at interior test values")

  # a predictor of the guessing rate can identify the model across cells, so
  # the per-cell deficit is announced, not warned about
  dat$bias <- rep(c("low", "high"), length.out = nrow(dat))
  expect_no_warning(expect_message(
    check_model(model, dat, bmf(Do ~ 1, Dn ~ 1, g ~ 0 + bias)),
    "within one design cell.*predictors on 'g' identify it across cells is not checked"
  ))

  # the rank does not read non-linear formulas, so a deficit under one is
  # announced, not warned about, even when it is real as here
  expect_no_warning(expect_message(
    check_model(model, dat, bmf(Do ~ inv_logit(phi), phi ~ 1, Dn ~ 1, g ~ 1)),
    "rank 2 for 3 free parameters.*formula\\(s\\) for 'Do' identify it is not checked"
  ))
})

test_that("a formula that ties parameters together is not reported as a rank deficit", {
  trees <- list(
    mpt_tree("old", list(yes = "Do + (1 - Do) * g", no = "(1 - Do) * (1 - g)")),
    mpt_tree("new", list(yes = "(1 - Dn) * g", no = "Dn + (1 - Dn) * (1 - g)"))
  )
  model <- mpt(trees, tree_id = "item_type")
  dat <- mpt_2htm_data()

  # Dn ~ Do identifies the model; the rank of the tree parameters cannot see it
  expect_no_warning(expect_message(
    tied <- check_model(model, dat, bmf(Do ~ 1, Dn ~ Do, g ~ 1)),
    "non-linear formula\\(s\\) for 'Dn' identify it is not checked"
  ))
  printed <- mpt_printed(tied)
  expect_match(printed, "1 combination\\(s\\) of all free parameters are not identified")
  expect_no_match(printed, "The model is not identified")
  expect_match(printed, "formula\\(s\\) for 'Dn' were not analysed")

  expect_silent(check_model(model, dat, bmf(Do ~ 1, Dn = 0.6, g ~ 1)))

  # u and v enter only as u + v, which the tree rank of 'a' cannot see
  single <- mpt(mpt_tree("x", list(A = "a", B = "1 - a")))
  sub_pars <- suppressMessages(check_model(
    single, data.frame(A = 5L, B = 5L), bmf(a ~ inv_logit(u + v), u ~ 1, v ~ 1)
  ))
  printed <- mpt_printed(sub_pars)
  expect_no_match(printed, "locally identified")
  expect_match(
    printed,
    "Jacobian rank 1 of 1 in the tree parameters .* formula\\(s\\) for 'a' were not analysed"
  )

  # the deficit text names a simplex group by its members, not its sticks
  guessing <- mpt(mpt_tree("t", list(
    A = "D + (1 - D) * gA", B = "(1 - D) * gB", C = "(1 - D) * gC"
  )), simplex = list(c("gA", "gB", "gC")))
  tied_simplex <- suppressMessages(check_model(
    guessing, data.frame(A = 5L, B = 3L, C = 2L),
    bmf(D ~ inv_logit(phi), phi ~ 1, gA ~ 1, gB ~ 1)
  ))
  expect_match(
    mpt_printed(tied_simplex),
    "1 combination\\(s\\) of .*the simplex group 'gA', 'gB', 'gC'.* 'D' were not analysed"
  )
})

test_that("deep chains of identified parameters keep full rank", {
  # category i is reached after i - 1 failures, so the columns of the last
  # parameters are products of many probabilities but stay identified
  chain_model <- function(k) {
    reach <- function(i) if (i > 1) paste0("(1 - a", seq_len(i - 1), ")")
    branches <- lapply(seq_len(k), function(i) {
      paste(c(reach(i), paste0("a", i)), collapse = " * ")
    })
    branches[[k + 1]] <- paste(reach(k + 1), collapse = " * ")
    mpt(mpt_tree("t", setNames(branches, paste0("c", seq_len(k + 1)))))
  }
  for (k in c(20, 30)) {
    printed <- mpt_printed(chain_model(k))
    expect_match(printed, glue("Jacobian rank {k} of {k} at interior test values"))
    expect_no_match(printed, "zero up to rounding")
  }

  # each p_i moves 1 - p1 * ... * p30 by 1e-10 only, yet it does move it
  path <- paste(paste0("p", 1:30), collapse = " * ")
  product <- mpt(mpt_tree("t", list(
    x = glue("{path} * q"), y = glue("{path} * (1 - q)"), z = glue("1 - {path}")
  )))
  printed <- mpt_printed(product)
  expect_match(printed, "has rank 2 for 31 free parameters")
  expect_no_match(printed, "zero up to rounding")
})

test_that("a parameter that cancels with a rounding residue is not counted as identified", {
  # D(x, "q") is a * b * c - c * b * a, a residue of 1e-17 at some test points
  model <- mpt(mpt_tree("t", list(
    x = "a * b * c * q + c * b * a * (1 - q)",
    y = "a * (1 - b)", z = "(1 - a) * c",
    w = "1 - a * b * c - a * (1 - b) - (1 - a) * c"
  )))
  printed <- mpt_printed(model)
  expect_match(printed, "has rank 3 for 4 free parameters")
  expect_match(printed, "derivative with respect to 'q' is zero up to rounding")
  expect_no_match(printed, "combination\\(s\\)")
  expect_warning(
    check_model(model, data.frame(x = 5, y = 5, z = 5, w = 5), bmf(a ~ 1, b ~ 1, c ~ 1, q ~ 1)),
    "'q' is zero up to rounding"
  )
})

test_that("a residue column stays a zero column when its partners are fixed", {
  # with a, b and c fixed, q is the only free column and a pure residue: at
  # these values a * b * c and c * b * a differ in the last bit
  model <- mpt(mpt_tree("t", list(
    x = "a * b * c * q + c * b * a * (1 - q)", y = "1 - a * b * c"
  )))
  formula <- bmf(q ~ 1, a = 0.11, b = 0.13, c = 0.69)
  dat <- data.frame(x = 5, y = 5)
  expect_warning(
    checked <- check_model(model, dat, formula),
    "derivative with respect to 'q' is zero up to rounding"
  )
  expect_match(mpt_printed(checked), "has rank 0 for 1 free parameters")
})

test_that("a model whose every free column is exactly zero reports rank 0", {
  model <- mpt(mpt_tree("t", list(x = "q * 0.25 + (1 - q) * 0.25", y = "0.75")))
  printed <- mpt_printed(model)
  expect_match(printed, "has rank 0 for 1 free parameters")
  expect_match(printed, "derivative with respect to 'q' is zero up to rounding")
})

test_that("a long list of entangled parameters is printed as its complement", {
  model <- mpt(list(
    mpt_tree("a", list(x = "D * r * s", y = "(1 - D * r * s) * g", z = "(1 - D * r * s) * (1 - g)")),
    mpt_tree("b", list(x = "g", y = "(1 - g) * h", z = "(1 - g) * (1 - h)"))
  ), tree_id = "tree")
  expect_match(
    mpt_printed(model),
    "2 combination\\(s\\) of all free parameters except 'g', 'h' cannot"
  )
})

test_that("the rank check works in the stick-breaking components of a simplex group", {
  guessing <- mpt(list(
    mpt_tree("srcA", list(A = "D + (1 - D) * gA", B = "(1 - D) * gB", N = "(1 - D) * gN")),
    mpt_tree("new", list(A = "gA", B = "gB", N = "gN"))
  ), tree_id = "tree", simplex = c("gA", "gB", "gN"))
  printed <- mpt_printed(guessing)
  expect_match(printed, "3 free parameter\\(s\\), 4 degrees of freedom")
  expect_match(printed, "Jacobian rank 3 of 3 at interior test values")
  dat <- data.frame(tree = c("srcA", "new"), A = c(5, 3), B = c(2, 3), N = c(3, 4))
  expect_no_warning(expect_no_message(
    check_model(guessing, dat, bmf(D ~ 1, gA ~ 1, gB ~ 1))
  ))

  # the members only enter through their sum, which is 1 whatever the sticks
  sum_only <- mpt(mpt_tree("t", list(
    A = "D * gA + D * gB + D * gN", B = "(1 - D) * h", N = "(1 - D) * (1 - h)"
  )), simplex = c("gA", "gB", "gN"))
  printed <- mpt_printed(sum_only)
  expect_match(printed, "has rank 2 for 4 free parameters")
  expect_match(
    printed,
    "derivative with respect to the simplex group 'gA', 'gB', 'gN' is zero"
  )
  expect_no_match(printed, "combination\\(s\\)")
  expect_no_match(sub(".*not identified", "", printed), "raw")
})

test_that("the rank check stacks rows over covariate values and never frees a covariate", {
  scaled <- mpt(list(
    mpt_tree("old", list(yes = "D + (1 - D) * g * x", no = "(1 - D) * (1 - g * x)")),
    mpt_tree("new", list(yes = "(1 - D) * g * x", no = "D + (1 - D) * (1 - g * x)"))
  ), tree_id = "item_type", covariates = "x")
  expect_match(
    mpt_printed(scaled),
    "Jacobian rank 2 of 2 at interior test values and 5 test values of the covariates"
  )
  dat <- data.frame(item_type = c("old", "new"), x = c(0.5, 1), yes = 5, no = 5)
  expect_no_warning(expect_no_message(check_model(scaled, dat, bmf(D ~ 1, g ~ 1))))

  # one tree with a 0/1 covariate: a and b are identified only across rows
  switch_tree <- mpt(
    mpt_tree("t", list(yes = "x * a + (1 - x) * b", no = "1 - x * a - (1 - x) * b")),
    covariates = "x"
  )
  printed <- mpt_printed(switch_tree)
  expect_match(printed, "in one design cell; covariate values that differ")
  expect_match(printed, "Jacobian rank 2 of 2")
  both <- data.frame(x = c(0, 1, 1), yes = 5, no = 5)
  expect_no_warning(expect_no_message(check_model(switch_tree, both, bmf(a ~ 1, b ~ 1))))
  only_a <- data.frame(x = c(1, 1), yes = 5, no = 5)
  expect_warning(
    check_model(switch_tree, only_a, bmf(a ~ 1, b ~ 1)),
    "rank 1 for 2 free parameters at .*covariate values in the data.*'b' is zero"
  )

  # a covariate is never listed among the parameters
  product <- mpt(
    mpt_tree("t", list(yes = "x * a * b", no = "1 - x * a * b")),
    covariates = "x"
  )
  printed <- mpt_printed(product)
  expect_match(printed, "rank 1 for 2 free parameters")
  expect_no_match(printed, "'x'")

  # a residue that cancels stays a zero column when the rows are stacked
  residue <- mpt(mpt_tree("t", list(
    yes = "x * (a * b * c * q + c * b * a * (1 - q))",
    no = "1 - x * (a * b * c * q + c * b * a * (1 - q))"
  )), covariates = "x")
  expect_match(mpt_printed(residue), "derivative with respect to 'q' is zero")
})

test_that("the rank check ranks the design over every distinct covariate setting", {
  # z adds settings but no information about a, b, c; only the single row at
  # x = 0.5 separates b from the other two, so a spread over a subset loses it
  bernstein <- mpt(
    mpt_tree("t", list(
      yes = "z * (a * (1 - x)^2 + b * 2 * x * (1 - x) + c * x^2) + (1 - z) * 0.5",
      no = "1 - z * (a * (1 - x)^2 + b * 2 * x * (1 - x) + c * x^2) - (1 - z) * 0.5"
    )),
    covariates = c("x", "z")
  )
  dat <- rbind(
    data.frame(x = 0, z = seq(0.2, 0.8, length.out = 30)),
    data.frame(x = 0.5, z = 0.5),
    data.frame(x = 1, z = seq(0.2, 0.8, length.out = 30))
  )
  dat$yes <- 5
  dat$no <- 5
  expect_equal(nrow(bmm:::.mpt_covariate_settings(bernstein, dat)$values[[1]]), 61L)
  expect_equal(bmm:::.mpt_identifiability(bernstein, dat)$rank, 3L)
  expect_no_warning(expect_no_message(
    check_model(bernstein, dat, bmf(a ~ 1, b ~ 1, c ~ 1))
  ))
})

test_that("the rank check uses only finite covariate values and the trees with rows", {
  switch_tree <- mpt(
    mpt_tree("t", list(yes = "x * a + (1 - x) * b", no = "1 - x * a - (1 - x) * b")),
    covariates = "x"
  )
  # no finite value in the data: check_data() reports it, the rank check
  # falls back to test values and stays silent
  expect_no_warning(check_model(
    switch_tree, data.frame(x = c(Inf, NA), yes = 5, no = 5), bmf(a ~ 1, b ~ 1)
  ))

  scaled <- mpt(list(
    mpt_tree("old", list(yes = "D + (1 - D) * g * x", no = "(1 - D) * (1 - g * x)")),
    mpt_tree("new", list(yes = "(1 - D) * g * x", no = "D + (1 - D) * (1 - g * x)"))
  ), tree_id = "item_type", covariates = "x")
  missing_old <- data.frame(item_type = c("old", "new"), x = c(NA, 1), yes = 5, no = 5)
  expect_no_warning(check_model(scaled, missing_old, bmf(D ~ 1, g ~ 1)))

  # a tree without rows adds nothing, whether or not it uses covariates
  mixed <- mpt(list(
    mpt_tree("old", list(yes = "Do + (1 - Do) * g * x", no = "(1 - Do) * (1 - g * x)")),
    mpt_tree("new", list(yes = "(1 - Dn) * g", no = "Dn + (1 - Dn) * (1 - g)"))
  ), tree_id = "item_type", covariates = "x")
  only_old <- data.frame(item_type = "old", x = c(1, 1), yes = 5, no = 5)
  expect_warning(
    check_model(mixed, only_old, bmf(Do ~ 1, Dn ~ 1, g ~ 1)),
    "rank 1 for 3 free parameters"
  )

  # with predictors, the message names the parameters itself, because print()
  # takes the covariates at test values rather than at the data
  only_a <- data.frame(x = c(1, 1), cond = c("p", "q"), yes = 5, no = 5)
  msg <- expect_message(
    check_model(switch_tree, only_a, bmf(a ~ 1, b ~ 0 + cond)),
    "the covariate values in the data"
  )
  expect_no_match(conditionMessage(msg), "print\\(model\\)")
  expect_match(
    conditionMessage(msg),
    "'b' is zero up to rounding at interior test values and the covariate values in the data"
  )
})

test_that("the Jacobian rank is reported as not computed when D() cannot differentiate", {
  model <- mpt(mpt_tree("t", list(x = "plogis(a)", y = "1 - plogis(a)")))
  printed <- mpt_printed(model)
  expect_match(printed, "Jacobian rank not computed: .*'plogis'")
  expect_no_warning(check_model(model, data.frame(x = 1, y = 1), bmf(a ~ 1)))
})

test_that("fixed parameter values stay probabilities and reach the prior on the latent scale", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  formula <- bmf(D ~ 1, g = 0.7)
  model_checked <- check_model(model, dat, formula)
  expect_equal(model_checked$fixed_parameters$g, 0.7)
  expect_equal(check_model(model_checked, dat, formula), model_checked)

  prior <- default_prior(formula, dat, model)
  constant_row <- prior[grepl("constant", prior$prior), ]
  expect_equal(constant_row$nlpar, "g")
  expect_equal(constant_row$prior, glue("constant({qlogis(0.7)})"))

  model_probit <- mpt(mpt_2htm_trees(), tree_id = "item_type", links = "probit")
  prior_probit <- default_prior(formula, dat, model_probit)
  constant_row <- prior_probit[grepl("constant", prior_probit$prior), ]
  expect_equal(constant_row$prior, glue("constant({qnorm(0.7)})"))

  expect_error(check_model(model, dat, bmf(D ~ 1, g = 1.5)), "strictly")
})

test_that("simplex parameters cannot be fixed to constants", {
  trees <- list(
    mpt_tree("A", list(a = "gA", b = "gB", c = "gC")),
    mpt_tree("B", list(a = "gB", b = "gC", c = "gA"))
  )
  model <- mpt(trees, tree_id = "tree", simplex = c("gA", "gB", "gC"))
  expect_error(
    check_model(model, formula = bmf(gA = 0.3)),
    "Fixing simplex parameters"
  )
  expect_error(
    check_model(model, formula = bmf(gAraw = 0.3)),
    "Fixing simplex parameters.*gAraw"
  )
})

test_that("mpt compiles with a simplex group via stick-breaking", {
  trees <- list(
    mpt_tree("sourceA", list(
      A = "dA + (1 - dA) * gA",
      B = "(1 - dA) * gB",
      New = "(1 - dA) * gNew"
    )),
    mpt_tree("sourceB", list(
      A = "(1 - dB) * gA",
      B = "dB + (1 - dB) * gB",
      New = "(1 - dB) * gNew"
    )),
    mpt_tree("new", list(A = "gA", B = "gB", New = "gNew"))
  )
  model <- mpt(trees, tree_id = "source", simplex = c("gA", "gB", "gNew"))
  expect_setequal(
    names(model$parameters),
    c("dA", "dB", "gA", "gB", "gNew", "gAraw", "gBraw")
  )
  expect_equal(model$links$gA, "identity")
  expect_equal(model$default_priors$gAraw$main, "logistic(0, 1)")
  expect_null(model$default_priors[["gA"]])

  dat <- expand.grid(
    id = factor(1:10), source = c("sourceA", "sourceB", "new"),
    stringsAsFactors = FALSE
  )
  counts <- t(rmultinom(nrow(dat), 30, c(0.4, 0.3, 0.3)))
  colnames(counts) <- c("A", "B", "New")
  dat <- cbind(dat, counts)
  formula <- bmf(dA ~ 1, dB ~ 1, gA ~ 1 + (1 | id), gB ~ 1)

  model_checked <- check_model(model, dat, formula)
  dat_checked <- check_data(model_checked, dat, formula)
  formula_checked <- suppressMessages(
    check_formula(model_checked, dat_checked, formula)
  )
  gA_rhs <- paste(deparse(formula_checked$gA[[3]]), collapse = " ")
  gNew_rhs <- paste(deparse(formula_checked$gNew[[3]]), collapse = " ")
  expect_true(grepl("inv_logit(gAraw)", gA_rhs, fixed = TRUE))
  expect_equal(gNew_rhs, "1 - (gA + gB)")
  gAraw_rhs <- paste(deparse(formula_checked$gAraw[[3]]), collapse = " ")
  expect_true(grepl("(1 | id)", gAraw_rhs, fixed = TRUE))

  expect_no_warning(
    suppressMessages(bmm(
      formula, dat, model,
      backend = "mock", mock_fit = 1, rename = FALSE
    )),
    message = "Non-linear"
  )
})

test_that("stick-breaking members form a simplex for any raw values", {
  tree <- mpt_tree("t", list(A = "gA", B = "gB", C = "gC", D = "gD"))
  model <- mpt(tree, simplex = c("gA", "gB", "gC", "gD"))
  dat <- data.frame(A = 5, B = 5, C = 5, D = 5)
  formula <- suppressMessages(check_formula(model, dat, bmf(gA ~ 1, gB ~ 1, gC ~ 1)))
  raw <- list(gAraw = 1.3, gBraw = -0.4, gCraw = 0.7)
  env <- list2env(c(raw, list(inv_logit = stats::plogis)))
  for (par in c("gA", "gB", "gC", "gD")) assign(par, eval(formula[[par]][[3]], env), env)
  s <- unname(stats::plogis(unlist(raw)))
  expect_equal(
    unlist(mget(c("gA", "gB", "gC", "gD"), env), use.names = FALSE),
    c(s[1], (1 - s[1]) * s[2], (1 - s[1]) * (1 - s[2]) * s[3], prod(1 - s))
  )
})

test_that("predictors on the derived simplex parameter are rejected", {
  trees <- list(mpt_tree("t", list(A = "gA", B = "gB", New = "gNew")))
  model <- mpt(trees, simplex = c("gA", "gB", "gNew"))
  dat <- data.frame(id = factor(1:5), A = 10, B = 10, New = 10)
  formula <- bmf(gA ~ 1, gB ~ 1, gNew ~ 1 + (1 | id))
  expect_error(
    suppressMessages(bmm(
      formula, dat, model,
      backend = "mock", mock_fit = 1, rename = FALSE
    )),
    "derived"
  )
})

test_that("check_data errors are informative", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()

  dat_missing <- dat[setdiff(names(dat), "new")]
  expect_error(
    check_data(model, dat_missing, bmf(D ~ 1, g ~ 1)),
    "Missing columns: 'new'"
  )

  dat_bad_cond <- dat
  dat_bad_cond$item_type[1] <- "unknown"
  expect_error(
    check_data(model, dat_bad_cond, bmf(D ~ 1, g ~ 1)),
    "Unmatched values: 'unknown'"
  )

  tree <- mpt_tree("main", list(
    correct = "D + (1 - D) * Gcorr",
    incorrect = "(1 - D) * (1 - Gcorr)"
  ))
  model_cov <- mpt(tree, covariates = "Gcorr")
  dat_cov <- data.frame(correct = 10, incorrect = 10)
  expect_error(
    check_data(model_cov, dat_cov, bmf(D ~ 1)),
    "covariates 'Gcorr' are missing"
  )
})

test_that("check_data requires a column for every declared covariate, used or not", {
  unused <- mpt(
    mpt_tree("main", list(correct = "D", incorrect = "1 - D")),
    covariates = "z"
  )
  expect_error(
    check_data(unused, data.frame(correct = 10, incorrect = 10), bmf(D ~ 1)),
    "covariates 'z' are missing"
  )
  expect_no_error(check_data(
    unused, data.frame(correct = 10, incorrect = 10, z = 1), bmf(D ~ 1)
  ))
})

test_that("check_data warns on missing counts and refuses the columns it builds", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()

  dat_na <- dat
  dat_na$old[c(1, 3)] <- NA
  expect_warning(
    checked <- check_data(model, dat_na, bmf(D ~ 1, g ~ 1)),
    "2 missing value\\(s\\), which are counted as 0"
  )
  expect_equal(unname(checked$Y[1, ]), c(0, dat$new[1]))

  dat_reserved <- dat
  dat_reserved$nTrials <- 50
  expect_error(
    check_data(model, dat_reserved, bmf(D ~ 1, g ~ 1)),
    "'nTrials' would be overwritten"
  )
  dat_reserved$Y <- 1
  expect_error(
    check_data(model, dat_reserved, bmf(D ~ 1, g ~ 1)),
    "'Y', 'nTrials'"
  )
})

test_that("mpt category probabilities match the production m3 likelihood", {
  # for the simple choice rule with fixed b, list size NL and response-set
  # size N, the simple-span M3 is an MPT with design-fixed guessing rates;
  # the two probability vectors must be identical through the analytic
  # bijection a = 2*b*Pm*(1-Pb)/(1-Pm), c = N*b*Pm*Pb/(1-Pm)
  b <- 0.1
  NL <- 4
  N <- 8

  mpt_model <- mpt(mpt_tree("main", list(
    correct = "Pm*Pb + Pm*(1 - Pb)*0.25 + (1 - Pm)*0.125",
    other = "Pm*(1 - Pb)*0.75 + (1 - Pm)*0.375",
    npl = "(1 - Pm)*0.5"
  )))

  m3_model <- m3(
    resp_cats = c("correct", "other", "npl"),
    num_options = c(1, NL - 1, N - NL),
    choice_rule = "simple",
    version = "custom"
  )
  act_funs <- bmf(correct ~ b + a + c, other ~ b + a, npl ~ b)

  grid <- expand.grid(Pm = seq(0.15, 0.9, 0.15), Pb = seq(0.15, 0.9, 0.15))
  max_diff <- max(vapply(seq_len(nrow(grid)), function(i) {
    Pm <- grid$Pm[i]
    Pb <- grid$Pb[i]
    p_mpt <- .mpt_probability_vector(
      pars = c(Pm = Pm, Pb = Pb), mpt_model = mpt_model
    )
    a <- 2 * b * Pm * (1 - Pb) / (1 - Pm)
    c_par <- N * b * Pm * Pb / (1 - Pm)
    p_m3 <- .compute_m3_probability_vector(
      pars = c(a = a, c = c_par, b = b), m3_model = m3_model,
      act_funs = act_funs
    )
    max(abs(p_mpt - p_m3))
  }, numeric(1)))

  expect_lt(max_diff, 1e-10)
})

test_that("simplex predictors may be given on the parameter or its raw component", {
  trees <- list(mpt_tree("t", list(A = "gA", B = "gB", C = "gC")))
  model <- mpt(trees, simplex = c("gA", "gB", "gC"))
  dat <- data.frame(id = factor(1:5), A = 10, B = 10, C = 10)

  # an explicit gA ~ 1 does not override predictors on gAraw
  on_raw <- suppressMessages(check_formula(
    model, dat, bmf(gA ~ 1, gAraw ~ 1 + (1 | id), gB ~ 1, gBraw ~ 1)
  ))
  expect_equal(deparse1(on_raw$gAraw[[3]]), "1 + (1 | id)")
  expect_equal(deparse1(on_raw$gA[[3]]), "inv_logit(gAraw)")

  expect_error(
    suppressMessages(check_formula(
      model, dat, bmf(gA ~ 0 + id, gAraw ~ 1 + (1 | id), gB ~ 1, gBraw ~ 1)
    )),
    "Conflicting"
  )
  expect_no_error(suppressMessages(check_formula(
    model, dat, bmf(gA ~ 0 + id, gAraw ~ 0 + id, gB ~ 1, gBraw ~ 1)
  )))

  dat$gAraw <- 1
  expect_error(check_data(model, dat, bmf(gA ~ 1, gB ~ 1)), "stick-breaking")
})

test_that("mpt supports multiple simplex groups", {
  tree <- mpt_tree("t", list(
    A = "m * gA + (1 - m) * hA",
    B = "m * gB + (1 - m) * hB",
    C = "m * gC + (1 - m) * hC"
  ))
  model <- mpt(
    tree,
    simplex = list(c("gA", "gB", "gC"), c("hA", "hB", "hC"))
  )
  expect_setequal(
    names(model$parameters),
    c("m", "gA", "gB", "gC", "hA", "hB", "hC", "gAraw", "gBraw", "hAraw", "hBraw")
  )

  dat <- data.frame(id = factor(1:8), A = 10, B = 10, C = 10)
  formula <- bmf(m ~ 1, gA ~ 1, gB ~ 1, hA ~ 1, hB ~ 1)
  # one tree has two degrees of freedom for m and four sticks
  expect_warning(
    expect_no_warning(
      suppressMessages(bmm(
        formula, dat, model,
        backend = "mock", mock_fit = 1, rename = FALSE
      )),
      message = "Non-linear"
    ),
    "rank 2 for 5 free parameters"
  )
})

test_that("factor tree identifier columns are matched to tree names", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  dat$item_type <- factor(dat$item_type, levels = c("old", "new"))
  checked <- check_data(model, dat, bmf(D ~ 1, g ~ 1))
  expect_equal(checked$Idx_old, as.integer(dat$item_type == "old"))
  expect_equal(checked$Idx_new, as.integer(dat$item_type == "new"))
})

test_that("mpt_tree validates impossible response categories", {
  branches <- list(a = "D", b = "1 - D")
  expect_error(
    mpt_tree("t", branches, impossible = "a"),
    "cannot be both impossible and have a branch"
  )
  expect_error(
    mpt_tree("t", branches, impossible = c("c", "c")),
    "must be unique"
  )
  expect_error(mpt_tree("t", branches, impossible = 1), "character vector")

  tree <- mpt_tree("t", branches, impossible = "c")
  expect_equal(tree$impossible, "c")
  expect_output(print(tree), "P\\(c\\) = 0 \\(structurally impossible\\)")
})

test_that("mpt requires impossible categories to exist in some other tree", {
  trees <- list(
    mpt_tree("one", list(a = "D", b = "1 - D"), impossible = "c"),
    mpt_tree("two", list(a = "1 - D", b = "D"), impossible = "c")
  )
  expect_error(
    mpt(trees, tree_id = "cond"),
    "impossible in every tree"
  )
  expect_error(
    mpt(mpt_tree("one", list(a = "D", b = "1 - D"), impossible = "c")),
    "impossible in every tree"
  )
})

test_that("mpt counts impossible categories when comparing trees", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  expect_equal(model$resp_vars$resp_cats, c("corr", "dist", "npl"))

  mismatched <- list(
    mpt_tree("one", list(a = "D", b = "1 - D"), impossible = "c"),
    mpt_tree("two", list(a = "D", b = "1 - D", d = "0"))
  )
  expect_error(mpt(mismatched, tree_id = "cond"), "same response categories")
})

test_that("check_data builds possibility indicators for impossible categories", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  checked <- check_data(model, dat, bmf(Pm ~ 1, Pb ~ 1))
  expect_equal(checked$Poss_dist, as.integer(dat$cond == "withdist"))

  dat_collide <- dat
  dat_collide$Poss_dist <- 1
  expect_error(check_data(model, dat_collide, bmf(Pm ~ 1, Pb ~ 1)), "reserved")

  dat_observed <- dat
  dat_observed$dist[dat_observed$cond == "same"][1] <- 3L
  expect_error(
    check_data(model, dat_observed, bmf(Pm ~ 1, Pb ~ 1)),
    "declared impossible for 1 observation"
  )
})

test_that("missing counts of an impossible category do not warn", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  dat$dist[dat$tree == "nodist"] <- NA
  expect_silent(checked <- check_data(model, dat, bmf(Pm ~ 1, Pb ~ 1)))
  expect_true(all(checked$Y[dat$tree == "nodist", "dist"] == 0))

  dat$dist[dat$tree == "withdist"][1] <- NA
  dat$corr[2] <- NA
  expect_warning(
    check_data(model, dat, bmf(Pm ~ 1, Pb ~ 1)),
    "2 missing value\\(s\\)"
  )
})

test_that("impossible categories are switched off in the linear predictor", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  formula <- bmf(Pm ~ 1, Pb ~ 1)
  checked_data <- check_data(model, dat, formula)
  checked_formula <- check_formula(model, checked_data, formula)

  # the tree that cannot produce the category contributes a placeholder, so
  # log() stays defined for its rows
  expect_match(deparse1(checked_formula$dist[[3]]), "Idx_nodist * (1)", fixed = TRUE)

  brms_formula <- configure_model(model, checked_data, checked_formula)$formula
  expect_match(
    deparse1(brms_formula$pforms$mudist[[3]]),
    "Poss_dist * log(dist) + (1 - Poss_dist) * (-100)",
    fixed = TRUE
  )
  expect_match(deparse1(brms_formula$formula[[3]]), "log(corr)", fixed = TRUE)
})

test_that("the brms formula reproduces the tree probabilities with impossible categories", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  formula <- bmf(Pm ~ 1, Pb ~ 1)
  checked_data <- check_data(model, dat, formula)
  checked_formula <- check_formula(model, checked_data, formula)
  brms_formula <- configure_model(model, checked_data, checked_formula)$formula
  resp_cats <- model$resp_vars$resp_cats

  # intercept-only parameters on the logit scale, evaluated as brms does:
  # category probabilities first, then the linear predictor of each category
  pars <- c(Pm = 0.4, Pb = 1.1)
  env <- c(as.list(checked_data), as.list(pars), list(inv_logit = stats::plogis))
  probs <- lapply(resp_cats, function(cat) eval(checked_formula[[cat]][[3]], env))
  names(probs) <- resp_cats
  env <- c(env, probs)
  eta <- vapply(resp_cats, function(cat) {
    rhs <- if (cat == resp_cats[1]) {
      brms_formula$formula[[3]]
    } else {
      brms_formula$pforms[[paste0("mu", cat)]][[3]]
    }
    eval(rhs, env)
  }, numeric(nrow(dat)))
  softmax <- exp(eta) / rowSums(exp(eta))

  expected <- t(vapply(seq_len(nrow(dat)), function(i) {
    .mpt_probability_vector(plogis(pars), model, tree = dat$tree[i])[resp_cats]
  }, numeric(length(resp_cats))))
  is_impossible <- dat$tree == "nodist"
  expect_true(all(softmax[is_impossible, "dist"] < 1e-40))
  expect_equal(
    unname(softmax[, c("corr", "npl")]), unname(expected[, c("corr", "npl")]),
    tolerance = 1e-12
  )
  expect_equal(unname(softmax[!is_impossible, ]), unname(expected[!is_impossible, ]),
    tolerance = 1e-12
  )
})

test_that("mpt compiles with structurally impossible categories", {
  expect_silent(bmm(
    bmf(Pm ~ 1 + (1 | id), Pb ~ 1),
    mpt_impossible_data(), mpt(mpt_impossible_trees(), tree_id = "tree"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
})

test_that("several levels of a factor can share one tree", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  checked <- check_data(model, dat, bmf(Pm ~ 0 + cond, Pb ~ 1))
  expect_equal(checked$Idx_withdist, as.integer(dat$cond == "withdist"))
  expect_equal(checked$Idx_nodist, as.integer(dat$cond %in% c("reord", "same")))

  # the experimental factor survives untouched for the parameter formulas
  expect_equal(checked$cond, dat$cond)
})

test_that("tree identifier values must match tree names", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  dat$tree[dat$cond == "reord"] <- "reord"
  expect_error(
    check_data(model, dat, bmf(Pm ~ 1, Pb ~ 1)),
    "Unmatched values: 'reord'"
  )

  missing_col <- mpt_impossible_data()
  missing_col$tree <- NULL
  expect_error(
    check_data(model, missing_col, bmf(Pm ~ 1, Pb ~ 1)),
    "not present in the data"
  )
})

test_that("check_data validates branch sums with observed covariate values", {
  # the tree sums to 1 only when Gcorr + Gother = 1, which synthetic test
  # values at construction cannot verify (mpt() stays silent) but the observed
  # covariate columns can (check_data() errors there)
  tree <- mpt_tree("main", list(
    correct = "D + (1 - D) * Gcorr",
    incorrect = "(1 - D) * Gother"
  ))
  model <- expect_silent(mpt(tree, covariates = c("Gcorr", "Gother")))
  dat <- data.frame(
    id = factor(1:6), Gcorr = 0.25, Gother = 0.75,
    correct = 10, incorrect = 30
  )
  expect_silent(check_data(model, dat, bmf(D ~ 1)))

  dat_bad <- dat
  dat_bad$Gother[3] <- 0.9
  expect_error(
    check_data(model, dat_bad, bmf(D ~ 1)),
    "do not sum to 1 for 1 row"
  )

  dat_na <- dat
  dat_na$Gcorr[c(2, 5)] <- NA
  expect_error(
    check_data(model, dat_na, bmf(D ~ 1)),
    "'Gcorr' is missing in row 2 of tree 'main'.*cannot be computed"
  )

  # swapped complements agree at the symmetric point D = 0.5, so the branch
  # expressions are checked at several parameter values
  swapped <- expect_silent(mpt(mpt_tree("main", list(
    correct = "D + (1 - D) * Gcorr", incorrect = "D * (1 - Gcorr)"
  )), covariates = "Gcorr"))
  dat_swap <- data.frame(Gcorr = 0.25, correct = 10, incorrect = 30)
  expect_error(check_data(swapped, dat_swap, bmf(D ~ 1)), "do not sum to 1")

  stray <- expect_silent(mpt(mpt_tree("main", list(
    correct = "D + (1 - D) * Gcorr", incorrect = "2 * D * (1 - D) * (1 - Gcorr)"
  )), covariates = "Gcorr"))
  expect_error(check_data(stray, dat_swap, bmf(D ~ 1)), "do not sum to 1")
})

test_that("check_data rejects covariate values that push a branch outside (0, 1]", {
  model <- mpt(mpt_tree("main", list(
    correct = "D + (1 - D) * G", incorrect = "(1 - D) * (1 - G)"
  )), covariates = "G")
  dat <- data.frame(G = c(0.25, 0.5, 0.75), correct = 10, incorrect = 30)
  expect_silent(check_data(model, dat, bmf(D ~ 1)))

  expect_silent(check_data(model, transform(dat, G = 0), bmf(D ~ 1)))

  dat_over <- dat
  dat_over$G[2] <- 1.2
  expect_error(
    check_data(model, dat_over, bmf(D ~ 1)),
    "'correct' in tree 'main' is [0-9.]+ in row 2, outside \\(0, 1\\]"
  )
  # no trailing blank when there is no hint about impossible categories
  expect_error(check_data(model, dat_over, bmf(D ~ 1)), "column\\(s\\): 'G'\\.$")

  # the branch D + (1 - D) * G is negative only for small D, below every
  # interior test point, so a slightly negative covariate needs the
  # near-boundary points
  for (g in c(-0.05, -0.1)) {
    dat_neg <- dat
    dat_neg$G[1] <- g
    expect_error(
      check_data(model, dat_neg, bmf(D ~ 1)),
      "'correct' in tree 'main' is -[0-9.]+ in row 1, outside \\(0, 1\\]"
    )
  }

  # a branch of exactly 0 is a structurally impossible category
  dat_one <- dat
  dat_one$G[3] <- 1
  expect_error(
    check_data(model, dat_one, bmf(D ~ 1)),
    "is 0 in row 3.*mpt_tree\\(impossible = \\)"
  )

  zero_int <- mpt(mpt_tree("main", list(
    hit = "D * G", miss = "1 - D * G"
  )), covariates = "G")
  dat_int <- data.frame(G = c(1L, 0L), hit = 10, miss = 30)
  expect_error(
    check_data(zero_int, dat_int, bmf(D ~ 1)),
    "'hit' in tree 'main' is 0 in row 2"
  )
})

test_that("the near-boundary points accept valid simplex trees with covariates", {
  tree <- mpt_tree("main", list(
    hit = "G * D + (1 - G * D) * gA",
    other = "(1 - G * D) * gB",
    miss = "(1 - G * D) * gC"
  ))
  model <- mpt(
    tree, covariates = "G", simplex = c("gA", "gB", "gC")
  )
  dat <- data.frame(G = c(0.2, 0.5, 0.9), hit = 10, other = 10, miss = 10)
  expect_silent(check_data(model, dat, bmf(D ~ 1, gA ~ 1, gB ~ 1)))
  expect_error(
    check_data(model, transform(dat, G = c(0.2, -0.1, 0.9)), bmf(D ~ 1, gA ~ 1, gB ~ 1)),
    "in row 2, outside \\(0, 1\\]"
  )
})

test_that("an NA covariate in a tree that does not use it is an error, not a dropped row", {
  trees <- list(
    mpt_tree("withdist", list(
      corr = "Pm * Pb + (1 - Pm) * (1 - Gd)",
      dist = "(1 - Pm) * Gd",
      npl = "Pm * (1 - Pb)"
    )),
    mpt_tree("nodist", list(
      corr = "Pm * Pb + (1 - Pm) * 0.2",
      npl = "Pm * (1 - Pb) + (1 - Pm) * 0.8"
    ), impossible = "dist")
  )
  model <- mpt(trees, tree_id = "tree", covariates = "Gd")
  dat <- mpt_impossible_data()
  dat$Gd <- ifelse(dat$tree == "withdist", 0.3, 0)
  expect_silent(check_data(model, dat, bmf(Pm ~ 1, Pb ~ 1)))

  dat_na <- dat
  dat_na$Gd[dat_na$tree == "nodist"] <- NA
  first_row <- which(dat_na$tree == "nodist")[1]
  expect_error(
    check_data(model, dat_na, bmf(Pm ~ 1, Pb ~ 1)),
    glue::glue("'Gd' is missing in row {first_row} of tree 'nodist'.*does not use 'Gd'")
  )

  # in a tree that uses the covariate, the branches cannot be computed
  dat_used <- dat
  dat_used$Gd[dat_used$tree == "withdist"][1] <- NA
  first_used <- which(dat_used$tree == "withdist")[1]
  expect_error(
    check_data(model, dat_used, bmf(Pm ~ 1, Pb ~ 1)),
    glue::glue("'Gd' is missing in row {first_used} of tree 'withdist'.*cannot be computed")
  )
})

test_that("the NA advice for an unused covariate does not promise that 0 is safe", {
  model <- mpt(list(
    mpt_tree("sstree", list(hit = "D + (1 - D) / ss", miss = "(1 - D) * (1 - 1 / ss)")),
    mpt_tree("plain", list(hit = "D + (1 - D) * g", miss = "(1 - D) * (1 - g)"))
  ), tree_id = "tree", covariates = "ss")
  dat <- data.frame(
    tree = rep(c("sstree", "plain"), each = 2), ss = c(4, 4, NA, NA),
    hit = 5, miss = 5
  )
  msg <- tryCatch(check_data(model, dat, bmf(D ~ 1, g ~ 1)), error = conditionMessage)
  expect_match(msg, "does not use 'ss'")
  expect_match(msg, "0 works unless one of them divides by 'ss'", fixed = TRUE)
  expect_no_match(msg, "cannot be computed")
})

test_that("a branch that is undefined on the rows of another tree is an error", {
  # brms evaluates every tree's branches on every row, so a 0 filled in for a
  # tree that does not use the covariate breaks a tree that divides by it
  model <- mpt(list(
    mpt_tree("sstree", list(hit = "D + (1 - D) / ss", miss = "(1 - D) * (1 - 1 / ss)")),
    mpt_tree("plain", list(hit = "D + (1 - D) * g", miss = "(1 - D) * (1 - g)"))
  ), tree_id = "tree", covariates = "ss")
  formula <- bmf(D ~ 1, g ~ 1)
  dat <- data.frame(
    tree = rep(c("sstree", "plain"), each = 2), ss = c(4, 4, 0, 0),
    hit = 5, miss = 5
  )
  expect_error(
    check_data(model, dat, formula),
    "'hit' in tree 'sstree' is not finite in row 3.*belongs to tree 'plain'"
  )
  dat$ss[3:4] <- 1
  expect_silent(check_data(model, dat, formula))
})

test_that("a branch that divides by a covariate equal to 0 in its own tree gets the package message", {
  model <- mpt(mpt_tree("main", list(
    hit = "D + (1 - D) / ss", miss = "(1 - D) * (1 - 1 / ss)"
  )), covariates = "ss")
  dat <- data.frame(ss = c(4, 0), hit = 5, miss = 5)
  expect_error(
    check_data(model, dat, bmf(D ~ 1)),
    "do not sum to 1.*first: row 2, sum = NaN"
  )
})

test_that("a valid branch that underflows to 0 near the boundary is accepted", {
  model <- mpt(mpt_tree("main", list(
    hit = "1 - (1 - D)^n", miss = "(1 - D)^n"
  )), covariates = "n")
  expect_silent(check_data(model, data.frame(hit = 5L, miss = 5L, n = 200), bmf(D ~ 1)))
  expect_silent(check_data(model, data.frame(hit = 5L, miss = 5L, n = 2), bmf(D ~ 1)))

  neg <- mpt(mpt_tree("main", list(
    correct = "D + (1 - D) * G", incorrect = "(1 - D) * (1 - G)"
  )), covariates = "G")
  expect_error(
    check_data(neg, data.frame(G = -0.05, correct = 10, incorrect = 30), bmf(D ~ 1)),
    "-[0-9.]+ in row 1"
  )
})

test_that("the sum and range messages name only the covariates the tree uses", {
  sums <- mpt(mpt_tree("main", list(
    a = "D", b = "(1 - D) * H"
  )), covariates = c("G", "H"))
  msg <- tryCatch(
    check_data(sums, data.frame(G = 0.5, H = 0.5, a = 5, b = 5), bmf(D ~ 1)),
    error = conditionMessage
  )
  expect_match(msg, "do not sum to 1 for 1 row.*at the test parameter values")
  expect_match(msg, "column(s): 'H'", fixed = TRUE)
  expect_no_match(msg, "'G'")

  # a tree without covariates beside a covariate tree is refused at mpt(),
  # before any covariate value exists to blame
  trees <- list(
    mpt_tree("t1", list(x = "2 * a", y = "1 - 2 * a")),
    mpt_tree("t2", list(x = "D + (1 - D) * G", y = "(1 - D) * (1 - G)"))
  )
  msg <- tryCatch(mpt(trees, tree_id = "tree", covariates = "G"), error = conditionMessage)
  expect_match(msg, "category 'x' in tree 't1' is [0-9.]+ at the test values a = .*outside \\(0, 1\\]")
  expect_no_match(msg, "'G'")
})

test_that("mpt() refuses a tree that leaves (0, 1] only at a boundary corner, with or without covariates", {
  tree <- mpt_tree("t", list(yes = "1.2 * a - 0.2", no = "1.2 - 1.2 * a"))
  corner_msg <- "category 'yes' in tree 't' is -0.1988 at the test values a = 0.001, outside \\(0, 1\\]"
  expect_error(mpt(tree), corner_msg)
  expect_error(mpt(tree, covariates = "G"), corner_msg)

  # 0.001^200 underflows to exactly 0 at the corner D = 0.999
  expect_silent(mpt(mpt_tree("u", list(hit = "1 - (1 - D)^200", miss = "(1 - D)^200"))))
})

test_that("mpt() alone decides a tree without covariates in a model with covariates and a simplex", {
  covariate_tree <- mpt_tree("t1", list(
    A = "gA * D + (1 - D) * x", B = "gB * D", C = "gC * D + (1 - D) * (1 - x)"
  ))
  build <- function(free_tree) {
    mpt(
      list(covariate_tree, free_tree),
      tree_id = "tree", covariates = "x", simplex = c("gA", "gB", "gC")
    )
  }
  expect_error(
    build(mpt_tree("t2", list(
      A = "1.25 * D", B = "(1 - 1.25 * D) * 0.5", C = "(1 - 1.25 * D) * 0.5"
    ))),
    "category 'A' in tree 't2' is 1.24875 at the test values D = 0.999"
  )

  # the branches of t2 sum to 1 at every value of D that mpt() tries, and to
  # something else at the values a data check in another symbol order would try
  # (t2 adds no symbol, so t1 alone sets the order)
  construction <- .mpt_tree_parameters(list(covariate_tree), "x")
  checked_d <- c(
    vapply(.mpt_test_points(construction, list(c("gA", "gB", "gC"))), `[[`, numeric(1), "D"),
    0.001, 0.999
  )
  bump <- paste0("(D - ", sprintf("%.10f", checked_d), ")", collapse = " * ")
  model <- build(mpt_tree("t2", list(A = "D", B = glue("1 - D + {bump}")), impossible = "C"))
  dat <- data.frame(tree = c("t1", "t2"), x = 0.5, A = 3, B = 3, C = c(3, 0))
  expect_silent(check_data(model, dat, bmf(D ~ 1)))
})

test_that("the data check catches a tree that is right at the first interior and boundary points", {
  at_first <- sprintf("%.8f", .mpt_test_points("D", list())[[1]][["D"]])
  model <- mpt(mpt_tree("main", list(
    a = "D",
    b = glue("1 - D + (D - {at_first}) * (D - 0.001) * (D - 0.999) * G")
  )), covariates = "G")
  expect_error(
    check_data(model, data.frame(G = 0.5, a = 5, b = 5), bmf(D ~ 1)),
    "do not sum to 1"
  )
})

test_that("an NA in a declared covariate that no branch uses keeps every row", {
  model <- mpt(mpt_tree("main", list(
    correct = "D + (1 - D) * G", incorrect = "(1 - D) * (1 - G)"
  )), covariates = c("G", "H"))
  dat <- data.frame(G = 0.25, H = c(NA, 1, 2, NA), correct = 10, incorrect = 30)
  expect_silent(checked <- check_data(model, dat, bmf(D ~ 1)))
  expect_equal(nrow(checked), nrow(dat))
})

test_that("check_data catches a typo between parameters four places apart", {
  branches <- list(
    r1 = "A * B", r2 = "A * (1 - B) * C", r3 = "A * (1 - B) * (1 - C) * E",
    r4 = "A * (1 - B) * (1 - C) * (1 - E) * G",
    r5 = "A * (1 - B) * (1 - C) * (1 - E) * (1 - G) + (1 - A) * F"
  )
  dat <- data.frame(G = c(0.2, 0.6), r1 = 5, r2 = 5, r3 = 5, r4 = 5, r5 = 5, r6 = 5)
  formula <- bmf(A ~ 1, B ~ 1, C ~ 1, E ~ 1, F ~ 1)

  good <- mpt(
    mpt_tree("main", c(branches, r6 = "(1 - A) * (1 - F)")), covariates = "G"
  )
  expect_silent(check_data(good, dat, formula))

  typo <- mpt(
    mpt_tree("main", c(branches, r6 = "(1 - A) * (1 - A)")), covariates = "G"
  )
  expect_error(check_data(typo, dat, formula), "do not sum to 1")
})

test_that("covariate sum check respects tree membership", {
  trees <- list(
    mpt_tree("cued", list(
      correct = "D + (1 - D) * Gcorr",
      incorrect = "(1 - D) * Gother"
    )),
    mpt_tree("free", list(
      correct = "D",
      incorrect = "1 - D"
    ))
  )
  model <- mpt(trees, tree_id = "cond", covariates = c("Gcorr", "Gother"))
  dat <- data.frame(
    cond = rep(c("cued", "free"), each = 3),
    Gcorr = c(0.25, 0.25, 0.25, 99, 99, 99),
    Gother = c(0.75, 0.75, 0.75, 99, 99, 99),
    correct = 10, incorrect = 30
  )
  # the invalid covariate values sit in rows of the tree that does not use
  # the covariates, so no warning should be raised
  expect_silent(check_data(model, dat, bmf(D ~ 1)))
})

test_that("a tree without rows is still checked on the rows of the other trees", {
  trees <- list(
    mpt_tree("sstree", list(hit = "D + (1 - D) / ss", miss = "(1 - D) * (1 - 1 / ss)")),
    mpt_tree("plain", list(hit = "D + (1 - D) * g", miss = "(1 - D) * (1 - g)"))
  )
  model <- mpt(trees, tree_id = "tt", covariates = "ss")
  # brms evaluates the sstree branches on the plain rows, where 1 / 0 is Inf
  dat <- data.frame(tt = "plain", ss = 0, hit = c(60, 70, 65), miss = c(40, 30, 35))
  expect_error(
    suppressWarnings(check_data(model, dat, bmf(D ~ 1, g ~ 1))),
    "category 'hit' in tree 'sstree' is not finite in row 1"
  )
})

test_that("the item-memory-first MPT matches the simple-rule m3 with a distractor category", {
  # bijection for act_funs corr ~ b+a+c, other ~ b+a, dist ~ b+d, npl ~ b
  # with candidate counts (1, 4, 5, 5) and S = 15b + 5a + c + 5d:
  #   Pi = (5a + c + 5d)/S, Pb = c/(5a + c + 5d), Pd = d/(a + d)
  b <- 0.1
  mpt_dist <- mpt(mpt_tree("newdist", list(
    corr = "Pi * Pb + Pi * (1 - Pb) * (1 - Pd) * (1/5) + (1 - Pi) * (1/15)",
    other = "Pi * (1 - Pb) * (1 - Pd) * (4/5) + (1 - Pi) * (4/15)",
    dist = "Pi * (1 - Pb) * Pd + (1 - Pi) * (5/15)",
    npl = "(1 - Pi) * (5/15)"
  )))
  m3_dist <- m3(
    resp_cats = c("corr", "other", "dist", "npl"),
    num_options = c(1, 4, 5, 5), choice_rule = "simple", version = "custom"
  )
  acts_dist <- bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b)

  grid <- expand.grid(a = c(0.2, 1, 3), c = c(0.5, 2, 6), d = c(0.1, 0.8, 2))
  max_diff <- max(vapply(seq_len(nrow(grid)), function(i) {
    a <- grid$a[i]
    c_par <- grid$c[i]
    d <- grid$d[i]
    denom <- 5 * a + c_par + 5 * d
    p_mpt <- .mpt_probability_vector(
      pars = c(
        Pi = denom / (15 * b + denom),
        Pb = c_par / denom,
        Pd = d / (a + d)
      ),
      mpt_model = mpt_dist
    )
    p_m3 <- .compute_m3_probability_vector(
      pars = c(a = a, c = c_par, d = d, b = b),
      m3_model = m3_dist, act_funs = acts_dist
    )
    max(abs(p_mpt - p_m3))
  }, numeric(1)))
  expect_lt(max_diff, 1e-10)

  # no-distractor condition: counts (1, 4, 10), S = 15b + 5a + c
  mpt_nodist <- mpt(mpt_tree("nodist", list(
    corr = "Pi * Pb + Pi * (1 - Pb) * (1/5) + (1 - Pi) * (1/15)",
    other = "Pi * (1 - Pb) * (4/5) + (1 - Pi) * (4/15)",
    npl = "(1 - Pi) * (10/15)"
  )))
  m3_nodist <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 10), choice_rule = "simple", version = "custom"
  )
  acts_nodist <- bmf(corr ~ b + a + c, other ~ b + a, npl ~ b)
  max_diff_nodist <- max(vapply(seq_len(nrow(grid)), function(i) {
    a <- grid$a[i]
    c_par <- grid$c[i]
    p_mpt <- .mpt_probability_vector(
      pars = c(
        Pi = (5 * a + c_par) / (15 * b + 5 * a + c_par),
        Pb = c_par / (5 * a + c_par)
      ),
      mpt_model = mpt_nodist
    )
    p_m3 <- .compute_m3_probability_vector(
      pars = c(a = a, c = c_par, b = b),
      m3_model = m3_nodist, act_funs = acts_nodist
    )
    max(abs(p_mpt - p_m3))
  }, numeric(1)))
  expect_lt(max_diff_nodist, 1e-10)
})

test_that("conditional_effects() shows mpt parameters on the native scale", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  # unequal trials per row, with a mean that is not an integer
  dat <- data.frame(
    item_type = rep(c("old", "new"), 2), cond = rep(c("x", "y"), each = 2),
    old = c(70, 20, 45, 30), new = c(30, 61, 15, 70)
  )
  fit <- bmm(
    bmf(D ~ cond, g ~ 1), dat, mpt(mpt_2htm_trees(), "item_type"),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, silent = 2
  )

  ce <- conditional_effects(fit, par = "D", robust = TRUE)
  np <- native_parameters(fit, pars = "D")
  medians <- c(tapply(np$value, np$cond, stats::median))
  expect_equal(ce$cond$estimate__, unname(medians[as.character(ce$cond$cond)]), tolerance = 1e-3)
  expect_named(conditional_effects(fit), "D.cond")
})

test_that("conditional_effects() shows simplex members on the probability scale", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  trees <- list(
    mpt_tree("x", list(A = "D + (1 - D) * a", B = "(1 - D) * b", C = "(1 - D) * c")),
    mpt_tree("y", list(A = "(1 - D) * a", B = "D + (1 - D) * b", C = "(1 - D) * c"))
  )
  dat <- data.frame(
    tt = rep(c("x", "y"), 2), cond = rep(c("p", "q"), each = 2),
    A = c(60, 15, 50, 30), B = c(25, 55, 20, 40), C = c(15, 21, 30, 35)
  )
  fit <- bmm(
    bmf(D ~ 1, a ~ cond, b ~ cond), dat, mpt(trees, "tt", simplex = c("a", "b", "c")),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, silent = 2
  )

  np <- native_parameters(fit, pars = c("a", "b"))
  # b = (1 - a) * stick, so its panel depends on both stick-breaking components
  for (p in c("a", "b")) {
    ce <- conditional_effects(fit, par = p, robust = TRUE)$cond
    medians <- c(tapply(np$value[np$parameter == p], np$cond[np$parameter == p], stats::median))
    expect_true(all(ce$lower__ > 0 & ce$upper__ < 1))
    expect_equal(ce$estimate__, unname(medians[as.character(ce$cond)]), tolerance = 1e-3)
  }
})

test_that("mpt keeps the trees as given before the restrictions", {
  trees <- setNames(mpt_2htm_trees(), c("old", "new"))
  restricted <- mpt(trees, tree_id = "item_type", restrictions = "g = 0.5")
  expect_equal(restricted$other_vars$unrestricted_trees, trees)
  expect_equal(deparse(restricted$other_vars$trees$old$branches$old), "D + (1 - D) * 0.5")
  expect_equal(deparse(restricted$other_vars$trees$new$branches$new), "D + (1 - D) * 0.5")
  expect_null(mpt(trees, tree_id = "item_type")$other_vars$unrestricted_trees)
})
