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

test_that("mpt stores its derived state once and can rebuild itself", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  expect_equal(model$other_vars$link, "logit")
  expect_equal(model$other_vars$indicators$tree, c(withdist = "Idx_withdist", nodist = "Idx_nodist"))
  expect_equal(model$other_vars$indicators$possible, c(dist = "Poss_dist"))
  expect_null(model$other_vars$simplex_raw)
  expect_setequal(names(model$parameters), c("Pm", "Pb"))
  expect_equal(model$links$Pm, "logit")
  expect_equal(model$default_priors$Pm$main, "logistic(0, 1)")
  expect_equal(model$default_priors$Pm$effects, "logistic(0, 1)")

  # the recorded call differs by construction; every other field must match
  without_call <- function(m) {
    attr(m, "call") <- NULL
    m
  }
  rebuilt <- do.call("mpt", .mpt_constructor_args(model))
  expect_equal(without_call(rebuilt), without_call(model))

  single <- mpt(mpt_tree("t", list(A = "gA", B = "gB", C = "gC")),
    simplex = c("gA", "gB", "gC"), links = "probit"
  )
  expect_null(single$other_vars$indicators$tree)
  expect_null(single$other_vars$indicators$possible)
  expect_equal(single$other_vars$simplex_raw, c(gA = "gAraw", gB = "gBraw"))
  expect_equal(single$links$gA, "identity")
  expect_equal(single$default_priors$gAraw$main, "normal(0, 1)")
  expect_equal(single$default_priors$gAraw$effects, "normal(0, 1)")
  rebuilt_single <- do.call("mpt", .mpt_constructor_args(single))
  expect_equal(without_call(rebuilt_single), without_call(single))

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
  zero_product <- mpt_tree("d", list(x = "0 * D", y = "1 - 0 * D"))
  expect_s3_class(mpt(zero_product), "mpt")

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
    deviations <- .mpt_tree_sum_deviations(
      restricted$other_vars$trees, names(restricted$parameters), character(0), list()
    )
    expect_true(all(is.na(deviations)))
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
  expect_no_warning(
    suppressMessages(bmm(
      formula, dat, model,
      backend = "mock", mock_fit = 1, rename = FALSE
    )),
    message = "Non-linear"
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
  # values at construction cannot verify (mpt() warns there) but the observed
  # covariate columns can (check_data() errors there)
  tree <- mpt_tree("main", list(
    correct = "D + (1 - D) * Gcorr",
    incorrect = "(1 - D) * Gother"
  ))
  model <- suppressWarnings(mpt(tree, covariates = c("Gcorr", "Gother")))
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
    "do not sum to 1 for 2 row"
  )

  # swapped complements agree at the symmetric point D = 0.5, so the branch
  # expressions are checked at several parameter values
  swapped <- suppressWarnings(mpt(mpt_tree("main", list(
    correct = "D + (1 - D) * Gcorr", incorrect = "D * (1 - Gcorr)"
  )), covariates = "Gcorr"))
  dat_swap <- data.frame(Gcorr = 0.25, correct = 10, incorrect = 30)
  expect_error(check_data(swapped, dat_swap, bmf(D ~ 1)), "do not sum to 1")

  stray <- suppressWarnings(mpt(mpt_tree("main", list(
    correct = "D + (1 - D) * Gcorr", incorrect = "2 * D * (1 - D) * (1 - Gcorr)"
  )), covariates = "Gcorr"))
  expect_error(check_data(stray, dat_swap, bmf(D ~ 1)), "do not sum to 1")
})

test_that("check_data rejects covariate values that push a branch outside (0, 1]", {
  model <- suppressWarnings(mpt(mpt_tree("main", list(
    correct = "D + (1 - D) * G", incorrect = "(1 - D) * (1 - G)"
  )), covariates = "G"))
  dat <- data.frame(G = c(0.25, 0.5, 0.75), correct = 10, incorrect = 30)
  expect_silent(check_data(model, dat, bmf(D ~ 1)))

  dat_over <- dat
  dat_over$G[2] <- 1.2
  expect_error(
    check_data(model, dat_over, bmf(D ~ 1)),
    "'correct' in tree 'main' is [0-9.]+ in row 2, outside \\(0, 1\\]"
  )

  # a branch of exactly 0 is a structurally impossible category
  dat_one <- dat
  dat_one$G[3] <- 1
  expect_error(
    check_data(model, dat_one, bmf(D ~ 1)),
    "is 0 in row 3.*mpt_tree\\(impossible = \\)"
  )

  zero_int <- suppressWarnings(mpt(mpt_tree("main", list(
    hit = "D * G", miss = "1 - D * G"
  )), covariates = "G"))
  dat_int <- data.frame(G = c(1L, 0L), hit = 10, miss = 30)
  expect_error(
    check_data(zero_int, dat_int, bmf(D ~ 1)),
    "'hit' in tree 'main' is 0 in row 2"
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
  model <- suppressWarnings(mpt(trees, tree_id = "tree", covariates = "Gd"))
  dat <- mpt_impossible_data()
  dat$Gd <- ifelse(dat$tree == "withdist", 0.3, 0)
  expect_silent(check_data(model, dat, bmf(Pm ~ 1, Pb ~ 1)))

  dat_na <- dat
  dat_na$Gd[dat_na$tree == "nodist"] <- NA
  first_row <- which(dat_na$tree == "nodist")[1]
  expect_error(
    check_data(model, dat_na, bmf(Pm ~ 1, Pb ~ 1)),
    glue::glue("'Gd' is missing in row {first_row} of tree 'nodist'")
  )

  # in a tree that uses the covariate, the sum check reports it first
  dat_used <- dat
  dat_used$Gd[dat_used$tree == "withdist"][1] <- NA
  expect_error(
    check_data(model, dat_used, bmf(Pm ~ 1, Pb ~ 1)),
    "do not sum to 1 for 1 row"
  )
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
  model <- suppressWarnings(
    mpt(trees, tree_id = "cond", covariates = c("Gcorr", "Gother"))
  )
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
