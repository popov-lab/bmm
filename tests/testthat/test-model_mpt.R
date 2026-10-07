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
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  expect_equal(model$other_vars$link, "logit")
  expect_equal(model$other_vars$indicators$tree, c(old = "Idx_old", new = "Idx_new"))
  expect_setequal(names(model$parameters), c("D", "g"))
  expect_equal(model$links$D, "logit")
  expect_equal(model$default_priors$D$main, "logistic(0, 1)")
  expect_equal(model$default_priors$D$effects, "logistic(0, 1)")

  single <- mpt(mpt_tree("t", list(A = "p", B = "1 - p")), links = "probit")
  expect_null(single$other_vars$indicators$tree)
  expect_equal(single$links$p, "probit")
  expect_equal(single$default_priors$p$main, "normal(0, 1)")
  expect_equal(single$default_priors$p$effects, "normal(0, 1)")
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

test_that("mpt errors on a branch that is the constant 0", {
  zero_padded <- list(
    mpt_tree("a", list(x = "D", y = "1 - D", z = "0")),
    mpt_tree("b", list(x = "1 - 1", y = "g", z = "1 - g"))
  )
  expect_error(mpt(zero_padded, tree_id = "t"), "constant 0")
  expect_error(mpt(zero_padded, tree_id = "t"), "'z' in tree 'a', 'x' in tree 'b'")
  expect_error(mpt(mpt_tree("c", list(x = "(0)", y = "D + (1 - D)"))), "constant 0")
  # not the literal 0 this guard reads, but 0 at every interior test point,
  # which the range check refuses
  zero_product <- mpt_tree("d", list(x = "0 * D", y = "1 - 0 * D"))
  expect_error(mpt(zero_product), "category 'x' in tree 'd' is 0 .*outside \\(0, 1\\]")
})

test_that("mpt errors on name collisions and reserved names", {
  tree_collision <- mpt_tree("t", list(D = "D + g", other = "1 - D - g"))
  expect_error(mpt(tree_collision), "both a parameter and a response category")
  tree_reserved <- mpt_tree("t", list(Y = "p", other = "1 - p"))
  expect_error(mpt(tree_reserved), "reserved")
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

test_that("response categories are refused as predictors before brms sees them", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  expect_error(
    bmm(bmf(D ~ old, g ~ 1), dat, model, backend = "mock", mock_fit = 1, rename = FALSE),
    "response counts 'old' .* cannot be predictors.*formula\\(s\\) for: 'D'"
  )
  expect_error(
    check_model(model, dat, bmf(D ~ inv_logit(k * new), k ~ 1, g ~ 1)),
    "response counts 'new' .* cannot be predictors.*formula\\(s\\) for: 'D'"
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

test_that("generated indicator columns do not trigger the clash warning", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  expect_no_warning(bmm(
    bmf(D ~ 1 + (1 | id), g ~ 1), dat, model,
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
  # the rank covers the tree parameters only, so u + v would pass unnoticed
  expect_message(
    check_model(model, dat, formula),
    "'a', 'b' .*Whether the data identify these parameters is not checked"
  )
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

test_that("printing an mpt model lists trees and the identifiability bound", {
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
  expect_message(
    check_model(model, dat, bmf(Do ~ 1, Dn ~ Do, g ~ 0 + bias)),
    "within one design cell.*'g' identify it across cells or the .* for 'Dn'"
  )

  # the rank does not read non-linear formulas, so a deficit under one is
  # announced, not warned about, even when it is real as here
  expect_no_warning(expect_message(
    check_model(model, dat, bmf(Do ~ inv_logit(phi), phi ~ 1, Dn ~ 1, g ~ 1)),
    "by the branch expressions alone \\(Jacobian rank 2 for 3 free.*'Do' identify it is not checked"
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
    "by the branch expressions alone .*formula\\(s\\) for 'Dn' identify it is not checked"
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
})

test_that("a formula that reaches no parameter of the deficit keeps the warning", {
  model <- mpt(list(
    mpt_tree("a", list(x = "D * r", y = "(1 - D * r) * g", z = "(1 - D * r) * (1 - g)")),
    mpt_tree("b", list(x = "g", y = "(1 - g) * h", z = "(1 - g) * (1 - h)"))
  ), tree_id = "tree")
  dat <- data.frame(tree = rep(c("a", "b"), each = 3), x = 5L, y = 5L, z = 5L)

  # D and r enter only as D * r; formulas for h cannot separate them
  expect_warning(
    suppressMessages(check_model(model, dat, bmf(D ~ 1, r ~ 1, g ~ 1, h ~ inv_logit(k), k ~ 1))),
    "combination\\(s\\) of 'D', 'r' cannot be estimated"
  )
  expect_warning(
    check_model(model, dat, bmf(D ~ 1, r ~ 1, g ~ 1, h ~ g)),
    "combination\\(s\\) of 'D', 'r' cannot be estimated"
  )
  expect_no_warning(expect_message(
    check_model(model, dat, bmf(D ~ 1, r ~ D, g ~ 1, h ~ 1)),
    "non-linear formula\\(s\\) for 'r' identify it is not checked"
  ))
  # h ~ D reads D, which tree b identifies through h
  expect_no_warning(expect_message(
    check_model(model, dat, bmf(D ~ 1, r ~ 1, g ~ 1, h ~ D)),
    "non-linear formula\\(s\\) for 'h' identify it is not checked"
  ))
})

test_that("a predictor on a parameter outside the deficit is not met with a warning", {
  # h is identified by tree t2; within one cell D and r enter only through
  # h * r + (1 - h) * D, but two values of h separate them
  model <- mpt(list(
    mpt_tree("t1", list(yes = "h * r + (1 - h) * D", no = "h * (1 - r) + (1 - h) * (1 - D)")),
    mpt_tree("t2", list(yes = "h", no = "1 - h"))
  ), tree_id = "tree")
  dat <- expand.grid(tree = c("t1", "t2"), cond = c("A", "B"), stringsAsFactors = FALSE)
  dat$yes <- 20L
  dat$no <- 20L

  expect_no_warning(expect_message(
    check_model(model, dat, bmf(D ~ 1, r ~ 1, h ~ 0 + cond)),
    "within one design cell.*predictors on 'h' identify it across cells is not checked"
  ))
  expect_no_warning(suppressMessages(
    check_model(model, dat, bmf(D ~ 1, r ~ 1, h ~ inv_logit(phi), phi ~ cond))
  ))
  expect_warning(
    check_model(model, dat, bmf(D ~ 1, r ~ 1, h ~ 1)),
    "all free parameters except 'h' cannot be estimated"
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

test_that("factor tree identifier columns are matched to tree names", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  dat <- mpt_2htm_data()
  dat$item_type <- factor(dat$item_type, levels = c("old", "new"))
  checked <- check_data(model, dat, bmf(D ~ 1, g ~ 1))
  expect_equal(checked$Idx_old, as.integer(dat$item_type == "old"))
  expect_equal(checked$Idx_new, as.integer(dat$item_type == "new"))
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
