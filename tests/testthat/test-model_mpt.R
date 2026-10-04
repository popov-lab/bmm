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
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  expect_equal(model$other_vars$link, "logit")
  expect_equal(model$other_vars$indicators$tree, c(old = "Idx_old", new = "Idx_new"))
  expect_setequal(names(model$parameters), c("D", "g"))
  expect_equal(model$links$D, "logit")
  expect_equal(model$default_priors$D$main, "logistic(0, 1)")
  expect_equal(model$default_priors$D$effects, "logistic(0, 1)")

  # the recorded call differs by construction; every other field must match
  without_call <- function(m) {
    attr(m, "call") <- NULL
    m
  }
  rebuilt <- do.call("mpt", .mpt_constructor_args(model))
  expect_equal(without_call(rebuilt), without_call(model))

  single <- mpt(mpt_tree("t", list(A = "p", B = "1 - p")), links = "probit")
  expect_null(single$other_vars$indicators$tree)
  expect_equal(single$links$p, "probit")
  expect_equal(single$default_priors$p$main, "normal(0, 1)")
  expect_equal(single$default_priors$p$effects, "normal(0, 1)")
  rebuilt_single <- do.call("mpt", .mpt_constructor_args(single))
  expect_equal(without_call(rebuilt_single), without_call(single))
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
