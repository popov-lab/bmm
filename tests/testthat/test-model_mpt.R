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

test_that("mpt identifies parameters and excludes covariates", {
  tree <- mpt_tree("main", list(
    correct = "Pb + (1 - Pb) * Pi * GcorrPi",
    other = "(1 - Pb) * Pi * (1 - GcorrPi) + (1 - Pb) * (1 - Pi)"
  ))
  model <- mpt(tree, covariates = "GcorrPi")
  expect_setequal(names(model$parameters), c("Pb", "Pi"))
  expect_equal(model$other_vars$covariates, "GcorrPi")
  expect_equal(model$links$Pb, "logit")
  expect_equal(model$default_priors$Pb$main, "logistic(0, 1)")
  expect_equal(model$default_priors$Pb$effects, "logistic(0, 0.5)")
})

test_that("mpt stores its derived state once and can rebuild itself", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  expect_equal(model$link, "logit")
  expect_null(model$other_vars$link)
  expect_equal(model$indicators$tree, c(withdist = "Idx_withdist", nodist = "Idx_nodist"))
  expect_equal(model$indicators$possible, c(dist = "Poss_dist"))
  expect_null(model$simplex_raw)

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
  expect_null(single$indicators$tree)
  expect_null(single$indicators$possible)
  expect_equal(single$simplex_raw, c(gA = "gAraw", gB = "gBraw"))
  rebuilt_single <- do.call("mpt", .mpt_constructor_args(single))
  expect_equal(without_call(rebuilt_single), without_call(single))
})

test_that("importers record their own call", {
  model <- mpt_from_string(
    "D + (1 - D) * g # old\n(1 - D) * (1 - g) # new", tree_names = "old"
  )
  expect_equal(deparse(attr(model, "call")[[1]]), "mpt_from_string")
  expect_output(print(model), "mpt_from_string")
})

test_that("mpt uses matched priors for the probit link", {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type", links = "probit")
  expect_equal(model$links$D, "probit")
  expect_equal(model$default_priors$D$main, "normal(0, 1)")
  expect_equal(model$default_priors$D$effects, "normal(0, 0.5)")
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

  expect_warning(
    suppressMessages(bmm(
      formula, dat, model,
      backend = "mock", mock_fit = 1, rename = FALSE
    )),
    "Non-linear transformations"
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

test_that("simplex formula handling tolerates omitted parameter formulas", {
  trees <- list(mpt_tree("t", list(A = "gA", B = "gB", C = "gC")))
  model <- mpt(trees, simplex = c("gA", "gB", "gC"))
  moved <- .mpt_move_simplex_formulas(model, bmf(gA ~ 1))
  expect_false(any(c("gA", "gB", "gC") %in% names(moved)))
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
  expect_warning(
    suppressMessages(bmm(
      formula, dat, model,
      backend = "mock", mock_fit = 1, rename = FALSE
    )),
    "Non-linear transformations"
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

test_that("impossible categories are switched off in the linear predictor", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  formula <- bmf(Pm ~ 1, Pb ~ 1)
  checked_data <- check_data(model, dat, formula)
  checked_formula <- check_formula(model, checked_data, formula)

  # the tree that cannot produce the category contributes a placeholder, so
  # log() stays defined for its rows
  expect_match(.mpt_rhs_chr(checked_formula$dist), "Idx_nodist * (1)", fixed = TRUE)

  brms_formula <- configure_model(model, checked_data, checked_formula)$formula
  expect_match(
    .mpt_rhs_chr(brms_formula$pforms$mudist),
    "Poss_dist * log(dist) + (1 - Poss_dist) * (-100)",
    fixed = TRUE
  )
  expect_match(.mpt_rhs_chr(brms_formula$formula), "log(corr)", fixed = TRUE)
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
  # covariate columns can
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
  expect_warning(
    check_data(model, dat_bad, bmf(D ~ 1)),
    "do not sum to 1 for 1 row"
  )

  dat_na <- dat
  dat_na$Gcorr[c(2, 5)] <- NA
  expect_warning(
    check_data(model, dat_na, bmf(D ~ 1)),
    "do not sum to 1 for 2 row"
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
