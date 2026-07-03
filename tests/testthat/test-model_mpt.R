mpt_2htm_trees <- function() {
  list(
    mpt_tree("old", list(
      old = "D + (1 - D) * g",
      new = "(1 - D) * (1 - g)"
    )),
    mpt_tree("new", list(
      old = "(1 - D) * g",
      new = "D + (1 - D) * (1 - g)"
    ))
  )
}

mpt_2htm_data <- function(n_id = 10, n_items = 50) {
  dat <- expand.grid(
    id = factor(seq_len(n_id)), item_type = c("old", "new"),
    stringsAsFactors = FALSE
  )
  p_old <- ifelse(dat$item_type == "old", 0.85, 0.15)
  dat$old <- rbinom(nrow(dat), n_items, p_old)
  dat$new <- n_items - dat$old
  dat
}

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

test_that("mpt_tree folds integer fractions into decimal literals", {
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1/4)",
    b = "(1 - p) * (3/4)"
  ))
  expect_false(grepl("/", tree$branches$a, fixed = TRUE))
  expect_true(grepl("0.25", tree$branches$a, fixed = TRUE))
  expect_true(grepl("0.75", tree$branches$b, fixed = TRUE))
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

test_that("mpt uses matched priors for the probit link", {
  model <- mpt(mpt_2htm_trees(), condition = "item_type", links = "probit")
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

test_that("mpt errors on inconsistent trees and missing condition", {
  trees <- list(
    mpt_tree("t1", list(a = "p", b = "1 - p")),
    mpt_tree("t2", list(a = "p", c = "1 - p"))
  )
  expect_error(mpt(trees, condition = "cond"), "same response categories")
  expect_error(mpt(mpt_2htm_trees()), "require the condition argument")
  dup_trees <- list(
    mpt_tree("t1", list(a = "p", b = "1 - p")),
    mpt_tree("t1", list(a = "p", b = "1 - p"))
  )
  expect_error(mpt(dup_trees, condition = "cond"), "unique")
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
    mpt(trees, condition = "item_type", simplex = c("g", "x")),
    "Unknown"
  )
  expect_error(
    mpt(trees, condition = "item_type", simplex = "g"),
    "at least two parameters"
  )
})

test_that("mpt warns when branch probabilities do not sum to 1", {
  bad_tree <- mpt_tree("t", list(a = "D * g", b = "(1 - D) * g"))
  expect_warning(mpt(bad_tree), "sum to")
})

test_that("check_formula generates linked category formulas", {
  model <- mpt(mpt_2htm_trees(), condition = "item_type")
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
  model <- mpt(mpt_2htm_trees(), condition = "item_type", links = "probit")
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
  model <- mpt(mpt_2htm_trees(), condition = "item_type")
  dat <- mpt_2htm_data()
  formula <- bmf(old ~ 1, D ~ 1, g ~ 1)
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE),
    "cannot be predicted directly"
  )
})

test_that("mpt compiles for a multi-tree binary-category model", {
  model <- mpt(mpt_2htm_trees(), condition = "item_type")
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

test_that("fixed parameter values are transformed to the latent scale", {
  model <- mpt(mpt_2htm_trees(), condition = "item_type")
  dat <- mpt_2htm_data()
  formula <- bmf(D ~ 1, g = 0.5)
  model_checked <- suppressMessages(check_model(model, dat, formula))
  expect_equal(model_checked$fixed_parameters$g, qlogis(0.5))

  prior <- suppressMessages(default_prior(formula, dat, model))
  constant_row <- prior[grepl("constant", prior$prior), ]
  expect_equal(constant_row$nlpar, "g")
  expect_equal(constant_row$prior, "constant(0)")

  formula_bad <- bmf(D ~ 1, g = 1.5)
  expect_error(
    suppressMessages(check_model(model, dat, formula_bad)),
    "strictly"
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
  model <- mpt(trees, condition = "source", simplex = c("gA", "gB", "gNew"))
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
  model <- mpt(mpt_2htm_trees(), condition = "item_type")
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
    p_mpt <- .compute_mpt_probability_vector(
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

test_that("mpt_from_string parses MPTinR-style model definitions", {
  model_2htm <- "
  D + (1 - D) * g        # old
  (1 - D) * (1 - g)      # new

  (1 - D) * g            # old
  D + (1 - D) * (1 - g)  # new
  "
  model <- mpt_from_string(
    model_2htm, tree_names = c("old", "new"), condition = "item_type"
  )
  manual <- mpt(mpt_2htm_trees(), condition = "item_type")
  expect_equal(model$trees, manual$trees)
  expect_equal(names(model$parameters), names(manual$parameters))

  summed <- mpt_from_string(
    "D # hit\n(1 - D) * g # hit\n(1 - D) * (1 - g) # miss",
    tree_names = "old"
  )
  expect_equal(summed$trees$old$branches$hit, "(D) + ((1 - D) * g)")

  no_comments <- mpt_from_string(
    "D + (1 - D) * g\n(1 - D) * (1 - g)",
    tree_names = "old", categories = c("hit", "miss")
  )
  expect_equal(names(no_comments$trees$old$branches), c("hit", "miss"))

  expect_error(
    mpt_from_string("D # a\n1 - D # b", tree_names = c("t1", "t2")),
    "tree block"
  )
  expect_error(
    mpt_from_string("D + (1 - D) * g\n(1 - D) * (1 - g)", tree_names = "old"),
    "categories"
  )
})

test_that("mpt_from_eqn imports EQN files with restrictions and renaming", {
  eqn_file <- tempfile(fileext = ".eqn")
  writeLines(c(
    "6",
    "old  old_hit   D_o",
    "old  old_hit   (1-D_o)*g_uess*G_fix",
    "old  old_miss  (1-D_o)*(1-g_uess*G_fix)",
    "new  new_fa    (1-D_n)*g_uess*G_fix",
    "new  new_cr    D_n",
    "new  new_cr    (1-D_n)*(1-g_uess*G_fix)"
  ), eqn_file)

  category_map <- c(
    old_hit = "yes", new_fa = "yes", old_miss = "no", new_cr = "no"
  )
  model <- suppressMessages(mpt_from_eqn(
    eqn_file,
    restrictions = c(G_fix = 1 / 4),
    categories = category_map,
    condition = "item_type"
  ))
  expect_setequal(names(model$parameters), c("Do", "Dn", "guess"))
  expect_setequal(model$resp_vars$resp_cats, c("yes", "no"))
  expect_false(any(grepl("G_fix", unlist(model$trees$old$branches))))
  expect_true(grepl("0.25", model$trees$old$branches$yes, fixed = TRUE))

  renaming <- attr(model, "mpt_renaming")
  expect_equal(renaming[["D_o"]], "Do")
  expect_equal(renaming[["old_hit"]], "yes")

  expect_error(
    suppressMessages(mpt_from_eqn(
      eqn_file, restrictions = c(G_fix = "D_o"), categories = category_map
    )),
    "numeric constants"
  )

  # per-tree category labels without a mapping cannot be combined
  expect_error(
    suppressMessages(mpt_from_eqn(eqn_file, restrictions = c(G_fix = 0.25))),
    "categories argument"
  )
})

test_that("mpt_from_eqn errors on names that clash after sanitizing", {
  eqn_file <- tempfile(fileext = ".eqn")
  writeLines(c(
    "t  a  d_A + dA",
    "t  b  1 - d_A - dA"
  ), eqn_file)
  expect_error(mpt_from_eqn(eqn_file), "duplicated names")
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

test_that("factor condition columns are matched to tree names", {
  model <- mpt(mpt_2htm_trees(), condition = "item_type")
  dat <- mpt_2htm_data()
  dat$item_type <- factor(dat$item_type, levels = c("old", "new"))
  checked <- check_data(model, dat, bmf(D ~ 1, g ~ 1))
  expect_equal(checked$Idx_old, as.integer(dat$item_type == "old"))
  expect_equal(checked$Idx_new, as.integer(dat$item_type == "new"))
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
    mpt(trees, condition = "cond", covariates = c("Gcorr", "Gother"))
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
