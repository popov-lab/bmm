test_that("impossible categories reach the Stan code as a -100 predictor", {
  model <- mpt(mpt_impossible_trees(), tree_id = "tree")
  dat <- mpt_impossible_data()
  formula <- bmf(Pm ~ 1, Pb ~ 1)

  # brms replaces the indicator columns by C_<parameter>_<i>[n] in the code
  code <- stancode(formula, data = dat, model = model)
  expect_match(
    code,
    "mudist[n] = (C_mudist_1[n] * log(nlp_dist[n]) + (1 - C_mudist_1[n]) * ( - 100))",
    fixed = TRUE
  )
  # the placeholder keeps log() defined for the tree without the category
  expect_match(code, "nlp_dist[n] = (C_dist_1[n] * ((1 - inv_logit(nlp_Pm[n])) * 0.2) + C_dist_2[n] * (1))", fixed = TRUE)
  expect_match(code, "mucorr[n] = (log(nlp_corr[n]))", fixed = TRUE)
  expect_match(code, "multinomial_logit2_lpmf(Y[n] | mu[n])", fixed = TRUE)

  # brms passes the generated indicator columns as non-linear covariates,
  # named per category formula in writing order: Idx_withdist, Idx_nodist
  # for each category and Poss_dist for the guarded linear predictor
  sdata <- standata(formula, data = dat, model = model)
  expect_equal(sdata$N, nrow(dat))
  expect_equal(as.integer(sdata$C_corr_1), as.integer(dat$tree == "withdist"))
  expect_equal(as.integer(sdata$C_corr_2), as.integer(dat$tree == "nodist"))
  expect_equal(as.integer(sdata$C_mudist_1), as.integer(dat$cond == "withdist"))
})

test_that("decimal constants are emitted literally, never in scientific notation", {
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1 - 0.001)",
    b = "(1 - p) * 0.001"
  ))
  dat <- data.frame(a = c(30L, 28L), b = c(0L, 2L))
  code <- stancode(bmf(p ~ 1), data = dat, model = mpt(tree))
  expect_match(code, "0.001", fixed = TRUE)
  expect_false(grepl("e-", code, fixed = TRUE))
})

test_that("simplex sticks get the beta priors of a symmetric Dirichlet(1)", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  tree <- mpt_tree("t", list(A = "gA", B = "gB", C = "gC", D = "gD"))
  dat <- data.frame(cond = factor(c("x", "y")), A = 5L, B = 3L, C = 2L, D = 1L)
  formula <- bmf(gA ~ 1, gB ~ 1, gC ~ 1)
  for (link in c("logit", "probit")) {
    model <- mpt(tree, simplex = c("gA", "gB", "gC", "gD"), links = link)
    code <- suppressMessages(stancode(formula, dat, model))
    expect_match(code, paste0("real ", link, "beta_lpdf(real x"), fixed = TRUE)
    expect_match(code, paste0("real ", link, "beta_rng("), fixed = TRUE)
    expect_match(code, paste0(link, "beta_lpdf(b_gAraw[1] | 1, 3)"), fixed = TRUE)
    expect_match(code, paste0(link, "beta_lpdf(b_gBraw[1] | 1, 2)"), fixed = TRUE)
    expect_match(code, "_lpdf(b_gCraw[1] | 0, 1)", fixed = TRUE)
  }
  # under cell-means coding every cell of a stick gets its beta prior
  code <- suppressMessages(stancode(
    bmf(gA ~ 0 + cond, gB ~ 0 + cond, gC ~ 1), dat,
    mpt(tree, simplex = c("gA", "gB", "gC", "gD"))
  ))
  expect_match(code, "logitbeta_lpdf(b_gAraw | 1, 3)", fixed = TRUE)
})

test_that("the per-stick uniform prior of ?mpt replaces the beta priors in every cell", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  tree <- mpt_tree("t", list(A = "gA", B = "gB", C = "gC", D = "gD"))
  dat <- data.frame(cond = factor(c("x", "y")), A = 5L, B = 3L, C = 2L, D = 1L)
  latent <- c(logit = "logistic(0, 1)", probit = "normal(0, 1)")
  formulas <- list(bmf(gA ~ 1, gB ~ 1, gC ~ 1), bmf(gA ~ 0 + cond, gB ~ 0 + cond, gC ~ 1))
  for (link in names(latent)) {
    model <- mpt(tree, simplex = c("gA", "gB", "gC", "gD"), links = link)
    treebugs <- brms::set_prior(latent[[link]], nlpar = c("gAraw", "gBraw"))
    for (formula in formulas) {
      code <- suppressMessages(stancode(formula, dat, model, prior = treebugs))
      expect_no_match(code, "beta_lpdf\\(b_")
    }
  }
})

test_that("models without a simplex group carry no simplex prior functions", {
  withr::local_seed(445)
  pcm_tree <- mpt_tree("study", list(
    C = "cp + (1 - cp) * rp * rp",
    E = "2 * (1 - cp) * rp * (1 - rp)",
    U = "(1 - cp) * (1 - rp) * (1 - rp)"
  ))
  pcm_dat <- data.frame(id = factor(1:2), C = 20L, E = 12L, U = 8L)
  codes <- list(
    suppressMessages(stancode(
      bmf(D ~ 1 + (1 | id), g ~ 1), mpt_2htm_data(n_id = 2),
      mpt(mpt_2htm_trees(), tree_id = "item_type")
    )),
    suppressMessages(stancode(bmf(cp ~ 1 + (1 | id), rp ~ 1), pcm_dat, mpt(pcm_tree)))
  )
  for (code in codes) {
    expect_no_match(code, "beta_lpdf|beta_rng")
  }
})

# Golden Stan code: any change to the program of these three models fails the
# snapshots in _snaps/model_mpt-stancode.md. After an intended change, read the
# diff that the failure prints, then accept it with
# testthat::snapshot_accept("model_mpt-stancode") (or
# snapshot_review("model_mpt-stancode") to see both versions side by side).
# The header line naming the brms and bmm versions is left out.
# brms writes part of the program, so the snapshots skip under any brms other
# than golden_brms: a brms release would otherwise fail every open PR at once,
# and tests/local/test-golden-stancode.R fails weekly until they are re-recorded.
# After a brms release, set golden_brms to the new version, run this file under
# it, read the diff and accept it as above.
golden_brms <- "2.23.0"

skip_unless_golden_brms <- function() {
  installed <- as.character(utils::packageVersion("brms"))
  skip_if(
    installed != golden_brms,
    paste0("golden Stan code recorded under brms ", golden_brms, ", installed brms ", installed)
  )
}

golden_stancode <- function(formula, data, model) {
  code <- suppressMessages(stancode(formula, data, model))
  cat(sub("^// generated with [^\n]*\n", "", code))
}

test_that("the Stan code of the two-high-threshold model is unchanged", {
  skip_on_cran() # expect_snapshot() skips there anyway; this skips building the code
  skip_unless_golden_brms()
  withr::local_seed(255)
  expect_snapshot(golden_stancode(
    bmf(D ~ 1 + (1 | id), g ~ 1), mpt_2htm_data(n_id = 2),
    mpt(mpt_2htm_trees(), tree_id = "item_type")
  ))
})

test_that("the Stan code of the pair-clustering model is unchanged", {
  skip_on_cran() # expect_snapshot() skips there anyway; this skips building the code
  skip_unless_golden_brms()
  pcm_tree <- mpt_tree("study", list(
    C = "cp + (1 - cp) * rp * rp",
    E = "2 * (1 - cp) * rp * (1 - rp)",
    U = "(1 - cp) * (1 - rp) * (1 - rp)"
  ))
  expect_snapshot(golden_stancode(
    bmf(cp ~ 1 + (1 | id), rp ~ 1), data.frame(id = factor(1:2), C = 20L, E = 12L, U = 8L),
    mpt(pcm_tree)
  ))
})

test_that("the Stan code of the Oberauer and Lewandowsky (2019) model is unchanged", {
  skip_on_cran() # expect_snapshot() skips there anyway; this skips building the code
  skip_unless_golden_brms()
  # the model and formula of the real-data section of the MPT article
  dat <- oberauer_lewandowsky_2019_e1
  dat$cond <- factor(dat$cond,
    levels = c("new distractors", "old reordered", "old same"),
    labels = c("newdist", "oldreord", "oldsame")
  )
  dat$tree <- ifelse(dat$cond == "newdist", "newdist", "nodist")
  trees <- list(
    mpt_tree("newdist", list(
      corr  = "Pi * Pb + Pi * (1 - Pb) * (1 - Pd) * (1/5) + (1 - Pi) * (1/15)",
      other = "Pi * (1 - Pb) * (1 - Pd) * (4/5) + (1 - Pi) * (4/15)",
      dist  = "Pi * (1 - Pb) * Pd + (1 - Pi) * (5/15)",
      npl   = "(1 - Pi) * (5/15)"
    )),
    mpt_tree("nodist", list(
      corr  = "Pi * Pb + Pi * (1 - Pb) * (1/5) + (1 - Pi) * (1/15)",
      other = "Pi * (1 - Pb) * (4/5) + (1 - Pi) * (4/15)",
      npl   = "(1 - Pi) * (10/15)"
    ), impossible = "dist")
  )
  expect_snapshot(golden_stancode(
    bmf(
      Pi ~ 0 + cond + (0 + cond || ID),
      Pb ~ 0 + cond + (0 + cond || ID),
      Pd ~ 1 + (1 | ID)
    ),
    dat, mpt(trees, tree_id = "tree")
  ))
})
