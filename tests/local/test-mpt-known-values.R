# Known-value fits of mpt(): one per model feature. The data come from
# rmultinom() on branch probabilities written out below (draw_counts() in
# helper-mpt.R), with large counts so that each posterior is tight, and
# expect_known_values() reads the estimates on the probability scale. Each fit
# takes about a second once its Stan program is compiled; compiling it takes
# 10 to 40 seconds.

fit_known <- function(formula, data, model, seed) {
  bmm(
    formula, data, model,
    backend = "cmdstanr", chains = 2, iter = 2000, warmup = 1000,
    refresh = 0, seed = seed
  )
}

# n: the trials of both trees, or of the old-item and the new-item tree
two_htm_counts <- function(D_old, D_new, g, n) {
  n <- rep_len(n, 2)
  rbind(
    draw_counts(c(old = D_old + (1 - D_old) * g, new = (1 - D_old) * (1 - g)), n[1]),
    draw_counts(c(old = (1 - D_new) * g, new = D_new + (1 - D_new) * (1 - g)), n[2])
  )
}

test_that("an mpt fit with an impossible category recovers its parameters", {
  withr::local_seed(11)
  Pm <- 0.6
  Pb <- 0.7
  with_dist <- c(corr = Pm * Pb + (1 - Pm) * 0.2, dist = (1 - Pm) * 0.2,
                 npl = Pm * (1 - Pb) + (1 - Pm) * 0.6)
  no_dist <- c(corr = Pm * Pb + (1 - Pm) * 0.2, dist = 0, npl = Pm * (1 - Pb) + (1 - Pm) * 0.8)
  dat <- data.frame(tree = c("withdist", "nodist", "nodist"))
  counts <- rbind(draw_counts(with_dist, 4000), draw_counts(no_dist, 4000), draw_counts(no_dist, 4000))
  dat <- cbind(dat, counts)

  fit <- fit_known(bmf(Pm ~ 1, Pb ~ 1), dat, mpt(mpt_impossible_trees(), "tree"), seed = 11)
  expect_known_values(fit, list(Pm = Pm, Pb = Pb))
})

test_that("an mpt fit with a guessing rate given as data recovers its parameter", {
  withr::local_seed(12)
  D <- 0.55
  guessing <- c(0.2, 0.25, 0.5)
  # unequal trials per row: with equal ones, covariate values read in the
  # wrong row order shift D by about one posterior SD and pass
  counts <- t(vapply(seq_along(guessing), function(i) {
    Gcorr <- guessing[i]
    draw_counts(c(correct = D + (1 - D) * Gcorr, incorrect = (1 - D) * (1 - Gcorr)), c(3000, 1500, 500)[i])
  }, integer(2)))
  dat <- data.frame(Gcorr = guessing, counts)
  model <- mpt(
    mpt_tree("main", list(correct = "D + (1 - D) * Gcorr", incorrect = "(1 - D) * (1 - Gcorr)")),
    covariates = "Gcorr"
  )

  fit <- fit_known(bmf(D ~ 1), dat, model, seed = 12)
  expect_known_values(fit, list(D = D))
})

test_that("an mpt fit with a fixed and an equality restriction recovers the free parameter", {
  withr::local_seed(13)
  D <- 0.65
  # with equal trials per tree a wrong guessing constant moves the old-item
  # and the new-item estimate of D in opposite directions and almost cancels;
  # at 12000/2000 a constant off by 0.1 moves D by about eight posterior SDs
  dat <- data.frame(item_type = c("old", "new"), two_htm_counts(D, D, g = 0.5, n = c(12000, 2000)))
  trees <- list(
    mpt_tree("old", list(old = "Do + (1 - Do) * g", new = "(1 - Do) * (1 - g)")),
    mpt_tree("new", list(old = "(1 - Dn) * g", new = "Dn + (1 - Dn) * (1 - g)"))
  )
  model <- mpt(trees, tree_id = "item_type", restrictions = c("g = 0.5", "Dn = Do"))
  expect_named(model$parameters, "Do")

  fit <- fit_known(bmf(Do ~ 1), dat, model, seed = 13)
  expect_known_values(fit, list(Do = D))
})

test_that("an mpt fit on the probit link recovers an intercept and an effect", {
  withr::local_seed(14)
  D <- c(x = 0.55, y = 0.8)
  g <- 0.4
  dat <- data.frame(
    item_type = rep(c("old", "new"), 2), cond = rep(c("x", "y"), each = 2),
    rbind(two_htm_counts(D[["x"]], D[["x"]], g, 1500), two_htm_counts(D[["y"]], D[["y"]], g, 1500))
  )

  fit <- fit_known(bmf(D ~ 1 + cond, g ~ 1), dat, mpt(mpt_2htm_trees(), "item_type", links = "probit"), seed = 14)
  expect_known_values(fit, list(D = unname(D), g = rep(g, 2)), newdata = data.frame(cond = c("x", "y")))
})

test_that("an mpt fit with a random intercept recovers the population and each participant", {
  withr::local_seed(15)
  n_id <- 8
  D_mean <- 0.7
  g <- 0.4
  D_id <- stats::plogis(stats::qlogis(D_mean) + stats::rnorm(n_id, 0, 0.5))
  dat <- data.frame(
    id = factor(rep(seq_len(n_id), each = 2)), item_type = rep(c("old", "new"), n_id),
    do.call(rbind, lapply(D_id, function(D) two_htm_counts(D, D, g, 1000)))
  )

  fit <- fit_known(bmf(D ~ 1 + (1 | id), g ~ 1), dat, mpt(mpt_2htm_trees(), "item_type"), seed = 15)
  # the population median is plogis() of the mean the participants were drawn
  # from; with eight participants its posterior is wide, hence max_rel_sd
  expect_known_values(fit, list(D = D_mean), max_rel_sd = 0.1)
  expect_known_values(fit, list(g = g))
  expect_known_values(
    fit, list(D = D_id), newdata = data.frame(id = factor(seq_len(n_id))), re_formula = NULL
  )
})

test_that("an mpt fit with a non-linear formula for a tree parameter recovers it", {
  withr::local_seed(16)
  a <- 0.3
  b <- 1.2
  g <- 0.4
  D <- stats::plogis(a + b * c(0, 1))
  dat <- data.frame(
    item_type = rep(c("old", "new"), 2), x = rep(c(0, 1), each = 2),
    rbind(two_htm_counts(D[1], D[1], g, 1500), two_htm_counts(D[2], D[2], g, 1500))
  )
  formula <- bmf(D ~ inv_logit(a + b * x), a ~ 1, b ~ 1, g ~ 1)

  fit <- suppressMessages(fit_known(formula, dat, mpt(mpt_2htm_trees(), "item_type"), seed = 16))
  expect_known_values(fit, list(D = D, g = rep(g, 2)), newdata = data.frame(x = c(0, 1)))
})

test_that("a single-tree mpt fit recovers its parameters", {
  withr::local_seed(17)
  cp <- 0.4
  rp <- 0.6
  probs <- c(C = cp + (1 - cp) * rp^2, E = 2 * (1 - cp) * rp * (1 - rp), U = (1 - cp) * (1 - rp)^2)
  dat <- data.frame(rbind(draw_counts(probs, 4000), draw_counts(probs, 4000)))
  model <- mpt(mpt_tree("study", list(
    C = "cp + (1 - cp) * rp * rp",
    E = "2 * (1 - cp) * rp * (1 - rp)",
    U = "(1 - cp) * (1 - rp) * (1 - rp)"
  )))

  fit <- fit_known(bmf(cp ~ 1, rp ~ 1), dat, model, seed = 17)
  expect_known_values(fit, list(cp = cp, rp = rp))
})

test_that("an mpt model imported from an EQN file recovers its parameters", {
  withr::local_seed(18)
  D <- 0.7
  g <- 0.4
  eqn_file <- withr::local_tempfile(fileext = ".eqn")
  writeLines(c(
    "6",
    "old hit D",
    "old hit (1-D)*g",
    "old miss (1-D)*(1-g)",
    "new fa (1-D)*g",
    "new cr D",
    "new cr (1-D)*(1-g)"
  ), eqn_file)
  model <- suppressMessages(mpt_from_eqn(
    eqn_file, categories = c(hit = "old", fa = "old", miss = "new", cr = "new"),
    tree_id = "item_type"
  ))
  dat <- data.frame(item_type = c("old", "new"), two_htm_counts(D, D, g, 1500))

  fit <- fit_known(bmf(D ~ 1, g ~ 1), dat, model, seed = 18)
  expect_known_values(fit, list(D = D, g = g))
})
