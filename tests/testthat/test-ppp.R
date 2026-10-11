# A bmmfit stand-in for ppp() with fixed draws: .ppp_counts() is mocked, so
# only the model class, the formula and the data are read from the fit
ppp_toy_fit <- function(data, formula = bmf(D ~ 1, g ~ 1)) {
  structure(
    list(bmm = list(model = mpt(mpt_2htm_trees(), "item_type"), user_formula = formula),
         data = data),
    class = c("bmmfit", "brmsfit")
  )
}

# draws given as a list of rows x categories matrices, one per draw
ppp_toy_draws <- function(...) {
  aperm(simplify2array(list(...)), c(3L, 1L, 2L))
}

local_ppp_counts <- function(observed, expected, replicated, data,
                             env = parent.frame()) {
  counts <- list(observed = observed, expected = expected,
                 replicated = replicated, data = data)
  local_mocked_bindings(.ppp_counts = function(model, fit, draw_ids) counts,
                        .env = env)
}

# Klauer (2010), Eq. 18, written out: cells is a list of participants x
# categories matrices of expected counts, one per cell, with trials n per
# participant and cell
klauer_t2 <- function(observed, cells, n) {
  n_people <- nrow(observed)
  blocks <- lapply(seq_along(cells), function(j) {
    Reduce(`+`, lapply(seq_len(n_people), function(t) {
      e <- cells[[j]][t, ]
      diag(e, length(e)) - tcrossprod(e) / n[t, j]
    })) / n_people
  })
  within <- matrix(0, ncol(observed), ncol(observed))
  offset <- 0
  for (b in blocks) {
    idx <- offset + seq_len(ncol(b))
    within[idx, idx] <- b
    offset <- offset + ncol(b)
  }
  sigma <- stats::cov(do.call(cbind, cells)) + (n_people - 1) / n_people * within
  sum((stats::cov(observed) - sigma)^2 / sqrt(outer(diag(sigma), diag(sigma))))
}

test_that("ppp() computes T1 by hand on a two-row, two-category example", {
  data <- data.frame(row = c("a", "b"))
  local_ppp_counts(
    observed = matrix(c(3, 6, 7, 4), 2, dimnames = list(NULL, c("old", "new"))),
    expected = ppp_toy_draws(matrix(c(4, 5, 6, 5), 2), matrix(c(6, 5, 4, 5), 2)),
    replicated = ppp_toy_draws(matrix(c(5, 5, 5, 5), 2), matrix(c(3, 7, 7, 3), 2)),
    data = data
  )
  fit <- ppp_toy_fit(data)

  # pooled: observed (9, 11); expected (9, 11) and (11, 9); replicated (10, 10)
  pooled <- ppp(fit, statistic = "T1", by = character(0))
  t1_obs <- c(0, 2^2 / 11 + 2^2 / 9)
  t1_rep <- c(1 / 9 + 1 / 11, 1 / 11 + 1 / 9)
  expect_equal(attr(pooled, "draws")$T1[, "observed"], t1_obs)
  expect_equal(attr(pooled, "draws")$T1[, "replicated"], t1_rep)
  expect_equal(pooled$observed, mean(t1_obs))
  expect_equal(pooled$replicated, mean(t1_rep))
  expect_equal(pooled$ppp, 0.5)
  expect_equal(pooled$ndraws, 2L)

  # one cell per row
  per_row <- ppp(fit, statistic = "T1", by = "row")
  expect_equal(
    attr(per_row, "draws")$T1[, "observed"],
    c(1 / 4 + 1 / 6 + 1 / 5 + 1 / 5, 3^2 / 6 + 3^2 / 4 + 1 / 5 + 1 / 5)
  )
})

test_that("ppp() computes T2 as in Klauer (2010), across cells", {
  # three participants with two cells each; trials 10 per cell
  data <- data.frame(id = factor(rep(1:3, 2)), item_type = rep(c("old", "new"), each = 3))
  observed <- matrix(c(7, 5, 9, 2, 4, 3, 3, 5, 1, 8, 6, 7), 6)
  expected <- matrix(c(6, 6, 8, 3, 3, 2, 4, 4, 2, 7, 7, 8), 6)
  replicated <- matrix(c(5, 7, 8, 3, 2, 3, 5, 3, 2, 7, 8, 7), 6)
  local_ppp_counts(observed, ppp_toy_draws(expected), ppp_toy_draws(replicated), data)
  fit <- ppp_toy_fit(data, bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id)))

  res <- ppp(fit, statistic = "T2", by = "item_type")
  # cells are sorted: "new" (rows 4:6) before "old" (rows 1:3)
  wide <- function(m) cbind(m[4:6, ], m[1:3, ])
  cells <- list(expected[4:6, ], expected[1:3, ])
  n <- matrix(10, 3, 2)
  expect_equal(unname(attr(res, "draws")$T2[, "observed"]), klauer_t2(wide(observed), cells, n))
  expect_equal(unname(attr(res, "draws")$T2[, "replicated"]), klauer_t2(wide(replicated), cells, n))
})

test_that("ppp() leaves participants without a matching set of cells out of T2", {
  data <- data.frame(id = factor(c(1, 2, 3, 1, 2)), item_type = c(rep("old", 3), "new", "new"))
  observed <- matrix(c(7, 5, 9, 2, 4, 3, 5, 1, 8, 6), 5)
  local_ppp_counts(observed, ppp_toy_draws(observed + 0.5), ppp_toy_draws(observed), data)
  fit <- ppp_toy_fit(data, bmf(D ~ 1 + (1 | id), g ~ 1))

  expect_message(res <- ppp(fit, statistic = "T2"), "leaves out 1 participant")
  expect_true(is.finite(res$observed))
})

test_that("ppp() skips impossible categories without Inf or NaN", {
  data <- data.frame(id = factor(1:3), item_type = "old")
  observed <- cbind(a = c(4, 5, 6), b = c(6, 5, 4), c = 0)
  expected <- cbind(c(5, 5, 5), c(5, 5, 5), 1e-43)
  replicated <- cbind(c(5, 4, 6), c(5, 6, 4), 0)
  local_ppp_counts(observed, ppp_toy_draws(expected), ppp_toy_draws(replicated), data)
  fit <- ppp_toy_fit(data, bmf(D ~ 1 + (1 | id), g ~ 1))

  res <- ppp(fit)
  expect_true(all(is.finite(res$observed)) && all(is.finite(res$replicated)))
  expect_equal(unname(attr(res, "draws")$T1[, "observed"]), 0)
  # the same counts without the impossible category give the same statistics
  local_ppp_counts(observed[, 1:2], ppp_toy_draws(expected[, 1:2]),
                   ppp_toy_draws(replicated[, 1:2]), data)
  expect_equal(unclass(ppp(fit))$observed, res$observed)
})

test_that("ppp() refuses observed counts in an impossible category", {
  data <- data.frame(item_type = "old")
  observed <- cbind(a = 4, b = 5, c = 1)
  expected <- cbind(5, 5, 1e-43)
  local_ppp_counts(observed, ppp_toy_draws(expected), ppp_toy_draws(expected), data)
  expect_error(ppp(ppp_toy_fit(data), statistic = "T1"), "probability 0 to observed responses in categories 'c'")
})

test_that("ppp() asks for the participant column when it cannot find one", {
  data <- data.frame(id = factor(1:2), item = factor(1:2), item_type = "old")
  observed <- matrix(c(4, 6, 6, 4), 2)
  local_ppp_counts(observed, ppp_toy_draws(observed), ppp_toy_draws(observed), data)

  no_re <- ppp_toy_fit(data)
  expect_message(res <- ppp(no_re), "computes T1 only")
  expect_equal(res$statistic, "T1")
  expect_error(ppp(no_re, statistic = "T2"), "Name the column with 'participant'")
  expect_error(ppp(no_re, level = "participant"), "level = 'participant' needs")
  expect_no_message(ppp(no_re, statistic = "T2", participant = "id"))

  crossed <- ppp_toy_fit(data, bmf(D ~ 1 + (1 | id) + (1 | item), g ~ 1))
  expect_error(ppp(crossed), "several grouping variables \\('id', 'item'\\)")
  expect_no_error(ppp(crossed, statistic = "T1"))
})

test_that("ppp() checks its arguments", {
  data <- data.frame(id = factor(1:2), item_type = "old")
  observed <- matrix(c(4, 6, 6, 4), 2)
  local_ppp_counts(observed, ppp_toy_draws(observed), ppp_toy_draws(observed), data)
  fit <- ppp_toy_fit(data, bmf(D ~ 1 + (1 | id), g ~ 1))

  expect_error(ppp(fit, level = "participant", statistic = "T2"), "group-level statistic")
  expect_error(ppp(fit, participant = "subject"), "'participant' must name a column")
  expect_error(ppp(fit, by = "cond"), "'by' must name columns")
  expect_error(ppp(fit, statistic = "T3"), "should be one of")

  not_counts <- fit
  not_counts$bmm$model <- mixture2p(resp_error = "y")
  expect_error(ppp(not_counts), "available for the count models")
})

test_that("ppp() returns one T1 per participant at the participant level", {
  data <- data.frame(id = factor(c(1, 2, 1, 2)), item_type = rep(c("old", "new"), each = 2))
  observed <- matrix(c(7, 5, 2, 4, 3, 5, 8, 6), 4)
  expected <- matrix(c(6, 6, 3, 3, 4, 4, 7, 7), 4)
  local_ppp_counts(observed, ppp_toy_draws(expected), ppp_toy_draws(expected), data)
  fit <- ppp_toy_fit(data, bmf(D ~ 1 + (1 | id), g ~ 1))

  res <- ppp(fit, level = "participant")
  expect_equal(as.character(res$participant), c("1", "2"))
  # participant 1: rows 1 and 3, each its own cell
  expect_equal(res$observed[1], 1 / 6 + 1 / 4 + 1 / 3 + 1 / 7)
  expect_equal(res$ppp, c(0, 0))
  expect_output(print(res), "'item_type' \\(2 cells\\)")
})

test_that(".ppp_complement() adds the second category of a binomial count", {
  draws <- matrix(c(3, 4, 7, 8), 2)
  out <- .ppp_complement(draws, trials = c(10, 20))
  expect_equal(dim(out), c(2L, 2L, 2L))
  expect_equal(out[, , 1], draws)
  expect_equal(out[, , 2], matrix(c(7, 6, 13, 12), 2))
})

test_that("ppp() works on an m3 fit", {
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "M3 fixture not available (excluded by .Rbuildignore)")
  withr::local_seed(1)
  res <- ppp(readRDS(path), ndraws = 20)
  expect_s3_class(res, "bmm_ppp")
  expect_equal(res$statistic, c("T1", "T2"))
  expect_true(all(res$ppp >= 0 & res$ppp <= 1))
  expect_equal(res$ndraws, c(20L, 20L))
})

test_that("ppp() does not flag a well-specified mpt model and flags a misspecified one", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  withr::local_seed(2)
  withr::local_options("bmm.silent" = 2)
  n_id <- 15
  dat <- expand.grid(id = factor(seq_len(n_id)), item_type = c("old", "new"),
                     stringsAsFactors = FALSE)
  D <- stats::plogis(stats::qlogis(0.6) + stats::rnorm(n_id, 0, 0.5))[dat$id]
  p_old <- ifelse(dat$item_type == "old", D + (1 - D) * 0.5, (1 - D) * 0.5)
  dat$old <- stats::rbinom(nrow(dat), 60, p_old)
  dat$new <- 60L - dat$old

  fit <- bmm(
    bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id)), dat, mpt(mpt_2htm_trees(), "item_type"),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, seed = 2
  )
  res <- ppp(fit)
  expect_equal(attr(res, "by"), "item_type")
  expect_true(all(res$ppp > 0.05))
  expect_equal(nrow(ppp(fit, level = "participant")), n_id)

  misfit <- bmm(
    bmf(D ~ 1, g = 0.9), dat, mpt(mpt_2htm_trees(), "item_type"),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, seed = 2
  )
  withr::local_options("bmm.silent" = 1)
  expect_message(res_misfit <- ppp(misfit), "computes T1 only")
  expect_lt(res_misfit$ppp, 0.01)
})

test_that("ppp() handles the impossible categories of a fitted mpt model", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  withr::local_seed(3)
  withr::local_options("bmm.silent" = 2)
  dat <- mpt_impossible_data()
  dat$corr <- dat$corr + stats::rbinom(nrow(dat), 5, 0.5)

  fit <- bmm(
    bmf(Pm ~ 1 + (1 | id), Pb ~ 1), dat, mpt(mpt_impossible_trees(), "tree"),
    backend = "cmdstanr", chains = 2, iter = 1000, refresh = 0, seed = 3
  )
  res <- ppp(fit)
  expect_equal(attr(res, "by"), "tree")
  expect_true(all(is.finite(c(res$observed, res$replicated, res$ppp))))
})
