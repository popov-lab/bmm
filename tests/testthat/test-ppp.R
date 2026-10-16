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
  # under cov()'s divisor T - 1 the expected sample covariance is
  # cov(expected) + within; Klauer's factor (T - 1) / T goes with the divisor T
  sigma <- stats::cov(do.call(cbind, cells)) + within
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

test_that("T2 is 0 when the observed covariance is the one the model implies", {
  # E[cov(counts)] by enumeration: 3 participants, 3 binomial trials each
  n_trials <- 3
  p <- c(0.2, 0.5, 0.7)
  outcomes <- as.matrix(expand.grid(0:n_trials, 0:n_trials, 0:n_trials))
  prob <- apply(outcomes, 1, function(k) prod(stats::dbinom(k, n_trials, p)))
  expected_cov <- Reduce(`+`, lapply(seq_len(nrow(outcomes)), function(o) {
    prob[o] * stats::cov(cbind(outcomes[o, ], n_trials - outcomes[o, ]))
  }))
  # observed counts whose sample covariance is exactly expected_cov, which has
  # rank 1 because the two counts of a participant sum to n_trials
  eig <- eigen(expected_cov, symmetric = TRUE)
  centred <- outer(c(1, -1, 0), sqrt(eig$values[1]) * eig$vectors[, 1])
  observed <- sweep(centred, 2, n_trials * c(0.5, 0.5), "+")
  colnames(observed) <- c("old", "new")

  data <- data.frame(id = factor(1:3), item_type = "old")
  expected <- cbind(n_trials * p, n_trials * (1 - p))
  local_ppp_counts(observed, ppp_toy_draws(expected), ppp_toy_draws(expected), data)
  res <- ppp(ppp_toy_fit(data, bmf(D ~ 1 + (1 | id), g ~ 1)), statistic = "T2")
  expect_equal(unname(attr(res, "draws")$T2[, "observed"]), 0)
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

test_that("ppp() leaves out categories whose expected count is exactly 0", {
  data <- data.frame(id = factor(1:3), item_type = "old")
  observed <- cbind(a = c(4, 5, 6), b = c(6, 5, 4), c = 0)
  expected <- cbind(c(5, 4, 6), c(5, 6, 4), 0)
  replicated <- cbind(c(6, 3, 6), c(4, 7, 4), 0)
  fit <- ppp_toy_fit(data, bmf(D ~ 1 + (1 | id), g ~ 1))
  local_ppp_counts(observed, ppp_toy_draws(expected), ppp_toy_draws(replicated), data)
  group <- attr(ppp(fit), "draws")
  person <- attr(ppp(fit, level = "participant"), "draws")

  local_ppp_counts(observed[, 1:2], ppp_toy_draws(expected[, 1:2]),
                   ppp_toy_draws(replicated[, 1:2]), data)
  expect_equal(group, attr(ppp(fit), "draws"))
  expect_equal(person, attr(ppp(fit, level = "participant"), "draws"))
})

test_that("T2 leaves out a cell in which a single category is possible", {
  # e.g. an m3 set size without "other" candidates: every response is correct
  data <- data.frame(id = factor(rep(1:3, 2)), item_type = rep(c("old", "new"), each = 3))
  observed <- cbind(old = c(10, 10, 10, 4, 6, 5), new = c(0, 0, 0, 6, 4, 5))
  expected <- cbind(old = c(10, 10, 10, 5, 4, 6), new = c(1e-43, 1e-43, 1e-43, 5, 6, 4))
  replicated <- cbind(old = c(10, 10, 10, 5, 3, 6), new = c(0, 0, 0, 5, 7, 4))
  fit <- ppp_toy_fit(data, bmf(D ~ 1 + (1 | id), g ~ 1))
  local_ppp_counts(observed, ppp_toy_draws(expected), ppp_toy_draws(replicated), data)
  with_forced_cell <- attr(ppp(fit, statistic = "T2"), "draws")$T2

  rows <- 4:6
  local_ppp_counts(observed[rows, ], ppp_toy_draws(expected[rows, ]),
                   ppp_toy_draws(replicated[rows, ]), data[rows, ])
  expect_equal(with_forced_cell, attr(ppp(fit, statistic = "T2"), "draws")$T2)
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

test_that("ppp() says where the participant column has to be", {
  # brms keeps only formula variables, so a model without person effects has
  # no participant column in its data
  data <- data.frame(item_type = c("old", "old"))
  observed <- matrix(c(4, 6, 6, 4), 2)
  local_ppp_counts(observed, ppp_toy_draws(observed), ppp_toy_draws(observed), data)
  no_re <- ppp_toy_fit(data)
  expect_match(tryCatch(ppp(no_re), message = conditionMessage),
               "contain only the variables of the formula")
  expect_match(tryCatch(ppp(no_re, participant = "id"), error = conditionMessage),
               "contain only the variables of the formula")
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

test_that("ppp() compares each replicated data set with the expectation of its own draw", {
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  withr::local_seed(1)
  res <- ppp(fit, statistic = "T1", by = character(0), ndraws = 5)
  expect_s3_class(res, "bmm_ppp")

  # the same random numbers as ppp(): the draw ids first, then the simulation
  withr::local_seed(1)
  draw_ids <- sample.int(brms::ndraws(fit), 5)
  yrep <- brms::posterior_predict(fit, draw_ids = draw_ids)
  mu <- brms::posterior_epred(fit, draw_ids = draw_ids)
  pooled <- function(x) apply(x, c(1L, 3L), sum)
  x2 <- function(n) rowSums((n - pooled(mu))^2 / pooled(mu))
  observed <- matrix(colSums(unclass(fit$data$Y)), 5, dim(mu)[3], byrow = TRUE)
  expect_equal(unname(attr(res, "draws")$T1[, "observed"]), x2(observed))
  expect_equal(unname(attr(res, "draws")$T1[, "replicated"]), x2(pooled(yrep)))
})

test_that("ppp() takes ndraws as a whole number of draws", {
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  expect_error(ppp(fit, ndraws = 2.5), "'ndraws' must be a whole number")
})

test_that("ppp() counts the second category of a yes/no model as trials minus responses", {
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  fit <- list(data = data.frame(n_old = c(30, 70), n_trials = c(100, 80)))
  draws <- matrix(c(28, 31, 69, 72), 2, 2)
  local_mocked_bindings(posterior_predict = function(...) draws,
                        posterior_epred = function(...) draws + 0.5, .package = "brms")
  counts <- .ppp_counts(model, fit, draw_ids = 1:2)
  expect_equal(counts$observed, cbind(n_old = c(30, 70), other = c(70, 10)))
  expect_equal(dim(counts$replicated), c(2, 2, 2))
  expect_equal(counts$replicated[, , 1], draws)
  expect_equal(counts$replicated[, , 2], cbind(100 - draws[, 1], 80 - draws[, 2]))
  expect_equal(counts$expected[, , 2], cbind(100 - draws[, 1], 80 - draws[, 2]) - 0.5)
})

test_that("ppp() rebuilds the tree column of an mpt fit that no formula kept", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  local_ppd <- function(fit) {
    draws <- array(1, c(2, nrow(fit$data), 2))
    local_mocked_bindings(posterior_predict = function(...) draws,
                          posterior_epred = function(...) draws, .package = "brms",
                          .env = parent.frame())
  }
  cases <- model_test_cases("stored_frame")
  fit <- stored_frame_fit(cases$mpt)
  expect_false("item_type" %in% names(fit$data))
  local_ppd(fit)
  data <- .ppp_counts(fit$bmm$model, fit, draw_ids = 1:2)$data
  expect_identical(levels(data$item_type), c("old", "new"))
  expect_identical(as.character(data$item_type), ifelse(fit$data$Idx_old == 1, "old", "new"))

  # a formula that names the tree column keeps it, and it is left as it is
  fit <- stored_frame_fit(cases$mpt_tree_predictor)
  local_ppd(fit)
  expect_identical(.ppp_counts(fit$bmm$model, fit, draw_ids = 1:2)$data$item_type,
                   fit$data$item_type)
})
