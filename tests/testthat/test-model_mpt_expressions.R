test_that("mpt_tree folds integer fractions into decimal literals", {
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1/4)",
    b = "(1 - p) * (3/4)"
  ))
  expect_equal(deparse1(tree$branches$a), "p + (1 - p) * 0.25")
  expect_equal(deparse1(tree$branches$b), "(1 - p) * 0.75")
})

test_that("mpt_tree folds compound integer constants into decimal literals", {
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1/(2*2))",
    b = "(1 - p) * (-1/4)",
    c = "p * (1/2^2)"
  ))
  expect_equal(deparse1(tree$branches$a), "p + (1 - p) * 0.25")
  expect_equal(deparse1(tree$branches$b), "(1 - p) * -0.25")
  expect_equal(deparse1(tree$branches$c), "p * 0.25")
  dat <- data.frame(a = c(30L, 28L), b = c(0L, 2L))
  code <- stancode(
    bmf(p ~ 1), data = dat,
    model = mpt(mpt_tree("t", list(a = "p + (1 - p) * (1/(2*2))", b = "(1 - p) * (1 - 1/(2*2))")))
  )
  expect_match(code, "* 0.25", fixed = TRUE)
  expect_false(grepl("1 / (2 * 2)", code, fixed = TRUE))
})

test_that("mpt_tree stores branch expressions as parsed calls", {
  tree <- mpt_tree("t", list(a = "D + (1 - D) * g", b = "(1 - D) * (1 - g)"))
  expect_identical(tree$branches$a, quote(D + (1 - D) * g))
  expect_equal(.mpt_expr_vars(tree), c("D", "g"))
  expect_equal(
    .mpt_eval_branches(tree, list(D = 0.7, g = 0.5)),
    c(a = 0.85, b = 0.15)
  )
})

test_that("mpt errors when branch probabilities do not sum to 1", {
  bad_tree <- mpt_tree("t", list(a = "D * g", b = "(1 - D) * g"))
  expect_error(mpt(bad_tree), "sum to")
  expect_error(mpt(bad_tree), "or tie them in the formula \\(e.g. Dn ~ Do\\)")

  good_tree <- mpt_tree("u", list(a = "D + (1 - D) * g", b = "(1 - D) * (1 - g)"))
  branch_errors <- .mpt_tree_branch_errors(
    list(t = bad_tree, u = good_tree), c("D", "g")
  )
  # the branches reduce to g, so the first test point reports g's value there
  expect_match(
    branch_errors[["t"]],
    glue("sum to {signif(.mpt_test_points(c('D', 'g'))[[1]][['g']], 6)} instead of 1")
  )
  expect_true(is.na(branch_errors[["u"]]))
})

test_that("mpt errors when a branch probability leaves (0, 1]", {
  # the branches sum to 1 for every value of a
  doubled <- mpt_tree("t", list(yes = "2 * a", no = "1 - 2 * a"))
  expect_error(
    mpt(doubled),
    "category 'yes' in tree 't' is 1\\.7.* at the test values a = 0\\.85.*, outside \\(0, 1\\]"
  )
  zero <- mpt_tree("z", list(yes = "a", no = "(1 - a) * 0", maybe = "1 - a"))
  expect_error(
    mpt(zero),
    "category 'no' in tree 'z' is 0 .*outside \\(0, 1\\].*makes the likelihood undefined"
  )
  expect_no_error(mpt(mpt_tree("t", list(yes = "a", no = "1 - a"))))
})

test_that("test points give every symbol its own interior value at every point", {
  for (n in c(2, 4, 5, 8, 12, 30)) {
    symbols <- paste0("p", seq_len(n))
    points <- .mpt_test_points(symbols)
    expect_length(points, 5)
    for (vals in points) {
      expect_named(vals, symbols)
      expect_true(all(vals > 0 & vals < 1))
      expect_equal(anyDuplicated(round(vals, 6)), 0L)
    }
    expect_equal(anyDuplicated(round(unlist(lapply(points, `[[`, 1)), 6)), 0L)
  }
})

test_that("no sum of two test values equals another such sum at every point", {
  # a linear sequence in the symbol index gives v_i + v_j == v_k + v_l at
  # every point whenever i + j == k + l, which the rank check would read as a
  # property of the model
  n <- 12
  values <- do.call(rbind, .mpt_test_points(paste0("p", seq_len(n))))
  pairs <- t(utils::combn(n, 2))
  sums <- apply(pairs, 1, function(pair) values[, pair[1]] + values[, pair[2]])
  coincide <- outer(seq_len(nrow(pairs)), seq_len(nrow(pairs)), Vectorize(
    function(a, b) a < b && all(abs(sums[, a] - sums[, b]) < 1e-9)
  ))
  expect_false(any(coincide))
})

test_that("a typo between parameters four places apart in the symbol order is caught", {
  # a short cycle of test values would give the first and the fifth parameter
  # the same value at every point, so (1 - A) * (1 - A) would pass for
  # (1 - A) * (1 - F)
  branches <- list(
    r1 = "A * B", r2 = "A * (1 - B) * C", r3 = "A * (1 - B) * (1 - C) * E",
    r4 = "A * (1 - B) * (1 - C) * (1 - E)", r5 = "(1 - A) * F"
  )
  good <- mpt_tree("main", c(branches, r6 = "(1 - A) * (1 - F)"))
  typo <- mpt_tree("main", c(branches, r6 = "(1 - A) * (1 - A)"))
  expect_no_error(mpt(good))
  expect_error(mpt(typo), "sum to")
})

test_that("a tree that sums to 1 at only one test point is still caught", {
  at_first <- sprintf("%.8f", .mpt_test_points("D")[[1]][["D"]])
  tree <- mpt_tree("t", list(a = "D", b = glue("1 - D + (D - {at_first})")))
  expect_true(abs(
    sum(.mpt_eval_branches(tree, as.list(.mpt_test_points("D")[[1]]))) - 1
  ) < 1e-6)
  expect_error(mpt(tree), "sum to")
})

test_that("constants that Stan would receive in scientific notation error", {
  # brms deparses formula constants into the Stan code and spaces out
  # operators, so 6.7e-05 would become the subtraction '6.7e - 05'
  expect_error(
    mpt_tree("t", list(
      a = "p + (1 - p) * (1 - 0.001/15.001)",
      b = "(1 - p) * (0.001/15.001)"
    )),
    "scientific"
  )
  # constants that deparse in fixed notation are fine
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1 - 0.001)",
    b = "(1 - p) * 0.001"
  ))
  expect_equal(deparse1(tree$branches$b), "(1 - p) * 0.001")
})
