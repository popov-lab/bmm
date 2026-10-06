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

test_that("restrictions parse from MPTinR strings and named lists alike", {
  from_strings <- .mpt_parse_restrictions(c("Dn = Do", "g = 1/4", "G1 = G2 = G3"))
  expect_equal(
    from_strings,
    list(Dn = quote(Do), g = 0.25, G1 = quote(G3), G2 = quote(G3))
  )
  expect_equal(
    .mpt_parse_restrictions(list(Dn = "Do", g = 0.25)),
    list(Dn = quote(Do), g = 0.25)
  )
  expect_equal(.mpt_parse_restrictions(c(g = 1 / 4)), list(g = 0.25))
  expect_equal(.mpt_parse_restrictions(from_strings), from_strings)
  expect_equal(.mpt_parse_restrictions(NULL), list())

  expect_error(.mpt_parse_restrictions("Do < Dn"), "Order constraints")
  expect_error(.mpt_parse_restrictions("Do >= Dn"), "Order constraints")
  expect_error(.mpt_parse_restrictions("Do = Dn * 2"), "does neither")
  expect_error(.mpt_parse_restrictions("Do"), "form")
  expect_error(.mpt_parse_restrictions("Do = "), "parse")
  expect_error(.mpt_parse_restrictions(list("Do")), "the form")
})

test_that("an unnamed list of restrictions parses like the character vector", {
  expect_equal(
    .mpt_parse_restrictions(list("Dn = Do", "g = 0.5")),
    .mpt_parse_restrictions(c("Dn = Do", "g = 0.5"))
  )
  expect_equal(
    .mpt_parse_restrictions(list("Dn = Do", "g = 0.5")),
    list(Dn = quote(Do), g = 0.5)
  )
})

test_that("restriction constants must be strictly between 0 and 1", {
  expect_error(.mpt_parse_restrictions("g = 1.5"), "strictly between 0 and 1")
  expect_error(.mpt_parse_restrictions("g = -0.2"), "strictly between 0 and 1")
  expect_error(.mpt_parse_restrictions("g = 0"), "strictly between 0 and 1")
  expect_error(.mpt_parse_restrictions("g = 1"), "strictly between 0 and 1")
  expect_error(.mpt_parse_restrictions(list(g = 1)), "strictly between")
  expect_equal(.mpt_parse_restrictions("g = 0.5"), list(g = 0.5))
  expect_equal(.mpt_parse_restrictions("g = 1/4"), list(g = 0.25))
})

test_that("the message for a constant 0 or 1 holds for every model", {
  expect_error(.mpt_parse_restrictions("D = 0"), "reduced tree")
  expect_error(.mpt_parse_restrictions("D = 1"), "mpt_tree\\(impossible = \\)")
})

test_that("an out-of-range constant is not explained with 0/1 branches", {
  for (r in c("g = 1.5", "g = -0.2")) {
    msg <- tryCatch(.mpt_parse_restrictions(r), error = conditionMessage)
    expect_match(msg, "strictly between 0 and 1")
    expect_no_match(msg, "branch")
  }
})

test_that("a restrictions file path gets a pointer to readLines()", {
  path <- withr::local_tempfile(fileext = ".restr", lines = "Dn = Do")
  expect_error(.mpt_parse_restrictions(path), "does not read restriction")
  expect_error(.mpt_parse_restrictions("models/2htm.restr"), "readLines")
  expect_error(.mpt_parse_restrictions("models/2htm.txt"), "readLines")
  expect_error(.mpt_parse_restrictions(list("models/2htm.restr")), "readLines")
  expect_equal(.mpt_parse_restrictions("Dn = Do"), list(Dn = quote(Do)))
})

test_that("an existing file without a known extension is recognised", {
  path <- withr::local_tempfile(lines = "Dn = Do")
  expect_error(.mpt_parse_restrictions(path), "does not read restriction")
})

test_that("the path hint quotes the path as valid R and names it once", {
  msg <- tryCatch(
    .mpt_parse_restrictions("C:\\models\\2htm.restr"),
    error = conditionMessage
  )
  expect_match(msg, 'readLines("C:\\\\models\\\\2htm.restr")', fixed = TRUE)
  expect_equal(lengths(regmatches(msg, gregexpr("2htm.restr", msg, fixed = TRUE))), 1L)
})

test_that("restrictions with a file name in a comment or an existing file are valid", {
  expect_equal(.mpt_parse_restrictions("g = 0.5 # half.txt"), list(g = 0.5))
  expect_equal(
    .mpt_parse_restrictions("Dn = Do # from broeder.2htm.restr"),
    list(Dn = quote(Do))
  )
  expect_equal(
    .mpt_parse_restrictions("Dn = Do # see notes.TXT "),
    list(Dn = quote(Do))
  )
  withr::with_dir(withr::local_tempdir(), {
    file.create("g = 0.5")
    expect_equal(.mpt_parse_restrictions("g = 0.5"), list(g = 0.5))
  })
  expect_error(
    .mpt_parse_restrictions("D1 < D2 # notes.txt"), "Order constraints"
  )
})

test_that("blank and comment lines of readLines() are skipped", {
  path <- withr::local_tempfile(
    lines = c("# 2HTM", "Dn = Do # equal detection", "", "  ", "g = 0.5", "  # end")
  )
  expect_equal(
    .mpt_parse_restrictions(readLines(path)),
    list(Dn = quote(Do), g = 0.5)
  )
  expect_equal(.mpt_parse_restrictions(c("", "# only a comment")), list())
})

test_that("chained order constraints get the order-constraint message", {
  expect_error(.mpt_parse_restrictions("G1 < G2 < G3"), "Order constraints")
  expect_error(.mpt_parse_restrictions("G1 > G2 > G3"), "Order constraints")
  expect_error(.mpt_parse_restrictions("D1 <= D2 <= D3"), "Order constraints")
  expect_error(.mpt_parse_restrictions("D1 < D2 < D3"), "'D1 < D2 < D3'")
  expect_error(.mpt_parse_restrictions("Do = "), "Cannot parse")
})

test_that("restriction chains resolve to their final target", {
  resolved <- .mpt_resolve_restrictions(list(A = quote(B), B = quote(C), g = 0.5))
  expect_equal(resolved, list(A = quote(C), B = quote(C), g = 0.5))
})

test_that("mpt errors when branch probabilities do not sum to 1", {
  bad_tree <- mpt_tree("t", list(a = "D * g", b = "(1 - D) * g"))
  expect_error(mpt(bad_tree), "sum to")

  good_tree <- mpt_tree("u", list(a = "D + (1 - D) * g", b = "(1 - D) * (1 - g)"))
  deviations <- .mpt_tree_sum_deviations(
    list(t = bad_tree, u = good_tree), c("D", "g"), list()
  )
  # the branches reduce to g, so the first test point reports g's value there
  expect_equal(
    deviations[["t"]], .mpt_test_points(c("D", "g"), list())[[1]][["g"]]
  )
  expect_true(is.na(deviations[["u"]]))
})

test_that("test points give every symbol its own interior value at every point", {
  for (n in c(2, 4, 5, 8, 12, 30)) {
    symbols <- paste0("p", seq_len(n))
    points <- .mpt_test_points(symbols, list())
    expect_length(points, 5)
    for (vals in points) {
      expect_named(vals, symbols)
      expect_true(all(vals > 0 & vals < 1))
      expect_equal(anyDuplicated(round(vals, 6)), 0L)
    }
    expect_equal(anyDuplicated(round(unlist(lapply(points, `[[`, 1)), 6)), 0L)
  }

  points <- .mpt_test_points(c("a", "b", "c", "d"), list(c("a", "b", "c")))
  for (vals in points) {
    expect_equal(sum(vals[c("a", "b", "c")]), 1)
  }
})

test_that("no sum of two test values equals another such sum at every point", {
  # a linear sequence in the symbol index gives v_i + v_j == v_k + v_l at
  # every point whenever i + j == k + l, which the rank check would read as a
  # property of the model
  n <- 12
  values <- do.call(rbind, .mpt_test_points(paste0("p", seq_len(n)), list()))
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
  at_first <- sprintf("%.8f", .mpt_test_points("D", list())[[1]][["D"]])
  tree <- mpt_tree("t", list(a = "D", b = glue("1 - D + (D - {at_first})")))
  expect_true(abs(
    sum(.mpt_eval_branches(tree, as.list(.mpt_test_points("D", list())[[1]]))) - 1
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
