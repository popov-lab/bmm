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

test_that("mpt_tree says why a constant is refused and which rewrite works", {
  expect_error(
    mpt_tree("t", list(a = "0.00010 + 0.9999 * p", b = "0.9999 * (1 - p)")),
    "'0.00010' \\(1e-04\\) in tree 't' in scientific notation.*\\(1 \\+ 9999 \\* p\\) / 10000"
  )
  # a folded constant has no written form, so brms's form is not repeated
  expect_error(
    mpt_tree("t", list(a = "1/10000 + 0.9999 * p", b = "0.9999 * (1 - p)")),
    "constant\\(s\\) '1e-04' in tree 't' in scientific notation, which"
  )
  # a large constant gets a rewrite for large numbers
  for (large in c("100000", "1e5")) {
    expect_error(
      mpt_tree("t", list(a = glue("p / {large}"), b = glue("1 - p / {large}"))),
      glue(
        "'{large}' \\(1e\\+05\\) in tree 't'.*write p / 100000 as p / 1000 / \\
        100, or provide the value as a data column declared in the covariates \\
        argument\\.$"
      )
    )
  }
  expect_no_error(mpt_tree("t", list(a = "p / 1000 / 100", b = "1 - p / 1000 / 100")))
  rescaled <- mpt(mpt_tree("t", list(
    a = "(1 + 9999 * p) / 10000", b = "9999 * (1 - p) / 10000"
  )))
  code <- stancode(
    bmf(p ~ 1), data = data.frame(a = c(30L, 28L), b = c(0L, 2L)), model = rescaled
  )
  expect_match(code, "(1 + 9999 * inv_logit(nlp_p[n])) / 10000", fixed = TRUE)
  expect_no_match(code, "e - 0")
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
  expect_error(.mpt_parse_restrictions("g = 0"), "mpt_tree(impossible = )", fixed = TRUE)
  expect_error(.mpt_parse_restrictions("g = 1"), "strictly between 0 and 1")
  expect_error(.mpt_parse_restrictions(list(g = 1)), "strictly between")
  expect_equal(.mpt_parse_restrictions("g = 0.5"), list(g = 0.5))
  expect_equal(.mpt_parse_restrictions("g = 1/4"), list(g = 0.25))
})

test_that("the message for a constant 0 or 1 holds for every model", {
  for (r in c("D = 0", "D = 1")) {
    msg <- tryCatch(.mpt_parse_restrictions(r), error = conditionMessage)
    expect_match(msg, "reduced tree")
    expect_no_match(msg, "turns a branch into the constant 0")
  }
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
  expect_error(mpt(bad_tree), "or tie them in the formula \\(e.g. Dn ~ Do\\)")
  expect_error(mpt(bad_tree), "use the restrictions argument")

  good_tree <- mpt_tree("u", list(a = "D + (1 - D) * g", b = "(1 - D) * (1 - g)"))
  branch_errors <- .mpt_tree_branch_errors(
    list(t = bad_tree, u = good_tree), c("D", "g"), list()
  )
  # the branches reduce to g, so the first test point reports g's value there
  expect_match(
    branch_errors[["t"]],
    glue("sum to {signif(.mpt_test_points(c('D', 'g'), list())[[1]][['g']], 6)} instead of 1")
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

test_that("the stick-breaking Jacobian matches central differences of the map", {
  sticks_to_members <- function(sticks) cumprod(c(1, 1 - sticks)) * c(sticks, 1)
  for (sticks in list(0.3, c(0.3, 0.6), c(0.3, 0.6, 0.2), c(0.7, 0.2, 0.5, 0.4))) {
    members <- sticks_to_members(sticks)
    numeric_map <- vapply(seq_along(sticks), function(m) {
      step <- replace(numeric(length(sticks)), m, 1e-6)
      (sticks_to_members(sticks + step) - sticks_to_members(sticks - step)) / 2e-6
    }, numeric(length(members)))
    names <- paste0("s", seq_along(sticks))
    exact <- .mpt_stick_jacobian(diag(length(members)), members, names)
    expect_equal(
      unname(exact), matrix(numeric_map, nrow = length(members)), tolerance = 1e-8
    )
    expect_equal(colnames(exact), names)
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

test_that(".mpt_branch_vertices visits the 2^k vertices of a multilinear branch", {
  vertices <- bmm:::.mpt_branch_vertices(quote(a * (1 - b) * c * G), c("a", "b", "c"), list())
  expect_true(attr(vertices, "exact"))
  expect_equal(dim(vertices), c(8L, 3L))
  expect_equal(colnames(vertices), c("a", "b", "c"))
  expect_setequal(unique(as.vector(vertices)), c(1e-9, 1 - 1e-9))
  expect_equal(nrow(unique(vertices)), 8L)

  # a parameter in several terms is still of degree 1
  twice <- bmm:::.mpt_branch_vertices(quote(a * (1 - b) + (1 - a) * b), c("a", "b"), list())
  expect_equal(dim(twice), c(4L, 2L))

  constant <- bmm:::.mpt_branch_vertices(quote(0.5 * G), c("a", "b"), list())
  expect_equal(dim(constant), c(1L, 0L))
})

test_that(".mpt_branch_vertices adds 0.5 for a parameter of higher degree", {
  vertices <- bmm:::.mpt_branch_vertices(quote(4 * u * (1 - u) * a), c("u", "a"), list())
  expect_equal(sort(unique(vertices[, "u"])), c(1e-9, 0.5, 1 - 1e-9))
  expect_equal(sort(unique(vertices[, "a"])), c(1e-9, 1 - 1e-9))
  expect_equal(nrow(vertices), 6L)

  # a derivative that stats::D() cannot form counts as higher degree
  unknown <- bmm:::.mpt_branch_vertices(quote(besselJ(a, 1)), "a", list())
  expect_equal(sort(unique(unknown[, "a"])), c(1e-9, 0.5, 1 - 1e-9))
  expect_false(attr(unknown, "differentiable"))
  expect_true(attr(vertices, "differentiable"))
})

test_that(".mpt_form bounds an expression by its form alone (#560)", {
  bound <- function(expr) bmm:::.mpt_form(expr, c("a", "b"), c("a", "b"))$bound
  # 2: in [0, 1]; 1: non-negative; 0: neither by form
  expect_equal(bound(quote(a * (1 - b)^2)), 2L)
  expect_equal(bound(quote(1 - a * (1 - b))), 2L)
  expect_equal(bound(quote(a * (1/5))), 2L)
  expect_equal(bound(quote(2 * a * (1 - a) + b)), 1L)
  expect_equal(bound(quote(1 - 2 * a)), 0L)
  expect_equal(bound(quote(a - b)), 0L)
  expect_equal(bound(quote(a / 2)), 0L)
  expect_equal(bound(quote(a^0.5)), 0L)
  expect_equal(bound(quote(exp(a))), 0L)
  expect_equal(bound(quote(G * a)), 0L)
  expect_equal(bmm:::.mpt_form(quote(G * a), c("a", "G"), "a")$bound, 2L)
  expect_equal(bmm:::.mpt_form(quote((1 - a)^n), "a", "a", c(n = 3))$bound, 2L)

  degree <- function(expr) bmm:::.mpt_form(expr, c("a", "b"), c("a", "b"))$degree
  expect_equal(degree(quote(2 * a * (1 - a)^2 + b)), c(3, 1))
  expect_equal(bmm:::.mpt_form(quote((1 - a)^n), "a", "a", c(n = 3))$degree, 3)
})

test_that(".mpt_certified_tree needs the branches to sum to 1 everywhere, not at the test points (#560)", {
  tree <- function(...) list(name = "t", branches = lapply(list(...), str2lang))
  expect_true(bmm:::.mpt_certified_tree(tree("a * b", "a * (1 - b)", "1 - a"), c("a", "b"), list()))
  expect_false(bmm:::.mpt_certified_tree(tree("a * b", "1 - a"), c("a", "b"), list()))

  # a^2 + 1 - a is 1 at a = 0 and a = 1 only; the grid of degree + 1 values
  # per parameter lies between them
  expect_false(bmm:::.mpt_certified_tree(tree("a^2", "1 - a"), "a", list()))

  # a simplex group sums to 1 only on the simplex
  group <- list(c("gA", "gB", "gC"))
  expect_true(bmm:::.mpt_certified_tree(tree("gA", "gB", "gC"), c("gA", "gB", "gC"), group))
  expect_false(bmm:::.mpt_certified_tree(tree("gA", "gB", "gC"), c("gA", "gB", "gC"), list()))
  # (gA + gB + gC)^2, expanded: 1 on the simplex, where the last member is
  # replaced by 1 minus the others and its degree adds to theirs
  squared <- tree("gA^2 + gB^2 + gC^2", "2 * gA * gB + 2 * gA * gC", "2 * gB * gC")
  expect_true(bmm:::.mpt_certified_tree(squared, c("gA", "gB", "gC"), group))
  expect_false(bmm:::.mpt_certified_tree(squared, c("gA", "gB", "gC"), list()))
  # with gB = 1 - gA these sum to 1 + 0.9 (gA - 1/6) (gA - 1/3): degree 2 in
  # gA, not the degree 1 that gA has in the branches, and exactly 1 at the two
  # grid values that degree 1 would take
  quadratic <- tree("0.5", "0.45 * gB", "0.9 * gA * gB")
  expect_false(bmm:::.mpt_certified_tree(quadratic, c("gA", "gB"), list(c("gA", "gB"))))
  expect_true(bmm:::.mpt_certified_tree(
    tree("(gA + gB + gC)^2 * D", "1 - D"), c("D", "gA", "gB", "gC"), group
  ))

  # with covariates, the sum must be 1 at each covariate row
  expect_true(bmm:::.mpt_certified_tree(
    tree("D + (1 - D) * G", "(1 - D) * H"), "D", list(), c("D", "G", "H"),
    list(G = c(0.2, 0.7), H = c(0.8, 0.3))
  ))
  expect_false(bmm:::.mpt_certified_tree(
    tree("D + (1 - D) * G", "(1 - D) * H"), "D", list(), c("D", "G", "H"),
    list(G = c(0.2, 0.7), H = c(0.8, 0.4))
  ))
})

test_that("the range search visits only points in the box and on the simplex (#560)", {
  # a branch undefined outside [0, 1] stays defined during the search
  vertices <- bmm:::.mpt_branch_vertices(quote(sqrt(u * (1 - u))), "u", list())
  expect_null(bmm:::.mpt_search_range_violation(
    quote(sqrt(u * (1 - u))), vertices, "u", list(), list(), 1e-6
  )$violation)
  # the reported witness is the extreme found, not the first value past 1
  found <- bmm:::.mpt_search_range_violation(
    quote(7 * u * (1 - u)^2),
    bmm:::.mpt_branch_vertices(quote(7 * u * (1 - u)^2), "u", list()),
    "u", list(), list(), 1e-6
  )
  expect_equal(found$violation$value, 28 / 27, tolerance = 1e-6)
})

test_that("the range search finds a narrow dip below a constant baseline (#560)", {
  # the minimum -0.00999 at u = 0.05 lies between the grid values 1e-9 and 0.5,
  # where the branch is 0.06 to within 4e-6
  dip <- "0.06 - 3.7095 * u * (1 - u)^19"
  expect_error(
    mpt(mpt_tree("t", list(y = dip, n = glue("1 - ({dip})")))),
    "category 'y' in tree 't' is -0\\.00998966 at the test values u = 0\\.05, outside"
  )
})

test_that("the range search reaches the later starts of every cell (#560)", {
  # 32 vertices of a1, ..., a5; only at all ai = 1 does the peak 1.0937 near
  # u = 0.049 appear, and the search finds it only from that cell's third grid
  # point u = 1 - 1e-9, which ranks below the third points of the 31 other
  # cells, all at 0.9
  product <- "a1 * a2 * a3 * a4 * a5"
  y <- glue(
    "{product} * (8 * u^2 * (1 - u)^2 + 60 * u * (1 - u)^20) + (1 - {product}) * 0.9"
  )
  expect_error(
    mpt(mpt_tree("t", list(y = y, n = glue("1 - ({y})")))),
    paste0(
      "category 'y' in tree 't' is 1\\.09373 at the test values ",
      "a1 = 1, a2 = 1, a3 = 1, a4 = 1, a5 = 1, u = 0\\.049"
    )
  )
})

test_that("the range search says which cells its budget left out (#560)", {
  branch <- quote((1 - G) * H + G * 0.9 * 27 / 4 * u * (1 - u)^2)
  vertices <- bmm:::.mpt_branch_vertices(branch, "u", list())
  rows <- list(G = c(rep(0, 9), 1), H = c(seq(0.80, 0.85, length.out = 9), 0.5))
  with_budget <- function(max_starts) {
    bmm:::.mpt_search_range_violation(
      branch, vertices, "u", list(), rows, 1e-6, max_starts = max_starts
    )
  }
  # each direction starts in 4 rows, the two together in 8
  four <- with_budget(4L)
  expect_null(four$violation)
  expect_equal(
    as.character(four$skipped),
    "left unsearched 2 of the 10 distinct covariate settings and capped its further starts"
  )
  expect_equal(
    as.character(with_budget(10L)$skipped),
    "searched all 10 distinct covariate settings and capped its further starts"
  )
  expect_null(with_budget(20L)$skipped)

  msgs <- testthat::capture_messages(bmm:::.mpt_message_sampled(c(
    "form|left unsearched 2 of the 10 distinct covariate settings" = "category 'y' in tree 't'",
    "form|left unsearched 2 of the 10 distinct covariate settings" = "category 'n' in tree 't'",
    "form|capped its further starts" = "category 'y' in tree 'u'",
    "grid" = "category 'a' in tree 's'",
    "grid" = "category 'b' in tree 'v'"
  ), "search"))
  expect_length(msgs, 2L)
  expect_match(
    msgs[1],
    paste0(
      "category 'y' in tree 't', category 'n' in tree 't', category 'y' in tree 'u' ",
      "is approximate: .*their trees are not written.*local search.*",
      "The search left unsearched 2 of the 10 distinct covariate settings for ",
      "category 'y' in tree 't', category 'n' in tree 't'; it capped its further ",
      "starts for category 'y' in tree 'u'\\."
    )
  )
  expect_match(msgs[2], "is approximate: their trees have the form")
  expect_no_match(msgs[2], "The search (left|searched|capped)")
})

test_that("the range search counts a cell searched in either direction (#560)", {
  # each direction has 256 first starts, towards the largest G and the
  # smallest; together they reach all 300 settings
  branch <- quote(0.4 * u * (1 - u) + 0.2 * G)
  search <- bmm:::.mpt_search_range_violation(
    branch, bmm:::.mpt_branch_vertices(branch, "u", list()), "u", list(),
    list(G = seq(0, 1, length.out = 300)), 1e-6
  )
  expect_null(search$violation)
  expect_equal(
    as.character(search$skipped),
    "searched all 300 distinct covariate settings and capped its further starts"
  )
})

test_that(".mpt_search_starts puts each row's best first and caps the later ranks apart (#560)", {
  # six cells of one grid point, three per covariate row, for 3 starts: the
  # second row's best comes before the first row's second best
  expect_equal(
    bmm:::.mpt_search_starts(
      rbind(c(9, 8, 7), c(1, 2, 3)), list(1L, 2L, 3L), max_starts = 3L, max_further = 0L
    ),
    list(row = c(1L, 2L, 1L), vertex = c(1L, 3L, 2L))
  )
  # four cells of three grid points: the first and second points take 8 of the
  # 10 starts, the third points only max_further = 1 of the 2 left
  expect_equal(
    bmm:::.mpt_search_starts(
      rbind(c(6, 5, 4, 3, 2, 1), c(1, 2, 3, 4, 5, 6)), list(1:3, 4:6),
      max_starts = 10L, max_further = 1L
    ),
    list(
      row = c(1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L),
      vertex = c(1L, 6L, 3L, 4L, 2L, 5L, 2L, 5L, 3L)
    )
  )
})

test_that(".mpt_branch_vertices takes a simplex group as one dimension", {
  group <- list(c("gA", "gB", "gC"))
  vertices <- bmm:::.mpt_branch_vertices(quote(D * gA + (1 - D) * gB), c("D", "gA", "gB", "gC"), group)
  expect_equal(colnames(vertices), c("D", "gA", "gB"))
  # gC takes the mass in the third vertex, where gA and gB are both small
  expect_equal(nrow(vertices), 2L * 3L)
  expect_true(all(rowSums(vertices[, c("gA", "gB")]) < 1.0001))

  all_members <- bmm:::.mpt_branch_vertices(quote(gA + 2 * gB + 3 * gC), c("gA", "gB", "gC"), group)
  expect_equal(nrow(all_members), 3L)
  expect_equal(unname(rowSums(all_members)), rep(1, 3))

  # a product of members peaks inside the simplex
  crossed <- bmm:::.mpt_branch_vertices(quote(gA * gB), c("gA", "gB", "gC"), group)
  expect_equal(nrow(crossed), 3L + 3L + 1L)
  expect_equal(unname(colMeans(crossed[7, , drop = FALSE])), rep(1 / 3, 2))
})

test_that(".mpt_branch_vertices falls back to a fixed sample above max_points", {
  parameters <- paste0("p", 1:11)
  branch <- str2lang(paste(parameters, collapse = " * "))
  withr::local_seed(1)
  seed_before <- .Random.seed
  vertices <- bmm:::.mpt_branch_vertices(branch, parameters, list())
  expect_identical(.Random.seed, seed_before)

  expect_false(attr(vertices, "exact"))
  expect_equal(dim(vertices), c(2L + 64L, 11L))
  expect_setequal(unique(as.vector(vertices)), c(1e-9, 1 - 1e-9))
  expect_equal(unname(vertices[1:2, 1]), c(1e-9, 1 - 1e-9))
  expect_identical(vertices, bmm:::.mpt_branch_vertices(branch, parameters, list()))
  # the sample uses both ends of every parameter
  expect_true(all(apply(vertices, 2, function(column) length(unique(column)) == 2L)))

  ten <- bmm:::.mpt_branch_vertices(str2lang(paste(parameters[1:10], collapse = " * ")), parameters, list())
  expect_true(attr(ten, "exact"))
  expect_equal(nrow(ten), 1024L)
})

test_that(".mpt_certified_tree checks the whole simplex when a group member is unused (#560)", {
  tree <- function(...) list(name = "t", branches = lapply(list(...), str2lang))
  group <- list(c("a", "b", "c"))
  # a + c is 1 - b, which is 1 only on the face b = 0
  unused <- bmm:::.mpt_certified_tree(tree("a", "c"), c("a", "b", "c"), group)
  expect_false(unused)
  expect_equal(attr(unused, "reason"), "identity")
  expect_true(bmm:::.mpt_certified_tree(tree("a", "b", "c"), c("a", "b", "c"), group))
  expect_false(bmm:::.mpt_certified_tree(tree("a^2", "c * (1 + a)"), c("a", "b", "c"), group))
  # the last member alone is computed from the members the tree does not use
  expect_true(bmm:::.mpt_certified_tree(
    tree("c^2", "c * (1 - c)", "1 - c"), c("a", "b", "c"), group
  ))
})

test_that(".mpt_certified_tree checks the sum at the ends of the parameter range (#560)", {
  tree <- function(...) list(name = "t", branches = lapply(list(...), str2lang))
  # off by 1e-8 only near D = 1, beyond the last of the equally spaced
  # interior grid values the check used before
  expect_false(bmm:::.mpt_certified_tree(
    tree("1 - (1 - D)^200", "(1 - D)^200 + 1e-8 * D^200"), "D", list()
  ))
  expect_true(bmm:::.mpt_certified_tree(tree("1 - (1 - D)^200", "(1 - D)^200"), "D", list()))
  expect_true(bmm:::.mpt_certified_tree(
    tree("1 - (1 - D)^n", "(1 - D)^n"), "D", list(), "D", list(n = c(2, 200)),
    c(n = 200)
  ))
})

test_that(".mpt_certified_tree says why a tree is not certified (#560)", {
  tree <- function(...) list(name = "t", branches = lapply(list(...), str2lang))
  outside <- bmm:::.mpt_certified_tree(tree("a / 2", "1 - a / 2"), "a", list())
  expect_equal(attr(outside, "reason"), "form")
  too_many <- bmm:::.mpt_certified_tree(
    tree("1 - (1 - D)^n", "(1 - D)^n"), "D", list(), "D", list(n = 1:300), c(n = 300)
  )
  expect_false(too_many)
  expect_equal(attr(too_many, "reason"), "grid_covariates")
  without_covariates <- bmm:::.mpt_certified_tree(
    tree("1 - (1 - D)^70000", "(1 - D)^70000"), "D", list()
  )
  expect_equal(attr(without_covariates, "reason"), "grid")
})

test_that("mpt() builds a model whose tree uses only the last member of a simplex group (#560)", {
  expect_silent(mpt(list(
    mpt_tree("t1", list(x = "gA", y = "gB")),
    mpt_tree("t2", list(x = "gB^2", y = "1 - gB^2"))
  ), tree_id = "tree", simplex = c("gA", "gB")))
  expect_silent(mpt(list(
    mpt_tree("t1", list(x = "gA", y = "gB", z = "gC")),
    mpt_tree("t2", list(x = "gC^2", y = "gC * (1 - gC)", z = "1 - gC"))
  ), tree_id = "tree", simplex = c("gA", "gB", "gC")))
})

test_that("the range search starts in every covariate row (#560)", {
  # the rows with G = 0 are flat at H; the row with G = 1 peaks at 1.0768 near
  # u = 0.048, while its grid values are all below 3e-5, below every flat row's
  model <- mpt(mpt_tree("t", list(
    y = "(1 - G) * H + G * 60 * u * (1 - u)^20",
    n = "1 - (1 - G) * H - G * 60 * u * (1 - u)^20"
  )), covariates = c("G", "H"))
  data <- data.frame(
    G = c(rep(0, 6), 1), H = c(0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.5), y = 5, n = 5
  )
  expect_error(
    check_data(model, data, bmf(u ~ 1)),
    paste0(
      "category 'y' in tree 't' is 1\\.07683 in row 7, outside \\(0, 1\\] ",
      "at the test parameter values u = 0\\.0476"
    )
  )
  # beyond the 64 starts of an earlier budget
  many_flat <- data.frame(
    G = c(rep(0, 100), 1), H = c(seq(0.80, 0.85, length.out = 100), 0.5), y = 5, n = 5
  )
  expect_error(
    check_data(model, many_flat, bmf(u ~ 1)),
    "category 'y' in tree 't' is 1\\.07683 in row 101"
  )
})

test_that("the range search starts at every vertex of the multilinear parameters (#560)", {
  # the peak 1.0768 sits at the one vertex P = 1, whose grid values rank
  # below those of the 31 vertices at 0.9
  product <- "a1 * a2 * a3 * a4 * a5"
  y <- glue("{product} * 60 * u * (1 - u)^20 + (1 - {product}) * 0.9")
  expect_error(
    mpt(mpt_tree("t", list(y = y, n = glue("1 - ({y})")))),
    paste0(
      "category 'y' in tree 't' is 1\\.07683 at the test values ",
      "a1 = 1, a2 = 1, a3 = 1, a4 = 1, a5 = 1, u = 0\\.0476"
    )
  )
})

test_that("the range search stays on the simplex at its corners (#560)", {
  # off the simplex by 2e-9, (gA + gB + gC)^2000 would be 1 + 4e-6
  power <- "(gA + gB + gC)^2000 * gA * D"
  expect_message(
    mpt(mpt_tree("t", list(A = power, B = glue("1 - {power}"))), simplex = c("gA", "gB", "gC")),
    "category 'A' in tree 't', category 'B' in tree 't' is approximate.*local search"
  )
})

test_that("the search message states why the branches are searched (#560)", {
  expect_message(
    mpt(mpt_tree("t", list(y = "plogis(a)", n = "1 - plogis(a)"))),
    paste0(
      "category 'y' in tree 't', category 'n' in tree 't' is approximate: ",
      "stats::D\\(\\) cannot differentiate them.*local search"
    )
  )
  repeated <- mpt(mpt_tree("t", list(hit = "1 - (1 - D)^n", miss = "(1 - D)^n")), covariates = "n")
  expect_message(
    check_data(repeated, data.frame(n = 1:300, hit = 5, miss = 5), bmf(D ~ 1)),
    paste0(
      "category 'hit' in tree 't', category 'miss' in tree 't' is approximate: ",
      "their tree has the form for which the check is exact, but the exact ",
      "check would need more than 65,536 points \\(a grid over the parameters ",
      "for each distinct covariate setting in the tree's rows\\).*local search"
    )
  )
  expect_message(
    mpt(mpt_tree("t", list(y = "0.5 + 0.5 * (a - 0.3)^2", n = "0.5 - 0.5 * (a - 0.3)^2"))),
    paste0(
      "is approximate: a parameter enters them more than once.*",
      "not written with parameters.*local search"
    )
  )
  # one message per reason
  msgs <- testthat::capture_messages(bmm:::.mpt_message_sampled(
    c(grid = "category 'a' in tree 's'", derivative = "category 'b' in tree 't'"), "search"
  ))
  expect_length(msgs, 2L)
  expect_match(msgs[1], "category 'a' in tree 's' is approximate: their tree has the form")
  expect_match(msgs[2], "category 'b' in tree 't' is approximate: stats::D\\(\\)")
})

test_that("the identity and grid reasons read as limits of the check (#560)", {
  # a valid tree whose high power leaves the sum of its branches 1 only to
  # within the rounding of the check off the simplex
  expect_message(
    mpt(list(
      mpt_tree("t", list(x = "(1 - gC)^60", y = "1 - (1 - gC)^60"), impossible = "z"),
      mpt_tree("g", list(x = "gA", y = "gB", z = "gC"))
    ), tree_id = "tr", simplex = c("gA", "gB", "gC")),
    paste0(
      "category 'x' in tree 't', category 'y' in tree 't' is approximate: their ",
      "tree has the form for which the check is exact, but the sum of its ",
      "branches could not be confirmed as 1 within the numerical precision of ",
      "the check, as happens with high powers such as \\(1 - g\\)\\^60 of a ",
      "member of a simplex group\\. A local search"
    )
  )
  # without covariates the grid is over the parameters alone, here 41^3 points
  expect_message(
    mpt(mpt_tree("t", list(y = "(u * v * w)^40", n = "1 - (u * v * w)^40"))),
    "more than 65,536 points \\(a grid over the parameters\\)\\. A local search"
  )
})

test_that("the search message names the reason for pair clustering over the grid limit (#560)", {
  skip_on_cran()
  # in the certified form, but 12 grid points for each of 6000 covariate rows
  pairs <- mpt(mpt_tree("pairs", list(
    E1 = "c * r * G", E2 = "(1 - c * G) * u^2",
    E3 = "(1 - c * G) * 2 * u * (1 - u)", E4 = "c * G * (1 - r) + (1 - c * G) * (1 - u)^2"
  )), covariates = "G")
  data <- data.frame(G = seq(0.2, 0.9, length.out = 6000), E1 = 1, E2 = 2, E3 = 3, E4 = 4)
  expect_message(
    check_data(pairs, data, bmf(c ~ 1, r ~ 1, u ~ 1)),
    paste0(
      "category 'E2' in tree 'pairs', category 'E3' in tree 'pairs', category ",
      "'E4' in tree 'pairs' is approximate: their tree has the form for which ",
      "the check is exact.*more than 65,536 points.*local search.*",
      "The search left unsearched 11,488 of the 12,000 combinations of a ",
      "distinct covariate setting with a vertex of the parameters it holds ",
      "fixed and capped its further starts for category 'E2' in tree 'pairs', ",
      "category 'E3' in tree 'pairs'; it left unsearched 23,488 of the 24,000 .*",
      "for category 'E4' in tree 'pairs'\\."
    )
  )
})
