test_that("mpt_tree folds integer fractions into decimal literals", {
  tree <- mpt_tree("t", list(
    a = "p + (1 - p) * (1/4)",
    b = "(1 - p) * (3/4)"
  ))
  expect_equal(deparse1(tree$branches$a), "p + (1 - p) * 0.25")
  expect_equal(deparse1(tree$branches$b), "(1 - p) * 0.75")
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
  expect_error(.mpt_parse_restrictions(list("Do")), "named")
})

test_that("restriction chains resolve to their final target", {
  resolved <- .mpt_resolve_restrictions(list(A = quote(B), B = quote(C), g = 0.5))
  expect_equal(resolved, list(A = quote(C), B = quote(C), g = 0.5))
})

test_that("mpt warns when branch probabilities do not sum to 1", {
  bad_tree <- mpt_tree("t", list(a = "D * g", b = "(1 - D) * g"))
  expect_warning(mpt(bad_tree), "sum to")
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
