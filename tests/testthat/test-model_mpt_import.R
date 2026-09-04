test_that("mpt_from_string parses MPTinR-style model definitions", {
  model_2htm <- "
  D + (1 - D) * g        # old
  (1 - D) * (1 - g)      # new

  (1 - D) * g            # old
  D + (1 - D) * (1 - g)  # new
  "
  model <- mpt_from_string(
    model_2htm, tree_names = c("old", "new"), tree_id = "item_type"
  )
  manual <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  expect_equal(model$trees, manual$trees)
  expect_equal(names(model$parameters), names(manual$parameters))

  summed <- mpt_from_string(
    "D # hit\n(1 - D) * g # hit\n(1 - D) * (1 - g) # miss",
    tree_names = "old"
  )
  expect_equal(deparse1(summed$trees$old$branches$hit), "(D) + ((1 - D) * g)")

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
    tree_id = "item_type"
  ))
  expect_setequal(names(model$parameters), c("Do", "Dn", "guess"))
  expect_setequal(model$resp_vars$resp_cats, c("yes", "no"))
  expect_false("Gfix" %in% .mpt_expr_vars(model$trees$old))
  expect_true(grepl("0.25", deparse1(model$trees$old$branches$yes), fixed = TRUE))

  renaming <- attr(model, "mpt_renaming")
  expect_equal(renaming[["D_o"]], "Do")
  expect_equal(renaming[["old_hit"]], "yes")

  # restrictions are written with the EQN names and renamed with the branches
  string_form <- suppressMessages(mpt_from_eqn(
    eqn_file, restrictions = "G_fix = 1/4", categories = category_map,
    tree_id = "item_type"
  ))
  expect_equal(string_form$trees, model$trees)
  expect_equal(string_form$restrictions, list(Gfix = 0.25))

  equated <- suppressMessages(mpt_from_eqn(
    eqn_file, restrictions = c("G_fix = 1/4", "D_n = D_o"),
    categories = category_map, tree_id = "item_type"
  ))
  expect_setequal(names(equated$parameters), c("Do", "guess"))
  expect_equal(deparse1(equated$trees$new$branches$no), "(Do) + ((1 - Do) * (1 - guess * 0.25))")

  expect_error(
    suppressMessages(mpt_from_eqn(
      eqn_file, restrictions = "D_o > D_n", categories = category_map
    )),
    "Order constraints"
  )

  # per-tree category labels without a mapping cannot be combined
  expect_error(
    suppressMessages(mpt_from_eqn(eqn_file, restrictions = c(G_fix = 0.25))),
    "categories argument"
  )
})

test_that("mpt_from_string passes restrictions on to mpt()", {
  model <- mpt_from_string(
    "D + (1 - D) * g # old\n(1 - D) * (1 - g) # new",
    tree_names = "old", restrictions = "g = 0.5"
  )
  expect_equal(names(model$parameters), "D")
  expect_equal(model$restrictions, list(g = 0.5))
})

test_that("mpt_from_eqn renames correctly when one name prefixes another", {
  # renaming substitutes whole symbols, so 'd_A' never touches 'd_A.x'
  eqn_file <- tempfile(fileext = ".eqn")
  writeLines(c(
    "t  a  d_A*d_A.x",
    "t  b  1-d_A*d_A.x"
  ), eqn_file)
  model <- suppressMessages(mpt_from_eqn(eqn_file))
  expect_setequal(names(model$parameters), c("dA", "dAx"))
  expect_equal(attr(model, "mpt_renaming")[["d_A.x"]], "dAx")
})

test_that("mpt_from_eqn errors on names that clash after sanitizing", {
  eqn_file <- tempfile(fileext = ".eqn")
  writeLines(c(
    "t  a  d_A + dA",
    "t  b  1 - d_A - dA"
  ), eqn_file)
  expect_error(mpt_from_eqn(eqn_file), "duplicated names")
})
