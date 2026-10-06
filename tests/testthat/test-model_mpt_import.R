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
  expect_equal(model$other_vars$trees, manual$other_vars$trees)
  expect_equal(names(model$parameters), names(manual$parameters))

  summed <- suppressMessages(mpt_from_string(
    "D # hit\n(1 - D) * g # hit\n(1 - D) * (1 - g) # miss",
    tree_names = "old"
  ))
  expect_equal(deparse1(summed$other_vars$trees$old$branches$hit), "(D) + ((1 - D) * g)")

  no_comments <- mpt_from_string(
    "D + (1 - D) * g\n(1 - D) * (1 - g)",
    tree_names = "old", categories = c("hit", "miss")
  )
  expect_equal(names(no_comments$other_vars$trees$old$branches), c("hit", "miss"))

  expect_error(
    mpt_from_string("D # a\n1 - D # b", tree_names = c("t1", "t2")),
    "tree block"
  )
  expect_error(
    mpt_from_string("D + (1 - D) * g\n(1 - D) * (1 - g)", tree_names = "old"),
    "categories"
  )
})

test_that("mpt_from_string says when it sums lines that share a category", {
  expect_message(
    mpt_from_string(
      "D # hit\n(1 - D) * g # hit\n(1 - D) * (1 - g) # miss",
      tree_names = "old"
    ),
    "tree old: 'line 1', 'line 2' -> hit"
  )
  # the label may come from the categories argument, and each tree is named
  expect_message(
    mpt_from_string(
      c("D", "(1 - D) * g", "(1 - D) * (1 - g)", "", "g", "1 - g"),
      tree_names = c("old", "new"), tree_id = "item_type",
      categories = list(old = c("hit", "hit", "miss"), new = c("hit", "miss"))
    ),
    "tree old: 'line 1', 'line 2' -> hit"
  )
  expect_no_message(
    mpt_from_string(
      "D # hit\n1 - D # miss\n\ng # hit\n1 - g # miss",
      tree_names = c("old", "new"), tree_id = "item_type"
    ),
    message = "share a response category"
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
  expect_false("Gfix" %in% .mpt_expr_vars(model$other_vars$trees$old))
  expect_true(grepl("0.25", deparse1(model$other_vars$trees$old$branches$yes), fixed = TRUE))

  renaming <- attr(model, "mpt_renaming")
  expect_equal(renaming$bmm_name[renaming$file_name == "D_o"], "Do")
  expect_equal(
    renaming[renaming$file_name == "old_hit", c("kind", "tree", "bmm_name")],
    data.frame(kind = "category", tree = "old", bmm_name = "yes"),
    ignore_attr = "row.names"
  )

  # restrictions are written with the EQN names and renamed with the branches
  string_form <- suppressMessages(mpt_from_eqn(
    eqn_file, restrictions = "G_fix = 1/4", categories = category_map,
    tree_id = "item_type"
  ))
  expect_equal(string_form$other_vars$trees, model$other_vars$trees)
  expect_equal(string_form$other_vars$restrictions, list(Gfix = 0.25))

  equated <- suppressMessages(mpt_from_eqn(
    eqn_file, restrictions = c("G_fix = 1/4", "D_n = D_o"),
    categories = category_map, tree_id = "item_type"
  ))
  expect_setequal(names(equated$parameters), c("Do", "guess"))
  expect_equal(deparse1(equated$other_vars$trees$new$branches$no), "(Do) + ((1 - Do) * (1 - guess * 0.25))")

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
  expect_equal(model$other_vars$restrictions, list(g = 0.5))
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
  renaming <- attr(model, "mpt_renaming")
  expect_equal(renaming$bmm_name[renaming$file_name == "d_A.x"], "dAx")
})

test_that("sanitized names map onto brms-safe names", {
  expect_equal(
    .mpt_sanitized_names(c("d_A", "d.B", "g")),
    list(d_A = "dA", d.B = "dB", g = "g")
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

write_eqn <- function(lines) {
  eqn_file <- tempfile(fileext = ".eqn")
  writeLines(lines, eqn_file)
  eqn_file
}

eqn_2htm_lines <- c(
  "old  hit   D_o",
  "old  hit   (1-D_o)*g",
  "old  miss  (1-D_o)*(1-g)",
  "new  fa    (1-D_n)*g",
  "new  cr    D_n",
  "new  cr    (1-D_n)*(1-g)"
)
eqn_2htm_map <- c(hit = "yes", fa = "yes", miss = "no", cr = "no")

test_that("mpt_from_eqn strips the tree prefix and sanitizes unmapped categories", {
  model <- suppressMessages(mpt_from_eqn(write_eqn(c(
    "old  old_yes  D_o",
    "old  old_yes  (1-D_o)*g",
    "old  old_no   (1-D_o)*(1-g)",
    "new  new_yes  (1-D_n)*g",
    "new  new_no   D_n",
    "new  new_no   (1-D_n)*(1-g)"
  )), tree_id = "item_type"))
  expect_setequal(model$resp_vars$resp_cats, c("yes", "no"))
  renaming <- attr(model, "mpt_renaming")
  expect_equal(
    renaming$bmm_name[renaming$tree %in% "old" & renaming$file_name == "old_yes"],
    "yes"
  )

  underscored <- suppressMessages(mpt_from_eqn(write_eqn(c(
    "t  my_hit   D",
    "t  my_hit   (1-D)*g",
    "t  my_miss  (1-D)*(1-g)"
  ))))
  expect_setequal(underscored$resp_vars$resp_cats, c("myhit", "mymiss"))
})

test_that("mpt_from_eqn errors when the tree prefix strip merges categories", {
  clash_file <- write_eqn(c(
    "A  A_x  0.2*g",
    "A  x    0.3*g",
    "A  y    1-0.5*g"
  ))
  expect_error(
    suppressMessages(mpt_from_eqn(clash_file)),
    "A_x.*x.*would all be named 'x'"
  )

  # the strip alone, without a clash, is fine
  expect_no_error(suppressMessages(mpt_from_eqn(write_eqn(c(
    "A  A_x  0.2*g",
    "A  y    1-0.2*g"
  )))))

  # categories the user maps on purpose may be merged
  merged <- suppressMessages(mpt_from_eqn(
    clash_file, categories = c(A_x = "z", x = "z")
  ))
  expect_setequal(merged$resp_vars$resp_cats, c("z", "y"))
})

test_that("mpt_from_eqn keeps the tree prefix of a mapped category", {
  # the prefix strip applies to file labels only, never to a name the user
  # chose in categories
  model <- suppressMessages(mpt_from_eqn(
    write_eqn(c(
      "old  hit   D",
      "old  hit   (1-D)*g",
      "old  miss  (1-D)*(1-g)"
    )),
    categories = c(hit = "old_yes", miss = "old_no")
  ))
  expect_equal(model$resp_vars$resp_cats, c("oldyes", "oldno"))
})

test_that("mpt_from_eqn errors when sanitizing merges categories across trees", {
  expect_error(
    mpt_from_eqn(write_eqn(c(
      "A  d_o  g",
      "A  y    1-g",
      "B  do   g",
      "B  y    1-g"
    )), tree_id = "t"),
    "duplicated names"
  )
})

test_that("mpt_from_eqn errors when a sanitized parameter name equals a covariate", {
  clash_file <- write_eqn(c("t yes a_b", "t no 1-a_b", "u yes ab*g", "u no 1-ab*g"))
  expect_error(
    mpt_from_eqn(clash_file, covariates = "ab", tree_id = "tr"),
    "duplicated names or clashes with a covariate: 'ab'"
  )
  # without the covariate the same file is caught by the duplicate check
  expect_error(mpt_from_eqn(clash_file, tree_id = "tr"), "duplicated names")

  # a covariate that the sanitizing leaves alone is no clash
  model <- suppressMessages(mpt_from_eqn(
    write_eqn(c("t yes a_b*x_1", "t no 1-a_b*x_1")), covariates = "x_1"
  ))
  expect_equal(names(model$parameters), "ab")
})

test_that("mpt_from_eqn skips # comment lines", {
  commented <- suppressMessages(mpt_from_eqn(
    write_eqn(c(
      "# 2HTM recognition model, see the paper",
      "6",
      eqn_2htm_lines[1:3],
      "# the new-item tree",
      eqn_2htm_lines[4:6]
    )),
    categories = eqn_2htm_map, tree_id = "item_type"
  ))
  plain <- suppressMessages(mpt_from_eqn(
    write_eqn(c("6", eqn_2htm_lines)),
    categories = eqn_2htm_map, tree_id = "item_type"
  ))
  expect_equal(commented$other_vars$trees, plain$other_vars$trees)
})

test_that("mpt_from_string takes a vector of lines like a single string", {
  model_text <- c(
    "D + (1 - D) * g        # old",
    "(1 - D) * (1 - g)      # new",
    "",
    "(1 - D) * g            # old",
    "D + (1 - D) * (1 - g)  # new"
  )
  from_vector <- mpt_from_string(
    model_text, tree_names = c("old", "new"), tree_id = "item_type"
  )
  from_string <- mpt_from_string(
    paste(model_text, collapse = "\n"),
    tree_names = c("old", "new"), tree_id = "item_type"
  )
  expect_equal(from_vector$other_vars$trees, from_string$other_vars$trees)
  expect_equal(names(from_vector$other_vars$trees), c("old", "new"))

  from_file <- tempfile()
  writeLines(model_text, from_file)
  expect_equal(
    mpt_from_string(
      readLines(from_file), tree_names = c("old", "new"), tree_id = "item_type"
    )$other_vars$trees,
    from_string$other_vars$trees
  )
})

test_that("mpt_from_eqn errors on categories that are not in the file", {
  eqn_file <- write_eqn(eqn_2htm_lines)
  expect_error(
    suppressMessages(mpt_from_eqn(
      eqn_file, categories = c(eqn_2htm_map, hitt = "yes"), tree_id = "item_type"
    )),
    "hitt.*The file uses: 'hit', 'miss', 'fa', 'cr'"
  )
})

test_that("mpt_from_eqn errors on restrictions that use names not in the file", {
  eqn_file <- write_eqn(c(
    "old  hit   D_o",
    "old  hit   (1-D_o)*g",
    "old  miss  (1-D_o)*(1-g)",
    "new  fa    (1-D_n)*g",
    "new  cr    D_n",
    "new  cr    (1-D_n)*(1-g)"
  ))
  expect_error(
    suppressMessages(mpt_from_eqn(
      eqn_file, restrictions = "Dn = Do", categories = eqn_2htm_map,
      tree_id = "item_type"
    )),
    "not in the EQN file: 'Dn', 'Do'.*as the file does"
  )
  expect_error(
    suppressMessages(mpt_from_eqn(
      eqn_file, restrictions = list(D_n = "Do"), categories = eqn_2htm_map,
      tree_id = "item_type"
    )),
    "not in the EQN file: 'Do'"
  )
  expect_no_error(suppressMessages(mpt_from_eqn(
    eqn_file, restrictions = "D_n = D_o", categories = eqn_2htm_map,
    tree_id = "item_type"
  )))
  expect_error(
    suppressMessages(mpt_from_eqn(
      eqn_file, restrictions = "g = FE", categories = eqn_2htm_map,
      tree_id = "item_type"
    )),
    "TreeBUGS"
  )
})

test_that("mpt_from_eqn does not fuse whitespace-separated tokens", {
  expect_error(
    suppressMessages(mpt_from_eqn(write_eqn(c("t  a  g h", "t  b  1 - g")))),
    "unexpected symbol"
  )
})

test_that("mpt_from_eqn reads a file without a final newline silently", {
  eqn_file <- tempfile(fileext = ".eqn")
  cat("t  a  g\nt  b  1 - g", file = eqn_file)
  expect_no_warning(mpt_from_eqn(eqn_file))
})

test_that("mpt_from_eqn skips a lone integer header and checks its count", {
  plain <- suppressMessages(mpt_from_eqn(
    write_eqn(eqn_2htm_lines), categories = eqn_2htm_map, tree_id = "item_type"
  ))
  expect_no_warning(with_header <- suppressMessages(mpt_from_eqn(
    write_eqn(c("  6  ", eqn_2htm_lines)),
    categories = eqn_2htm_map, tree_id = "item_type"
  )))
  expect_equal(with_header$other_vars$trees, plain$other_vars$trees)

  expect_warning(
    suppressMessages(mpt_from_eqn(
      write_eqn(c("8", eqn_2htm_lines)),
      categories = eqn_2htm_map, tree_id = "item_type"
    )),
    "8 equations.*6 equation lines"
  )
})

test_that("mpt_from_eqn reports short lines and unparsable lines by number", {
  expect_error(
    mpt_from_eqn(write_eqn(c("t a", "t  a  g", "t  b  1 - g"))),
    "line 1: 't a'"
  )
  expect_error(
    mpt_from_eqn(write_eqn(c("t  a  g", "", "t b", "t  b  1 - g"))),
    "line 3: 't b'"
  )
  # a prose header is read as an equation, as every line after it
  expect_error(
    mpt_from_eqn(write_eqn(c("2HTM recognition model", eqn_2htm_lines))),
    "'2HTM' \\(first on line 1\\)"
  )
  expect_error(
    mpt_from_eqn(write_eqn(c("Model for recognition data", eqn_2htm_lines))),
    "line 1 .*'Model for recognition data'"
  )
  expect_error(
    mpt_from_eqn(write_eqn(c("t  a  g", "t  b  1 - g h"))),
    "line 2.*'t  b  1 - g h'"
  )
})

pc_eqn_lines <- c(
  "pairs1  E1  c1*r1",
  "pairs1  E2  (1-c1)*u1*u1",
  "pairs1  E3  2*(1-c1)*u1*(1-u1)",
  "pairs1  E4  c1*(1-r1)+(1-c1)*(1-u1)*(1-u1)",
  "single1 F1  a1",
  "single1 F2  1-a1",
  "pairs2  E1  c2*r2",
  "pairs2  E2  (1-c2)*u2*u2",
  "pairs2  E3  2*(1-c2)*u2*(1-u2)",
  "pairs2  E4  c2*(1-r2)+(1-c2)*(1-u2)*(1-u2)",
  "single2 F1  a2",
  "single2 F2  1-a2"
)
pc_string <- "
c1 * r1                                         # E1
(1 - c1) * u1 * u1                              # E2
2 * (1 - c1) * u1 * (1 - u1)                    # E3
c1 * (1 - r1) + (1 - c1) * (1 - u1) * (1 - u1)  # E4

a1      # F1
1 - a1  # F2

c2 * r2                                         # E1
(1 - c2) * u2 * u2                              # E2
2 * (1 - c2) * u2 * (1 - u2)                    # E3
c2 * (1 - r2) + (1 - c2) * (1 - u2) * (1 - u2)  # E4

a2      # F1
1 - a2  # F2
"
pc_trees <- c("pairs1", "single1", "pairs2", "single2")
pc_impossible <- list(
  pairs1 = c("F1", "F2"), single1 = c("E1", "E2", "E3", "E4"),
  pairs2 = c("F1", "F2"), single2 = c("E1", "E2", "E3", "E4")
)

# the pair-clustering model built by hand, with a1 = u1 and a2 = u2 as in
# MPTinR's restriction file for it
pc_manual <- function() {
  pair_tree <- function(k) {
    mpt_tree(paste0("pairs", k), lapply(list(
      E1 = "c# * r#", E2 = "(1 - c#) * u# * u#",
      E3 = "2 * (1 - c#) * u# * (1 - u#)",
      E4 = "c# * (1 - r#) + (1 - c#) * (1 - u#) * (1 - u#)"
    ), gsub, pattern = "#", replacement = k), impossible = c("F1", "F2"))
  }
  single_tree <- function(k) {
    mpt_tree(
      paste0("single", k), list(F1 = paste0("a", k), F2 = paste0("1 - a", k)),
      impossible = c("E1", "E2", "E3", "E4")
    )
  }
  mpt(
    list(pair_tree(1), single_tree(1), pair_tree(2), single_tree(2)),
    tree_id = "tree", restrictions = c("a1 = u1", "a2 = u2")
  )
}

test_that("both importers build the pair-clustering model with impossible", {
  manual <- pc_manual()
  from_eqn <- mpt_from_eqn(
    write_eqn(pc_eqn_lines), restrictions = c("a1 = u1", "a2 = u2"),
    tree_id = "tree", impossible = pc_impossible
  )
  from_string <- mpt_from_string(
    pc_string, tree_names = pc_trees, tree_id = "tree",
    restrictions = c("a1 = u1", "a2 = u2"), impossible = pc_impossible
  )
  ignored <- c("call", "mpt_source", "mpt_renaming")
  expect_equal(from_eqn, manual, ignore_attr = ignored)
  expect_equal(from_string, manual, ignore_attr = ignored)
})

test_that("importers name the categories each tree lacks", {
  # line 1 opens a tree that later lines continue, so it is no title
  expect_no_match(
    tryCatch(
      mpt_from_eqn(write_eqn(pc_eqn_lines), tree_id = "tree"),
      error = conditionMessage
    ),
    "title"
  )
  expect_error(
    mpt_from_eqn(write_eqn(pc_eqn_lines), tree_id = "tree"),
    paste0(
      "pairs1: 'F1', 'F2'.*single1: 'E1', 'E2', 'E3', 'E4'.*",
      "impossible = list\\(pairs1 = c\\('F1', 'F2'\\).*categories"
    )
  )
  expect_error(
    mpt_from_string(pc_string, tree_names = pc_trees, tree_id = "tree"),
    "pairs1: 'F1', 'F2'.*impossible.*same category label"
  )
  # a partial impossible list names only the trees still lacking categories
  expect_error(
    mpt_from_eqn(
      write_eqn(pc_eqn_lines), tree_id = "tree",
      impossible = pc_impossible[c("pairs1", "single1")]
    ),
    "trees use:\n *pairs2: 'F1', 'F2'\n *single2"
  )
})

test_that("importers reject impossible entries the model does not have", {
  eqn_file <- write_eqn(pc_eqn_lines)
  expect_error(
    mpt_from_eqn(
      eqn_file, tree_id = "tree", impossible = c(pairs1 = "F1")
    ),
    "named list"
  )
  expect_error(
    mpt_from_eqn(
      eqn_file, tree_id = "tree",
      impossible = c(pc_impossible, list(pairs3 = "F1"))
    ),
    "trees that are not in the model: 'pairs3'"
  )
  expect_error(
    mpt_from_string(
      pc_string, tree_names = pc_trees, tree_id = "tree",
      impossible = modifyList(pc_impossible, list(pairs1 = c("F1", "F_2")))
    ),
    "categories that are not in the model: 'F_2'"
  )
})

test_that("mpt_from_string skips # comment lines without splitting trees", {
  commented <- mpt_from_string(
    c(
      "# 2HTM, Snodgrass & Corwin (1988)",
      "  # old-item tree",
      "D + (1 - D) * g        # old",
      "# the miss branch",
      "(1 - D) * (1 - g)      # new",
      "",
      "# new-item tree",
      "(1 - D) * g            # old",
      "D + (1 - D) * (1 - g)  # new"
    ),
    tree_names = c("old", "new"), tree_id = "item_type"
  )
  manual <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  expect_equal(commented$other_vars$trees, manual$other_vars$trees)
})

test_that("mpt_from_eqn announces categories that its map merges in a tree", {
  clash_file <- write_eqn(c(
    "A  A_x  0.2*g",
    "A  x    0.3*g",
    "A  y    1-0.5*g"
  ))
  suppressMessages(expect_message(
    mpt_from_eqn(clash_file, categories = c(A_x = "z", x = "z")),
    "tree A: 'A_x', 'x' -> z"
  ))
  suppressMessages(expect_message(
    mpt_from_eqn(clash_file, categories = c(A_x = "y", x = "x")),
    "tree A: 'A_x', 'y' -> y"
  ))
  # the same label in different trees is no merge
  suppressMessages(expect_no_message(
    mpt_from_eqn(
      write_eqn(eqn_2htm_lines), categories = eqn_2htm_map,
      tree_id = "item_type"
    ),
    message = "several categories"
  ))
})

test_that("mpt_from_eqn errors on an unnamed categories vector", {
  eqn_file <- write_eqn(c("t  a  g", "t  b  1 - g"))
  expect_error(
    mpt_from_eqn(eqn_file, categories = c("x", "y")),
    "named character vector"
  )
  expect_error(
    mpt_from_eqn(eqn_file, categories = c(a = "x", "y")),
    "named character vector"
  )
})

test_that("mpt_from_eqn refuses a categories or tree_names map that is not one name per entry", {
  eqn_file <- write_eqn(c(
    "old hit D", "old hit (1-D)*g", "old miss (1-D)*(1-g)",
    "new fa (1-D)*g", "new cr D", "new cr (1-D)*(1-g)"
  ))
  expect_error(
    mpt_from_eqn(eqn_file, categories = list(
      hit = c("yes", "no"), miss = "no", fa = "yes", cr = "no"
    )),
    "categories.*'hit' has 2 values"
  )
  expect_error(
    mpt_from_eqn(eqn_file, categories = as.list(eqn_2htm_map)),
    "categories.*named character vector.*It is a list"
  )
  expect_error(
    mpt_from_eqn(
      eqn_file, categories = c(hit = "yes", hit = "no", miss = "no", fa = "yes")
    ),
    "categories.*'hit' is given more than once"
  )
  expect_error(
    mpt_from_eqn(eqn_file, categories = c(hit = "", miss = "no")),
    "categories.*'hit' has no value"
  )
  expect_error(
    mpt_from_eqn(eqn_file, categories = c(hit = "yes", "no")),
    "categories.*named character vector"
  )

  numbered <- write_eqn(c("1 yes g", "1 no 1-g", "2 yes h", "2 no 1-h"))
  expect_error(
    mpt_from_eqn(numbered, tree_names = list("1" = c("old", "new"))),
    "tree_names.*'1' has 2 values"
  )
  expect_error(
    mpt_from_eqn(numbered, tree_names = c("1" = "old", "1" = "new")),
    "tree_names.*'1' is given more than once"
  )
  expect_error(
    mpt_from_eqn(numbered, tree_names = list("1" = "old", "2" = "new")),
    "tree_names.*named character vector.*It is a list"
  )
})

test_that("mpt_from_string checks the entries of a categories list by tree", {
  model_text <- c("D", "1 - D", "", "g", "1 - g")
  by_tree <- function(categories) {
    mpt_from_string(
      model_text, tree_names = c("old", "new"), tree_id = "item_type",
      categories = categories
    )
  }
  expect_no_error(by_tree(list(old = c("yes", "no"), new = c("yes", "no"))))
  expect_error(
    by_tree(list(old = c("yes", "no"), old = c("no", "yes"), new = c("yes", "no"))),
    "names a tree more than once: 'old'"
  )
  expect_error(
    by_tree(list(old = c("yes", NA), new = c("yes", "no"))),
    "character vectors without missing or empty values.*'old'"
  )
  expect_error(
    by_tree(list(old = c("yes", "no"), new = 1:2)),
    "character vectors without missing or empty values.*'new'"
  )
})

test_that("mpt_from_eqn maps tree labels through tree_names", {
  numbered <- write_eqn(c(
    "1  hit   D_o",
    "1  hit   (1-D_o)*g",
    "1  miss  (1-D_o)*(1-g)",
    "2  fa    (1-D_n)*g",
    "2  cr    D_n",
    "2  cr    (1-D_n)*(1-g)"
  ))
  named <- suppressMessages(mpt_from_eqn(
    write_eqn(eqn_2htm_lines), categories = eqn_2htm_map, tree_id = "item_type"
  ))
  renamed <- suppressMessages(mpt_from_eqn(
    numbered, categories = eqn_2htm_map, tree_id = "item_type",
    tree_names = c("1" = "old", "2" = "new")
  ))
  expect_equal(renamed$other_vars$trees, named$other_vars$trees)
  renaming <- attr(renamed, "mpt_renaming")
  expect_equal(
    renaming[renaming$kind == "tree", c("file_name", "bmm_name")],
    data.frame(file_name = c("1", "2"), bmm_name = c("old", "new")),
    ignore_attr = "row.names"
  )

  expect_error(
    suppressMessages(mpt_from_eqn(
      numbered, categories = eqn_2htm_map, tree_id = "item_type"
    )),
    "rename '1' \\(first on line 1\\), '2' \\(first on line 4\\) with"
  )
  # an unmapped label keeps its name
  expect_error(
    suppressMessages(mpt_from_eqn(
      numbered, categories = eqn_2htm_map, tree_id = "item_type",
      tree_names = c("1" = "old")
    )),
    "'2'.*tree_names"
  )
  expect_error(
    suppressMessages(mpt_from_eqn(
      numbered, categories = eqn_2htm_map, tree_id = "item_type",
      tree_names = c("1" = "old", "3" = "new")
    )),
    "not in the EQN file: '3'"
  )
  expect_error(
    suppressMessages(mpt_from_eqn(
      numbered, categories = eqn_2htm_map, tree_id = "item_type",
      tree_names = c("old", "new")
    )),
    "named character vector"
  )
  expect_error(
    suppressMessages(mpt_from_eqn(
      numbered, categories = eqn_2htm_map, tree_id = "item_type",
      tree_names = c("1" = "same", "2" = "same")
    )),
    "same name.*'same'"
  )
})

test_that("importers record where the model came from", {
  # the lines are stored as read, untrimmed
  eqn_file <- write_eqn(c("# 2HTM", "  6  ", eqn_2htm_lines))
  model <- suppressMessages(mpt_from_eqn(
    eqn_file, categories = eqn_2htm_map, tree_id = "item_type"
  ))
  expect_equal(attr(model, "mpt_source"), list(
    file = normalizePath(eqn_file),
    md5 = unname(tools::md5sum(eqn_file)),
    lines = c("# 2HTM", "  6  ", eqn_2htm_lines)
  ))

  model_text <- c("D + (1 - D) * g # old", "(1 - D) * (1 - g) # new")
  from_string <- mpt_from_string(model_text, tree_names = "old")
  expect_equal(
    attr(from_string, "mpt_source"),
    list(text = paste(model_text, collapse = "\n"))
  )
})

test_that("mpt_renaming is a data frame of every name and whether it changed", {
  model <- suppressMessages(mpt_from_eqn(write_eqn(c(
    "old  old_yes  D_o",
    "old  old_yes  (1-D_o)*g",
    "old  old_no   (1-D_o)*(1-g)",
    "new  new_yes  (1-D_n)*g",
    "new  new_no   D_n",
    "new  new_no   (1-D_n)*(1-g)"
  )), tree_id = "item_type"))
  expect_equal(
    attr(model, "mpt_renaming"),
    data.frame(
      kind = c(rep("parameter", 3), rep("tree", 2), rep("category", 4)),
      tree = c(rep(NA, 5), "old", "old", "new", "new"),
      file_name = c(
        "D_o", "g", "D_n", "old", "new", "old_yes", "old_no", "new_yes", "new_no"
      ),
      bmm_name = c("Do", "g", "Dn", "old", "new", "yes", "no", "yes", "no"),
      renamed = c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
    )
  )

  # the message lists one parameter per line and the categories per tree
  messages <- testthat::capture_messages(mpt_from_eqn(
    write_eqn(eqn_2htm_lines), categories = eqn_2htm_map, tree_id = "item_type"
  ))
  expect_match(
    messages,
    paste0(
      "\n +parameter +D_o -> Do\n +parameter +D_n -> Dn\n",
      " +categories in tree old: hit -> yes, miss -> no\n",
      " +categories in tree new: fa -> yes, cr -> no"
    ),
    all = FALSE
  )

  # no empty parameter line when only categories were renamed
  expect_match(
    testthat::capture_messages(mpt_from_eqn(
      write_eqn(c("t  t_a  g", "t  t_b  1 - g"))
    )),
    "^[^\n]*\n  categories in tree t: t_a -> a, t_b -> b\n$"
  )

  # nothing renamed: every name is still listed, and no message is shown
  expect_no_message(
    unchanged <- mpt_from_eqn(write_eqn(c("t  a  g", "t  b  1 - g")))
  )
  expect_equal(
    attr(unchanged, "mpt_renaming"),
    data.frame(
      kind = c("parameter", "tree", "category", "category"),
      tree = c(NA, NA, "t", "t"), file_name = c("g", "t", "a", "b"),
      bmm_name = c("g", "t", "a", "b"), renamed = FALSE
    )
  )
})

# MPTinR's pair-clustering model file: no category labels, '#' header lines
pc_mptinr_model <- c(
  "# Pair-clustering model (Batchelder & Riefer, 1986), two list lengths",
  "# pairs tree: E1 both adjacent, E2 both non-adjacent, E3 one, E4 none",
  "# singleton tree: F1 recalled, F2 not recalled",
  "",
  "c1 * r1",
  "(1-c1) * u1 * u1",
  "2 * (1-c1) * u1 * (1-u1)",
  "c1 * (1-r1) + (1-c1) * (1-u1) * (1-u1)",
  "",
  "a1",
  "1 - a1",
  "",
  "c2 * r2",
  "(1-c2) * u2 * u2",
  "2 * (1-c2) * u2 * (1-u2)",
  "c2 * (1-r2) + (1-c2) * (1-u2) * (1-u2)",
  "",
  "a2",
  "1 - a2"
)
pc_categories <- list(
  pairs1 = c("E1", "E2", "E3", "E4"), single1 = c("F1", "F2"),
  pairs2 = c("E1", "E2", "E3", "E4"), single2 = c("F1", "F2")
)

test_that("mpt_from_string labels unlabelled MPTinR trees with a list by tree", {
  from_mptinr <- mpt_from_string(
    pc_mptinr_model, tree_names = pc_trees, categories = pc_categories,
    tree_id = "tree", restrictions = c("a1 = u1", "a2 = u2"),
    impossible = pc_impossible
  )
  expect_equal(
    from_mptinr, pc_manual(),
    ignore_attr = c("call", "mpt_source", "mpt_renaming")
  )

  # an inline label takes precedence over the entry for its line
  relabelled <- mpt_from_string(
    c("D + (1 - D) * g # yes", "(1 - D) * (1 - g)", "", "(1 - D) * g",
      "D + (1 - D) * (1 - g)"),
    tree_names = c("old", "new"), tree_id = "item_type",
    categories = list(old = c("hit", "no"), new = c("yes", "no"))
  )
  expect_equal(
    names(relabelled$other_vars$trees$old$branches), c("yes", "no")
  )
})

test_that("mpt_from_string checks a categories list against the trees", {
  pc_call <- function(categories) {
    mpt_from_string(
      pc_mptinr_model, tree_names = pc_trees, categories = categories,
      tree_id = "tree", impossible = pc_impossible
    )
  }
  expect_error(pc_call(unname(pc_categories)), "list by tree")
  expect_error(
    pc_call(c(pc_categories, list(pairs3 = "E1"))),
    "trees that are not in the model: 'pairs3'"
  )
  expect_error(
    pc_call(pc_categories[-2]),
    "Tree 'single1'.*positions 1, 2.*line i, so the tree needs 2 entries, but categories gives 0"
  )
  expect_error(
    pc_call(modifyList(pc_categories, list(pairs2 = c("E1", "E2", "E3")))),
    "Tree 'pairs2'.*needs 4 entries, but categories gives 3"
  )
  expect_error(
    pc_call(modifyList(pc_categories, list(single2 = c("F1", "F2", "F3")))),
    "single2.*2 lines.*3 entries"
  )
})

test_that("mpt_from_string checks a shared categories vector", {
  two_trees <- c("D + (1 - D) * g", "(1 - D) * (1 - g)", "",
                 "(1 - D) * g", "D + (1 - D) * (1 - g)")
  expect_error(
    mpt_from_string(
      two_trees, tree_names = c("old", "new"), categories = c("yes", "no", "x")
    ),
    "3 entries.*longest tree has 2 lines"
  )
  expect_error(
    mpt_from_string(
      two_trees, tree_names = c("old", "new"),
      categories = c(hit = "yes", miss = "no")
    ),
    "by position.*list by tree"
  )
  expect_error(
    mpt_from_string(two_trees, tree_names = c("old", "new"), categories = "yes"),
    "Tree 'old'.*needs 2 entries, but categories gives 1"
  )
  # a shorter tree takes the leading entries of a shared vector: here the
  # two-option tree lacks the third option
  three_options <- mpt_from_string(
    c("D + (1 - D) / 3", "(1 - D) / 3", "(1 - D) / 3", "",
      "D + (1 - D) / 2", "(1 - D) / 2"),
    tree_names = c("three", "two"), tree_id = "tree",
    categories = c("target", "lure", "other"), impossible = list(two = "other")
  )
  expect_equal(
    names(three_options$other_vars$trees$two$branches), c("target", "lure")
  )

  # trees of different lengths under one vector: the list form comes first
  expect_error(
    mpt_from_string(
      pc_mptinr_model, tree_names = pc_trees, tree_id = "tree",
      categories = c("E1", "E2", "E3", "E4")
    ),
    paste0(
      "single1: 'E3', 'E4'.*\n[^\n]*different numbers of lines.*",
      "categories = list\\(pairs1 = .*\n.*impossible"
    )
  )
})

test_that("mpt_renaming keeps the labels that a by-map reshape needs", {
  # TreeBUGS file: the pair tree prefixes its labels, the singleton tree not
  model <- suppressMessages(mpt_from_eqn(
    write_eqn(c(
      "P  P_E1  c*r",
      "P  P_E2  (1-c)*u*u",
      "P  P_E3  2*(1-c)*u*(1-u)",
      "P  P_E4  c*(1-r)+(1-c)*(1-u)*(1-u)",
      "S  F1    u",
      "S  F2    1-u"
    )),
    tree_id = "tree",
    impossible = list(P = c("F1", "F2"), S = c("E1", "E2", "E3", "E4"))
  ))
  map <- attr(model, "mpt_renaming")
  map <- map[map$kind == "category", ]
  expect_equal(
    map[map$tree == "S", c("file_name", "bmm_name", "renamed")],
    data.frame(file_name = c("F1", "F2"), bmm_name = c("F1", "F2"), renamed = FALSE),
    ignore_attr = "row.names"
  )
  wide <- data.frame(
    P_E1 = 5, P_E2 = 3, P_E3 = 4, P_E4 = 8, F1 = 6, F2 = 4
  )
  long <- do.call(rbind, lapply(split(map, map$tree), function(tree_map) {
    counts <- wide[tree_map$file_name]
    names(counts) <- tree_map$bmm_name
    # categories impossible in this tree get a count of 0
    counts[setdiff(model$resp_vars$resp_cats, names(counts))] <- 0
    data.frame(tree = tree_map$tree[1], counts[model$resp_vars$resp_cats])
  }))
  expect_equal(long$tree, c("P", "S"))
  expect_equal(unname(rowSums(long[-1])), c(20, 10))
})

test_that("mpt_from_eqn rejects categories maps with missing values", {
  eqn_file <- write_eqn(eqn_2htm_lines)
  expect_error(
    mpt_from_eqn(
      eqn_file, categories = c(hit = NA, fa = "yes", miss = "no", cr = "no")
    ),
    "named character vector"
  )
  expect_error(
    mpt_from_eqn(
      eqn_file, categories = setNames(c("yes", "no"), c("hit", NA))
    ),
    "named character vector"
  )
})

test_that("mpt_from_eqn wants covariates spelled as in the file", {
  eqn_file <- write_eqn(c("t  a  g*G_c", "t  b  1-g*G_c"))
  expect_error(
    mpt_from_eqn(eqn_file, covariates = "Gc"),
    "covariates.*not in the EQN file: 'Gc'.*The file has 'G_c'"
  )
  expect_error(
    mpt_from_eqn(eqn_file, covariates = c("G_c", "Z")),
    "not in the EQN file: 'Z'"
  )
  model <- mpt_from_eqn(eqn_file, covariates = "G_c")
  expect_equal(names(model$parameters), "g")
})

test_that("mpt_from_eqn points numeric category labels to categories", {
  numeric_cats <- write_eqn(c("t  1  g", "t  2  1-g"))
  expect_error(
    mpt_from_eqn(numeric_cats),
    "'1' \\(first on line 1\\), '2' \\(first on line 2\\).*categories argument"
  )
  # every line uses the labels, so line 1 is no title
  expect_no_match(
    tryCatch(mpt_from_eqn(numeric_cats), error = conditionMessage), "title"
  )
  model <- suppressMessages(
    mpt_from_eqn(numeric_cats, categories = c("1" = "a", "2" = "b"))
  )
  expect_equal(model$resp_vars$resp_cats, c("a", "b"))
})

test_that("mpt_from_eqn suggests removing a title on the first line", {
  hint <- "If line 1 is a title, delete it or start it with #"
  expect_error(
    mpt_from_eqn(write_eqn(c("2HTM recognition model", eqn_2htm_lines))),
    hint, fixed = TRUE
  )
  expect_error(
    mpt_from_eqn(write_eqn(c("Model for recognition data", eqn_2htm_lines))),
    hint, fixed = TRUE
  )
  expect_error(
    mpt_from_eqn(write_eqn(c("2HTM model", eqn_2htm_lines))),
    hint, fixed = TRUE
  )
  # a title with valid names forms a tree that lacks the other categories
  expect_error(
    suppressMessages(mpt_from_eqn(
      write_eqn(c("Two high threshold", eqn_2htm_lines)),
      categories = eqn_2htm_map, tree_id = "item_type"
    )),
    hint, fixed = TRUE
  )
  # no hint when the problem is on a later line
  later <- tryCatch(
    mpt_from_eqn(write_eqn(c("t  a  g", "t  b  1 - g h"))),
    error = conditionMessage
  )
  expect_no_match(later, "title")
})

test_that("mpt_from_string says how many categories entries a tree needs", {
  expect_error(
    mpt_from_string(
      c("a", "b # y", "1 - a - b"), tree_names = "t1",
      categories = list(t1 = c("x", "z"))
    ),
    "positions 1, 3\\. Entry i of categories labels line i, so the tree needs 3 entries, but categories gives 2"
  )
})

test_that("mpt_from_eqn names the prefix strip in a merge message", {
  merge_file <- write_eqn(c(
    "old  old_yes  D",
    "old  hit      (1-D)*g",
    "old  old_no   (1-D)*(1-g)"
  ))
  suppressMessages(expect_message(
    mpt_from_eqn(merge_file, categories = c(hit = "yes")),
    "tree old: 'old_yes', 'hit' -> yes \\(after removing the prefix 'old_' from 'old_yes'\\)"
  ))
})
