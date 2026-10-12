eqn_2htm_model <- function(...) {
  eqn_file <- withr::local_tempfile(fileext = ".eqn", .local_envir = parent.frame())
  writeLines(c(
    "6",
    "old  hit   D_o",
    "old  hit   (1-D_o)*g",
    "old  miss  (1-D_o)*(1-g)",
    "new  fa    (1-D_n)*g",
    "new  cr    D_n",
    "new  cr    (1-D_n)*(1-g)"
  ), eqn_file)
  suppressMessages(mpt_from_eqn(
    eqn_file, categories = c(hit = "yes", fa = "yes", miss = "no", cr = "no"),
    tree_id = "item_type", restrictions = "D_n = D_o", ...
  ))
}

wide_2htm <- function() {
  data.frame(hit = c(30, 25), miss = c(10, 15), fa = c(8, 12), cr = c(32, 28))
}

test_that("mpt_data_reshape gives the same long data for TreeBUGS and MPTinR columns", {
  model <- eqn_2htm_model()
  by_name <- mpt_data_reshape(wide_2htm(), model)
  expect_equal(
    by_name,
    data.frame(
      id = c(1L, 1L, 2L, 2L), item_type = c("old", "new", "old", "new"),
      yes = c(30L, 8L, 25L, 12L), no = c(10L, 32L, 15L, 28L)
    )
  )
  # MPTinR sorts the trees (new, old) and then the labels (cr, fa, hit, miss)
  mptinr <- unname(as.matrix(wide_2htm()[c("cr", "fa", "hit", "miss")]))
  expect_identical(mpt_data_reshape(mptinr, model, columns = "mptinr"), by_name)
  expect_no_error(suppressMessages(
    bmm_data_check(bmf(Do ~ 1, g ~ 1), by_name, model)
  ))
})

test_that("mpt_data_reshape sorts MPTinR columns like check.mpt()", {
  eqn_text <- function(first_label) {
    c(
      "8",
      "10 12 a", "10 2 (1-a)*b", "10 1 (1-a)*(1-b)",
      paste("9", first_label, "a"), "9 4 1-a",
      "2 3 a", "2 11 (1-a)*b", "2 20 (1-a)*(1-b)"
    )
  }
  build_model <- function(first_label) {
    eqn_file <- withr::local_tempfile(fileext = ".eqn", .local_envir = parent.frame())
    writeLines(eqn_text(first_label), eqn_file)
    suppressMessages(mpt_from_eqn(
      eqn_file,
      categories = setNames(
        c("x", "y", "z", "x", "y", "x", "y", "z"),
        c("12", "2", "1", first_label, "4", "3", "11", "20")
      ),
      tree_names = c("10" = "t10", "9" = "t9", "2" = "t2"),
      impossible = list(t9 = "z"), tree_id = "tree"
    ))
  }
  # one distinct count per label, so that a wrong order changes the long data
  counts <- function(labels) {
    as.data.frame(
      as.list(setNames(seq_along(labels) * 10, labels)), check.names = FALSE
    )
  }

  # MPTinR 1.14.1, check.mpt(): numbers sort numerically if all labels of the
  # file are numbers (3 11 20 | 4 5 | 1 2 12) ...
  numeric_model <- build_model("5")
  numeric_order <- c("3", "11", "20", "4", "5", "1", "2", "12")
  expect_identical(
    mpt_data_reshape(
      unname(as.matrix(counts(numeric_order))), numeric_model, columns = "mptinr"
    ),
    mpt_data_reshape(counts(numeric_order), numeric_model)
  )
  # ... and as character strings if one label is not (11 20 3 | 4 a | 1 12 2)
  mixed_model <- build_model("a")
  mixed_order <- c("11", "20", "3", "4", "a", "1", "12", "2")
  expect_identical(
    mpt_data_reshape(
      unname(as.matrix(counts(mixed_order))), mixed_model, columns = "mptinr"
    ),
    mpt_data_reshape(counts(mixed_order), mixed_model)
  )
})

test_that("mpt_data_reshape sets the categories a tree cannot produce to 0", {
  model <- mpt(
    list(
      mpt_tree("pairs", list(E1 = "c * r", E2 = "1 - c * r"), impossible = "F1"),
      mpt_tree("singles", list(F1 = "u", E2 = "1 - u"), impossible = "E1")
    ),
    tree_id = "item_type"
  )
  wide <- data.frame(
    pairs.E1 = c(20, 22), pairs.E2 = c(20, 18),
    singles.F1 = c(25, 24), singles.E2 = c(15, 16)
  )
  long <- mpt_data_reshape(wide, model)
  expect_equal(
    long,
    data.frame(
      id = c(1L, 1L, 2L, 2L), item_type = rep(c("pairs", "singles"), 2),
      E1 = c(20L, 0L, 22L, 0L), E2 = c(20L, 15L, 18L, 16L),
      F1 = c(0L, 25L, 0L, 24L)
    )
  )
  expect_error(
    mpt_data_reshape(unname(as.matrix(wide)), model, columns = "mptinr"),
    "needs a model from mpt_from_eqn"
  )
  expect_no_error(suppressMessages(
    bmm_data_check(bmf(c ~ 1, r ~ 1, u ~ 1), long, model)
  ))
})

test_that("mpt_data_reshape fills the impossible categories of an EQN model with 0", {
  eqn_file <- withr::local_tempfile(fileext = ".eqn")
  writeLines(c(
    "P  P_E1  c*r",
    "P  P_E2  (1-c)*u*u",
    "P  P_E3  2*(1-c)*u*(1-u)",
    "P  P_E4  c*(1-r)+(1-c)*(1-u)*(1-u)",
    "S  F1    u",
    "S  F2    1-u"
  ), eqn_file)
  model <- suppressMessages(mpt_from_eqn(
    eqn_file, tree_id = "tree",
    impossible = list(P = c("F1", "F2"), S = c("E1", "E2", "E3", "E4"))
  ))
  # the labels keep the file's prefix; the long data use the bmm names
  wide <- data.frame(
    P_E1 = c(5, 1), P_E2 = c(3, 2), P_E3 = c(4, 3), P_E4 = c(8, 4),
    F1 = c(6, 5), F2 = c(4, 6)
  )
  long <- mpt_data_reshape(wide, model)
  expect_equal(
    long,
    data.frame(
      id = c(1L, 1L, 2L, 2L), tree = c("P", "S", "P", "S"),
      E1 = c(5L, 0L, 1L, 0L), E2 = c(3L, 0L, 2L, 0L),
      E3 = c(4L, 0L, 3L, 0L), E4 = c(8L, 0L, 4L, 0L),
      F1 = c(0L, 6L, 0L, 5L), F2 = c(0L, 4L, 0L, 6L)
    )
  )
  # MPTinR sorts the trees (P, S) and the labels within a tree
  expect_identical(
    mpt_data_reshape(unname(as.matrix(wide)), model, columns = "mptinr"), long
  )
  expect_no_error(suppressMessages(
    bmm_data_check(bmf(c ~ 1, r ~ 1, u ~ 1), long, model)
  ))
})

test_that("mpt_data_reshape carries non-count columns and uses the id column", {
  model <- eqn_2htm_model()
  wide <- cbind(
    subject = c("s1", "s2"), group = c("a", "b"), age = c(21, 34), wide_2htm()
  )
  long <- mpt_data_reshape(wide, model, id = "subject")
  expect_named(long, c("subject", "item_type", "group", "age", "yes", "no"))
  expect_equal(long$group, c("a", "a", "b", "b"))
  expect_equal(long$age, c(21, 21, 34, 34))
  expect_equal(long$subject, c("s1", "s1", "s2", "s2"))
})

test_that("mpt_data_reshape reads a label shared by trees from the qualified name", {
  eqn_file <- withr::local_tempfile(fileext = ".eqn")
  writeLines(c("4", "old 1 D", "old 2 1-D", "new 1 g", "new 2 1-g"), eqn_file)
  model <- suppressMessages(mpt_from_eqn(
    eqn_file, categories = c("1" = "yes", "2" = "no"), tree_id = "tree"
  ))
  wide <- data.frame(old.1 = 30, old.2 = 10, new.1 = 8, new.2 = 32)
  expect_equal(mpt_data_reshape(wide, model)$yes, c(30L, 8L))
  expect_error(
    mpt_data_reshape(data.frame(`1` = 1, `2` = 2, check.names = FALSE), model),
    "old.1"
  )
})

test_that("mpt_data_reshape stops when one wide column would feed two categories", {
  eqn_file <- withr::local_tempfile(fileext = ".eqn")
  writeLines(c("4", "old hit p", "old miss 1-p", "new old.hit q", "new old.miss 1-q"), eqn_file)
  model <- suppressMessages(mpt_from_eqn(
    eqn_file,
    categories = c(hit = "yes", miss = "no", old.hit = "yes", old.miss = "no"),
    tree_id = "tree"
  ))
  # 'old.hit' is the label of tree new and the qualified name of hit in tree old
  wide <- data.frame(hit = 1, miss = 2, old.hit = 30, old.miss = 40)
  expect_error(mpt_data_reshape(wide, model), "more than one response category.*old.hit")
  qualified <- data.frame(old.hit = 1, old.miss = 2, new.old.hit = 30, new.old.miss = 40)
  expect_equal(mpt_data_reshape(qualified, model)$yes, c(1L, 30L))
})

test_that("mpt_data_reshape sums the labels that categories merged", {
  eqn_file <- withr::local_tempfile(fileext = ".eqn")
  writeLines(c(
    "5", "old hit D", "old miss (1-D)*g", "old miss2 (1-D)*(1-g)",
    "new fa g", "new cr 1-g"
  ), eqn_file)
  model <- suppressMessages(mpt_from_eqn(
    eqn_file,
    categories = c(hit = "yes", fa = "yes", miss = "no", miss2 = "no", cr = "no"),
    tree_id = "tree"
  ))
  wide <- data.frame(hit = 30, miss = 4, miss2 = 6, fa = 8, cr = 32)
  expect_equal(mpt_data_reshape(wide, model)$no, c(10L, 32L))
})

test_that("mpt_data_reshape returns integer counts and keeps the types of carried columns", {
  model <- eqn_2htm_model()
  wide <- data.frame(
    subject = factor(c("s1", "s2")), when = as.Date("2026-01-01") + 0:1,
    flag = c(TRUE, NA), group = factor(c("a", "b"), levels = c("b", "a", "unused")),
    wide_2htm()
  )
  long <- mpt_data_reshape(wide, model, id = "subject")
  expect_type(long$yes, "integer")
  expect_type(long$no, "integer")
  expect_identical(long$subject, factor(c("s1", "s1", "s2", "s2")))
  expect_identical(long$when, rep(wide$when, each = 2))
  expect_identical(long$flag, c(TRUE, TRUE, NA, NA))
  expect_identical(long$group, factor(c("a", "a", "b", "b"), levels = levels(wide$group)))
})

test_that("mpt_data_reshape takes the qualified name where a bare label also exists", {
  model <- eqn_2htm_model()
  wide <- cbind(wide_2htm(), old.hit = c(5, 6))
  long <- mpt_data_reshape(wide, model)
  expect_equal(long$yes[long$item_type == "old"], c(5L, 6L))
})

test_that("mpt_data_reshape without a tree_id column in a single-tree model", {
  model <- mpt(mpt_tree("t", list(a = "p", b = "1 - p")))
  long <- mpt_data_reshape(data.frame(t.a = c(3, 4), t.b = c(7, 6)), model)
  expect_equal(long, data.frame(id = 1:2, a = c(3L, 4L), b = c(7L, 6L)))
})

test_that("mpt_data_reshape reports unusable input", {
  model <- eqn_2htm_model()
  expect_error(mpt_data_reshape(wide_2htm()), "model")
  expect_error(mpt_data_reshape(model = model), "data")
  expect_error(mpt_data_reshape(wide_2htm(), "model"), "mpt model")
  expect_error(mpt_data_reshape(list(1), model), "data frame or a matrix")
  expect_error(mpt_data_reshape(wide_2htm()[0, ], model), "no rows")
  expect_error(mpt_data_reshape(wide_2htm(), model, columns = "other"), "should be one of")

  expect_error(mpt_data_reshape(wide_2htm(), model, id = "subject"), "id column 'subject'")
  expect_error(mpt_data_reshape(wide_2htm(), model, id = c("a", "b")), "id argument")
  repeated <- cbind(subject = c("s1", "s1"), wide_2htm())
  expect_error(mpt_data_reshape(repeated, model, id = "subject"), "repeated values")
})

test_that("mpt_data_reshape refuses duplicated column names", {
  model <- eqn_2htm_model()
  repeated <- cbind(wide_2htm(), hit = c(1, 1))
  expect_error(mpt_data_reshape(repeated, model), "duplicated column names.*hit")
  positional <- unname(as.matrix(wide_2htm()))
  colnames(positional) <- c("a", "a", "b", "b")
  expect_error(
    mpt_data_reshape(positional, model, columns = "mptinr"),
    "duplicated column names.*a.*b"
  )
})

test_that("mpt_data_reshape names every missing count column", {
  model <- eqn_2htm_model()
  err <- expect_error(
    mpt_data_reshape(wide_2htm()[c("hit", "fa")], model), "lack count columns"
  )
  expect_match(conditionMessage(err), "old.miss or miss", fixed = TRUE)
  expect_match(conditionMessage(err), "new.cr or cr", fixed = TRUE)
})

test_that("mpt_data_reshape rejects counts that are not counts", {
  model <- eqn_2htm_model()
  with_hit <- function(value) {
    wide <- wide_2htm()
    wide$hit[1] <- value
    wide
  }
  expect_error(mpt_data_reshape(with_hit(-1), model), "non-negative whole numbers.*hit")
  expect_error(mpt_data_reshape(with_hit(2.5), model), "non-negative whole numbers.*hit")
  expect_error(mpt_data_reshape(with_hit(Inf), model), "non-negative whole numbers.*hit")
  expect_error(mpt_data_reshape(with_hit(NA), model), "missing values.*hit")
  wide <- wide_2htm()
  wide$hit <- as.character(wide$hit)
  expect_error(mpt_data_reshape(wide, model), "must be numeric.*hit")
})

test_that("mpt_data_reshape checks the column count in MPTinR mode", {
  model <- eqn_2htm_model()
  expect_error(
    mpt_data_reshape(matrix(1:6, 2), model, columns = "mptinr"),
    "expects 4 columns.*have 3"
  )
  expect_error(
    mpt_data_reshape(
      cbind(subject = 1:2, group = 1:2, unname(as.matrix(wide_2htm()))),
      model, id = "subject", columns = "mptinr"
    ),
    "expects 4 columns.*have 5"
  )
})

test_that("mpt_data_reshape refuses carried columns that the long data need", {
  model <- eqn_2htm_model()
  expect_error(
    mpt_data_reshape(cbind(id = 1:2, wide_2htm()), model),
    "id = 'name'"
  )
  expect_error(
    mpt_data_reshape(cbind(item_type = "x", wide_2htm()), model),
    "item_type"
  )
  expect_error(
    mpt_data_reshape(cbind(yes = 1, wide_2htm()), model),
    "yes"
  )
  expect_no_error(mpt_data_reshape(cbind(id = 1:2, wide_2htm()), model, id = "id"))
})
