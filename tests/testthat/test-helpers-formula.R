test_that("+.bmmformula method works", {
  withr::local_options("bmm.silent" = 2)
  f1 <- bmf(y ~ 1)
  f2 <- bmf(kappa ~ 1)
  f3 <- bmf(kappa ~ 1, m ~ 1)
  f4 <- bmf(kappa ~ 1, m ~ A + B + (A | ID))
  f5 <- bmf(c ~ set_size)
  f6 <- formula(c ~ 1)
  f7 <- formula(m ~ A + B + (A | ID))

  # adding two bmmformulas works with one formula in each
  expect_equal(f1 + f2, bmf(y ~ 1, kappa ~ 1))

  # adding two bmmformulas works with different number of formulas in each
  expect_equal(f1 + f3, bmf(y ~ 1, kappa ~ 1, m ~ 1))

  # adding two more complex bmmformulas works
  expect_equal(f1 + f4, bmf(y ~ 1, kappa ~ 1, m ~ A + B + (A | ID)))

  # adding three bmmformulas work
  expect_equal(f1 + f2 + f5, bmf(y ~ 1, kappa ~ 1, c ~ set_size))

  # adding a formula to a bmmformula works
  expect_equal(f1 + f6, bmf(y ~ 1, c ~ 1))
  expect_equal(f1 + f7, bmf(y ~ 1, m ~ A + B + (A | ID)))

  # adding a formula to a bmmformula overwrites shared parameters
  suppressMessages(expect_equal(f3 + f7, bmf(kappa ~ 1, m ~ A + B + (A | ID))))

  expect_error(f6 + f1, "The first argument must be a bmmformula.")
  expect_error(f1 + 1, "The second argument must be a formula or a bmmformula.")

  # adding a null formula to a bmmformula returns the same bmmformula (#264)
  base_f <- bmf(y ~ 1)
  expect_equal(base_f, base_f + formula(NULL))
  expect_equal(base_f, bmf(base_f[[1]], as.formula(NULL)))

  # the second formula must have a lhs variable
  expect_error(base_f + bmf(~ 1), "Formulas must have a left-hand-side variable")
  expect_error(base_f + formula(~ 1), "Formulas must have a left-hand-side variable")
})

describe("subsetting a bmmformula with [", {
  f <- bmf(y ~ exp(a) + b, a ~ 1 + (1 | id), b ~ 1, c = 3)

  it("returns a bmmformula", {
    expect_s3_class(f["y"], "bmmformula")
    expect_s3_class(f[c("y", "a")], "bmmformula")
    expect_s3_class(f[c("b", "c")], "bmmformula")
  })

  it("correctly resets the `nl` attribute of each element", {
    attributes_exist <- function(formula) {
      ats <- c("nl")
      all(ats %in% names(attributes(formula)))
    }

    expect_true(all(vapply(f["y"], attributes_exist, logical(1))))
    expect_true(all(vapply(f[c("y", "b")], attributes_exist, logical(1))))

    expect_false(attr(f[c("y", "c")]$y, "nl"))
    expect_true(attr(f[c("y", "a")]$y, "nl"))
  })

  it("returns the same object when subset with []", {
    expect_identical(f, f[])
  })

  it("can reassign elements of a bmmformula", {
    f_new <- f
    f_new["y"] <- y ~ 1
    expect_s3_class(f_new, "bmmformula")
    expect_length(f_new, length(f))
    expect_equal(names(f_new), names(f))
    expect_equal(f_new["y"], bmf(y ~ 1))
    expect_false(is_nl(f_new)["y"])

    f_new <- f
    f_new["y"] <- bmf(y ~ exp(b))
    expect_s3_class(f_new, "bmmformula")
    expect_length(f_new, length(f))
    expect_equal(names(f_new), names(f))
    expect_equal(f_new["y"], bmf(y ~ exp(b)))
    expect_true(is_nl(f_new)["y"])

    f_new <- f
    f_rep <- bmf(y ~ 1, c = 1)
    f_new[names(f_rep)] <- f_rep
    expect_s3_class(f_new, "bmmformula")
    expect_length(f_new, length(f))
    expect_equal(names(f_new), names(f))
    expect_equal(f_new["y"], bmf(y ~ 1))

    f_new <- f
    f_rep <- bmf(y ~ 1, c = 1, new ~ 1 + (1|id))
    f_new[names(f_rep)] <- f_rep
    expect_s3_class(f_new, "bmmformula")
    expect_length(f_new, length(f)+1)
    expect_true(all(names(f) %in% names(f_new)))
    expect_equal(f_new["new"], f_rep["new"])
  })

  it("reconstructs correctly when adding all elements", {
    f <- bmf(y ~ exp(a) + b, a ~ 1 + (1 | id), b ~ 1, c = 3)
    expect_equal(f, f[1]+f[2]+f[3]+f[4])
  })
})

test_that("unrecognized_parameters works", {
  f <- bmf(c ~ 1, a ~ 1, s ~ 1, kappa ~ 1)
  expect_equal(length(unrecognized_parameters(imm(NA, NA, NA, NA), f)), 0)
  expect_equal(unrecognized_parameters(imm(NA, NA, NA, NA, version = "bsc"), f), "a")
  expect_equal(unrecognized_parameters(sdm(NA), f), c("a", "s"))
})

test_that("unrecognized_parameters doesnt reject non-linear transformations", {
  f <- bmf(c ~ 1, a ~ 1, s ~ 1, kappa ~ exp(logkappa), logkappa ~ 1)
  expect_equal(length(unrecognized_parameters(imm(NA, NA, NA, NA), f)), 0)
  expect_equal(unrecognized_parameters(imm(NA, NA, NA, NA, version = "bsc"), f), "a")
  expect_equal(unrecognized_parameters(sdm(NA), f), c("a", "s"))
})

test_that("add_missing_parameters works", {
  f <- bmf(c ~ 1)
  model_pars <- names(imm(NA, NA, NA, NA)$parameters)
  expect_equal(
    names(suppressMessages(
      add_missing_parameters(imm(NA, NA, NA, NA), f)
    )),
    model_pars
  )

  f <- bmf(c ~ 1, s ~ 1, a ~ 1, kappa ~ 1)
  model_pars <- names(imm(NA, NA, NA, NA)$parameters)
  expect_equal(
    names(suppressMessages(add_missing_parameters(
      imm(NA, NA, NA, NA), f
    ))),
    model_pars
  )
})

test_that("check_formula gives expected errors", {
  expect_error(
    check_formula(sdm("dev_rad"),
      data = NULL,
      formula = brmsf <- brms::bf(dev_rad ~ 1, c ~ 1)
    ),
    "The provided formula is not a bmm formula"
  )

  expect_error(
    check_formula(sdm("dev_rad"),
      data = NULL,
      formula = bmf(c ~ 1, kappa1 ~ 1)
    ),
    "Unrecognized model parameters"
  )
})

test_that("check_formula works", {
  withr::local_options(bmm.silent = 2)
  expect_equal(
    names(check_formula(sdm("dev_rad"),
      data = NULL,
      formula = bmf(c ~ 1, kappa ~ 1)
    )),
    c("mu", "c", "kappa")
  )
  expect_equal(
    names(check_formula(sdm("dev_rad"),
      data = NULL,
      formula = bmf(c ~ 1)
    )),
    c("mu", "c", "kappa")
  )
})

test_that("has_intercept works", {
  expect_true(has_intercept(y ~ 1))
  expect_true(has_intercept(y ~ A))
  expect_true(has_intercept(y ~ A + B))
  expect_true(has_intercept(y ~ A * B))
  expect_true(has_intercept(y ~ 1 + A))
  expect_true(has_intercept(y ~ A + (A | ID)))

  expect_false(has_intercept(y ~ 0 + A))
  expect_false(has_intercept(y ~ 0 + A + B))
  expect_false(has_intercept(y ~ 0 + A * B))
  expect_false(has_intercept(y ~ 0 + A + (A | ID)))
})

test_that("rhs_vars works with bmmformulas", {
  f <- bmf(y ~ 1)
  expect_equal(rhs_vars(f), character(0))

  f <- bmf(y ~ a + b, x ~ c)
  expect_equal(rhs_vars(f), c("a", "b", "c"))

  f <- bmf(y ~ a + b + a:b, x ~ c)
  expect_equal(rhs_vars(f), c("a", "b", "c"))

  f <- bmf(y ~ a + b + a:b + (a | d), x ~ c)
  expect_equal(rhs_vars(f), c("a", "b", "d", "c"))

  f <- bmf(y ~ a + b + a:b + (a | d), x ~ c + d)
  expect_equal(rhs_vars(f), c("a", "b", "d", "c"))

  f <- bmf(y ~ a + b + a:b + (a | d), x ~ c + d, d ~ m)
  expect_equal(rhs_vars(f), c("a", "b", "d", "c", "m"))

  # test with non-linear transformations
  f <- bmf(y ~ a + b + a:b + (a | d), x ~ c + d, d ~ exp(m + j))
  expect_equal(rhs_vars(f), c("a", "b", "d", "c", "m", "j"))
})

test_that("lhs_vars works with regular formula", {
  f1 <- formula(y ~ x)
  f2 <- formula(~x)
  expect_equal(lhs_vars(f1), "y")
  expect_equal(lhs_vars(f2), character(0))
})

test_that("lhs_vars works with bmmformulas", {
  f <- bmf(y ~ 1)
  expect_equal(lhs_vars(f), "y")

  f <- bmf(y ~ a + b, x ~ c)
  expect_equal(lhs_vars(f), c("y", "x"))

  f <- bmf(y ~ a + b + a:b, x = 3)
  expect_equal(lhs_vars(f), c("y", "x"))

  f <- bmf(y ~ a + b + a:b + (a | d), x ~ c, d = 3)
  expect_equal(lhs_vars(f), c("y", "x", "d"))

  f <- bmf(y ~ a + b + a:b + (a | d), x ~ c + d, d ~ exp(m + j), m ~ 1)
  expect_equal(lhs_vars(f), c("y", "x", "d", "m"))
})

test_that("lhs_vars works with brmsformulas", {
  bf <- brms::bf
  nlf <- brms::nlf
  lf <- brms::lf

  bf(y ~ 1) |>
    lhs_vars() |>
    expect_equal("mu")
  bf(y ~ a, a ~ 1, nl = TRUE) |>
    lhs_vars() |>
    expect_equal(c("mu", "a"))
  bf(y ~ a, a ~ 1, sigma ~ b, nl = TRUE) |>
    lhs_vars() |>
    expect_equal(c("mu", "sigma", "a"))
})

test_that("assign_nl_attr works", {
  x <- bmf(y ~ c, c ~ a + b, a ~ d, m ~ 1)
  x <- assign_nl_attr(x)
  types <- is_nl(x)
  expect_equal(types, c(y = TRUE, c = TRUE, a = FALSE, m = FALSE))

  x <- bmf(y ~ 1)
  x <- assign_nl_attr(x)
  types <- is_nl(x)
  expect_equal(types, c(y = FALSE))

  f1 <- bmf(y ~ a)
  f2 <- bmf(a ~ 1)
  f3 <- f1 + f2
  f4 <- bmf(y ~ a, a ~ 1)
  expect_equal(f3, f4)
  types3 <- is_nl(f3)
  expect_equal(types3, c(y = TRUE, a = FALSE))
})

test_that("print.bmmformula works", {
  res <- utils::capture.output(bmf(a ~ 1, b = 2))
  expect_equal(res, c("a ~ 1", "b = 2"))
})

test_that("apply_links matches a directly written formula", {
  form <- bmf(x ~ a + c, kappa ~ 1, a ~ 1, c ~ 1)
  links <- list(a = "log", c = "logit")
  reform <- apply_links(form, links)
  expect_equal(reform, bmf(x ~ exp(a) + inv_logit(c), kappa ~ 1, a ~ 1, c ~ 1))
})

test_that("apply_links works with different spacing formula formatting", {
  form <- bmf(x ~ a + car, kappa ~ 1, a ~ 1, car ~ 1)
  links <- list(a = "log", car = "logit")
  reform <- apply_links(form, links)
  expect_equal(reform, bmf(x ~ exp(a) + inv_logit(car), kappa ~ 1, a ~ 1, car ~ 1))
})

test_that("apply_links works with links for multiple predicted parameters", {
  form <- bmf(x ~ a + c, kappa ~ b + d, a ~ 1, c ~ 1, b ~ 1, d ~ 1)
  links <- list(a = "log", c = "identity", d = "probit")
  reform <- apply_links(form, links)
  expect_equal(reform, bmf(x ~ exp(a) + c, kappa ~ b + Phi(d), a ~ 1, c ~ 1, b ~ 1, d ~ 1))
})

test_that("apply_links works when parameter is already part of a transformation", {
  form <- bmf(x ~ log(a) + c^2, kappa ~ 1, a ~ 1, c ~ 1)
  links <- list(a = "probit", c = "log")
  reform <- apply_links(form, links)
  expect_equal(reform, bmf(x ~ log(Phi(a)) + exp(c)^2, kappa ~ 1, a ~ 1, c ~ 1))
})

test_that("apply_links works when parameter appears in to parts of a formula", {
  form <- bmf(x ~ log(a^c) + c^2, kappa ~ 1, a ~ 1, c ~ 1)
  links <- list(a = "probit", c = "log")
  reform <- apply_links(form, links)
  expect_equal(reform, bmf(x ~ log(Phi(a)^exp(c)) + exp(c)^2, kappa ~ 1, a ~ 1, c ~ 1))
})

test_that("apply_links gives error when unknown link type is given", {
  form <- bmf(x ~ log(a^c) + c^2, kappa ~ 1, a ~ 1, c ~ 1)
  links <- list(a = "probit", c = "logggg")
  expect_error(apply_links(form, links), "should be one of")
})

test_that("apply_links works with identity link", {
  form <- bmf(x ~ a + c, kappa ~ 1, a ~ 1, c ~ 1)
  links <- list(a = "identity", c = "log")
  reform <- apply_links(form, links)
  expect_equal(reform, bmf(x ~ a + exp(c), kappa ~ 1, a ~ 1, c ~ 1))
})

test_that("apply_links handles empty links list", {
  form <- bmf(x ~ a + c, kappa ~ 1, a ~ 1, c ~ 1)
  links <- list()
  reform <- apply_links(form, links)
  expect_equal(reform, form)
})

test_that("apply_links handles NULL links input", {
  form <- bmf(x ~ a + c, kappa ~ 1, a ~ 1, c ~ 1)
  links <- NULL
  reform <- apply_links(form, links)
  expect_equal(reform, form)
})

test_that("apply_links is case sensitive for link names", {
  form <- bmf(x ~ a + c, kappa ~ 1, a ~ 1, c ~ 1)
  links <- list(a = "LOG")
  expect_error(apply_links(form, links), "should be one of")
})


