# a race model as far as the shared helpers read it; the class is a stand-in
# for lba, lnr or rdm
race_stub <- function(version = "simple", n_choices = 3L, resp_cats = NULL,
                      accumulators = NULL) {
  structure(
    list(
      resp_vars = list(rt = "rt", response = "resp"),
      other_vars = list(n_choices = n_choices, resp_cats = resp_cats,
                        accumulators = accumulators),
      version = version
    ),
    class = c("bmmodel", "racestub", paste0("racestub_", version))
  )
}

# =============================================================================
# EXPECTED DECISION TIME
# =============================================================================

test_that("race_expected_time() integrates an exponential race to its mean", {
  # n_j accumulators of rate r_j finish first at rate sum(n_j r_j)
  total <- sum(c(2, 0.5, 4) * c(1, 3, 2))
  log_surv <- function(t) matrix(-total * t, nrow = 1)
  t_lo <- stats::qexp(1e-12, total)
  t_hi <- stats::qexp(1e-12, total, lower.tail = FALSE)
  expect_lt(abs(race_expected_time(log_surv, t_lo, t_hi) * total - 1), 1e-7)
})

test_that("race_expected_time() keeps the draws apart and returns the same twice", {
  meanlog <- c(-1, 0, 0.5)
  sdlog <- c(0.4, 0.8, 1.2)
  log_surv <- function(t) {
    matrix(stats::plnorm(rep(t, each = 3), meanlog, sdlog, lower.tail = FALSE, log.p = TRUE),
           nrow = 3)
  }
  t_hi <- max(stats::qlnorm(1e-12, meanlog, sdlog, lower.tail = FALSE))
  t_lo <- max(min(stats::qlnorm(1e-12, meanlog, sdlog)), t_hi * 1e-12)
  first <- race_expected_time(log_surv, t_lo, t_hi)
  expect_identical(first, race_expected_time(log_surv, t_lo, t_hi))
  expect_lt(max(abs(first / exp(meanlog + sdlog^2 / 2) - 1)), 1e-6)
})

test_that("race_time_range() brackets power-law tails for race_expected_time()", {
  # S(t) = (1 + t)^-k has mean 1 / (k - 1) and a tail that decays like a power
  # of t, so the range has to reach far beyond the mean
  k <- c(3, 4)
  log_surv <- function(t) matrix(-k * rep(log1p(t), each = 2), nrow = 2)
  range <- race_time_range(log_surv)
  expect_true(all(exp(log_surv(range[2])) * range[2] <= 1e-12))
  expect_true(all(-expm1(log_surv(range[1])) <= 1e-12) || range[1] <= range[2] * 1e-12)
  expect_equal(race_expected_time(log_surv, range[1], range[2]), 1 / (k - 1), tolerance = 1e-6)
})

test_that("race_time_range() stops at t_max for a survivor without a finite mean", {
  log_surv <- function(t) matrix(-log1p(t), nrow = 1)
  expect_equal(race_time_range(log_surv, t_max = 1e5)[2], 1e5)
  nan_surv <- function(t) matrix(NaN, nrow = 1)
  expect_equal(race_time_range(nan_surv, t_max = 1e3)[2], 1e3)
})

# =============================================================================
# LOG-SPACE HELPERS
# =============================================================================

test_that("log_sum_exp() keeps infinities and NA", {
  expect_equal(log_sum_exp(log(2), log(3)), log(5))
  expect_equal(log_sum_exp(-Inf, -Inf), -Inf)
  expect_equal(log_sum_exp(Inf, Inf), Inf)
  expect_equal(log_sum_exp(-Inf, 1), 1)
  expect_true(is.na(log_sum_exp(NA_real_, 1)))
})

test_that("log_diff_exp_floored() floors a difference rounded away and keeps NA", {
  expect_equal(log_diff_exp_floored(log(5), log(3)), log(2))
  expect_equal(log_diff_exp_floored(c(0, -1), c(0, 0)), c(-Inf, -Inf))
  expect_equal(log_diff_exp_floored(0, 0, floor = log(1e-300)), log(1e-300))
  expect_equal(log_diff_exp_floored(c(NA, 1, -Inf), c(0, NA, -Inf)), c(NA, NA, -Inf))
  # where the naive form loses everything, the floored one does not
  expect_equal(log_diff_exp_floored(-1e-20, -2e-20), log(1e-20), tolerance = 1e-12)
})

test_that("log_Phi_diff() is accurate in both tails and across zero", {
  lo <- c(-2, -0.5, 0.3)
  hi <- c(-1, 0.7, 2)
  expect_equal(log_Phi_diff(lo, hi), log(stats::pnorm(hi) - stats::pnorm(lo)), tolerance = 1e-12)
  # far tails, where pnorm() itself rounds to 0 or 1
  expect_equal(log_Phi_diff(-Inf, -40), stats::pnorm(-40, log.p = TRUE))
  expect_equal(log_Phi_diff(40, Inf), stats::pnorm(40, lower.tail = FALSE, log.p = TRUE))
  expect_equal(log_Phi_diff(-Inf, Inf), 0)
  expect_equal(log_Phi_diff(c(1, 2), c(1, 1)), c(-Inf, -Inf))
  expect_equal(log_Phi_diff(c(NA, 0), c(1, NA)), c(NA_real_, NA_real_))
})

test_that("recycle_args() recycles to the longest argument and keeps the names", {
  expect_equal(recycle_args(a = 1, b = 1:3), list(a = c(1, 1, 1), b = 1:3))
})

# =============================================================================
# CATEGORY NAMES AND DATA CODING
# =============================================================================

test_that("race_category_names() returns the categories and refuses unusable names", {
  shared <- c("ndt", "s")
  expect_equal(race_category_names(bmf(left ~ 1, right ~ 1, ndt ~ 1), shared), c("left", "right"))
  expect_error(race_category_names(bmf(ndt ~ 1), shared), "at least one accumulator")
  expect_error(race_category_names(bmf(Mu ~ 1), shared), "reserved internal")
  expect_error(race_category_names(bmf(Intercept ~ 1), shared), "reserved internal")
  expect_error(race_category_names(bmf(data ~ 1), shared), "Stan reserved")
  expect_error(race_category_names(bmf(Y ~ 1), shared), "Stan reserved")
  # Stan is case-sensitive, so only the exact keyword is refused
  expect_equal(race_category_names(bmf(Data ~ 1), shared), "Data")
  expect_error(race_category_names(bmf(opt2 ~ 1), shared), "end in a number")
  expect_error(race_category_names(bmf(opt_a ~ 1), shared), "underscores")
})

test_that("race_check_data() refuses response times that cannot exceed ndt", {
  model <- race_stub()
  dat <- data.frame(rt = c(0.5, 0.7), resp = c(1, 2))
  expect_silent(race_check_data(model, dat))
  expect_error(race_check_data(model, transform(dat, rt = c(0, 0.7))), "zero or negative")
  expect_error(race_check_data(model, transform(dat, rt = c(NA, 0.7))), "NA values")
  expect_error(race_check_data(model, transform(dat, rt = c("a", "b"))), "double or integer")
  expect_error(race_check_data(model, transform(dat, resp = c(1, NA))), "NA values")
  expect_error(race_check_data(model, dat[, "rt", drop = FALSE]), "not present")
  expect_warning(race_check_data(model, transform(dat, rt = c(0.5, 12))), "larger than 10")
})

test_that("race_code_simple_response() codes correct against any error", {
  model <- race_stub(n_choices = 4L)
  dat <- data.frame(rt = 0.5, resp = factor(c("1", "3", "4", "1")))
  coded <- race_code_simple_response(model, dat, "racestub")
  expect_identical(coded$resp, c(1L, 3L, 4L, 1L))
  expect_identical(coded$.racestub_cat, c(1L, 2L, 2L, 1L))
  expect_identical(coded$.racestub_n1, rep(1L, 4))
  expect_identical(coded$.racestub_n2, rep(3L, 4))
  expect_error(race_code_simple_response(model, transform(dat, resp = "a"), "racestub"),
               "non-numeric label")
  expect_error(race_code_simple_response(model, transform(dat, resp = 5), "racestub"),
               "integers in 1:4")
})

test_that("race_code_custom_response() codes the categories and their accumulators", {
  dat <- data.frame(rt = 0.5, resp = c("left", "right", "left"),
                    n_left = c(1, 0, 2), n_right = c(2, 1, 0))
  model <- race_stub("custom", NULL, c("left", "right"))
  coded <- race_code_custom_response(model, dat, "racestub", c("ndt", "s"))
  expect_identical(coded$.racestub_cat, c(1L, 2L, 1L))
  expect_identical(coded$.racestub_n1, rep(1L, 3))

  model$other_vars$accumulators <- c(right = 3, left = 2)
  coded <- race_code_custom_response(model, dat, "racestub", c("ndt", "s"))
  expect_identical(coded$.racestub_n1, rep(2L, 3))
  expect_identical(coded$.racestub_n2, rep(3L, 3))

  # a category may sit a trial out, but not one it wins
  model$other_vars$accumulators <- c(left = "n_left", right = "n_right")
  coded <- race_code_custom_response(model, dat, "racestub", c("ndt", "s"))
  expect_identical(coded$.racestub_n2, c(2L, 1L, 0L))
  expect_error(
    race_code_custom_response(model, transform(dat, n_right = c(2, 0, 0)), "racestub", c("ndt", "s")),
    "Category 'right' wins on trials"
  )

  expect_error(
    race_code_custom_response(model, transform(dat, resp = "up"), "racestub", c("ndt", "s")),
    "not specified in the formula"
  )
  expect_error(
    race_code_custom_response(model, transform(dat, resp = "ndt"), "racestub", c("ndt", "s")),
    "reserved internal"
  )
  expect_warning(
    race_code_custom_response(model, transform(dat, resp = "left", n_left = 1), "racestub", c("ndt", "s")),
    "not present"
  )
})

test_that("race_code_accumulators() refuses counts that do not name every category once", {
  dat <- data.frame(rt = 0.5, resp = c("left", "right"), .racestub_cat = 1:2,
                    n_left = c(1, 1.5))
  model <- race_stub("custom", NULL, c("left", "right"))
  refuses <- function(accumulators, pattern) {
    model$other_vars$accumulators <- accumulators
    expect_error(race_code_accumulators(model, dat, "racestub"), pattern)
  }
  refuses(c(left = 1, right = 1, up = 1), "exactly the formula categories")
  refuses(c(left = 1), "exactly the formula categories")
  refuses(c(1, 2), "uniquely named")
  refuses(c(left = 1, left = 2), "uniquely named")
  refuses(c(left = 1, right = 2.5), "right = 2.5")
  refuses(c(left = 0, right = 1), "left = 0")
  refuses(c(left = "n_left", right = "n_right"), "not found")
  refuses(c(left = "n_left", right = "n_left"), "integers >= 0")
  refuses(list(left = 1, right = 1), "NULL, a named numeric vector")
})

# =============================================================================
# VINT CODING THROUGH BRMS
# =============================================================================

test_that("race_bf() puts the coding columns into vint() in order", {
  formula <- race_bf(race_stub(), "racestub", 3)
  expect_equal(
    deparse1(formula$formula),
    "rt | vint(.racestub_cat, .racestub_n1, .racestub_n2, .racestub_n3) ~ 1"
  )
})

test_that("race_family_vars() slices the vint columns only where brms threads", {
  withr::local_options(brms.threads = NULL)
  expect_equal(race_family_vars(2), c("vint1", "vint2", "vint3"))
  withr::local_options(brms.threads = brms::threading(2))
  expect_equal(race_family_vars(2), paste0("vint", 1:3, "[start:end]"))
  withr::local_options(brms.threads = brms::threading(2, force = TRUE))
  expect_equal(race_family_vars(2), c("vint1", "vint2", "vint3"))
})

test_that("race_counts() reads the counts of one observation", {
  prep <- list(data = list(vint1 = c(1L, 2L), vint2 = c(1L, 3L), vint3 = c(2L, 0L)))
  expect_identical(race_counts(prep, 2, 2), c(3L, 0L))
})

test_that("race_revert_check_data() rebuilds what brms dropped and nothing else", {
  model <- race_stub("custom", NULL, c("left", "right"),
                     accumulators = c(left = "n_left", right = "n_right"))
  frame <- data.frame(rt = c(0.5, 0.6), .racestub_cat = c(2L, 1L),
                      .racestub_n1 = c(1L, 2L), .racestub_n2 = c(3L, 1L))
  reverted <- race_revert_check_data(model, frame, "racestub")
  expect_identical(reverted$resp, c("right", "left"))
  expect_identical(reverted$n_left, c(1L, 2L))
  expect_identical(reverted$n_right, c(3L, 1L))
  expect_setequal(attr(reverted, "rebuilt"), c("resp", "n_left", "n_right"))
  expect_false(any(grepl("^\\.racestub_", colnames(reverted))))

  frame$resp <- c("right", "left")
  reverted <- race_revert_check_data(model, frame, "racestub")
  expect_false("resp" %in% attr(reverted, "rebuilt"))

  simple <- race_revert_check_data(race_stub(), frame[1:4], "racestub")
  expect_identical(simple$resp, c(2L, 1L))
})

test_that("race_pp_observables() labels the category codes", {
  spec <- race_pp_observables(race_stub())
  expect_identical(spec$observed, c(rt = "Y", response = "vint1"))
  expect_match(spec$checks$response$label, "'1 = correct', '2 = error'", fixed = TRUE)
  spec <- race_pp_observables(race_stub("custom", NULL, c("left", "right")))
  expect_match(spec$checks$response$label, "'1 = left', '2 = right'", fixed = TRUE)
})
