# Tests for conditional_effects.bmmfit() and its internal helpers
# Tier 1: Unit tests (always run, no fitted model)
# Tier 2: Fixture-based integration tests (skip on CRAN)
# Tier 3: Model-fitting integration tests (skip on CRAN)

load_sdm_fit <- function() {
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

# ===========================================================================
# Tier 1: Unit tests — .extract_re_grouping_vars()
# ===========================================================================

test_that(".extract_re_grouping_vars extracts single-bar grouping var", {
  f <- y ~ x + (1 | id)
  expect_equal(.extract_re_grouping_vars(f), "id")
})

test_that(".extract_re_grouping_vars extracts double-bar grouping var", {
  f <- y ~ x + (1 || id)
  expect_equal(.extract_re_grouping_vars(f), "id")
})

test_that(".extract_re_grouping_vars extracts correlation-ID and grouping var", {
  f <- y ~ x + (1 |ID1| id)
  result <- .extract_re_grouping_vars(f)
  expect_true("id" %in% result)
  expect_true("ID1" %in% result)
})

test_that(".extract_re_grouping_vars extracts gr() grouping var", {
  f <- y ~ x + (1 | gr(id, by = exp))
  expect_equal(.extract_re_grouping_vars(f), "id")
})

test_that(".extract_re_grouping_vars extracts gr() with cor arg", {
  f <- y ~ x + (1 | gr(id, cor = FALSE))
  expect_equal(.extract_re_grouping_vars(f), "id")
})

test_that(".extract_re_grouping_vars extracts mm() grouping vars", {
  f <- y ~ x + (1 | mm(g1, g2))
  result <- .extract_re_grouping_vars(f)
  expect_true("g1" %in% result)
  expect_true("g2" %in% result)
  expect_length(result, 2)
})

test_that(".extract_re_grouping_vars extracts crossed grouping vars", {
  f <- y ~ x + (1 | id:group)
  result <- .extract_re_grouping_vars(f)
  expect_true("id" %in% result)
  expect_true("group" %in% result)
})

test_that(".extract_re_grouping_vars handles multiple RE terms", {
  f <- y ~ x + (1 | id) + (1 | group)
  result <- .extract_re_grouping_vars(f)
  expect_true("id" %in% result)
  expect_true("group" %in% result)
})

test_that(".extract_re_grouping_vars returns empty for no RE", {
  f <- y ~ x
  expect_equal(.extract_re_grouping_vars(f), character(0))
})

test_that(".extract_re_grouping_vars returns empty for intercept only", {
  f <- y ~ 1
  expect_equal(.extract_re_grouping_vars(f), character(0))
})


# ===========================================================================
# Tier 1: Unit tests — .ce_summarize_draws()
# ===========================================================================

test_that(".ce_summarize_draws computes mean/SD summary", {
  set.seed(42)
  draws <- matrix(rnorm(1000 * 3), nrow = 1000, ncol = 3)
  result <- .ce_summarize_draws(draws)

  expect_named(result, c("estimate", "lower", "upper", "se"))
  expect_length(result$estimate, 3)
  expect_length(result$lower, 3)
  expect_length(result$upper, 3)
  expect_length(result$se, 3)

  # Estimates should be close to column means
  expect_equal(result$estimate, colMeans(draws), tolerance = 1e-10)
  # SE should be close to column SDs
  expect_equal(result$se, apply(draws, 2, sd), tolerance = 1e-10)
})

test_that(".ce_summarize_draws uses median/MAD when robust = TRUE", {
  set.seed(42)
  draws <- matrix(rnorm(1000 * 2), nrow = 1000, ncol = 2)
  result <- .ce_summarize_draws(draws, robust = TRUE)

  expect_equal(result$estimate, apply(draws, 2, median), tolerance = 1e-10)
  expect_equal(result$se, apply(draws, 2, mad), tolerance = 1e-10)
})

test_that(".ce_summarize_draws handles single-row draws", {
  draws <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
  result <- .ce_summarize_draws(draws)

  expect_equal(result$estimate, c(1, 2, 3))
  expect_length(result$lower, 3)
  expect_length(result$upper, 3)
})

test_that(".ce_summarize_draws prob argument controls CI width", {
  set.seed(42)
  draws <- matrix(rnorm(5000 * 2), nrow = 5000, ncol = 2)

  wide <- .ce_summarize_draws(draws, prob = 0.95)
  narrow <- .ce_summarize_draws(draws, prob = 0.50)

  # Wider prob → wider interval
  expect_true(all(wide$upper - wide$lower > narrow$upper - narrow$lower))
})


# ===========================================================================
# Tier 1: Unit tests — .apply_link_transform()
# ===========================================================================

# Helper to create mock brms_conditional_effects objects
mock_ce <- function(...) {
  dfs <- list(...)
  class(dfs) <- c("brms_conditional_effects", "list")
  dfs
}

mock_ce_df <- function(estimate, lower, upper) {
  data.frame(
    x = seq_along(estimate),
    estimate__ = estimate,
    lower__ = lower,
    upper__ = upper
  )
}

test_that(".apply_link_transform is no-op for identity link", {
  ce <- mock_ce(
    eff1 = mock_ce_df(c(1, 2, 3), c(0.5, 1.5, 2.5), c(1.5, 2.5, 3.5))
  )
  result <- .apply_link_transform(ce, "identity", inverse = TRUE)
  expect_equal(result[[1]]$estimate__, c(1, 2, 3))
  expect_equal(result[[1]]$lower__, c(0.5, 1.5, 2.5))
  expect_equal(result[[1]]$upper__, c(1.5, 2.5, 3.5))
})

test_that(".apply_link_transform applies inverse log (exp)", {
  ce <- mock_ce(
    eff1 = mock_ce_df(c(0, 1, 2), c(-0.5, 0.5, 1.5), c(0.5, 1.5, 2.5))
  )
  result <- .apply_link_transform(ce, "log", inverse = TRUE)
  expect_equal(result[[1]]$estimate__, exp(c(0, 1, 2)), tolerance = 1e-10)
  expect_equal(result[[1]]$lower__, exp(c(-0.5, 0.5, 1.5)), tolerance = 1e-10)
  expect_equal(result[[1]]$upper__, exp(c(0.5, 1.5, 2.5)), tolerance = 1e-10)
})

test_that(".apply_link_transform applies forward log", {
  ce <- mock_ce(
    eff1 = mock_ce_df(c(1, 2, 3), c(0.5, 1.5, 2.5), c(1.5, 2.5, 3.5))
  )
  result <- .apply_link_transform(ce, "log", inverse = FALSE)
  expect_equal(result[[1]]$estimate__, log(c(1, 2, 3)), tolerance = 1e-10)
  expect_equal(result[[1]]$lower__, log(c(0.5, 1.5, 2.5)), tolerance = 1e-10)
  expect_equal(result[[1]]$upper__, log(c(1.5, 2.5, 3.5)), tolerance = 1e-10)
})

test_that(".apply_link_transform applies inverse logit (plogis)", {
  ce <- mock_ce(
    eff1 = mock_ce_df(c(-1, 0, 1), c(-2, -1, 0), c(0, 1, 2))
  )
  result <- .apply_link_transform(ce, "logit", inverse = TRUE)
  expect_equal(result[[1]]$estimate__, plogis(c(-1, 0, 1)), tolerance = 1e-10)
  expect_equal(result[[1]]$lower__, plogis(c(-2, -1, 0)), tolerance = 1e-10)
})

test_that(".apply_link_transform preserves class and names", {
  ce <- mock_ce(
    set_size = mock_ce_df(c(1, 2), c(0.5, 1.5), c(1.5, 2.5)),
    condition = mock_ce_df(c(3, 4), c(2.5, 3.5), c(3.5, 4.5))
  )
  result <- .apply_link_transform(ce, "log", inverse = TRUE)
  expect_s3_class(result, "brms_conditional_effects")
  expect_named(result, c("set_size", "condition"))
})

test_that(".apply_link_transform transforms all elements in list", {
  ce <- mock_ce(
    eff1 = mock_ce_df(c(0, 1), c(-0.5, 0.5), c(0.5, 1.5)),
    eff2 = mock_ce_df(c(2, 3), c(1.5, 2.5), c(2.5, 3.5))
  )
  result <- .apply_link_transform(ce, "log", inverse = TRUE)
  expect_equal(result[[1]]$estimate__, exp(c(0, 1)), tolerance = 1e-10)
  expect_equal(result[[2]]$estimate__, exp(c(2, 3)), tolerance = 1e-10)
})


# ===========================================================================
# Tier 1: Unit tests — .has_category_dpars()
# ===========================================================================

test_that(".has_category_dpars is TRUE for the multinomial models only (#512)", {
  mock <- function(formula, data, model) {
    bmm(formula, data, model, backend = "mock", mock_fit = 1, rename = FALSE)
  }
  rating <- cbind(data.frame(stimulus = c(0L, 1L)),
                  rsdt_rating(2, 50, c(0L, 1L), d = 1, thresholds = c(-0.5, 0, 0.5)))
  ranking <- as.data.frame(rsdt_ranking(2, 50, m = 4, d = 1))
  cdp <- cbind(data.frame(stimulus = c(0L, 1L)),
               rsdt_cdp(2, 50, c(0L, 1L), dfam = 0.8, drec = 1, rcrit = 0.5,
                        thresholds = .cdp_make_thresholds(0, -0.3, 3, 3, "parsimonious"),
                        n_new = 3))

  expect_true(.has_category_dpars(mock(
    bmf(c ~ 1, a ~ 1), oberauer_lewandowsky_2019_e1,
    m3(resp_cats = c("corr", "other", "npl"),
       num_options = c("n_corr", "n_other", "n_npl"),
       choice_rule = "simple", version = "ss")
  )))
  expect_true(.has_category_dpars(mock(
    bmf(d ~ 1, criterion ~ 1, spacing ~ 1), rating,
    sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  )))
  expect_true(.has_category_dpars(mock(
    bmf(d ~ 1), ranking, sdt_ranking(paste0("rank", 1:4), m = 4)
  )))
  expect_true(.has_category_dpars(mock(
    bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, spacing ~ 1, rcrit ~ 1), cdp,
    sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  )))
  expect_false(.has_category_dpars(mock(
    bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm(resp_error = "dev_rad")
  )))
})


# ===========================================================================
# Tier 1: Unit tests — .filter_internal_effects()
# ===========================================================================

test_that(".filter_internal_effects removes internal variables", {
  # Build a mock bmmfit with minimal structure
  mock_bmmfit <- list(
    bmm = list(
      model = structure(
        list(other_vars = list()),
        class = c("sdm", "bmmodel")
      )
    )
  )

  ce <- mock_ce(
    set_size = mock_ce_df(1:3, 0:2, 2:4),
    LureIdx1 = mock_ce_df(1:3, 0:2, 2:4),
    Idx_corr = mock_ce_df(1:3, 0:2, 2:4),
    inv_ss = mock_ce_df(1:3, 0:2, 2:4),
    Item1_Col_rad = mock_ce_df(1:3, 0:2, 2:4),
    expS = mock_ce_df(1:3, 0:2, 2:4)
  )

  result <- .filter_internal_effects(ce, mock_bmmfit)
  expect_named(result, "set_size")
  expect_s3_class(result, "brms_conditional_effects")
})

test_that(".filter_internal_effects keeps all user vars", {
  mock_bmmfit <- list(
    bmm = list(
      model = structure(
        list(other_vars = list()),
        class = c("sdm", "bmmodel")
      )
    )
  )

  ce <- mock_ce(
    set_size = mock_ce_df(1:3, 0:2, 2:4),
    condition = mock_ce_df(1:3, 0:2, 2:4)
  )

  result <- .filter_internal_effects(ce, mock_bmmfit)
  expect_named(result, c("set_size", "condition"))
})


# ===========================================================================
# Tier 2: Fixture-based integration tests
# ===========================================================================

test_that("conditional_effects returns correct class for par = 'c'", {
  skip_on_cran()
  fit <- load_sdm_fit()

  ce <- conditional_effects(fit, par = "c")
  expect_s3_class(ce, "brms_conditional_effects")
  expect_true(length(ce) > 0)
})

test_that("conditional_effects works for intercept-only par = 'kappa'", {
  skip_on_cran()
  fit <- load_sdm_fit()

  ce <- conditional_effects(fit, par = "kappa")
  expect_s3_class(ce, "brms_conditional_effects")
})

test_that("conditional_effects with par = NULL returns all estimated params", {
  skip_on_cran()
  fit <- load_sdm_fit()

  ce <- conditional_effects(fit)
  expect_s3_class(ce, "brms_conditional_effects")
  # SDM fixture has estimated params: c and kappa
  # Effect names are prefixed with par name: "c.set_size", "kappa.1"
  effect_names <- names(ce)
  expect_true(any(grepl("^c\\.", effect_names)))
  expect_true(any(grepl("^kappa\\.", effect_names)))
})

test_that("conditional_effects errors for invalid par name", {
  skip_on_cran()
  fit <- load_sdm_fit()

  expect_error(
    conditional_effects(fit, par = "nonexistent"),
    "not found in model"
  )
})

test_that("conditional_effects errors for non-character par", {
  skip_on_cran()
  fit <- load_sdm_fit()

  expect_error(
    conditional_effects(fit, par = 42),
    "must be a single character string"
  )
})

test_that("scale = 'native' gives positive values for log-linked par", {
  skip_on_cran()
  fit <- load_sdm_fit()

  ce <- conditional_effects(fit, par = "c", scale = "native")
  # c has log link, so native scale = exp(sampling) → all positive
  estimates <- ce[[1]]$estimate__
  expect_true(all(estimates > 0))
})

test_that("scale = 'sampling' can give negative values for log-linked par", {
  skip_on_cran()
  fit <- load_sdm_fit()

  ce <- conditional_effects(fit, par = "c", scale = "sampling")
  # On log scale, values can be any real number
  # Just verify it returns successfully and has different values from native
  ce_native <- conditional_effects(fit, par = "c", scale = "native")
  expect_false(
    isTRUE(all.equal(ce[[1]]$estimate__, ce_native[[1]]$estimate__))
  )
})

test_that("effects argument limits output to specified effect", {
  skip_on_cran()
  fit <- load_sdm_fit()

  ce <- conditional_effects(fit, par = "c", effects = "set_size")
  expect_length(ce, 1)
  expect_true("set_size" %in% names(ce))
})

test_that("plotting conditional_effects works", {
  skip_on_cran()
  fit <- load_sdm_fit()

  ce <- conditional_effects(fit, par = "c")
  p <- plot(ce, plot = FALSE)
  expect_true(length(p) > 0)
})


# ===========================================================================
# Tier 2: Multinomial families (m3 fixture; rows differ in nTrials)
# ===========================================================================

load_m3_fit <- function() {
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "m3 fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

test_that("m3 conditional_effects is the parameter evaluated on the brms grid (#510, #512)", {
  skip_on_cran()
  fit <- load_m3_fit()
  expect_gt(length(unique(fit$data$nTrials)), 1)

  ce <- conditional_effects(fit, par = "c", scale = "sampling")
  expect_named(ce, "cond")
  draws <- brms::posterior_linpred(fit, newdata = ce$cond, nlpar = "c",
                                   re_formula = NA)
  expect_equal(ce$cond$estimate__, apply(draws, 2, median))
  expect_equal(attr(ce$cond, "response"), "c")
})

test_that("m3 conditional_effects applies conditions and re_formula = NULL (#512)", {
  skip_on_cran()
  fit <- load_m3_fit()

  ce <- conditional_effects(fit, par = "c", scale = "sampling",
                            conditions = data.frame(ID = c(1, 2)),
                            re_formula = NULL)
  expect_setequal(ce$cond$ID, c(1, 2))
  draws <- brms::posterior_linpred(fit, newdata = ce$cond, nlpar = "c",
                                   re_formula = NULL)
  expect_equal(ce$cond$estimate__, apply(draws, 2, median))
  by_id <- split(ce$cond$estimate__, ce$cond$ID)
  expect_false(isTRUE(all.equal(by_id[[1]], by_id[[2]])))
})

test_that("m3 conditional_effects keeps only effects of the parameter's own predictors", {
  skip_on_cran()
  fit <- load_m3_fit()

  expect_setequal(names(conditional_effects(fit)), c("a.cond", "c.cond"))
  expect_length(conditional_effects(fit, par = "d"), 0)
})


# ===========================================================================
# Tier 3: Model-fitting integration tests
# ===========================================================================

test_that("conditional_effects works on sdt_rating fits (#512)", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)))
  withr::local_seed(512)

  dat <- expand.grid(id = 1:6, stimulus = c(0L, 1L), cond = c("A", "B"))
  dat <- cbind(dat, rsdt_rating(nrow(dat), 150, dat$stimulus,
                                d = ifelse(dat$cond == "A", 1, 2),
                                thresholds = c(-0.5, 0, 0.5)))
  dat$x <- stats::rnorm(nrow(dat))
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  fit <- bmm(bmf(d ~ 1 + cond + x, criterion ~ 1, spacing ~ 1), dat, model,
             backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  expect_setequal(names(suppressMessages(conditional_effects(fit))),
                  c("d.cond", "d.x"))

  ce <- suppressMessages(conditional_effects(
    fit, par = "d", effects = "x:cond", int_conditions = list(cond = "B")
  ))
  expect_true(all(ce[[1]]$cond == "B"))

  ce <- suppressMessages(conditional_effects(
    fit, par = "d", effects = "x", spaghetti = TRUE, draw_ids = 1:10
  ))
  draws <- brms::posterior_linpred(fit, newdata = ce$x, nlpar = "d",
                                   re_formula = NA, draw_ids = 1:10)
  expect_equal(ce$x$estimate__, apply(draws, 2, median))
  expect_equal(attr(ce$x, "spaghetti")$estimate__, as.numeric(t(draws)))
})
