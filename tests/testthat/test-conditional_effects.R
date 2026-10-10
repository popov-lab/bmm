# Tests for conditional_effects.bmmfit() and its internal helpers
# Tier 1: Unit tests (always run, no fitted model; the multinomial-family
#         routes run on mock fits with brms stubbed out)
# Tier 2: Fixture-based integration tests (skip on CRAN; the fixtures are not
#         in the built package, so they skip on CI too)
# Tier 3: Model-fitting integration tests (skip on CRAN and CI)

load_sdm_fit <- function() {
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

load_m3_fit <- function() {
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "m3 fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

mock_bmm <- function(formula, data, model) {
  bmm(formula, data, model, backend = "mock", mock_fit = 1, rename = FALSE)
}

mock_rating_fit <- function(formula = bmf(d ~ 1 + x, criterion ~ 1, spacing ~ 1)) {
  dat <- cbind(
    data.frame(stimulus = c(0L, 1L), x = c(0, 1), id = c(1L, 2L),
               cond = factor(c("A", "B")), cond2 = factor(c("u", "v"))),
    rsdt_rating(2, 50, c(0L, 1L), d = 1, thresholds = c(-0.5, 0, 0.5))
  )
  suppressWarnings(
    mock_bmm(formula, dat, sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus"))
  )
}

mock_sdm_fit <- function() {
  mock_bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm(resp_error = "dev_rad"))
}

mock_mixture3p_fit <- function() {
  dat <- oberauer_lin_2017
  dat$id <- factor(dat$ID)
  mock_bmm(
    bmf(kappa ~ 1, thetat ~ 0 + set_size + (1 | id), thetant ~ 1), dat,
    mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")
  )
}

# What brms returns for the conditions grid of a category dpar: one element per
# element of `effects` (a character vector, or a list of vectors of one or two
# variable names), `n` grid rows per level of a second variable, and the
# attributes brms attaches.
stub_grid <- function(effects, n = 3, spaghetti = FALSE) {
  one <- function(vars) {
    df <- if (length(vars) == 2) {
      expand.grid(seq_len(n), c("A", "B"), stringsAsFactors = FALSE)
    } else {
      data.frame(seq_len(n))
    }
    df <- stats::setNames(df, vars)
    df$cond__ <- factor(1)
    df$effect1__ <- df[[vars[1]]]
    df[c("estimate__", "se__", "lower__", "upper__")] <- 0
    points <- data.frame(effect1__ = 1:6, cond__ = factor(1))
    points$resp__ <- matrix(1:12, 6, 2)
    df <- structure(df, effects = vars, response = "mur1", surface = FALSE,
                    categorical = FALSE, ordinal = FALSE, points = points)
    if (spaghetti) {
      attr(df, "spaghetti") <- data.frame(estimate__ = NA_real_)
    }
    df
  }
  effects <- as.list(effects)
  structure(
    stats::setNames(lapply(effects, one),
                    vapply(effects, paste, "", collapse = ":")),
    class = "brms_conditional_effects"
  )
}

# Replaces the two brms calls of the multinomial-family routes: the conditions
# grid is `grid`, the parameter's draws are `draws` (a matrix, a function of
# the arguments of the call, or by default ones). The arguments of every call
# are recorded in `$grid` and `$linpred` of the returned environment. The fit
# has 100 draws.
stub_category_route <- function(grid, draws = NULL, .env = parent.frame()) {
  seen <- new.env()
  seen$grid <- list()
  seen$linpred <- list()
  local_mocked_bindings(
    .brms_conditional_effects = function(x, ...) {
      seen$grid[[length(seen$grid) + 1]] <- list(...)
      grid
    },
    .env = .env
  )
  local_mocked_bindings(
    posterior_linpred = function(object, ...) {
      args <- list(...)
      seen$linpred[[length(seen$linpred) + 1]] <- args
      if (is.function(draws)) {
        draws(args)
      } else {
        draws %||% matrix(1, 2, nrow(args$newdata))
      }
    },
    .package = "brms",
    .env = .env
  )
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms", .env = .env)
  seen
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

test_that(".extract_re_grouping_vars leaves out the correlation ID", {
  expect_equal(.extract_re_grouping_vars(y ~ x + (1 |ID1| id)), "id")
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 |p| id)), "id")
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 + x |p| id) + (1 |q| g2)), c("id", "g2"))
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 |1| id)), "id")
})

test_that(".extract_re_grouping_vars combines a correlation ID with other syntaxes", {
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 |p| gr(id, by = x))), "id")
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 |p| gr(id, cor = FALSE))), "id")
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 |p| mm(g1, g2))), c("g1", "g2"))
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 |p| g1:g2)), c("g1", "g2"))
})

test_that(".extract_re_grouping_vars takes the first unnamed gr() argument as the group", {
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 | gr(by = grp, ID))), "ID")
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 | gr(cor = FALSE, ID))), "ID")
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 |p| gr(by = grp, ID))), "ID")
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 | gr(by = grp, group = ID))), "ID")
})

test_that(".extract_re_grouping_vars extracts nested grouping vars", {
  expect_equal(.extract_re_grouping_vars(~ 1 + (1 | g1/g2)), c("g1", "g2"))
})

test_that(".extract_re_grouping_vars skips named arguments of mm()", {
  f <- ~ 1 + (1 | mm(g1, g2, weights = cbind(w1, w2), scale = FALSE))
  expect_equal(.extract_re_grouping_vars(f), c("g1", "g2"))
})

test_that(".extract_re_grouping_vars ignores bars outside random-effects terms", {
  expect_equal(.extract_re_grouping_vars(~ 1 + x + (1 | id) + z), "id")
  expect_equal(.extract_re_grouping_vars(y ~ x[, 1] + (1 | id)), "id")
})

test_that(".extract_re_cor_ids returns the IDs that tie random effects together", {
  expect_equal(.extract_re_cor_ids(~ 1 + (1 |p| id)), "p")
  expect_equal(.extract_re_cor_ids(~ 1 + (1 |p| gr(id, cor = FALSE)) + (x |q| id)), c("p", "q"))
  expect_equal(.extract_re_cor_ids(~ 1 + (1 | id) + (1 || g2)), character(0))
  expect_equal(.extract_re_cor_ids(~ 1 + x), character(0))
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
  withr::local_seed(42)
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
  withr::local_seed(42)
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
  withr::local_seed(42)
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

mock_spaghetti <- function(estimate) {
  data.frame(
    x = rep(1:2, times = length(estimate) / 2),
    estimate__ = estimate,
    sample__ = factor(rep(c("1", "2"), each = length(estimate) / 2))
  )
}

test_that(".apply_link_transform puts spaghetti lines on the same scale as the summary (#512)", {
  draws <- c(0.2, 0.3, 0.5, 0.6)
  for (link in c("log", "logit", "tan_half")) {
    for (inverse in c(TRUE, FALSE)) {
      df <- mock_ce_df(c(0.2, 0.4), c(0.1, 0.3), c(0.3, 0.5))
      attr(df, "spaghetti") <- mock_spaghetti(draws)
      result <- .apply_link_transform(mock_ce(eff = df), link, inverse = inverse)[[1]]
      expect_equal(
        attr(result, "spaghetti")$estimate__,
        link_transform(draws, link, inverse = inverse)
      )
      expect_equal(
        result$estimate__,
        link_transform(c(0.2, 0.4), link, inverse = inverse)
      )
    }
  }
})

test_that(".apply_link_transform leaves the other spaghetti columns alone (#512)", {
  df <- mock_ce_df(c(0, 1), c(-0.5, 0.5), c(0.5, 1.5))
  spaghetti <- mock_spaghetti(c(0.2, 0.3, 0.5, 0.6))
  attr(df, "spaghetti") <- spaghetti
  result <- .apply_link_transform(mock_ce(eff = df), "log", inverse = TRUE)
  expect_equal(attr(result[[1]], "spaghetti")[c("x", "sample__")], spaghetti[c("x", "sample__")])
})

test_that(".apply_link_transform transforms the spaghetti of every effect (#512)", {
  one <- mock_ce_df(c(0, 1), c(-0.5, 0.5), c(0.5, 1.5))
  two <- mock_ce_df(c(2, 3), c(1.5, 2.5), c(2.5, 3.5))
  attr(one, "spaghetti") <- mock_spaghetti(c(0, 1, 0.5, 1.5))
  attr(two, "spaghetti") <- mock_spaghetti(c(2, 3, 2.5, 3.5))
  result <- .apply_link_transform(mock_ce(a = one, b = two), "log", inverse = TRUE)
  expect_equal(attr(result$a, "spaghetti")$estimate__, exp(c(0, 1, 0.5, 1.5)))
  expect_equal(attr(result$b, "spaghetti")$estimate__, exp(c(2, 3, 2.5, 3.5)))
})

test_that(".apply_link_transform copes with effects that carry no spaghetti", {
  ce <- mock_ce(eff = mock_ce_df(c(0, 1), c(-0.5, 0.5), c(0.5, 1.5)))
  result <- .apply_link_transform(ce, "log", inverse = TRUE)
  expect_null(attr(result[[1]], "spaghetti"))
})


# ===========================================================================
# Tier 1: Unit tests — .has_category_dpars()
# ===========================================================================

test_that(".has_category_dpars is TRUE for the multinomial models only (#512)", {
  ranking <- as.data.frame(rsdt_ranking(2, 50, m = 4, d = 1))
  cdp <- cbind(data.frame(stimulus = c(0L, 1L)),
               rsdt_cdp(2, 50, c(0L, 1L), dfam = 0.8, drec = 1, rcrit = 0.5,
                        thresholds = .cdp_make_thresholds(0, -0.3, 3, 3, "parsimonious"),
                        n_new = 3))

  expect_true(.has_category_dpars(mock_bmm(
    bmf(c ~ 1, a ~ 1), oberauer_lewandowsky_2019_e1,
    m3(resp_cats = c("corr", "other", "npl"),
       num_options = c("n_corr", "n_other", "n_npl"),
       choice_rule = "simple", version = "ss")
  )))
  expect_true(.has_category_dpars(mock_rating_fit(
    bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  )))
  expect_true(.has_category_dpars(mock_bmm(
    bmf(d ~ 1), ranking, sdt_ranking(paste0("rank", 1:4), m = 4)
  )))
  expect_true(.has_category_dpars(mock_bmm(
    bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, spacing ~ 1, rcrit ~ 1), cdp,
    sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  )))
  expect_false(.has_category_dpars(mock_sdm_fit()))
})

test_that(".has_category_dpars agrees with brms's conv_cats_dpars for every brms family (#512)", {
  skip_on_cran()
  expect_true(exists("conv_cats_dpars", asNamespace("brms")),
              label = "brms still has conv_cats_dpars()")
  family_names <- sub("^[.]family_", "",
                      grep("^[.]family_", ls(asNamespace("brms"), all.names = TRUE),
                           value = TRUE))
  by_brms <- vapply(family_names, function(name) {
    tryCatch(
      brms:::conv_cats_dpars(structure(list(family = name), class = "brmsfamily")),
      error = function(e) NA
    )
  }, logical(1))
  by_brms <- by_brms[!is.na(by_brms)]

  expect_gt(length(by_brms), 40)
  expect_gt(sum(by_brms), 0)
  for (name in names(by_brms)) {
    expect_identical(.has_category_dpars(list(family = list(family = name))),
                     by_brms[[name]], label = name)
  }
})

test_that(".has_category_dpars agrees with brms's conv_cats_dpars for the families bmm uses (#512)", {
  skip_on_cran()
  expect_true(exists("conv_cats_dpars", asNamespace("brms")),
              label = "brms still has conv_cats_dpars()")
  families <- suppressMessages(list(
    brms::mixture(brms::von_mises, brms::von_mises),
    brms::custom_family("cf_int", dpars = "mu", links = "identity", type = "int"),
    brms::custom_family("cf_real", dpars = "mu", links = "identity", type = "real"),
    mock_sdm_fit()$family,
    mock_rating_fit()$family
  ))

  for (family in families) {
    expect_identical(.has_category_dpars(list(family = family)),
                     brms:::conv_cats_dpars(family),
                     label = family$family)
  }
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

test_that(".filter_internal_effects removes the Poss_ columns an MPT fit generates", {
  fit <- bmm(
    bmf(Pm ~ 1, Pb ~ 1),
    mpt_impossible_data(),
    mpt(mpt_impossible_trees(), tree_id = "tree"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  ce <- mock_ce(
    cond = mock_ce_df(1:3, 0:2, 2:4),
    Poss_dist = mock_ce_df(1:3, 0:2, 2:4),
    `cond:Poss_dist` = mock_ce_df(1:3, 0:2, 2:4)
  )

  expect_named(.filter_internal_effects(ce, fit), "cond")
})

test_that(".filter_internal_effects keeps a user column named Poss_*", {
  mock_bmmfit <- list(
    bmm = list(
      model = structure(
        list(other_vars = list()),
        class = c("sdm", "bmmodel")
      )
    )
  )
  mpt_fit <- bmm(
    bmf(Pm ~ 1, Pb ~ 1),
    mpt_impossible_data(),
    mpt(mpt_impossible_trees(), tree_id = "tree"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  ce <- mock_ce(
    set_size = mock_ce_df(1:3, 0:2, 2:4),
    Poss_load = mock_ce_df(1:3, 0:2, 2:4),
    `Poss_load:set_size` = mock_ce_df(1:3, 0:2, 2:4)
  )

  expect_named(
    .filter_internal_effects(ce, mock_bmmfit),
    c("set_size", "Poss_load", "Poss_load:set_size")
  )
  expect_named(
    .filter_internal_effects(ce, mpt_fit),
    c("set_size", "Poss_load", "Poss_load:set_size")
  )
})


# ===========================================================================
# Tier 1: Unit tests — conditional_effects() on multinomial families
# (mock fits, brms stubbed: these run on CI)
# ===========================================================================

test_that("conditional_effects rejects arguments a multinomial-family parameter cannot honour (#512)", {
  local_mocked_bindings(
    .ce_single_parameter = function(...) "reached",
    .ce_all_parameters = function(...) "reached"
  )
  fit <- mock_rating_fit()

  expect_error(conditional_effects(fit, par = "d", method = "posterior_predict"),
               "latent quantities.*'method'")
  expect_error(conditional_effects(fit, par = "d", categorical = TRUE),
               "'categorical'")
  expect_error(conditional_effects(fit, par = "d", ordinal = TRUE),
               "'ordinal'")
  expect_error(conditional_effects(fit, par = "d", select_points = 0.5),
               "'select_points'")
  expect_error(conditional_effects(fit, par = "d", transform = exp),
               "'transform'")
  expect_error(conditional_effects(fit, categorical = TRUE), "'categorical'")
  expect_error(
    conditional_effects(fit, par = "d", categorical = TRUE, transform = exp),
    "'categorical', 'transform'"
  )

  expect_equal(conditional_effects(fit, par = "d", method = "posterior_epred"),
               "reached")
  expect_equal(conditional_effects(fit, par = "d", method = "posterior_linpred",
                                   categorical = FALSE, select_points = 0),
               "reached")
})

test_that("the method is read as brms reads it: aliases of epred and linpred pass, predictive ones do not (#512)", {
  skip_on_cran()
  local_mocked_bindings(.ce_single_parameter = function(...) "reached")
  fit <- mock_rating_fit()
  methods <- c("posterior_epred", "fitted", "pp_expect", "posterior_linpred",
               "posterior_predict", "predict", "pp", "predictive_error", "residuals")

  for (method in methods) {
    by_brms <- brms:::validate_pp_method(method) %in%
      c("posterior_epred", "posterior_linpred")
    result <- tryCatch(conditional_effects(fit, par = "d", method = method),
                       error = function(e) conditionMessage(e))
    expect_identical(identical(result, "reached"), by_brms, label = method)
    if (!by_brms) {
      expect_match(result, "Set 'method' to \"posterior_epred\" or \"posterior_linpred\"",
                   label = method)
      expect_no_match(result, "Remove", label = method)
    }
  }
})

test_that("the error says what to remove and what to set, and names the bmm model (#512)", {
  fit <- mock_rating_fit()

  expect_error(conditional_effects(fit, par = "d", categorical = TRUE),
               "Remove 'categorical'.$")
  expect_error(conditional_effects(fit, par = "d", method = "predict"),
               "points. Set 'method' to .*which agree for a parameter.$")
  expect_error(conditional_effects(fit, par = "d", categorical = TRUE, method = "predict"),
               "Remove 'categorical'. Set 'method' to")
  expect_error(conditional_effects(fit, par = "d", categorical = TRUE),
               "parameters of the 'sdt_rating' model are latent quantities")
})

test_that("conditional_effects rejects the arguments brms accepts by a partial name (#512)", {
  local_mocked_bindings(.ce_single_parameter = function(...) "reached")
  fit <- mock_rating_fit()

  expect_error(conditional_effects(fit, par = "d", categ = TRUE), "'categorical'")
  expect_error(conditional_effects(fit, par = "d", ordi = TRUE), "'ordinal'")
  expect_error(conditional_effects(fit, par = "d", meth = "posterior_predict"),
               "'method'")
  expect_error(conditional_effects(fit, par = "d", select = 0.5), "'select_points'")
  expect_error(conditional_effects(fit, par = "d", trans = exp), "'transform'")
  expect_error(
    conditional_effects(fit, par = "d", method = c("posterior_epred", "posterior_linpred")),
    "'method'"
  )
})

test_that("a prefix of several brms arguments is an error on both routes, as in brms (#512)", {
  stub_category_route(stub_grid("x"))

  for (call in list(
    function(...) conditional_effects(mock_rating_fit(), par = "d", ...),
    function(...) conditional_effects(mock_mixture3p_fit(), par = "thetat", ...)
  )) {
    expect_error(call(re = NULL), "matches multiple formal arguments")
    expect_error(call(pro = 0.5), "matches multiple formal arguments")
    expect_no_error(call(re_formula = NULL, res = 5))
  }
})

test_that("conditional_effects leaves the arguments of other models to brms (#512)", {
  local_mocked_bindings(.ce_single_parameter = function(...) "reached")

  expect_equal(
    conditional_effects(mock_sdm_fit(), par = "c", method = "posterior_predict",
                        categorical = TRUE, transform = exp, select_points = 0.5),
    "reached"
  )
})

test_that(".brms_ce_default reads the defaults of the arguments the routes use (#512)", {
  expect_identical(.brms_ce_default("re_formula"), NA)
  expect_identical(.brms_ce_default("prob"), 0.95)
  expect_identical(.brms_ce_default("robust"), TRUE)
})

test_that(".ce_match_args matches a call as R does: exact, partial and positional (#512)", {
  expect_named(
    .ce_match_args(list(categ = 1, meth = 2, spaghetti = 3, sample_new_levels = 4,
                        re_formula = 5)),
    c("re_formula", "method", "spaghetti", "categorical", "sample_new_levels")
  )
  expect_named(.ce_match_args(list("x", prob = 0.5, "y")), c("effects", "conditions", "prob"))
  expect_equal(.ce_match_args(list("x", "y")), list(effects = "x", conditions = "y"))
  expect_equal(.ce_match_args(list()), list())
})

test_that(".ce_match_args reads the partial names of the draw arguments as brms does (#512)", {
  expect_equal(.ce_match_args(list(draw = 1:2)), list(draw_ids = 1:2))
  expect_equal(.ce_match_args(list(nd = 3)), list(ndraws = 3))
  expect_equal(.ce_match_args(list(nsamp = 3)), list(nsamples = 3))
  expect_equal(.ce_match_args(list(sub = 2)), list(subset = 2))
})

test_that(".ce_match_args refuses a name that fits several formals, in R's words (#512)", {
  for (name in c("re", "r", "pro", "s", "t")) {
    expect_error(.ce_match_args(stats::setNames(list(1), name)),
                 "argument 1 matches multiple formal arguments", label = name)
  }
  expect_named(.ce_match_args(list(prob = 1)), "prob")
  expect_named(.ce_match_args(list(zzz = 1)), "zzz")
})

test_that(".ce_match_args matches the draw arguments only after brms's own formals (#512)", {
  # brms matches its own formals first and the draw arguments where it hands the
  # rest on, so `su` is `surface` (not `subset`) and `n` is `ndraws`
  expect_equal(.ce_match_args(list(su = TRUE)), list(surface = TRUE))
  expect_equal(.ce_match_args(list(n = 3)), list(ndraws = 3))
  expect_equal(.ce_match_args(list(n = 3, su = TRUE)), list(surface = TRUE, ndraws = 3))
  expect_equal(.ce_match_args(list(foo = quote(stop("evaluated!")))),
               list(foo = quote(stop("evaluated!"))))
})

# ---------------------------------------------------------------------------
# .ce_split_args(): which argument reaches which brms call
# ---------------------------------------------------------------------------

test_that(".ce_split_args sends the arguments of conditional_effects to the grid call (#512)", {
  conditions <- data.frame(cond = "A")
  split <- .ce_split_args(
    effects = "x", conditions = conditions, resolution = 7, categ = FALSE,
    transform = NULL, method = "posterior_epred", spaghetti = TRUE, prob = 0.5,
    probs = c(0.1, 0.9), robust = FALSE, sample_new_levels = "gaussian"
  )

  expect_setequal(
    names(split$grid),
    c("effects", "conditions", "resolution", "categorical", "transform", "method",
      "spaghetti", "prob", "probs", "robust", "sample_new_levels", "re_formula",
      "draw_ids")
  )
  expect_equal(split$grid$effects, "x")
  expect_equal(split$grid$conditions, conditions)
  expect_equal(split$grid$probs, c(0.1, 0.9))
  expect_equal(split$grid$draw_ids, 1)
  expect_identical(split$grid$re_formula, NA)

  expect_setequal(names(split$linpred),
                  c("sample_new_levels", "re_formula", "allow_new_levels"))
  expect_equal(split$linpred$sample_new_levels, "gaussian")
  expect_identical(split$linpred$re_formula, NA)
  expect_true(split$linpred$allow_new_levels)

  expect_equal(split$effects, "x")
  expect_equal(split$summary, list(prob = 0.5, probs = c(0.1, 0.9), robust = FALSE))
})

test_that(".ce_split_args gives the summary brms's default prob and leaves the rest unset (#512)", {
  split <- .ce_split_args()
  expect_equal(split$summary, list(prob = 0.95, probs = NULL, robust = NULL))
  expect_null(split$effects)
  expect_length(split$draws, 0)
})

test_that(".ce_split_args keeps NULL as a re_formula (#512)", {
  split <- .ce_split_args(re_form = NULL)
  expect_true("re_formula" %in% names(split$grid))
  expect_true("re_formula" %in% names(split$linpred))
  expect_null(split$grid$re_formula)
  expect_null(split$linpred$re_formula)
})

test_that(".ce_split_args takes the draw arguments out of both calls (#512)", {
  split <- .ce_split_args(ndraws = 7, draw_ids = 2:3, nsamples = 4, subset = 5)

  expect_named(split$draws, c("ndraws", "draw_ids", "nsamples", "subset"))
  expect_equal(split$grid$draw_ids, 1)
  expect_false(any(c("ndraws", "nsamples", "subset") %in% names(split$grid)))
  expect_false(any(c("ndraws", "draw_ids", "nsamples", "subset") %in% names(split$linpred)))
})

test_that(".ce_resolve_draw_ids resolves ndraws into ids and lets draw_ids win (#512)", {
  fit <- mock_rating_fit()
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")

  expect_null(.ce_resolve_draw_ids(fit, list()))
  expect_equal(
    withr::with_seed(512, .ce_resolve_draw_ids(fit, list(ndraws = 7))),
    withr::with_seed(512, sample.int(100L, 7L))
  )
  expect_equal(.ce_resolve_draw_ids(fit, list(ndraws = 7, draw_ids = 4:5)), 4:5)
})

test_that(".ce_resolve_draw_ids lets a deprecated alias win and warns as brms does (#512)", {
  fit <- mock_rating_fit()
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")

  expect_warning(
    ids <- withr::with_seed(512, .ce_resolve_draw_ids(fit, list(ndraws = 3, nsamples = 7))),
    "Argument 'nsamples' is deprecated. Please use argument 'ndraws' instead.",
    fixed = TRUE
  )
  expect_equal(ids, withr::with_seed(512, sample.int(100L, 7L)))

  expect_warning(
    ids <- .ce_resolve_draw_ids(fit, list(draw_ids = 1:2, subset = c(2L, 5L))),
    "Argument 'subset' is deprecated. Please use argument 'draw_ids' instead.",
    fixed = TRUE
  )
  expect_equal(ids, c(2L, 5L))
})

test_that("the deprecation warnings of the draw arguments are brms's (#512)", {
  skip_on_cran()
  warning_of <- function(expr) tryCatch(expr, warning = conditionMessage)
  ndraws <- NULL
  nsamples <- 7
  draw_ids <- NULL
  subset <- 2

  expect_identical(
    warning_of(.ce_resolve_draw_ids(mock_rating_fit(), list(nsamples = 7))),
    warning_of(brms:::use_alias(ndraws, nsamples))
  )
  expect_identical(
    warning_of(.ce_resolve_draw_ids(mock_rating_fit(), list(subset = 2))),
    warning_of(brms:::use_alias(draw_ids, subset))
  )
})

test_that(".validate_draw_ids reads the draws as brms does: same ids, same messages (#512)", {
  skip_on_cran()
  fit <- mock_rating_fit()
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")
  # brms reports a call without draws when it summarises them
  by_brms <- function(...) {
    ids <- brms:::validate_draw_ids(fit, ...)
    if (!is.null(ids)) brms::posterior_summary(matrix(0, length(ids), 1))
    ids
  }
  outcome <- function(fun, args) {
    withr::with_seed(
      1,
      tryCatch(suppressWarnings(do.call(fun, args)), error = conditionMessage)
    )
  }
  cases <- list(
    list(ndraws = 0), list(ndraws = 101), list(ndraws = NA_real_), list(ndraws = c(2, 3)),
    list(ndraws = "5"), list(ndraws = "a"), list(ndraws = TRUE), list(ndraws = 2.5),
    list(ndraws = -1), list(ndraws = 5), list(ndraws = 100),
    list(draw_ids = 0), list(draw_ids = 101), list(draw_ids = integer(0)),
    list(draw_ids = "3"), list(draw_ids = NA), list(draw_ids = c(1, NA)),
    list(draw_ids = c(1, 1)), list(draw_ids = 3:1), list(draw_ids = 2.5),
    list(ndraws = 5, draw_ids = 3:4), list(ndraws = 0, draw_ids = 3:4),
    list(ndraws = 5, draw_ids = 0), list(ndraws = 5, draw_ids = integer(0)),
    list()
  )

  for (args in cases) {
    label <- paste(names(args), vapply(args, deparse1, ""), collapse = ", ")
    expect_identical(
      outcome(.validate_draw_ids, c(list(fit), args)),
      outcome(by_brms, args),
      label = label
    )
  }
})

# ---------------------------------------------------------------------------
# Both routes: the arguments of conditional_effects() reach brms's grid call,
# the draws are resolved after it, and posterior_linpred() sees no argument of
# conditional_effects()
# ---------------------------------------------------------------------------

ce_argument_values <- list(
  effects = "x", conditions = data.frame(x = 1), int_conditions = list(x = 1),
  re_formula = NULL, prob = 0.5, robust = FALSE, method = "posterior_epred",
  spaghetti = FALSE, surface = FALSE, categorical = FALSE, ordinal = FALSE,
  transform = NULL, resolution = 7, select_points = 0, too_far = 0,
  probs = c(0.1, 0.9)
)

test_that("brms's conditional_effects arguments reach the grid call, not posterior_linpred, on both routes (#512)", {
  expect_setequal(names(ce_argument_values),
                  setdiff(names(.brms_ce_formals()), c("x", "...")))
  rating <- mock_rating_fit()
  seen_rating <- stub_category_route(stub_grid("x"))
  do.call(conditional_effects, c(list(rating, par = "d"), ce_argument_values))

  mixture3p <- mock_mixture3p_fit()
  seen_softmax <- stub_category_route(stub_grid("x"))
  do.call(conditional_effects, c(list(mixture3p, par = "thetat"), ce_argument_values))

  for (seen in list(seen_rating, seen_softmax)) {
    expect_gt(length(seen$linpred), 0)
    for (args in seen$linpred) {
      expect_setequal(setdiff(names(args), c("object", "newdata", "nlpar")),
                      c("re_formula", "draw_ids", "allow_new_levels"))
    }
    for (name in names(ce_argument_values)) {
      expect_true(name %in% names(seen$grid[[1]]), label = name)
    }
    expect_equal(seen$grid[[1]]$prob, 0.5)
    expect_equal(seen$grid[[1]]$probs, c(0.1, 0.9))
  }
})

test_that("both routes evaluate every parameter of a call on the draws the call selected (#512)", {
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")
  rating <- mock_rating_fit()
  mixture3p <- mock_mixture3p_fit()

  for (draws in list(list(ndraws = 7), list(nsamples = 7))) {
    withr::local_seed(512)
    seen <- stub_category_route(stub_grid("x"))
    # the second route draws the same ids as the first only by the same seed
    suppressWarnings(do.call(conditional_effects, c(list(mixture3p, par = "thetat"), draws)))
    suppressWarnings(do.call(conditional_effects, c(list(rating, par = "d"), draws)))
    softmax_ids <- lapply(seen$linpred[1:2], `[[`, "draw_ids")
    expect_identical(softmax_ids[[1]], softmax_ids[[2]])
    expect_length(unique(softmax_ids[[1]]), 7)
    expect_true(all(softmax_ids[[1]] %in% 1:100))
    expect_equal(seen$grid[[1]]$draw_ids, 1)
    expect_false(any(c("ndraws", "nsamples") %in% names(seen$linpred[[3]])))
    expect_length(unique(seen$linpred[[3]]$draw_ids), 7)
  }

  seen <- stub_category_route(stub_grid("x"))
  suppressWarnings(conditional_effects(mixture3p, par = "thetat", subset = c(2L, 5L)))
  for (args in seen$linpred) {
    expect_equal(args$draw_ids, c(2L, 5L))
    expect_false("subset" %in% names(args))
  }
})

test_that("both routes refuse invalid draws with brms's message, after the grid call and before any draw (#512)", {
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")
  fits <- list(
    list(fit = mock_rating_fit(), par = "d"),
    list(fit = mock_mixture3p_fit(), par = "thetat")
  )

  for (route in fits) {
    seen <- stub_category_route(stub_grid("x"))
    expect_error(
      conditional_effects(route$fit, par = route$par, ndraws = 0),
      "Argument 'ndraws' should be between 1 and the maximum number of draws (100).",
      fixed = TRUE
    )
    expect_error(
      conditional_effects(route$fit, par = route$par, ndraws = 10^6),
      "Argument 'ndraws' should be between 1 and the maximum number of draws (100).",
      fixed = TRUE
    )
    expect_error(conditional_effects(route$fit, par = route$par, draw_ids = 0),
                 "Some 'draw_ids' indices are out of range.", fixed = TRUE)
    suppressWarnings(expect_error(
      conditional_effects(route$fit, par = route$par, nsamples = 0),
      "Argument 'ndraws' should be between 1"
    ))
    expect_length(seen$grid, 4)
    expect_length(seen$linpred, 0)
  }
})

test_that("the deprecated draw arguments warn on both routes, per parameter (#512)", {
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")
  stub_category_route(stub_grid("x"))
  expect_warning(
    conditional_effects(mock_mixture3p_fit(), par = "thetat", nsamples = 3),
    "Argument 'nsamples' is deprecated. Please use argument 'ndraws' instead.",
    fixed = TRUE
  )
  expect_warning(
    conditional_effects(mock_rating_fit(), par = "d", subset = 1:3),
    "Argument 'subset' is deprecated. Please use argument 'draw_ids' instead.",
    fixed = TRUE
  )
})

test_that("both routes read the arguments after `scale` as brms reads a call: by position and by partial name (#512)", {
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")
  routes <- list(
    list(fit = mock_rating_fit(), par = "d"),
    list(fit = mock_mixture3p_fit(), par = "thetat")
  )

  for (route in routes) {
    seen <- stub_category_route(stub_grid("x"))
    conditional_effects(route$fit, route$par, "native", "x", re_formula = NULL)
    expect_equal(seen$grid[[1]]$effects, "x")
    expect_null(seen$grid[[1]]$re_formula)
    for (args in seen$linpred) {
      expect_false(any(names(args) == ""))
      expect_false("transform" %in% names(args))
    }

    seen <- stub_category_route(stub_grid("x"))
    conditional_effects(route$fit, route$par, "native", "x", data.frame(x = 1))
    expect_equal(seen$grid[[1]]$conditions, data.frame(x = 1))

    seen <- stub_category_route(stub_grid("x"))
    conditional_effects(route$fit, route$par, draw = 4:5)
    for (args in seen$linpred) {
      expect_equal(args$draw_ids, 4:5)
      expect_false("draw" %in% names(args))
    }
  }
})

test_that("both routes refuse an argument that fits several formals in R's words (#512)", {
  stub_category_route(stub_grid("x"))
  expect_error(conditional_effects(mock_rating_fit(), "d", "native", "x", re = NULL),
               "argument 2 matches multiple formal arguments")
  expect_error(conditional_effects(mock_mixture3p_fit(), "thetat", "native", re = NULL),
               "matches multiple formal arguments")
})

test_that("both routes read `su` and `n` as brms does and leave quoted arguments unevaluated (#512)", {
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")
  routes <- list(
    list(fit = mock_rating_fit(), par = "d"),
    list(fit = mock_mixture3p_fit(), par = "thetat")
  )

  for (route in routes) {
    seen <- stub_category_route(stub_grid("x"))
    conditional_effects(route$fit, route$par, effects = "x", su = TRUE, n = 3,
                        foo = quote(stop("evaluated!")))
    expect_true(seen$grid[[1]]$surface)
    expect_identical(seen$grid[[1]]$foo, quote(stop("evaluated!")))
    for (args in seen$linpred) {
      expect_length(args$draw_ids, 3)
      expect_identical(args$foo, quote(stop("evaluated!")))
    }
  }
})

# ---------------------------------------------------------------------------
# Category-family route
# ---------------------------------------------------------------------------

test_that("the category-family route reports the summary of the draws it evaluates (#512)", {
  fit <- mock_rating_fit()
  draws <- rbind(c(1, 2, 3), c(2, 4, 6), c(3, 6, 9), c(10, 20, 30))
  stub_category_route(stub_grid("x"), draws)
  expect_summary <- function(ce, probs, robust) {
    summ <- brms::posterior_summary(draws, probs = probs, robust = robust)
    expect_equal(ce$estimate__, summ[, 1], ignore_attr = TRUE)
    expect_equal(ce$se__, summ[, 2], ignore_attr = TRUE)
    expect_equal(ce$lower__, summ[, 3], ignore_attr = TRUE)
    expect_equal(ce$upper__, summ[, 4], ignore_attr = TRUE)
  }

  ce <- .ce_nlpar_category_family(fit, "d")
  expect_s3_class(ce, "brms_conditional_effects")
  expect_equal(attr(ce$x, "response"), "d")
  expect_summary(ce$x, c(0.025, 0.975), robust = TRUE)
  expect_equal(ce$x$estimate__, apply(draws, 2, median))
  expect_equal(ce$x$se__, apply(draws, 2, mad))

  ce <- .ce_nlpar_category_family(fit, "d", robust = FALSE)
  expect_summary(ce$x, c(0.025, 0.975), robust = FALSE)
  expect_equal(ce$x$estimate__, colMeans(draws))
  expect_equal(ce$x$se__, apply(draws, 2, sd))

  ce <- .ce_nlpar_category_family(fit, "d", prob = 0.5)
  expect_summary(ce$x, c(0.25, 0.75), robust = TRUE)

  ce <- .ce_nlpar_category_family(fit, "d", prob = 0.5, probs = c(0.1, 0.9))
  expect_summary(ce$x, c(0.1, 0.9), robust = TRUE)
})

test_that("the category-family route asks brms for the grid at one draw and evaluates the parameter on it (#512)", {
  fit <- mock_rating_fit()
  grid <- stub_grid("x")
  seen <- stub_category_route(grid)

  .ce_nlpar_category_family(fit, "d")
  expect_named(seen$grid[[1]], c("dpar", "effects", "re_formula", "draw_ids"),
               ignore.order = TRUE)
  expect_equal(seen$grid[[1]]$effects, "x")
  dpars <- names(brms::brmsterms(fit$formula)$dpars)
  expect_gt(length(dpars), 1)
  expect_equal(seen$grid[[1]]$dpar, dpars[1])
  expect_equal(seen$grid[[1]]$draw_ids, 1)
  expect_identical(seen$grid[[1]]$re_formula, NA)
  expect_equal(
    names(seen$linpred[[1]]$newdata),
    setdiff(names(grid$x), c("estimate__", "se__", "lower__", "upper__"))
  )
  expect_equal(seen$linpred[[1]]$nlpar, "d")
  expect_identical(seen$linpred[[1]]$re_formula, NA)
  expect_true(seen$linpred[[1]]$allow_new_levels)
})

test_that("the category-family route hands the arguments of conditional_effects to brms untouched (#512)", {
  fit <- mock_rating_fit()
  seen <- stub_category_route(stub_grid("x"))
  conditions <- data.frame(cond = "A")
  int_conditions <- list(x = c(-1, 1))

  .ce_nlpar_category_family(
    fit, "d", effects = "x", conditions = conditions,
    int_conditions = int_conditions, resolution = 7, too_far = 0.1,
    re_formula = NULL, draw_ids = c(2L, 5L), sample_new_levels = "gaussian"
  )
  grid_args <- seen$grid[[1]]
  expect_equal(grid_args$effects, "x")
  expect_equal(grid_args$conditions, conditions)
  expect_equal(grid_args$int_conditions, int_conditions)
  expect_equal(grid_args$resolution, 7)
  expect_equal(grid_args$too_far, 0.1)
  expect_true("re_formula" %in% names(grid_args))
  expect_null(grid_args$re_formula)
  expect_equal(grid_args$draw_ids, 1)

  linpred_args <- seen$linpred[[1]]
  expect_true("re_formula" %in% names(linpred_args))
  expect_null(linpred_args$re_formula)
  expect_equal(linpred_args$draw_ids, c(2L, 5L))
  expect_equal(linpred_args$sample_new_levels, "gaussian")

  .ce_nlpar_category_family(fit, "d", re_formula = ~ (1 | id))
  expect_equal(seen$grid[[2]]$re_formula, ~ (1 | id))
  expect_equal(seen$linpred[[2]]$re_formula, ~ (1 | id))
  expect_equal(seen$grid[[2]]$draw_ids, 1)
})

test_that("the category-family route returns spaghetti lines of the draws, labelled as brms labels them (#512)", {
  fit <- mock_rating_fit()
  non_summary <- function(grid) {
    setdiff(names(grid), c("estimate__", "se__", "lower__", "upper__"))
  }

  grid <- stub_grid("x", spaghetti = TRUE)
  draws <- rbind(c(1, 2, 3), c(2, 4, 6), c(3, 6, 9), c(10, 20, 30))
  stub_category_route(grid, draws)
  spaghetti <- attr(.ce_nlpar_category_family(fit, "d")$x, "spaghetti")
  expect_named(spaghetti, c(non_summary(grid$x), "estimate__", "sample__"))
  expect_equal(spaghetti$estimate__, as.numeric(t(draws)))
  expect_equal(spaghetti$x, rep(1:3, times = 4))
  expect_equal(spaghetti$sample__, factor(rep(1:4, each = 3)))

  grid <- stub_grid(list(c("x", "cond")), spaghetti = TRUE)
  draws <- matrix(seq_len(24), 4, 6)
  stub_category_route(grid, draws)
  spaghetti <- attr(.ce_nlpar_category_family(fit, "d", effects = "x:cond")[[1]],
                    "spaghetti")
  expect_equal(spaghetti$estimate__, as.numeric(t(draws)))
  expect_equal(spaghetti$cond, rep(rep(c("A", "B"), each = 3), times = 4))
  expect_equal(
    spaghetti$sample__,
    factor(paste0(rep(1:4, each = 6), "_", rep(rep(c("A", "B"), each = 3), 4)))
  )

  stub_category_route(stub_grid("x"))
  expect_null(attr(.ce_nlpar_category_family(fit, "d")$x, "spaghetti"))
})

test_that("brms still says that it sets the trials variables (#512)", {
  skip_on_cran()
  expect_true(
    any(grepl("Setting all 'trials' variables", deparse(brms:::prepare_conditions),
              fixed = TRUE))
  )
})

test_that("brms's 'trials' message is muffled on the category-family route, other conditions are not (#512)", {
  fit <- mock_rating_fit()
  stub_category_route(stub_grid("x"))
  local_mocked_bindings(
    .brms_conditional_effects = function(x, ...) {
      message("Setting all 'trials' variables to 1 by default if not specified otherwise.")
      message("another message")
      warning("another warning")
      stub_grid("x")
    }
  )

  messages <- character()
  warnings <- character()
  withCallingHandlers(
    .ce_nlpar_category_family(fit, "d"),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(messages, "another message\n")
  expect_equal(warnings, "another warning")
})

test_that("the category-family route returns no response points and plots with points = TRUE (#512)", {
  skip_if_not_installed("ggplot2")
  stub_category_route(stub_grid("x"))

  ce <- .ce_nlpar_category_family(mock_rating_fit(), "d")
  expect_equal(nrow(attr(ce$x, "points")), 0)
  expect_no_warning(
    suppressMessages(
      ggplot2::ggplot_build(plot(ce, points = TRUE, plot = FALSE)[[1]])
    )
  )
})

test_that("every element of the grid is evaluated on its own conditions (#512)", {
  stub_category_route(
    stub_grid(list("x", c("x", "cond"))),
    function(args) matrix(seq_len(nrow(args$newdata)), 2, nrow(args$newdata), byrow = TRUE)
  )
  ce <- .ce_nlpar_category_family(mock_rating_fit(), "d", effects = c("x", "x:cond"))

  expect_equal(ce$x$estimate__, 1:3)
  expect_equal(ce[["x:cond"]]$estimate__, 1:6)
})

test_that("the category-family route keeps the effects the user asks for (#512)", {
  stub_category_route(stub_grid("stimulus"))
  ce <- .ce_nlpar_category_family(mock_rating_fit(), "d", effects = "stimulus")
  expect_named(ce, "stimulus")
})

# ---------------------------------------------------------------------------
# Default effects: those brms lists for the parameter's own formula
# ---------------------------------------------------------------------------

effect_labels <- function(effects) {
  sort(vapply(effects, function(vars) paste(sort(vars), collapse = ":"), ""))
}

test_that(".ce_own_effects lists the effects brms lists for the parameter's own formula (#512)", {
  own <- function(formula, par = "d") {
    effect_labels(.ce_own_effects(mock_rating_fit(formula), par))
  }

  expect_equal(own(bmf(d ~ 1 + x + cond, criterion ~ 1 + x * cond2, spacing ~ 1)),
               c("cond", "x"))
  expect_equal(own(bmf(d ~ 1 + x * cond, criterion ~ 1, spacing ~ 1)),
               c("cond", "cond:x", "x"))
  expect_equal(own(bmf(d ~ 1 + x:cond, criterion ~ 1, spacing ~ 1)), "cond:x")
  expect_equal(own(bmf(d ~ 1 + x:cond, criterion ~ 1 + x + cond, spacing ~ 1)),
               "cond:x")
  expect_equal(
    own(bmf(d ~ 1 + log(abs(x) + 1) + I(x^2), criterion ~ 1, spacing ~ 1)),
    "x"
  )
  expect_setequal(
    own(bmf(d ~ 1 + x * cond * cond2, criterion ~ 1, spacing ~ 1)),
    c("x", "cond", "cond2", "cond:x", "cond2:x", "cond:cond2")
  )
  expect_equal(own(bmf(d ~ 1 + (1 | id), criterion ~ 1 + x, spacing ~ 1)),
               character(0))
  expect_equal(own(bmf(d ~ 1, criterion ~ 1 + x, spacing ~ 1)), character(0))
})

test_that(".ce_own_effects follows the formulas of the parameters a non-linear formula names (#512)", {
  own <- function(formula, par = "d") {
    effect_labels(.ce_own_effects(mock_rating_fit(formula), par))
  }

  expect_equal(
    own(bmf(d ~ exp(nlc), nlc ~ 1 + x * cond, criterion ~ 1, spacing ~ 1)),
    c("cond", "cond:x", "x")
  )
  expect_equal(
    own(bmf(d ~ exp(nlc), nlc ~ 1 + x + cond, criterion ~ 1 + x * cond, spacing ~ 1)),
    c("cond", "x")
  )
  expect_equal(
    own(bmf(d ~ nl1 + nl2, nl1 ~ 1 + x, nl2 ~ 1 + cond, criterion ~ 1, spacing ~ 1)),
    c("cond", "x")
  )
  expect_equal(
    own(bmf(d ~ exp(n1), n1 ~ n2 + x, n2 ~ 1 + cond, criterion ~ 1 + x * cond,
            spacing ~ 1)),
    c("cond", "x")
  )
  # a non-linear formula can make its own data variables interact
  expect_setequal(
    own(bmf(d ~ exp(nlc) * x * cond, nlc ~ 1 + cond2, criterion ~ 1, spacing ~ 1)),
    c("cond", "cond:x", "cond2", "x")
  )
})

# 72 rows: enough for brms to build smooth, gp and monotonic terms
mock_rating_terms_fit <- function(formula) {
  withr::local_seed(512)
  withr::local_options(bmm.sort_data = FALSE)
  dat <- expand.grid(id = 1:18, stimulus = c(0L, 1L), cond = factor(c("A", "B")))
  dat <- cbind(
    dat, rsdt_rating(nrow(dat), 50, dat$stimulus, d = 1, thresholds = c(-0.5, 0, 0.5)),
    x = stats::rnorm(nrow(dat)), z = stats::rnorm(nrow(dat)),
    sdx = stats::runif(nrow(dat), 0.1, 0.5),
    ord = factor(rep(1:3, length.out = nrow(dat)), ordered = TRUE)
  )
  suppressWarnings(mock_bmm(
    formula, dat, sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  ))
}

# brms's effects of the parameter's own formula, as it lists them for a
# model with this formula alone
brms_own_effects <- function(fit, par) {
  nlpar <- brms::brmsterms(fit$formula)$nlpars[[par]]
  effects <- brms:::get_all_effects(nlpar)
  unique(effects[lengths(effects) <= 2])
}

test_that(".ce_own_effects lists the effects brms lists for smooth, gp, monotonic and plain terms (#512)", {
  skip_on_cran()
  d_formulas <- list(
    d ~ 1 + x * cond,
    d ~ 0 + cond:x,
    d ~ poly(x, 2):cond,
    d ~ x + I(x^2),
    d ~ 1 + offset(x) + z,
    d ~ x * z * cond,
    d ~ s(x),
    d ~ s(x, by = cond),
    d ~ 1 + cond + s(x, by = cond),
    d ~ s(x, k = 4, bs = "cr"),
    d ~ s(x, z),
    d ~ t2(x, z),
    d ~ gp(x),
    d ~ gp(x, by = cond),
    d ~ gp(x, z),
    d ~ mo(ord),
    d ~ mo(ord) * x,
    d ~ mo(ord):x
  )

  for (formula in d_formulas) {
    fit <- mock_rating_terms_fit(
      bmf(formula, criterion ~ 1 + z + cond, spacing ~ 1)
    )
    expect_equal(
      effect_labels(.ce_own_effects(fit, "d")),
      effect_labels(brms_own_effects(fit, "d")),
      label = deparse1(formula)
    )
  }
})

test_that(".ce_own_effects takes the variables of a term from the arguments that name variables (#512)", {
  own <- function(formula) {
    effect_labels(.ce_own_effects(
      mock_rating_terms_fit(bmf(formula, criterion ~ 1, spacing ~ 1)), "d"
    ))
  }

  expect_equal(own(d ~ s(x, by = cond)), c("cond", "cond:x", "x"))
  expect_equal(own(d ~ gp(x, by = cond)), c("cond", "cond:x", "x"))
  expect_equal(own(d ~ t2(x, z)), c("x", "x:z", "z"))
  expect_equal(own(d ~ mo(ord):x), c("ord", "ord:x", "x"))
  expect_equal(own(d ~ s(x, k = 4, bs = "cr")), "x")
  expect_equal(effect_labels(.ce_term_effects("me(x, sdx)", c("x", "sdx"))), "x")
  expect_equal(effect_labels(.ce_term_effects("mi(x)", "x")), "x")
  expect_equal(
    effect_labels(.ce_term_effects("s(x, k = nk, by = cond)", c("x", "nk", "cond"))),
    c("cond", "cond:x", "x")
  )
  expect_equal(
    effect_labels(.ce_term_effects("gp(x, scale = sc)", c("x", "sc"))),
    "x"
  )
})

# The labels of the effects .ce_select_own_effects() keeps of `brms_effects`,
# in its order
kept_effects <- function(brms_effects, formula, par) {
  vapply(.ce_select_own_effects(brms_effects, .ce_own_effects(mock_rating_fit(formula), par)),
         paste, "", collapse = ":")
}

test_that("default effects of the category-family route are those brms lists for the parameter's own formula (#512)", {
  kept <- function(formula, par = "d") {
    brms_effects <- list("stimulus", "x", "cond", "cond2", c("x", "cond2"),
                         c("x", "cond"))
    kept_effects(brms_effects, formula, par)
  }

  expect_equal(
    kept(bmf(d ~ 1 + x + cond + cond2, criterion ~ 1 + x * cond2, spacing ~ 1)),
    c("x", "cond", "cond2")
  )
  expect_equal(
    kept(bmf(d ~ 1 + x * cond, criterion ~ 1 + x * cond2, spacing ~ 1)),
    c("x", "cond", "x:cond")
  )
  expect_equal(
    kept(bmf(d ~ 1 + x:cond, criterion ~ 1 + x * cond2, spacing ~ 1)),
    "x:cond"
  )
  expect_equal(
    kept(bmf(d ~ 1 + x:cond, criterion ~ 1 + x + cond, spacing ~ 1)),
    "x:cond"
  )
  expect_equal(
    kept(bmf(d ~ exp(nlc), nlc ~ 1 + x * cond, criterion ~ 1 + cond2, spacing ~ 1)),
    c("x", "cond", "x:cond")
  )
  expect_equal(
    kept(bmf(d ~ exp(nlc), nlc ~ 1 + x + cond, criterion ~ 1 + x * cond, spacing ~ 1)),
    c("x", "cond")
  )
  expect_equal(
    kept(bmf(d ~ 1 + (1 | id), criterion ~ 1 + x * cond, spacing ~ 1)),
    character(0)
  )
})

test_that("default effects of the category-family route list an effect once, in any variable order (#512)", {
  brms_effects <- list("stimulus", "x", "cond", c("x", "cond"), c("cond", "x"))
  kept <- function(formula, par = "d") kept_effects(brms_effects, formula, par)
  formula <- bmf(d ~ 1 + x * cond, criterion ~ 1 + cond * x, spacing ~ 1)

  expect_equal(kept(formula, "d"), c("x", "cond", "x:cond"))
  expect_equal(kept(formula, "criterion"), c("x", "cond", "cond:x"))
  expect_equal(
    kept(bmf(d ~ 1 + x + cond, criterion ~ 1 + cond:x, spacing ~ 1), "criterion"),
    "cond:x"
  )
  expect_equal(
    kept(bmf(d ~ 1 + x:cond, criterion ~ 1, spacing ~ 1), "d"),
    "x:cond"
  )
  expect_equal(
    kept(bmf(d ~ exp(nlc), nlc ~ 1 + cond * x, criterion ~ 1, spacing ~ 1), "d"),
    c("x", "cond", "cond:x")
  )
})

test_that("default effects of the category-family route keep the order of brms's grid (#512)", {
  brms_effects <- list("stimulus", "x", "cond", c("x", "cond"), "cond2", c("cond", "x"))
  formula <- bmf(d ~ exp(nlc) * cond * x, nlc ~ 1 + cond2, criterion ~ 1 + x * cond,
                 spacing ~ 1)

  expect_equal(kept_effects(brms_effects, formula, "d"), c("x", "cond", "cond2", "cond:x"))
})

# The effects brms's own conditional_effects() plots for a fit without `effects`,
# with each variable in the place brms gives it: brms runs as it is, only the
# evaluation of each element is replaced by its effects
brms_default_effects <- function(fit, re_formula = NA) {
  local_mocked_bindings(
    contains_draws = function(x) invisible(TRUE),
    conditional_effects.brmsterms = function(x, fit, cond_data, ...) {
      list(attr(cond_data, "effects"))
    },
    .package = "brms"
  )
  suppressMessages(unclass(brms:::conditional_effects.brmsfit(
    fit, re_formula = re_formula, dpar = names(brms::brmsterms(fit$formula)$dpars)[1]
  )))
}

test_that(".ce_brms_default_effects lists the effects brms plots by default, in its order and orientation (#512)", {
  skip_on_cran()
  formulas <- list(
    bmf(d ~ 1 + cond * x, criterion ~ 1 + z:x, spacing ~ 1),
    bmf(d ~ 1 + id:x + cond + (1 | id), criterion ~ 1 + x * z, spacing ~ 1),
    bmf(d ~ exp(nlc) * cond * x, nlc ~ 1 + ord, criterion ~ 1 + cond:z, spacing ~ 1),
    bmf(d ~ 1 + s(x, by = cond) + mo(ord), criterion ~ 1, spacing ~ 1)
  )

  for (formula in formulas) {
    fit <- mock_rating_terms_fit(formula)
    for (re_formula in list(NA, NULL)) {
      expect_equal(.ce_brms_default_effects(fit, re_formula),
                   brms_default_effects(fit, re_formula),
                   label = deparse1(formula$d))
    }
  }
  # a grouping variable counts as a factor once the random effects are kept
  fit <- mock_rating_terms_fit(formulas[[2]])
  expect_contains(.ce_brms_default_effects(fit, NA), list(c("id", "x")))
  expect_contains(.ce_brms_default_effects(fit, NULL), list(c("x", "id")))
})

test_that("without effects, the category-family route asks brms only for the parameter's own effects (#512)", {
  fit <- mock_rating_fit(
    bmf(d ~ 1 + cond * x, criterion ~ 1 + cond2, spacing ~ exp(nls), nls ~ 1 + x)
  )
  seen <- stub_category_route(stub_grid("x"))

  conditional_effects(fit, par = "d")
  conditional_effects(fit, par = "criterion")
  conditional_effects(fit, par = "spacing", re_formula = NULL)
  # brms puts the numeric variable of a factor-by-numeric effect first
  expect_equal(seen$grid[[1]]$effects, c("cond", "x", "x:cond"))
  expect_equal(seen$grid[[2]]$effects, "cond2")
  expect_equal(seen$grid[[3]]$effects, "x")
  expect_null(seen$grid[[3]]$re_formula)

  seen <- stub_category_route(stub_grid("x"))
  conditional_effects(fit, par = "d", effects = "cond:x")
  expect_equal(seen$grid[[1]]$effects, "cond:x")
})

test_that("a parameter without effects of its own asks brms for no grid and returns no effects (#512)", {
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")
  fit <- mock_rating_fit(bmf(d ~ 1 + x, criterion ~ 1 + (1 | id), spacing ~ 1))
  seen <- stub_category_route(stub_grid(list("stimulus", "x")))

  ce <- conditional_effects(fit, par = "criterion")
  expect_identical(
    ce, structure(list(), names = character(0), class = "brms_conditional_effects")
  )
  expect_identical(conditional_effects(fit, par = "spacing", scale = "sampling"), ce)
  expect_length(seen$grid, 0)
  expect_length(seen$linpred, 0)
  expect_error(conditional_effects(fit, par = "criterion", ndraws = 0),
               "Argument 'ndraws' should be between 1")
  expect_error(conditional_effects(fit, par = "criterion", effects = "x"), NA)
  expect_length(seen$grid, 1)
})

test_that("a parameter without effects of its own still draws its ids, so the next parameter's draws are unchanged (#512)", {
  local_mocked_bindings(ndraws = function(x) 100L, .package = "brms")
  fit <- mock_rating_fit(bmf(d ~ 1, criterion ~ 1 + x, spacing ~ 1))
  seen <- stub_category_route(stub_grid("x"))

  withr::with_seed(512, conditional_effects(fit, ndraws = 5))
  expect_length(seen$linpred, 1)
  expect_equal(seen$linpred[[1]]$draw_ids,
               withr::with_seed(512, {
                 sample.int(100, 5)
                 sample.int(100, 5)
               }))
})

test_that(".formula_nodes lists the formulas of a parameter and of the parameters it names (#512)", {
  fit <- mock_rating_fit(
    bmf(d ~ exp(n1) + x, n1 ~ n2 + cond, n2 ~ 1 + cond2 + (1 | id),
        criterion ~ 1, spacing ~ 1)
  )

  nodes <- .formula_nodes(fit, "d")
  expect_length(nodes, 3)
  expect_equal(vapply(nodes, function(node) length(node$sub_pars), 1L), c(1L, 1L, 0L))
  expect_equal(nodes[[1]]$sub_pars, "n1")
  expect_equal(nodes[[1]]$term_vars, list(character(0), "x"))
  expect_equal(nodes[[2]]$term_vars, list(character(0), "cond"))
  expect_equal(nodes[[3]]$term_vars, list("cond2"))
  expect_equal(nodes[[2]]$labels, c("n2", "cond"))
  expect_equal(nodes[[3]]$labels, "cond2")
  expect_equal(.formula_nodes(fit, "criterion")[[1]]$term_vars, list())
  expect_length(.formula_nodes(fit, "sdratio"), 0)

  # a model parameter without a user formula is a sub-parameter without a node
  fit_fixed <- mock_rating_fit(bmf(d ~ exp(sdratio) + x, criterion ~ 1, spacing ~ 1))
  nodes_fixed <- .formula_nodes(fit_fixed, "d")
  expect_length(nodes_fixed, 1)
  expect_equal(nodes_fixed[[1]]$sub_pars, "sdratio")
  expect_equal(nodes_fixed[[1]]$term_vars, list(character(0), "x"))
})

test_that(".np_grid_vars lists each variable once, in the order of the formulas, and follows the re_formula (#512)", {
  fit <- mock_rating_fit(
    bmf(d ~ exp(n1) + x, n1 ~ n2 + cond, n2 ~ 1 + cond2 + (1 | id),
        criterion ~ 1, spacing ~ 1)
  )

  expect_equal(.np_grid_vars(fit, c("d", "d"), NULL), c("x", "cond", "cond2", "id"))
  expect_equal(.np_grid_vars(fit, "d", NA), c("x", "cond", "cond2"))
  expect_equal(.np_grid_vars(fit, "criterion", NULL), character(0))
})

test_that(".np_grid_vars follows the sub-parameters in the order the formula names them (#512)", {
  fit <- mock_rating_fit(
    bmf(d ~ nl1 + nl2, nl1 ~ 1 + x, nl2 ~ 1 + cond, criterion ~ 1, spacing ~ 1)
  )

  expect_equal(.np_grid_vars(fit, "d", NA), c("x", "cond"))
})

test_that(".np_grid_vars skips a parameter fixed to a constant (#512)", {
  fit <- mock_rating_fit(bmf(d ~ 1 + x, criterion ~ 1, spacing = 0.3))

  expect_equal(.np_grid_vars(fit, c("d", "spacing"), NA), "x")
  expect_length(.formula_nodes(fit, "spacing"), 0)
})

# ===========================================================================
# Tier 1: Unit tests — mixture3p softmax route
# ===========================================================================

# The softmax of the nlpars is a function of their draws on the brms grid, so
# every argument that shapes those draws must reach brms as the user gave it.
softmax_ce <- function(par, draws_t, draws_nt, spaghetti = FALSE, ...) {
  seen <- stub_category_route(
    stub_grid("set_size", n = ncol(draws_t), spaghetti = spaghetti),
    function(args) if (args$nlpar == "thetat") draws_t else draws_nt
  )
  out <- conditional_effects(mock_mixture3p_fit(), par = par,
                             effects = "set_size", spaghetti = spaghetti, ...)
  list(out = out[[1]], seen = seen)
}

softmax_of <- function(draws_t, draws_nt) {
  1 / (1 + exp(draws_nt - draws_t) + exp(-draws_t))
}

test_that("mixture3p softmax route asks brms for the grid of the parameter it returns (#512)", {
  draws <- matrix(0, 3, 2)
  for (par in c("thetat", "thetant")) {
    res <- softmax_ce(par, draws, draws)
    expect_equal(res$seen$grid[[1]]$nlpar, par)
    expect_equal(vapply(res$seen$linpred, `[[`, "", "nlpar"), c("thetat", "thetant"))
  }
})

test_that("mixture3p softmax route returns a brms_conditional_effects object (#512)", {
  stub_category_route(stub_grid("set_size", n = 2))
  expect_s3_class(
    .compute_softmax_ce_inner(mock_mixture3p_fit(), "thetat"),
    "brms_conditional_effects"
  )
})

test_that("mixture3p softmax route forwards re_formula = NULL to brms (#512)", {
  draws <- matrix(0, 3, 2)
  res <- softmax_ce("thetat", draws, draws, re_formula = NULL)
  for (args in c(res$seen$grid, res$seen$linpred)) {
    expect_true("re_formula" %in% names(args))
    expect_null(args$re_formula)
  }
})

test_that("mixture3p softmax route forwards a re_formula formula to brms (#512)", {
  draws <- matrix(0, 3, 2)
  res <- softmax_ce("thetat", draws, draws, re_formula = ~ (1 | id))
  for (args in c(res$seen$grid, res$seen$linpred)) {
    expect_equal(args$re_formula, ~ (1 | id))
  }
})

test_that("mixture3p softmax route ignores group-level effects by default (#512)", {
  draws <- matrix(0, 3, 2)
  res <- softmax_ce("thetat", draws, draws)
  for (args in c(res$seen$grid, res$seen$linpred)) {
    expect_identical(args$re_formula, NA)
  }
})

test_that("mixture3p softmax route evaluates thetat and thetant on the draws the user selected (#512)", {
  draws <- matrix(0, 3, 2)
  res <- softmax_ce("thetat", draws, draws, draw_ids = c(2L, 5L))
  expect_equal(res$seen$grid[[1]]$draw_ids, 1)
  expect_length(res$seen$linpred, 2)
  for (args in res$seen$linpred) {
    expect_equal(args$draw_ids, c(2L, 5L))
  }
})

test_that("mixture3p softmax route sends the grid arguments to brms and the rest to posterior_linpred (#512)", {
  draws <- matrix(0, 3, 2)
  conditions <- data.frame(set_size = 2)
  res <- softmax_ce("thetat", draws, draws, spaghetti = TRUE,
                    conditions = conditions, int_conditions = list(set_size = 1:2),
                    resolution = 7, too_far = 0.1, sample_new_levels = "gaussian")

  grid_args <- res$seen$grid[[1]]
  expect_true(grid_args$spaghetti)
  expect_equal(grid_args$conditions, conditions)
  expect_equal(grid_args$int_conditions, list(set_size = 1:2))
  expect_equal(grid_args$resolution, 7)
  expect_equal(grid_args$too_far, 0.1)
  for (args in res$seen$linpred) {
    expect_equal(args$sample_new_levels, "gaussian")
  }
})

test_that("mixture3p softmax route leaves transform and method to the grid call (#512)", {
  draws <- matrix(0, 3, 2)
  res <- softmax_ce("thetat", draws, draws, transform = NULL, method = "posterior_predict")

  expect_equal(res$seen$grid[[1]]$method, "posterior_predict")
  expect_true("transform" %in% names(res$seen$grid[[1]]))
  for (args in res$seen$linpred) {
    expect_false(any(c("transform", "method") %in% names(args)))
  }
})

test_that("mixture3p softmax route summarises the softmax of the draws (#512)", {
  # one column is far in the tails, where a naive exp() overflows
  draws_t <- cbind(c(0.2, 1.1, -0.4, 2.3), c(900, 850, 800, 880), c(-900, -850, -800, -880))
  draws_nt <- cbind(c(-0.3, 0.5, 0.9, -1.2), c(-900, -850, 800, -880), c(-900, -870, -810, -890))

  for (par in c("thetat", "thetant")) {
    p <- if (par == "thetat") {
      softmax_of(draws_t[, 1], draws_nt[, 1])
    } else {
      softmax_of(draws_nt[, 1], draws_t[, 1])
    }
    res <- softmax_ce(par, draws_t, draws_nt)$out
    expect_equal(res$estimate__[1], mean(p))
    expect_equal(res$se__[1], sd(p))
    expect_equal(res$lower__[1], unname(quantile(p, 0.025)))
    expect_equal(res$upper__[1], unname(quantile(p, 0.975)))
    expect_true(all(is.finite(res$estimate__)))
  }
  res <- softmax_ce("thetat", draws_t, draws_nt)$out
  expect_equal(res$estimate__[2], mean(c(1, 1, 0.5, 1)), tolerance = 1e-12)
  # both nlpars far below the reference category: the probability is 0
  expect_equal(res$estimate__[3], 0)
})

test_that("mixture3p softmax route honours prob, probs and robust (#512)", {
  draws_t <- cbind(c(0.2, 1.1, -0.4, 2.3, 4.1), c(0, 0, 0, 0, 3))
  draws_nt <- draws_t * 0
  p <- softmax_of(draws_t[, 1], draws_nt[, 1])

  res <- softmax_ce("thetat", draws_t, draws_nt, prob = 0.5)$out
  expect_equal(res$lower__[1], unname(quantile(p, 0.25)))
  expect_equal(res$upper__[1], unname(quantile(p, 0.75)))

  res <- softmax_ce("thetat", draws_t, draws_nt, probs = c(0.1, 0.9))$out
  expect_equal(res$lower__[1], unname(quantile(p, 0.1)))
  expect_equal(res$upper__[1], unname(quantile(p, 0.9)))

  res <- softmax_ce("thetat", draws_t, draws_nt, robust = TRUE)$out
  expect_equal(res$estimate__[1], median(p))
  expect_equal(res$se__[1], mad(p))
  expect_false(isTRUE(all.equal(median(p), mean(p))))
})

test_that("mixture3p softmax route puts the spaghetti lines on the probability scale (#512)", {
  draws_t <- cbind(c(0.2, 1.1, -0.4, 2.3), c(0.5, -1, 0.3, 0.8))
  draws_nt <- cbind(c(-0.3, 0.5, 0.9, -1.2), c(0.1, 0.2, -0.6, 1.4))

  for (par in c("thetat", "thetant")) {
    p <- if (par == "thetat") {
      softmax_of(draws_t, draws_nt)
    } else {
      softmax_of(draws_nt, draws_t)
    }
    res <- softmax_ce(par, draws_t, draws_nt, spaghetti = TRUE)$out
    spaghetti <- attr(res, "spaghetti")
    expect_equal(spaghetti$estimate__, as.numeric(t(p)))
    expect_equal(spaghetti$set_size, rep(1:2, times = 4))
    expect_equal(spaghetti$sample__, factor(rep(1:4, each = 2)))
  }

  res <- softmax_ce("thetat", draws_t, draws_nt)$out
  expect_null(attr(res, "spaghetti"))
})

# ===========================================================================
# Tier 1: Unit tests — scale of the summary columns and spaghetti lines
# ===========================================================================

test_that("conditional_effects puts the category-family route on the requested scale (#512)", {
  dat <- cbind(
    data.frame(stimulus = c(0L, 1L), x = c(0, 1)),
    rsdt_rating(2, 50, c(0L, 1L), d = 1, thresholds = c(-0.5, 0, 0.5))
  )
  fit <- suppressWarnings(mock_bmm(
    bmf(d ~ 1 + x, criterion ~ 1, spacing ~ 1), dat,
    sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", links = list(d = "log"))
  ))
  draws <- rbind(c(0, 1, 2), c(1, 3, 2), c(2, 2, 5))
  stub_category_route(stub_grid("x", spaghetti = TRUE), draws)

  sampling <- conditional_effects(fit, par = "d", effects = "x", spaghetti = TRUE,
                                  scale = "sampling")$x
  native <- conditional_effects(fit, par = "d", effects = "x", spaghetti = TRUE)$x
  expect_equal(sampling$estimate__, apply(draws, 2, median))
  expect_equal(attr(sampling, "spaghetti")$estimate__, as.numeric(t(draws)))
  expect_equal(native$estimate__, exp(sampling$estimate__))
  expect_equal(native$lower__, exp(sampling$lower__))
  expect_equal(native$upper__, exp(sampling$upper__))
  expect_equal(native$se__, sampling$se__)
  expect_equal(attr(native, "spaghetti")$estimate__, exp(as.numeric(t(draws))))
})

test_that("conditional_effects puts a distributional parameter on the requested scale (#512)", {
  df <- data.frame(x = 1:2, estimate__ = c(1, 2), se__ = c(0.5, 0.5),
                   lower__ = c(0.5, 1), upper__ = c(2, 4))
  attr(df, "effects") <- "x"
  attr(df, "spaghetti") <- data.frame(x = 1:2, estimate__ = c(1, 3),
                                      sample__ = factor(1))
  local_mocked_bindings(.brms_conditional_effects = function(x, ...) {
    structure(list(x = df), class = "brms_conditional_effects")
  })
  fit <- mock_sdm_fit()

  native <- conditional_effects(fit, par = "kappa")$x
  sampling <- conditional_effects(fit, par = "kappa", scale = "sampling")$x
  expect_equal(native$estimate__, df$estimate__)
  expect_equal(sampling$estimate__, log(df$estimate__))
  expect_equal(sampling$lower__, log(df$lower__))
  expect_equal(sampling$se__, df$se__)
  expect_equal(attr(sampling, "spaghetti")$estimate__, log(c(1, 3)))
})

# ===========================================================================
# Tier 1: Unit tests — the par argument
# ===========================================================================

test_that("conditional_effects errors for invalid par name", {
  expect_error(
    conditional_effects(mock_sdm_fit(), par = "nonexistent"),
    "not found in model"
  )
})

test_that("conditional_effects errors for a par that is not a single character string", {
  fit <- mock_sdm_fit()

  expect_error(conditional_effects(fit, par = 42), "must be a single character string")
  expect_error(conditional_effects(fit, par = c("c", "kappa")),
               "must be a single character string")
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
  p <- suppressMessages(plot(ce, plot = FALSE))
  expect_true(length(p) > 0)
})


# ===========================================================================
# Tier 2: Multinomial families (m3 fixture; rows differ in nTrials)
# ===========================================================================

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

test_that("m3 conditional_effects reports invalid draw and interval arguments as brms does (#512)", {
  skip_on_cran()
  fit <- load_m3_fit()

  expect_error(conditional_effects(fit, par = "c", ndraws = 0),
               "Argument 'ndraws' should be between 1 and the maximum number of draws")
  expect_error(conditional_effects(fit, par = "c", ndraws = 10^6),
               "Argument 'ndraws' should be between 1 and the maximum number of draws")
  expect_error(conditional_effects(fit, par = "c", draw_ids = 0),
               "Some 'draw_ids' indices are out of range.", fixed = TRUE)
  expect_error(conditional_effects(fit, par = "c", prob = 1.5),
               "'prob' must be a single numeric value in [0, 1].", fixed = TRUE)
  expect_error(
    suppressWarnings(conditional_effects(fit, par = "c", probs = c(0.1, 0.5, 0.9))),
    "Arguments 'probs' must be of length 2.", fixed = TRUE
  )
  expect_warning(conditional_effects(fit, par = "c", probs = c(0.1, 0.9)),
                 "Argument 'probs' is deprecated. Please use 'prob' instead.",
                 fixed = TRUE)
  expect_error(conditional_effects(fit, par = "c", draw_ids = integer(0)),
               "No posterior draws supplied.", fixed = TRUE)
  expect_error(conditional_effects(fit, par = "c", ndraws = "a"),
               "Cannot coerce 'ndraws' to a single integer value.", fixed = TRUE)
  expect_error(conditional_effects(fit, par = "c", re = NULL),
               "matches multiple formal arguments")
  expect_error(conditional_effects(fit, par = "c", transform = exp),
               "parameters of the 'm3' model are latent quantities")
})

test_that("m3 conditional_effects keeps only effects of the parameter's own predictors", {
  skip_on_cran()
  fit <- load_m3_fit()

  expect_setequal(names(conditional_effects(fit)), c("a.cond", "c.cond"))
  expect_length(conditional_effects(fit, par = "d"), 0)
})


test_that("m3 conditional_effects is quiet and plots with points = TRUE (#512)", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()

  expect_no_message(ce <- conditional_effects(fit))
  expect_no_warning(
    suppressMessages(
      ggplot2::ggplot_build(plot(ce, points = TRUE, plot = FALSE)[[1]])
    )
  )
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
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", links = list(d = "log"))
  fit <- bmm(bmf(d ~ 1 + cond + x, criterion ~ 1, spacing ~ 1), dat, model,
             backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  expect_no_message(ce <- conditional_effects(fit))
  expect_setequal(names(ce), c("d.cond", "d.x"))

  ce <- conditional_effects(
    fit, par = "d", effects = "x:cond", int_conditions = list(cond = "B")
  )
  expect_true(all(ce[[1]]$cond == "B"))

  draw_ids <- 1:11
  ce <- conditional_effects(fit, par = "d", effects = "x", spaghetti = TRUE,
                            draw_ids = draw_ids, scale = "sampling")
  draws <- brms::posterior_linpred(fit, newdata = ce$x, nlpar = "d",
                                   re_formula = NA, draw_ids = draw_ids)
  expect_equal(ce$x$estimate__, apply(draws, 2, median))
  expect_equal(attr(ce$x, "spaghetti")$estimate__, as.numeric(t(draws)))

  ce <- conditional_effects(fit, par = "d", effects = "x", spaghetti = TRUE,
                            draw_ids = draw_ids)
  expect_equal(ce$x$estimate__, exp(apply(draws, 2, median)))
  expect_equal(attr(ce$x, "spaghetti")$estimate__, exp(as.numeric(t(draws))))
})

test_that("m3 conditional_effects does not evaluate the category dpar on other formulas' effects (#512)", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)))

  # n_dist is 5, 0, 0 across conditions and enters the first category dpar
  # through log(), so brms's grid of Idx_dist:n_dist (n_dist at mean - sd < 0)
  # makes it NaN
  formula <- bmf(corr ~ b + a + c, dist ~ b + d, other ~ b + a, npl ~ b,
                 c ~ 0 + cond, a ~ 1, d ~ 1, b = 0)
  model <- m3(resp_cats = c("dist", "corr", "other", "npl"),
              num_options = c("n_dist", "n_corr", "n_other", "n_npl"),
              version = "custom", links = list(a = "log", c = "log", d = "log"))
  fit <- suppressWarnings(suppressMessages(bmm(
    formula, oberauer_lewandowsky_2019_e1, model, backend = "cmdstanr",
    chains = 1, iter = 300, warmup = 150, refresh = 0, silent = 2, seed = 1
  )))

  expect_no_warning(ce_c <- conditional_effects(fit, par = "c"))
  expect_no_warning(ce_a <- conditional_effects(fit, par = "a"))
  expect_no_warning(ce_all <- conditional_effects(fit))
  expect_equal(ce_c, conditional_effects(fit, par = "c", effects = "cond"))
  expect_length(ce_a, 0)
  expect_named(ce_all, "c.cond")
  expect_equal(
    .ce_brms_default_effects(fit, NA),
    unname(lapply(suppressWarnings(suppressMessages(.brms_conditional_effects(
      fit, dpar = names(brms::brmsterms(fit$formula)$dpars)[1], draw_ids = 1
    ))), attr, "effects"))
  )
})
