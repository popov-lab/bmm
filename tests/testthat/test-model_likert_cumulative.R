############################################################################# !
# CONSTRUCTOR                                                            ####
############################################################################# !

test_that("likert_cumulative() declares activation and the thresholds of its type", {
  expect_named(likert_cumulative("y", "item", n_categories = 5)$parameters,
               c("activation", "delta1", "delta2", "delta3"))
  expect_named(likert_cumulative("y", "item", n_categories = 5, threshold_type = "parsimonious")$parameters,
               c("activation", "spacing"))
  expect_named(likert_cumulative("y", "item", n_categories = 5, threshold_type = "softmax")$parameters,
               c("activation", "spacing", "delta1", "delta2"))
  expect_named(likert_cumulative("y", "item", n_categories = 4, threshold_type = "log_ratio")$parameters,
               c("activation", "delta1", "delta2"))
  # with K unknown until the data arrive, only the K-free parts are declared
  expect_named(likert_cumulative("y", "item")$parameters, "activation")
  expect_named(likert_cumulative("y", "item", threshold_type = "equidistant")$parameters,
               c("activation", "spacing"))
  model <- likert_cumulative("y", "item", n_categories = 5)
  expect_s3_class(model, c("bmmodel", "likert_cumulative"), exact = TRUE)
  expect_false(any(grepl("^thres|^criterion$|^disc", names(model$parameters))))
  expect_identical(model$domain, "Personality and individual differences")
  expect_identical(model$task, "Questionnaire and rating-scale responses")
})

test_that("likert_cumulative() refuses arguments it cannot use", {
  for (bad in list(2, 1, 2.5, c(5, 6), "5", NA, Inf)) {
    expect_error(likert_cumulative("y", "item", n_categories = bad),
                 "n_categories must be NULL or a single whole number of at least 3")
  }
  expect_error(likert_cumulative("y", "item", n_categories = 3, threshold_type = "log_ratio"),
               "log_ratio thresholds require at least 4")
  expect_error(likert_cumulative(c("y1", "y2"), "item"), "response must name one column")
  expect_error(likert_cumulative("y", c("a", "b")), "item must name one column")
  expect_error(likert_cumulative("y", "item", keying = c("k1", "k2")),
               "keying must be NULL or the name of one column")
  expect_error(likert_cumulative("y", "item", keying = 1),
               "keying must be NULL or the name of one column")
  expect_error(likert_cumulative("y"), "item")
  expect_error(likert_cumulative("y", "item", dist = "cauchy"), "should be one of")
})

test_that("the default priors follow D16, with a K-matched delta mean on log_distance", {
  priors <- function(K, type = "log_distance") {
    likert_cumulative("y", "item", n_categories = K, threshold_type = type)$default_priors
  }
  expect_identical(priors(5)$activation,
                   list(main = "normal(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"))
  # log(2 * log(K - 1) / (K - 2)) to two decimals
  expect_identical(priors(3)$delta1$main, "normal(0.33, 0.5)")
  expect_identical(priors(4)$delta2$main, "normal(0.09, 0.5)")
  expect_identical(priors(5)$delta3$main, "normal(-0.08, 0.5)")
  expect_identical(priors(7)$delta5$main, "normal(-0.33, 0.5)")
  expect_identical(priors(7)$delta1[c("effects", "sd")],
                   list(effects = "normal(0, 0.3)", sd = "exponential(2)"))
  expect_identical(priors(5, "parsimonious")$spacing,
                   list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"))
  expect_identical(priors(5, "log_ratio")$delta1$main, "normal(0, 1)")
  expect_identical(priors(5, "softmax")$delta2$main, "normal(0, 1)")
})

test_that("only activation takes another link", {
  expect_identical(settable_links(likert_cumulative("y", "item", n_categories = 5)), "activation")
  model <- likert_cumulative("y", "item", n_categories = 5, links = list(activation = "log"))
  expect_identical(model$links$activation, "log")
  expect_error(likert_cumulative("y", "item", n_categories = 5, links = list(delta1 = "log")),
               "delta1")
})

test_that("the printed name follows the noise distribution", {
  name <- function(dist) likert_cumulative("y", "item", dist = dist)$name
  expect_match(name("normal"), "probit")
  expect_match(name("logistic"), "proportional odds")
  expect_match(name("gumbel_min"), "complementary log-log")
  expect_match(name("gumbel_max"), "\\(log-log\\)")
  expect_message(likert_cumulative("y", "item", dist = "logistic", threshold_type = "parsimonious"),
                 "the equidistant type is the rating scale model")
  expect_no_message(likert_cumulative("y", "item", dist = "logistic", threshold_type = "equidistant"))
  expect_no_message(likert_cumulative("y", "item", dist = "logistic"))
})


############################################################################# !
# CHECK_MODEL, CHECK_DATA, CHECK_FORMULA                                 ####
############################################################################# !

lc_check <- function(model, data, formula = bmf(activation ~ 0 + item + (1 | id))) {
  model <- check_model(model, data, formula)
  check_data(model, data, formula)
}

all_categories <- function(dat, K = 5L) {
  dat$y[seq_len(K)] <- seq_len(K)
  dat
}

test_that("an ordered factor gives K from its levels and keeps them for update()", {
  dat <- all_categories(likert_sim_frame())
  levs <- c("sd", "d", "n", "a", "sa")
  dat$y <- factor(levs[dat$y], levels = levs, ordered = TRUE)
  model <- check_model(likert_cumulative("y", "item"), dat, bmf(activation ~ 0 + item))
  expect_identical(model$other_vars$n_categories, 5L)
  expect_identical(model$other_vars$response_levels, levs)
  expect_named(model$parameters, c("activation", "delta1", "delta2", "delta3"))
  expect_identical(model$default_priors$delta3$main, "normal(-0.08, 0.5)")
  checked <- check_data(model, dat, bmf(activation ~ 0 + item))
  expect_identical(checked$lc_resp, as.integer(dat$y))

  reverted <- revert_check_data(model, checked[c("lc_resp", "lc_k", "lc_dist", "lc_thresh", "item")])
  expect_identical(reverted$y, dat$y)
  expect_false(any(grepl("^lc_", names(reverted))))
})

test_that("the response encodings are refused where K cannot be read safely", {
  dat <- all_categories(likert_sim_frame())
  expect_error(lc_check(likert_cumulative("y", "item"), dat),
               "n_categories.*integer responses")
  unordered <- dat
  unordered$y <- factor(unordered$y)
  expect_error(lc_check(likert_cumulative("y", "item"), unordered), "ordered factor")
  ordered <- dat
  ordered$y <- factor(ordered$y, ordered = TRUE)
  expect_error(lc_check(likert_cumulative("y", "item", n_categories = 6), ordered),
               "n_categories = 6.*5 levels")
  two <- dat
  two$y <- factor(ifelse(dat$y > 2, "agree", "disagree"), levels = c("disagree", "agree"),
                  ordered = TRUE)
  expect_error(lc_check(likert_cumulative("y", "item"), two), "at least 3 categories")
  three <- dat
  three$y <- factor(pmin(dat$y, 3L), ordered = TRUE)
  expect_error(lc_check(likert_cumulative("y", "item", threshold_type = "log_ratio"), three),
               "log_ratio thresholds require at least 4")
  expect_error(lc_check(likert_cumulative("missing", "item", n_categories = 5), dat),
               "'missing' missing in the data")
})

test_that("integer responses must be whole numbers from 1 to K", {
  dat <- all_categories(likert_sim_frame())
  model <- likert_cumulative("y", "item", n_categories = 5)
  for (bad in list(0, 6, 2.5, -1)) {
    dat_bad <- dat
    dat_bad$y[7] <- bad
    expect_error(lc_check(model, dat_bad), "whole numbers from 1 to 5")
  }
  as_char <- dat
  as_char$y <- as.character(dat$y)
  expect_error(lc_check(model, as_char), "integer category from 1 to 5 or an ordered factor")
  checked <- lc_check(model, dat)
  expect_identical(checked$lc_resp, as.integer(dat$y))
  expect_identical(unique(checked$lc_k), 5L)
  expect_identical(unique(checked$lc_dist), 1L)
  expect_identical(unique(checked$lc_thresh), 3L)
  expect_s3_class(checked$item, "factor")
  logistic <- lc_check(likert_cumulative("y", "item", n_categories = 5, dist = "logistic",
                                         threshold_type = "softmax"), dat)
  expect_identical(unique(logistic$lc_dist), 4L)
  expect_identical(unique(logistic$lc_thresh), 5L)
})

test_that("NA responses are dropped with a warning that counts them", {
  dat <- all_categories(likert_sim_frame())
  dat$y[c(8, 9)] <- NA
  expect_warning(checked <- lc_check(likert_cumulative("y", "item", n_categories = 5), dat),
                 "2 rows with a missing response")
  expect_identical(nrow(checked), nrow(dat) - 2L)
})

test_that("a category nobody chose is named in a warning", {
  dat <- likert_sim_frame()
  dat$y[dat$y == 5L] <- 4L
  expect_warning(lc_check(likert_cumulative("y", "item", n_categories = 5), dat),
                 "Category 5 of 5 is never observed")
})

test_that("item and keying are checked", {
  dat <- all_categories(likert_sim_frame())
  na_item <- dat
  na_item$item[3] <- NA
  expect_error(lc_check(likert_cumulative("y", "item", n_categories = 5), na_item),
               "Item variable 'item' must not contain NA")
  expect_error(lc_check(likert_cumulative("y", "it", n_categories = 5), dat),
               "Item variable 'it' missing in the data")
  keyed <- likert_cumulative("y", "item", n_categories = 5, keying = "keying")
  expect_no_warning(lc_check(keyed, dat))
  for (bad in list(c(1, 0), c(1, NA), c(1, 2))) {
    dat_bad <- dat
    dat_bad$keying[1:2] <- bad
    expect_error(lc_check(keyed, dat_bad), "must be 1 \\(regular\\) or -1 \\(reversed\\)")
  }
  varies <- dat
  varies$keying[1] <- -varies$keying[1]
  expect_error(lc_check(keyed, varies), "constant within each item.*i1")
  expect_error(lc_check(likert_cumulative("y", "item", n_categories = 5, keying = "key"), dat),
               "Keying variable 'key' missing in the data")
})

test_that("check_data() warns before overwriting a column it builds", {
  dat <- all_categories(likert_sim_frame())
  dat$lc_k <- 1
  expect_warning(lc_check(likert_cumulative("y", "item", n_categories = 5), dat),
                 "lc_k.*reserved")
})

test_that("activation needs a formula and item terms belong in it", {
  dat <- all_categories(likert_sim_frame())
  model <- likert_cumulative("y", "item", n_categories = 5)
  expect_error(check_model(model, dat, bmf(delta1 ~ 1)), "formula for activation")
  checked <- lc_check(model, dat)
  expect_warning(check_formula(model, checked, bmf(activation ~ 1 + (1 | id))),
                 "confounded with the intercept")
  expect_no_warning(suppressMessages(check_formula(model, checked, bmf(activation ~ 0 + item + (1 | id)))))
  expect_no_warning(suppressMessages(check_formula(model, checked, bmf(activation ~ 1 + (1 | item) + (1 | id)))))
  checked$item_loc <- as.numeric(checked$item)
  expect_no_warning(suppressMessages(check_formula(model, checked,
                                                   bmf(activation ~ 1 + offset(item_loc) + (1 | id)))))
})

test_that("the threshold terms that change what the model is are named", {
  dat <- all_categories(likert_sim_frame())
  pars <- likert_cumulative("y", "item", n_categories = 5, threshold_type = "parsimonious")
  checked <- lc_check(pars, dat)
  expect_message(check_formula(pars, checked, bmf(activation ~ 0 + item, spacing ~ 0 + item)),
                 "item-specific threshold spread.*not fitted")
  expect_message(check_formula(pars, checked, bmf(activation ~ 0 + item, spacing ~ 1 + (1 | id))),
                 "scaled threshold model")
  ld <- likert_cumulative("y", "item", n_categories = 5)
  expect_message(check_formula(ld, lc_check(ld, dat),
                               bmf(activation ~ 0 + item, delta1 ~ 1 + (1 | id))),
                 "Johnson \\(2003\\)")
  # an item grouping term is no person term
  notes <- function(model, formula) {
    capture_messages(check_formula(model, lc_check(model, dat), formula))
  }
  item_spacing <- notes(pars, bmf(activation ~ 0 + item, spacing ~ 1 + (1 | item)))
  expect_false(any(grepl("person term", item_spacing)))
  expect_true(any(grepl("item-specific threshold spread", item_spacing)))
  expect_false(any(grepl("Johnson", notes(ld, bmf(activation ~ 0 + item, delta1 ~ 1 + (1 | item))))))
  expect_true(any(grepl("scaled threshold model",
                        notes(pars, bmf(activation ~ 0 + item, spacing ~ 1 + (1 | id:item))))))
  # softmax scales its whole threshold set by exp(spacing): with item-specific
  # widths this is Lubbe and Schuster's model with item-specific thresholds
  sm <- likert_cumulative("y", "item", n_categories = 5, threshold_type = "softmax")
  expect_true(any(grepl("scaled threshold model with shared thresholds",
                        notes(sm, bmf(activation ~ 0 + item, spacing ~ 1 + (1 | id))))))
  expect_true(any(grepl("scaled threshold model with item-specific thresholds",
                        notes(sm, bmf(activation ~ 0 + item, spacing ~ 0 + item + (1 | id),
                                      delta1 ~ 0 + item, delta2 ~ 0 + item)))))
  # without an item term on spacing the items share one span
  shared_span <- notes(sm, bmf(activation ~ 0 + item, spacing ~ 1 + (1 | id), delta1 ~ 0 + item))
  expect_false(any(grepl("with item-specific thresholds", shared_span)))
  expect_true(any(grepl("item-specific interval proportions", shared_span)))
  # keying is constant within item, so grouping by it is no person term
  keyed <- likert_cumulative("y", "item", n_categories = 5, keying = "keying")
  expect_false(any(grepl("Johnson", notes(keyed, bmf(activation ~ 0 + item,
                                                     delta1 ~ 1 + (1 | item:keying))))))
})


############################################################################# !
# CONFIGURE_MODEL AND STAN CODE                                          ####
############################################################################# !

test_that("configure_model() builds the custom family with one dpar per threshold part", {
  dat <- all_categories(likert_sim_frame(K = 7L), 7L)
  for (K in 4:7) {
    dat_k <- dat
    dat_k$y <- pmin(dat$y, K)
    dat_k$y[seq_len(K)] <- seq_len(K)
    model <- likert_cumulative("y", "item", n_categories = K)
    formula <- bmf(activation ~ 0 + item)
    checked <- lc_check(model, dat_k, formula)
    cfg <- suppressMessages(configure_model(model, checked, check_formula(model, checked, formula)))
    fam <- cfg$formula$family
    expect_identical(fam$name, "likert_cumulative")
    expect_identical(fam$dpars, c("mu", "activation", paste0("delta", as.roman(seq_len(K - 2)))))
    expect_identical(fam$type, "int")
    expect_true(fam$loop)
    expect_identical(fam$vars, c("vint1[n]", "vint2[n]", "vint3[n]"))
    scode <- paste(vapply(cfg$stanvars, `[[`, "", "scode"), collapse = "\n")
    expect_match(scode, paste0(
      "real likert_cumulative_lpmf\\(int y, real mu, real activation, ",
      paste0("real delta", as.roman(seq_len(K - 2)), collapse = ", "),
      ", int K, int dist_type, int thresh_type\\)"
    ))
    expect_match(scode, paste0("array\\[", K - 2, "\\] real deltas"))
  }
  pars <- likert_cumulative("y", "item", n_categories = 5, threshold_type = "softmax")
  cfg <- suppressMessages(configure_model(pars, lc_check(pars, all_categories(likert_sim_frame())),
                                          bmf(activation ~ 1)))
  expect_identical(cfg$formula$family$dpars, c("mu", "activation", "spacing", "deltaI", "deltaII"))
})

test_that("stancode() carries the generated wrapper and the four chunks", {
  dat <- all_categories(likert_sim_frame())
  code <- suppressMessages(stancode(bmf(activation ~ 0 + item + (1 | id)), dat,
                                    likert_cumulative("y", "item", n_categories = 5)))
  expect_match(code, "real likert_cumulative_lpmf(int y, real mu, real activation, real deltaI, real deltaII, real deltaIII, int K, int dist_type, int thresh_type)",
               fixed = TRUE)
  expect_match(code, "real likert_cumulative_log_prob(int y, real activation, real spacing, array[] real deltas,",
               fixed = TRUE)
  expect_match(code, "vector make_thresholds(real criterion", fixed = TRUE)
  expect_match(code, "real ordinal_log_prob_cat(int cat", fixed = TRUE)
  expect_match(code, "real sdt_log_cumprob(", fixed = TRUE)
  expect_match(code, "deltaI[n] = (nlp_delta1[n]);", fixed = TRUE)
  expect_match(code, "target += likert_cumulative_lpmf(Y[n] | mu[n], activation[n], deltaI[n], deltaII[n], deltaIII[n], vint1[n], vint2[n], vint3[n]);",
               fixed = TRUE)
})

test_that("a mock fit runs and keeps the user's parameter names", {
  dat <- all_categories(likert_sim_frame())
  fit <- suppressMessages(bmm(bmf(activation ~ 0 + item + (1 | id), delta2 ~ 1 + (1 | id)), dat,
                              likert_cumulative("y", "item", n_categories = 5),
                              backend = "mock", mock_fit = 1, rename = FALSE))
  expect_s3_class(fit, "bmmfit")
  prior <- brms::prior_summary(fit)
  expect_true(all(c("delta1", "delta2", "delta3") %in% prior$nlpar))
  expect_true("activation" %in% prior$dpar)
})


############################################################################# !
# DISTRIBUTION FUNCTIONS                                                 ####
############################################################################# !

lc_grid_thresholds <- function(K) {
  .make_thresholds(0, K, "log_distance", deltas = seq(-0.3, 0.3, length.out = K - 2))
}

test_that("dlikert_cumulative() is the ordinal kernel at scale 1 and sums to 1", {
  for (K in c(3L, 5L, 7L)) {
    thr <- lc_grid_thresholds(K)
    for (dist in .noise_dist_names) {
      for (act in c(-8, -1.3, 0, 0.7, 8)) {
        p <- dlikert_cumulative(seq_len(K), act, thr, dist = dist)
        expect_equal(sum(p), 1, tolerance = 1e-12)
        expect_equal(dlikert_cumulative(seq_len(K), act, thr, dist = dist, log = TRUE),
                     .ordinal_category_log_probs(thr, act, 1, dist))
      }
    }
  }
  # one row of thresholds per observation
  thr <- rbind(c(-1, 0, 1), c(-2, 0.5, 3))
  expect_equal(dlikert_cumulative(c(1, 4), c(0, 1), thr),
               c(pnorm(-1), 1 - pnorm(3 - 1)))
})

test_that("dlikert_cumulative() matches brms's cumulative family", {
  # brms's cloglog is P(Y <= k) = 1 - exp(-exp(thr_k - eta)), the gumbel_min cdf
  withr::local_seed(662)
  links <- c(probit = "normal", logit = "logistic", cloglog = "gumbel_min")
  activation <- stats::rnorm(50, 0, 1)
  thr <- t(apply(matrix(stats::rnorm(50 * 4, 0, 1.2), 50), 1, sort))
  for (link in names(links)) {
    for (k in 1:5) {
      theirs <- as.numeric(brms:::dcumulative(k, activation, thr, link = link))
      # brms differences the cdf on the probability scale, so its log loses
      # digits where a category is rare; the comparison is on its scale
      ours <- dlikert_cumulative(k, activation, thr, dist = links[[link]], log = TRUE)
      expect_lt(max(abs(exp(ours) - theirs)), 1e-12)
    }
  }
  # the other Gumbel is a different model
  theirs <- as.numeric(brms:::dcumulative(2, 0.3, rbind(c(-1, 0, 1, 2)), link = "cloglog"))
  expect_gt(abs(dlikert_cumulative(2, 0.3, c(-1, 0, 1, 2), dist = "gumbel_max") - theirs), 0.01)
})

test_that("dlikert_cumulative() refuses impossible responses and thresholds", {
  expect_error(dlikert_cumulative(5, 0, c(-1, 0, 1)), "x must be whole numbers from 1 to 4")
  expect_error(dlikert_cumulative(1.5, 0, c(-1, 0, 1)), "x must be whole numbers")
  expect_error(dlikert_cumulative(1, 0, c(0, -1, 1)), "thresholds must be increasing")
  expect_error(dlikert_cumulative(1, 0, 0), "at least 2 thresholds")
  expect_error(dlikert_cumulative(1, 0, c(-1, 1), dist = "cauchy"), "should be one of")
})

test_that("rlikert_cumulative() draws categories with dlikert_cumulative()'s probabilities", {
  withr::local_seed(17)
  thr <- c(-1.2, -0.3, 0.4, 1.5)
  y <- rlikert_cumulative(40000, 0.25, thr, dist = "logistic")
  expect_true(is.integer(y))
  expect_equal(as.numeric(table(factor(y, 1:5))) / 40000,
               dlikert_cumulative(1:5, 0.25, thr, dist = "logistic"), tolerance = 0.02)
  rows <- rlikert_cumulative(3, c(-20, 0, 20), rbind(thr, thr, thr))
  expect_identical(rows[c(1, 3)], c(1L, 5L))
})

test_that("the threshold builder places the thresholds around 0 by type", {
  thr <- .likert_cumulative_thresholds(5L, "parsimonious", spacing = 0)
  expect_equal(thr, log(1:4 / (4:1)))
  ld <- .likert_cumulative_thresholds(5L, "log_distance", deltas = c(0, 0, 0))
  expect_equal(ld, c(-1.5, -0.5, 0.5, 1.5))
  draws <- .likert_cumulative_thresholds(4L, "log_distance", deltas = cbind(c(0, log(2)), 0))
  expect_equal(draws, rbind(c(-1, 0, 1), c(-2, 0, 1)))
})


############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT, POSTERIOR_EPRED                            ####
############################################################################# !

lc_fake_prep <- function(dpars, y, K, dist, thresh) {
  n <- length(y)
  epred_prep(dpars, list(Y = y, vint1 = rep(K, n), vint2 = rep(dist, n), vint3 = rep(thresh, n)))
}

test_that("log_lik, posterior_predict and posterior_epred agree with the distribution functions", {
  withr::local_seed(5)
  n_draws <- 7
  draws <- function(lo, hi) matrix(stats::runif(n_draws * 3, lo, hi), n_draws, 3)
  dpars <- list(activation = draws(-1, 1), spacing = draws(-0.4, 0.2),
                deltaI = draws(-0.5, 0.5), deltaII = draws(-0.5, 0.5))
  y <- c(1L, 3L, 5L)
  prep <- lc_fake_prep(dpars, y, 5L, 2L, 5L)
  for (i in 1:3) {
    thr <- .likert_cumulative_thresholds(5L, "softmax", spacing = dpars$spacing[, i],
                                         deltas = cbind(dpars$deltaI[, i], dpars$deltaII[, i]))
    expect_equal(log_lik_likert_cumulative(i, prep),
                 dlikert_cumulative(y[i], dpars$activation[, i], thr, dist = "gumbel_min", log = TRUE))
    pp <- posterior_predict_likert_cumulative(i, prep)
    expect_length(pp, n_draws)
    expect_true(all(pp %in% 1:5))
    expected_index <- vapply(seq_len(n_draws), function(s) {
      sum(1:5 * dlikert_cumulative(1:5, dpars$activation[s, i], thr[s, ], dist = "gumbel_min"))
    }, numeric(1))
    expect_equal(posterior_epred_likert_cumulative(prep)[, i], expected_index)
  }
  expect_epred_by_cell(posterior_epred_likert_cumulative, dpars,
                       list(vint1 = rep(5L, 3), vint2 = rep(2L, 3), vint3 = rep(5L, 3)))
})

test_that("posterior_predict draws the categories with the model's probabilities", {
  withr::local_seed(8)
  n_draws <- 30000
  dpars <- list(activation = matrix(0.4, n_draws, 1), deltaI = matrix(-0.2, n_draws, 1),
                deltaII = matrix(0.1, n_draws, 1), deltaIII = matrix(0.3, n_draws, 1))
  prep <- lc_fake_prep(dpars, 1L, 5L, 1L, 3L)
  pp <- posterior_predict_likert_cumulative(1, prep)
  thr <- .likert_cumulative_thresholds(5L, "log_distance", deltas = c(-0.2, 0.1, 0.3))
  expect_equal(as.numeric(table(factor(pp, 1:5))) / n_draws,
               dlikert_cumulative(1:5, 0.4, thr), tolerance = 0.03)
})


############################################################################# !
# STAN PARITY                                                            ####
############################################################################# !

# Each (K, threshold type) model writes its own likert_cumulative_lpmf; the
# copies are renamed and called from one fixed_param generated-quantities
# program, so the wrapper bmm generates is what meets the R twin, not the
# shared kernel alone. sig_figs = 17 round-trips a double.
likert_stan_logp <- function(cases, grid) {
  sc <- system.file("stan_chunks", package = "bmm")
  chunks <- c("sdt_dist_functions.stan", "threshold_functions.stan", "ordinal_functions.stan",
              "likert_cumulative_functions.stan")
  funs <- paste(vapply(file.path(sc, chunks), read_lines2, character(1)), collapse = "\n")
  wrappers <- vapply(seq_along(cases), function(j) {
    sub("likert_cumulative_lpmf", paste0("lc_case", j, "_lpmf"),
        .likert_cumulative_lpmf_stan(cases[[j]]), fixed = TRUE)
  }, character(1))
  calls <- vapply(seq_along(cases), function(j) {
    np <- length(.likert_threshold_dpars(cases[[j]]))
    args <- paste0(", par[i, ", seq_len(np), "]", collapse = "")
    paste0("    if (case_id[i] == ", j, ") lp[i] = lc_case", j, "_lpmf(y[i] | 0", ", act[i]",
           if (np) args, ", K[i], dist[i], thresh[i]);")
  }, character(1))
  program <- paste0(
    "functions {\n", funs, "\n", paste(wrappers, collapse = "\n"), "\n}\n",
    "data {\n  int N; int P; array[N] int case_id; array[N] int y; array[N] int K;\n",
    "  array[N] int dist; array[N] int thresh; vector[N] act; array[N, P] real par;\n}\n",
    "generated quantities {\n  vector[N] lp;\n  for (i in 1:N) {\n",
    paste(calls, collapse = "\n"), "\n  }\n}\n"
  )
  data <- list(N = nrow(grid), P = ncol(grid$par), case_id = grid$case_id, y = grid$y,
               K = grid$K, dist = grid$dist, thresh = grid$thresh, act = grid$act, par = grid$par)
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1,
    iter_warmup = 0, refresh = 0, show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  as.numeric(csv[1, paste0("lp.", seq_len(nrow(grid)))])
}

test_that("the generated Stan likelihood matches its R twin in both tails", {
  skip_unless_stan_tests()
  skip_on_cran()

  # activation 8 puts every threshold in the lower tail of the latent
  # variable, -8 in the upper tail; K = 5, 6 and 7 cover the odd and even
  # placement around 0 and log_ratio intervals beyond the first on each side
  cases <- list()
  rows <- list()
  for (K in 5:7) {
    for (type in .threshold_types) {
      model <- .model_likert_cumulative("y", "item", n_categories = K, threshold_type = type)
      cases[[length(cases) + 1]] <- model
      np <- length(.likert_threshold_dpars(model))
      pars <- seq(-0.4, 0.3, length.out = max(np, 1))[seq_len(np)]
      grid <- expand.grid(y = seq_len(K), act = c(-8, -0.9, 0, 1.7, 8),
                          dist = seq_along(.noise_dist_names))
      grid$case_id <- length(cases)
      grid$K <- K
      grid$thresh <- match(type, .threshold_types)
      grid$par <- matrix(c(pars, rep(0, 5 - np)), nrow(grid), 5, byrow = TRUE)
      rows[[length(rows) + 1]] <- grid
    }
  }
  grid <- do.call(rbind, rows)
  stan <- likert_stan_logp(cases, grid)

  r <- vapply(seq_len(nrow(grid)), function(i) {
    model <- cases[[grid$case_id[i]]]
    np <- length(.likert_threshold_dpars(model))
    par <- grid$par[i, seq_len(np)]
    type <- model$other_vars$threshold_type
    has_spacing <- "spacing" %in% names(model$parameters)
    thr <- .likert_cumulative_thresholds(
      grid$K[i], type,
      spacing = if (has_spacing) par[1],
      deltas = if (np > has_spacing) par[(1 + has_spacing):np]
    )
    dlikert_cumulative(grid$y[i], grid$act[i], thr, dist = .noise_dist_names[grid$dist[i]], log = TRUE)
  }, numeric(1))

  expect_true(all(is.finite(stan)))
  # per element: expect_equal()'s tolerance is a mean relative difference,
  # which the tail values (|log p| up to 6e4) dilute
  expect_lt(max(abs(stan - r) / pmax(1, abs(r))), 1e-10)
})

test_that("update() rebuilds item and keying when the formula kept neither", {
  dat <- all_categories(likert_sim_frame())
  model <- likert_cumulative("y", "item", n_categories = 5, keying = "keying")
  fit <- suppressWarnings(suppressMessages(bmm(bmf(activation ~ 1 + (1 | id)), dat, model,
                                               backend = "mock", mock_fit = 1, rename = FALSE)))
  expect_false(any(c("item", "keying", "y") %in% names(fit$data)))
  reverted <- revert_check_data(fit$bmm$model, fit$data)
  expect_identical(reverted$y, dat$y)
  expect_identical(as.integer(reverted$item), as.integer(factor(dat$item)))
  expect_identical(reverted$keying, as.integer(dat$keying))
  expect_setequal(attr(reverted, "rebuilt"), c("item", "keying"))
  rechecked <- suppressWarnings(check_stored_data(fit$bmm$model, fit$data, fit$bmm$user_formula))
  expect_equal(brms::standata(fit, newdata = rechecked), brms::standata(fit))
})
