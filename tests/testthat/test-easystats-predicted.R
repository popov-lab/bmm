# insight::get_predicted() and marginaleffects' get_predict() return a model
# parameter on its native scale, from native_parameters(); modelbased's
# estimate_means() and estimate_contrasts() reach the fit through the second,
# estimate_relation() and the other estimate_*() functions through the first.

# get_predicted() centres the draws on their mean, as insight does for brms
# fits; marginaleffects takes the median
native_centres <- function(fit, newdata, par, re_formula = NULL, centre = stats::median) {
  np <- native_parameters(fit, newdata = newdata, pars = par, re_formula = re_formula)
  n_cells <- nrow(newdata)
  apply(matrix(np$value, nrow = n_cells, byrow = TRUE), 1, centre)
}

native_means <- function(fit, newdata, par, re_formula = NULL) {
  native_centres(fit, newdata, par, re_formula, centre = mean)
}

native_medians <- function(fit, newdata, par, re_formula = NULL) {
  native_centres(fit, newdata, par, re_formula)
}

test_that("get_predicted() returns the draws of a parameter on its native scale", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  grid <- data.frame(set_size = factor(1:4))
  pred <- insight::get_predicted(fit, data = grid, predict = "c", ci = 0.95)
  expect_s3_class(pred, "get_predicted")
  expect_identical(attr(pred, "predict"), "c")
  iterations <- attr(pred, "iterations")
  expect_equal(dim(iterations), c(4, brms::ndraws(fit)))
  expect_identical(names(iterations)[1], "iter_1")
  draws <- as.data.frame(fit, variable = paste0("b_c_set_size", 1:4))
  expect_equal(unname(as.numeric(pred)), unname(vapply(draws, function(d) mean(exp(d)), numeric(1))))
  expect_equal(unname(unlist(iterations[2, ])), exp(draws$b_c_set_size2))
  expect_named(attr(pred, "ci_data"), c("SE", "CI_low", "CI_high"))
})

test_that("get_predicted() predicts the first estimated parameter by default and names it", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  grid <- data.frame(set_size = factor(1:4))
  expect_message(pred <- insight::get_predicted(fit, data = grid, predict = NULL), "'c'")
  expect_equal(pred, insight::get_predicted(fit, data = grid, predict = "c"), ignore_attr = "call")
})

test_that("get_predicted() leaves insight's own prediction types to brms and refuses anything else", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  grid <- data.frame(set_size = factor(1:4))
  link <- insight::get_predicted(fit, data = grid, predict = "link")
  expect_equal(unname(unlist(attr(link, "iterations")[1, ])),
               unname(brms::posterior_linpred(fit, newdata = grid)[, 1]))
  expect_error(insight::get_predicted(fit, data = grid, predict = "thetat"), "'mu', 'c', 'kappa'")
})

test_that("the mixture weights are predicted after the joint softmax", {
  skip_if_not_installed("insight")
  skip_if_not_installed("modelbased")
  fit <- load_fixture_fit("bmmfit_mixture3p.rds")
  grid <- data.frame(set_size = factor(c(1, 2, 4)))
  thetant <- insight::get_predicted(fit, data = grid, predict = "thetant")
  # set size 1 has no non-target, so the model switches the component off
  expect_identical(unname(unlist(attr(thetant, "iterations")[1, ])), rep(0, brms::ndraws(fit)))
  relation <- suppressMessages(modelbased::estimate_relation(fit, by = "set_size", predict = "thetat"))
  expect_equal(relation$Predicted, native_means(fit, grid, "thetat"))
  expect_true(all(relation$Predicted > 0.5 & relation$Predicted < 1))
})

test_that("estimate_means() and estimate_contrasts() summarise the native draws", {
  skip_if_not_installed("modelbased")
  skip_if_not_installed("marginaleffects")
  fit <- load_fixture_fit("bmmfit_mixture3p.rds")
  grid <- data.frame(set_size = factor(c(1, 2, 4)))
  for (par in c("thetat", "thetant", "kappa")) {
    means <- suppressMessages(modelbased::estimate_means(fit, by = "set_size", predict = par))
    expect_identical(as.character(means$set_size), c("1", "2", "4"), info = par)
    # modelbased names the column of a parameter brms also has, such as kappa, after it
    expect_equal(means[[attr(means, "coef_name")]], native_medians(fit, grid, par), info = par)
  }
  contrasts <- suppressMessages(
    modelbased::estimate_contrasts(fit, contrast = "set_size", predict = "thetat")
  )
  np <- native_parameters(fit, newdata = grid, pars = "thetat")
  draws <- split(np$value, np$set_size)
  expect_equal(contrasts$Median[contrasts$Level1 == "4" & contrasts$Level2 == "1"],
               stats::median(draws[["4"]] - draws[["1"]]))
})

test_that("a grid without the grouping factors predicts the population level", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_mpt_ppcheck.rds")
  population <- data.frame(nTrials = 40)
  expect_equal(as.numeric(insight::get_predicted(fit, data = population, predict = "D")),
               native_means(fit, population, "D", re_formula = NA))
  participants <- data.frame(nTrials = 40, id = factor(1:2))
  expect_equal(as.numeric(insight::get_predicted(fit, data = participants, predict = "D")),
               native_means(fit, participants, "D"))
  expect_equal(as.numeric(insight::get_predicted(fit, data = participants, predict = "D",
                                                 include_random = FALSE)),
               native_means(fit, participants, "D", re_formula = NA))
})

test_that("random effects are included only for levels the fit knows", {
  skip_if_not_installed("insight")
  skip_if_not_installed("modelbased")
  fit <- load_fixture_fit("bmmfit_mpt_ppcheck.rds")
  # modelbased's grid holds the grouping factor at NA
  unknown <- data.frame(nTrials = 40, id = NA)
  expect_equal(as.numeric(insight::get_predicted(fit, data = unknown, predict = "D")),
               native_means(fit, unknown, "D", re_formula = NA))
  relation <- suppressMessages(modelbased::estimate_relation(fit, predict = "D"))
  # the grid crosses the response counts, which D does not depend on
  expect_equal(unique(relation$Predicted), native_means(fit, data.frame(nTrials = 40), "D", re_formula = NA))
  participants <- data.frame(nTrials = 40, id = factor(1:2))
  for (include_random in list(FALSE, NA)) {
    expect_equal(as.numeric(insight::get_predicted(fit, data = participants, predict = "D",
                                                   include_random = include_random)),
                 native_means(fit, participants, "D", re_formula = NA))
  }
})

test_that("marginaleffects' re_formula, ndraws and new levels reach the predictions", {
  skip_if_not_installed("marginaleffects")
  fit <- load_fixture_fit("bmmfit_mpt_ppcheck.rds")
  participant <- data.frame(nTrials = 40, id = factor(1))
  population <- suppressWarnings(
    marginaleffects::predictions(fit, newdata = participant, dpar = "D", re_formula = NA)
  )
  expect_equal(population$estimate, native_medians(fit, participant, "D", re_formula = NA))
  new <- suppressWarnings(marginaleffects::predictions(
    fit, newdata = data.frame(nTrials = 40, id = "new"), dpar = "D",
    allow_new_levels = TRUE, sample_new_levels = "gaussian"
  ))
  expect_true(new$estimate > 0 && new$estimate < 1)
  few <- suppressWarnings(marginaleffects::predictions(fit, newdata = participant, dpar = "D", ndraws = 10))
  expect_identical(ncol(marginaleffects::get_draws(few, "rvar")$rvar[[1]] |> posterior::draws_of() |> t()), 10L)
})

test_that("comparisons() contrasts the native draws", {
  skip_if_not_installed("marginaleffects")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  comparison <- suppressWarnings(marginaleffects::avg_comparisons(
    fit, variables = list(set_size = c("1", "4")), dpar = "c",
    newdata = data.frame(set_size = factor(1, levels = 1:4))
  ))
  np <- native_parameters(fit, newdata = data.frame(set_size = factor(c(1, 4))), pars = "c")
  draws <- split(np$value, np$set_size)
  expect_equal(comparison$estimate, stats::median(draws[["4"]] - draws[["1"]]))
})

test_that("get_predicted() predicts the first parameter of several and warns", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  grid <- data.frame(set_size = factor(1:4))
  expect_warning(pred <- insight::get_predicted(fit, data = grid, predict = c("c", "kappa")), "'c'")
  expect_identical(attr(pred, "predict"), "c")
})

test_that("get_predicted() takes newdata, type and centrality_function as insight does", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  draws <- as.data.frame(fit, variable = c("b_c_set_size1", "b_c_set_size4", "b_kappa_Intercept"))
  pred <- insight::get_predicted(fit, newdata = data.frame(set_size = factor(c(1, 4), levels = 1:4)),
                                 predict = "c")
  expect_equal(as.numeric(pred), c(mean(exp(draws$b_c_set_size1)), mean(exp(draws$b_c_set_size4))))
  grid <- data.frame(set_size = factor(1:4))
  kappa <- insight::get_predicted(fit, data = grid, predict = "c", type = "kappa")
  expect_identical(attr(kappa, "predict"), "kappa")
  expect_equal(as.numeric(kappa), rep(mean(exp(draws$b_kappa_Intercept)), 4))
  link <- insight::get_predicted(fit, data = grid, predict = NULL, type = "link")
  expect_equal(unname(unlist(attr(link, "iterations")[1, ])),
               unname(brms::posterior_linpred(fit, newdata = grid)[, 1]))
  centred <- insight::get_predicted(fit, data = grid[c(1, 4), , drop = FALSE], predict = "c",
                                    centrality_function = stats::median)
  expect_equal(as.numeric(centred), c(stats::median(exp(draws$b_c_set_size1)),
                                      stats::median(exp(draws$b_c_set_size4))))
})

test_that("get_predicted() predicts from the number of draws in iterations", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  draws <- exp(as.data.frame(fit, variable = "b_c_set_size2")$b_c_set_size2)
  withr::local_seed(11)
  pred <- insight::get_predicted(fit, data = data.frame(set_size = factor(2, levels = 1:4)),
                                 predict = "c", iterations = 10)
  withr::local_seed(11)
  ids <- sort(sample.int(length(draws), 10))
  expect_equal(unname(unlist(attr(pred, "iterations"))), draws[ids])
})

test_that("a grouping factor at NA predicts the population level even with new levels allowed", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_mpt_ppcheck.rds")
  population <- mean(stats::plogis(as.data.frame(fit, variable = "b_D_Intercept")$b_D_Intercept))
  pred <- insight::get_predicted(fit, data = data.frame(nTrials = 40, id = NA), predict = "D",
                                 allow_new_levels = TRUE)
  expect_equal(as.numeric(pred), population)
})

test_that("get_predicted() warns when it drops the random effects asked for", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_mpt_ppcheck.rds")
  expect_warning(insight::get_predicted(fit, data = data.frame(nTrials = 40), predict = "D",
                                        include_random = TRUE), "no column 'id'")
  expect_warning(insight::get_predicted(fit, data = data.frame(nTrials = 40, id = "new"), predict = "D",
                                        include_random = TRUE), "levels the fit does not know")
  expect_no_warning(insight::get_predicted(fit, data = data.frame(nTrials = 40), predict = "D"))
})
