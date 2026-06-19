# Shared fixtures for the SDT post-processing tests (test-helpers-analysis.R,
# test-plot.R, test-pp_check.R). Mock/structural only: fake bmmfit objects plus
# a brms::posterior_linpred stand-in, so ROC/AUC/pp_check logic is exercised
# without fitting a Stan model.

n_draws_mock <- 40L

# brms::posterior_linpred stand-in: deterministic per-parameter draws, each
# condition column constant at the supplied value.
mock_linpred_factory <- function(draws) {
  function(object, dpar = NULL, nlpar = NULL, newdata = NULL, ...) {
    par <- if (!is.null(dpar)) dpar else nlpar
    n_cond <- if (is.null(newdata)) 1L else nrow(newdata)
    val <- rep_len(draws[[par]], n_cond)
    matrix(rep(val, each = n_draws_mock), nrow = n_draws_mock, ncol = n_cond)
  }
}

fake_binary_fit <- function(uv = FALSE, multi = FALSE) {
  model <- sdt_binary(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  if (uv) model$fixed_parameters$sdratio <- NULL

  if (multi) {
    data <- data.frame(
      condition = factor(rep(paste0("br", 1:5), each = 2L), levels = paste0("br", 1:5)),
      stimulus  = rep(0:1, 5L),
      n_old     = c(5, 40, 10, 55, 20, 70, 35, 85, 55, 95),
      n_trials  = 100L,
      dist_type = 1L
    )
    uf <- bmf(dprime ~ 1 + (1 | id), criterion ~ 0 + condition + (1 | id), sdratio ~ 1)
  } else {
    data <- data.frame(stimulus = 0:1, n_old = c(20L, 80L),
                       n_trials = 100L, dist_type = 1L)
    uf <- bmf(dprime ~ 1, criterion ~ 1)
  }
  structure(list(data = data, bmm = list(model = model, user_formula = uf)),
            class = c("bmmfit", "brmsfit"))
}

fake_rating_fit <- function(uv = FALSE, threshold_type = "parsimonious",
                            n_ratings = 6L) {
  resp <- paste0("r", seq_len(n_ratings))
  model <- sdt_rating(response = resp, stimulus = "stimulus",
                      threshold_type = threshold_type)
  if (uv) model$fixed_parameters$sdratio <- NULL

  Yn <- c(30, 25, 20, 13, 8, 4)
  Ys <- c(4, 8, 13, 20, 25, 30)
  data <- data.frame(stimulus = c(0, 1), nTrials = c(sum(Yn), sum(Ys)))
  data$Y <- rbind(Yn, Ys)
  uf <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  structure(list(data = data, bmm = list(model = model, user_formula = uf)),
            class = c("bmmfit", "brmsfit"))
}

fake_mafc_fit <- function(dist = "normal", m = 4L) {
  model <- sdt_mafc(response = "n_correct", n_trials = "n_trials", m = m, dist = dist)
  data <- data.frame(n_correct = 80L, n_trials = 100L, dist_type = 1L)
  uf <- bmf(dprime ~ 1)
  structure(list(data = data, bmm = list(model = model, user_formula = uf)),
            class = c("bmmfit", "brmsfit"))
}

fake_ranking_fit <- function(dist = "gumbel_min", m = 3L) {
  resp <- paste0("rank", seq_len(m))
  model <- sdt_ranking(response = resp, m = m, dist = dist)
  data <- data.frame(rank1 = 10L, rank2 = 6L, rank3 = 4L, nTrials = 20L)
  uf <- bmf(dprime ~ 1)
  structure(list(data = data, bmm = list(model = model, user_formula = uf)),
            class = c("bmmfit", "brmsfit"))
}

# Minimal multinomial bmmfit for pp_check dispatch tests: family() returns
# "multinomial" and the model class drives the rating default-group logic.
fake_multinomial_fit <- function(model_class, stimulus = "stimulus") {
  model <- structure(list(other_vars = list(stimulus = stimulus)),
                     class = c("bmmodel", model_class))
  structure(list(family = brms::multinomial(refcat = NA),
                 bmm = list(model = model)),
            class = c("bmmfit", "brmsfit"))
}
