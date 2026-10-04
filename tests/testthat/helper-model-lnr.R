# the cases the shared tests run for lnr; see helper-model-cases.R
model_cases_lnr <- function() {
  lnr_data <- data.frame(
    rt = rep(c(0.6, 0.8, 1.1, 0.7, 0.9), 4),
    choice = rep(1:4, 5),
    label = rep(c("correct", "similar", "other", "correct", "similar"), 4),
    n_correct = 1L, n_similar = rep(c(2L, 3L), 10), n_other = 4L,
    id = factor(rep(1:5, each = 4))
  )
  custom_formula <- bmf(correct ~ 1, similar ~ 1, other ~ 1, ndt ~ 1, s ~ 1)
  sim <- withr::with_seed(1, rlnr(20, m = c(-1, 0, 0), s = 0.5, ndt = 0.2))
  sim$label <- c("correct", "lure", "lure")[sim$response]
  observables_custom <- bmf(correct ~ 1, lure ~ 1)

  list(
    stored_frame = list(
      # every option is chosen, so a rebuilt response that lost which error was
      # made must not report options 3 and 4 as never chosen
      lnr = list(
        model = lnr("rt", "choice", n_choices = 4),
        formula = bmf(correct ~ 1, error ~ 1, ndt ~ 1, s ~ 1), data = lnr_data
      ),
      lnr_predictor = list(
        model = lnr("rt", "choice", n_choices = 4),
        formula = bmf(correct ~ 1 + id, error ~ 1, ndt ~ 1, s ~ 1), data = lnr_data
      ),
      # check_data() reads the categories that check_model() takes from the
      # formula, and update() passes check_data() the fit's checked model
      lnr_custom = list(
        model = check_model(lnr("rt", "label", version = "custom", accumulators = c(
          correct = "n_correct", similar = "n_similar", other = "n_other"
        )), formula = custom_formula),
        formula = custom_formula, data = lnr_data
      )
    ),
    check_data = list(
      lnr = list(
        model = lnr("rt", "response", n_choices = 2),
        data = data.frame(rt = rep(c(0.6, 0.8), 25), response = rep(1:2, 25))
      )
    ),
    observables = list(
      lnr = list(
        model = lnr(rt = "rt", response = "response", n_choices = 3),
        data = sim, formula = bmf(correct ~ 1)
      ),
      lnr_custom = list(
        model = check_model(lnr(rt = "rt", response = "label", version = "custom"),
                            formula = observables_custom),
        data = sim, formula = observables_custom
      )
    )
  )
}
