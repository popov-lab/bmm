# the cases the shared tests run for lba; see helper-model-cases.R
model_cases_lba <- function() {
  lba_data <- data.frame(
    rt = c(0.6, 0.8, 1.1, 0.7, 0.9, 0.65), response = c(1, 3, 1, 2, 1, 1),
    choice = c("left", "right", "left", "right", "right", "left"),
    n_left = c(1L, 1L, 2L, 2L, 1L, 2L), n_right = c(2L, 1L, 1L, 2L, 1L, 1L),
    cond = factor(rep(c("a", "b"), 3))
  )
  custom_formula <- bmf(left ~ 1, right ~ 1 + cond, gap ~ 1, sp ~ 1, ndt ~ 1)
  sim <- withr::with_seed(1, rlba(20, drift = c(2, 1), gap = 0.5, sp = 0.5, ndt = 0.2))
  sim$label <- c("a", "z")[sim$response]
  observables_custom <- bmf(a ~ 1, z ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)

  list(
    # the frame keeps only the vint() columns: the response and the
    # accumulator columns come back from them
    stored_frame = list(
      lba = list(
        model = lba("rt", "response", n_choices = 3),
        formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1),
        data = lba_data
      ),
      # update() works on the fit's model, which check_model() has given the
      # category names of the formula
      lba_custom = list(
        model = check_model(
          lba("rt", "choice", version = "custom",
              accumulators = c(left = "n_left", right = "n_right")),
          lba_data, custom_formula
        ),
        formula = custom_formula, data = lba_data
      )
    ),
    check_data = list(
      lba = list(
        model = lba("rt", "response", n_choices = 2),
        data = data.frame(rt = rep(c(0.6, 0.8), 25), response = rep(1:2, 25))
      )
    ),
    observables = list(
      lba = list(
        model = lba(rt = "rt", response = "response", n_choices = 2),
        data = sim, formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
      ),
      lba_custom = list(
        model = check_model(lba(rt = "rt", response = "label", version = "custom"),
                            NULL, observables_custom),
        data = sim, formula = observables_custom
      )
    )
  )
}
