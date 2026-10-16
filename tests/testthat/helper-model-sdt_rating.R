# the cases the shared tests run for sdt_rating; see helper-model-cases.R
model_cases_sdt_rating <- function() {
  responses <- paste0("r", 1:5)
  model <- sdt_rating(responses, "stimulus")
  data <- toy_rating_data()
  list(
    stored_frame = list(
      sdt_rating = list(
        model = model,
        formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1, sdratio ~ 1),
        data = data
      ),
      sdt_rating_deltas = list(
        model = sdt_rating(responses, "stimulus", threshold_type = "log_distance"),
        formula = bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1, delta3 ~ 1),
        data = data
      ),
      sdt_rating_dpsdt = list(
        model = sdt_rating(responses, "stimulus", version = "dpsdt"),
        formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1, Ro ~ 1),
        data = data
      ),
      # meta-d' needs an even number of categories
      sdt_rating_metad = list(
        model = sdt_rating(paste0("r", 1:4), "stimulus", version = "metad"),
        formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1, logmratio ~ 1 + (1 | id)),
        data = data[setdiff(names(data), "r5")]
      )
    ),
    sd_default = list(
      sdt_rating = list(
        model = model, data = data,
        formula = bmf(d ~ 1 + (1 | id), criterion ~ 1 + (1 | id),
                      spacing ~ 1 + (1 | id), sdratio ~ 1 + (1 | id)),
        sd = c(d = "exponential(1)", criterion = "exponential(1)",
               spacing = "exponential(2)", sdratio = "exponential(2)")
      )
    ),
    init = list(
      sdt_rating = list(model = model, data = data,
                        formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1))
    ),
    # one count column per rating category
    check_data = list(sdt_rating = list(model = model, data = data)),
    links = list(sdt_rating = list(model = model, settable = c("d", "criterion", "bias")))
  )
}
