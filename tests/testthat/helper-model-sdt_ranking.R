# the cases the shared tests run for sdt_ranking; see helper-model-cases.R
model_cases_sdt_ranking <- function() {
  ranks <- paste0("rank", 1:5)
  data <- meyer_grant_jakob_2025[as.integer(meyer_grant_jakob_2025$id) <= 4, ]
  unequal <- sdt_ranking(ranks, m = "set_size", dist = "normal")
  list(
    stored_frame = list(
      sdt_ranking = list(model = sdt_ranking(ranks, m = 5), formula = bmf(d ~ 1),
                         data = data[data$set_size == 5, ]),
      sdt_ranking_m_column = list(model = unequal, formula = bmf(d ~ 1, sdratio ~ 1), data = data),
      sdt_ranking_m_predictor = list(
        model = sdt_ranking(ranks, m = "set_size"),
        formula = bmf(d ~ 1 + set_size), data = data
      )
    ),
    sd_default = list(
      sdt_ranking = list(
        model = unequal, formula = bmf(d ~ 1 + (1 | id), sdratio ~ 1 + (1 | id)), data = data,
        sd = c(d = "exponential(1)", sdratio = "exponential(2)")
      )
    ),
    init = list(
      sdt_ranking = list(model = unequal, formula = bmf(d ~ 1 + (1 | id), sdratio ~ 1),
                         data = data)
    ),
    # a wide rank-frequency response, which the shared columns cannot express
    check_data = list(sdt_ranking = list(
      model = sdt_ranking(c("rank1", "rank2"), m = 2),
      data = data.frame(rank1 = rep(30L, 50), rank2 = rep(20L, 50))
    )),
    links = list(sdt_ranking = list(model = sdt_ranking(ranks, m = 5), settable = "d"))
  )
}
