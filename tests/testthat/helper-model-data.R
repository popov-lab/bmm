# Small data sets shared by the model-case providers (helper-model-<model>.R)
# and the tests. A mock fit runs on at most two participants or 200 rows
# unless it is about data size (ADR 0004, R3).

# the first n participants of a package data set
small_data <- function(data = oberauer_lin_2017, n = 2, id = "ID") {
  data[data[[id]] %in% utils::head(unique(data[[id]]), n), ]
}

toy_rt_data <- function() {
  data.frame(
    rt = rep(c(0.6, 0.8, 1.1, 0.7), 5),
    response = rep(c(1, 0), 10),
    id = factor(rep(1:5, each = 4))
  )
}

toy_ezdm_data <- function() {
  data.frame(
    mean_rt = rep(c(0.5, 0.6), 10), var_rt = rep(c(0.02, 0.03), 10),
    n_upper = rep(c(60, 70), 10), n_trials = 100, id = factor(rep(1:10, each = 2))
  )
}

# with or without the Know/Guess split and with a column prefix, which
# check_data() infers from the columns alone
cdp_data <- function(guess, prefix = "rk") {
  cols <- .sdt_cdp_response_cols(2, 2, guess, prefix)
  counts <- matrix(c(30L, 20L, 5L, 7L, 9L, 11L, 4L, 14L), 12, 8, byrow = TRUE)
  cbind(
    data.frame(stimulus = rep(0:1, 6), id = factor(rep(1:6, each = 2))),
    stats::setNames(as.data.frame(counts[, seq_along(cols)]), cols)
  )
}

toy_rating_data <- function() {
  data.frame(
    stimulus = rep(c(0L, 1L), 4), id = factor(rep(1:4, each = 2)),
    r1 = c(30, 8, 26, 10, 33, 6, 28, 9), r2 = c(25, 12, 27, 14, 22, 11, 24, 13),
    r3 = c(20, 15, 21, 16, 19, 17, 22, 15), r4 = c(15, 25, 16, 22, 17, 28, 14, 26),
    r5 = c(10, 40, 10, 38, 9, 38, 12, 37)
  )
}

toy_mafc_data <- function() {
  data.frame(
    n_correct = c(80, 55, 78, 60, 85, 52, 81, 58), n_trials = 100,
    n_afc = rep(c(2, 4), 4), cond = factor(rep(c("a", "b"), each = 4))
  )
}

toy_mpt_data <- function() {
  data <- data.frame(
    id = factor(rep(1:5, 2)), item_type = factor(rep(c("old", "new"), each = 5)),
    old = c(40:44, 8:12)
  )
  data$new <- 50L - data$old
  data
}
