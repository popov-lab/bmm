mpt_2htm_trees <- function() {
  list(
    mpt_tree("old", list(
      old = "D + (1 - D) * g",
      new = "(1 - D) * (1 - g)"
    )),
    mpt_tree("new", list(
      old = "(1 - D) * g",
      new = "D + (1 - D) * (1 - g)"
    ))
  )
}

mpt_2htm_data <- function(n_id = 10, n_items = 50) {
  dat <- expand.grid(
    id = factor(seq_len(n_id)), item_type = c("old", "new"),
    stringsAsFactors = FALSE
  )
  p_old <- ifelse(dat$item_type == "old", 0.85, 0.15)
  dat$old <- rbinom(nrow(dat), n_items, p_old)
  dat$new <- n_items - dat$old
  dat
}
