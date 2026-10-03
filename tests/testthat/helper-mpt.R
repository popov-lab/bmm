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

mpt_impossible_trees <- function() {
  list(
    mpt_tree("withdist", list(
      corr = "Pm * Pb + (1 - Pm) * 0.2",
      dist = "(1 - Pm) * 0.2",
      npl = "Pm * (1 - Pb) + (1 - Pm) * 0.6"
    )),
    mpt_tree("nodist", list(
      corr = "Pm * Pb + (1 - Pm) * 0.2",
      npl = "Pm * (1 - Pb) + (1 - Pm) * 0.8"
    ), impossible = "dist")
  )
}

# two of the three conditions share the nodist branch structure, so the tree
# identifier column is coarser than the experimental factor
mpt_impossible_data <- function(n_id = 6) {
  dat <- expand.grid(
    id = factor(seq_len(n_id)), cond = c("withdist", "reord", "same"),
    stringsAsFactors = FALSE
  )
  dat$tree <- ifelse(dat$cond == "withdist", "withdist", "nodist")
  dat$corr <- 30L
  dat$npl <- 20L
  dat$dist <- ifelse(dat$cond == "withdist", 10L, 0L)
  dat
}
