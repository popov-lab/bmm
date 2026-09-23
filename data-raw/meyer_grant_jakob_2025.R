# Reshape the aggregated Meyer-Grant & Jakob (2025) ranking data into the wide
# format sdt_ranking() consumes: one rank-count column per rank position
# (rank1..rank5), with structural zeros where the rank exceeds the trial's set
# size. The set_size column is kept so the per-row set-size feature (`m =
# "set_size"`) can fit all set sizes jointly.
#
# Source: long aggregate (id, set_size, rank, count), one row per rank position.

library(dplyr)
library(tidyr)

long <- get(load("data/meyer_grant_jakob_2025.rda"))

rank_cols <- paste0("rank", seq_len(max(long$rank)))

meyer_grant_jakob_2025 <- long |>
  mutate(rank = paste0("rank", rank)) |>
  pivot_wider(names_from = rank, values_from = count, values_fill = 0) |>
  select(id, set_size, all_of(rank_cols)) |>
  arrange(id, set_size) |>
  as.data.frame()

usethis::use_data(meyer_grant_jakob_2025, overwrite = TRUE)
