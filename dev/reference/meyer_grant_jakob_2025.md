# Ranking recognition data from Meyer-Grant & Jakob (2025)

Ranking signal-detection data from 60 subjects. On each trial,
participants saw 3, 4, or 5 face images—one a studied target—and ranked
them by perceived oldness; the rank assigned to the target is recorded.
The set size varies across trials, so the data are aggregated to target
rank-frequency counts per subject and set size, in the wide format
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md)
consumes: one count column per rank position, with structural zeros
where the rank exceeds the trial's set size. Fitting all set sizes
jointly relies on the per-row set-size feature of
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md):
pass the set-size column to `m`.

## Usage

``` r
meyer_grant_jakob_2025
```

## Format

### `meyer_grant_jakob_2025`

A data frame with 180 rows (60 subjects x set sizes 3, 4, 5) and 7
columns:

- id:

  Factor with sequential codes `s01`–`s60` assigned by bmm; the original
  participant identifiers are not shipped

- set_size:

  Integer number of ranked items on the trial (3, 4, or 5)

- rank1, rank2, rank3, rank4, rank5:

  Integer number of trials in which the target received that rank (rank1
  = most likely target). Columns beyond `set_size` are structural zeros.
  The counts sum to 56 within each row.

## Source

Meyer-Grant, C. G., & Jakob, M. (2025). Ranking tasks in recognition
memory: A direct test of the two-high-threshold contrast model. *Journal
of Experimental Psychology: General*, 154(5), 1445–1455.
[doi:10.1037/xge0001700](https://doi.org/10.1037/xge0001700) . Data on
OSF: <https://osf.io/gtzu7/>.

## Examples

``` r
if (FALSE) { # \dontrun{
# Ranking SDT with set size varying per row: pass the set-size column to `m`
# so trials with 3, 4, and 5 alternatives are fit jointly.
model <- sdt_ranking(
  response = c("rank1", "rank2", "rank3", "rank4", "rank5"), m = "set_size"
)
fit <- bmm(
  formula = bmf(d ~ 1 + (1 | id)),
  data = meyer_grant_jakob_2025,
  model = model,
  backend = "cmdstanr"
)
} # }
```
