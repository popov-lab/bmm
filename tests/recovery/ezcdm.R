# Recovery plan for ezcdm() (the EZ circular diffusion model), three-parameter
# version: the drift angle is 0 because responses are coded as deviations from
# the target.
#
# The generator runs the model forward and shares no code with bmm. Each
# trial is a two-dimensional Wiener process with noise 1 per dimension that
# starts at the centre of a disk, drifts with speed `drift` towards angle 0
# and stops at the first step outside the circle of radius `bound`; the
# response angle is the angle of the exit point and the decision time is the
# time of the step. The non-decision time is added to the decision time. The
# trials of a cell are then summarised the way ezcdm() defines its data: the
# circular mean and the circular variance (1 - mean resultant length) of the
# angles, and the mean and the variance (denominator n - 1) of the response
# times, as ezcdm_summary_stats(method = "simple") returns them. bmm
# evaluates the closed-form moments of the same process (a / v) * I1 / I0 and
# a von Mises concentration a * v, plus a Gamma/normal approximation to the
# sampling distribution of the RT statistics, so an error in dezcdm(),
# rezcdm() or the Stan function does not cancel.
#
# Euler steps: dt is 0.001 s at the truth below, smaller when the drift is
# large and larger (up to 0.1 s) when the boundary is large, so that one step
# stays a small fraction of the boundary at any prior draw of the calibration.
# A path observed at discrete steps crosses the circle later than the
# continuous path; the boundary is therefore moved inward by 0.5826 * sqrt(dt)
# (the Broadie-Glasserman-Kou continuity correction, which holds for any
# diffusion with unit noise at the boundary).
ezcdm_euler <- function(drift, bound) {
  n <- length(drift)
  dt <- pmin(1e-3 * pmin(pmax(bound^2, 1), 100), 0.02 * bound / pmax(drift, 1e-3))
  reach <- bound - 0.5826 * sqrt(dt)
  x <- y <- elapsed <- numeric(n)
  angle <- time <- numeric(n)
  active <- seq_len(n)
  while (length(active) > 0) {
    m <- length(active)
    step <- sqrt(dt[active])
    x[active] <- x[active] + drift[active] * dt[active] + step * stats::rnorm(m)
    y[active] <- y[active] + step * stats::rnorm(m)
    elapsed[active] <- elapsed[active] + dt[active]
    out <- x[active]^2 + y[active]^2 >= reach[active]^2
    done <- active[out]
    angle[done] <- atan2(y[done], x[done])
    time[done] <- elapsed[done]
    active <- active[!out]
  }
  list(angle = angle, time = time)
}

# `group` numbers the cells, `n_per_cell` is the number of trials of each
ezcdm_cells <- function(drift, bound, ndt, n_per_cell) {
  group <- rep(seq_along(drift), n_per_cell)
  path <- ezcdm_euler(drift[group], bound[group])
  rt <- ndt[group] + path$time
  n <- as.vector(rowsum(rep(1, length(group)), group))
  mean_cos <- as.vector(rowsum(cos(path$angle), group)) / n
  mean_sin <- as.vector(rowsum(sin(path$angle), group)) / n
  mean_rt <- as.vector(rowsum(rt, group)) / n
  data.frame(
    mean_angle = atan2(mean_sin, mean_cos),
    var_angle = 1 - sqrt(mean_cos^2 + mean_sin^2),
    mean_rt = mean_rt,
    var_rt = as.vector(rowsum((rt - mean_rt[group])^2, group)) / (n - 1),
    n_trials = n_per_cell
  )
}

# one subject's row, for the recovery grid
ezcdm_generator <- function(pars, n_trials, model) {
  out <- ezcdm_cells(pars$drift, pars$bound, pars$ndt, n_trials)
  names(out) <- c(unlist(model$resp_vars), model$other_vars$n_trials)
  out
}

# the calibration's data sets: bmmtools hands over one prior draw under the
# fit's names (log links on drift, bound and ndt) and the design data, one row
# per subject; the statistics are drawn again from the subject values
ezcdm_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    exp(draws[[sprintf("b_%s_Intercept", par)]] +
          unname(draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))]))
  }
  out <- ezcdm_cells(subject("drift"), subject("bound"), subject("ndt"), data$n_trials)
  data[c("mean_angle", "var_angle", "mean_rt", "var_rt")] <-
    out[c("mean_angle", "var_angle", "mean_rt", "var_rt")]
  data
}

bmmtools::validation_plan(
  bmm::ezcdm(
    mean_angle = "mean_angle", var_angle = "var_angle",
    mean_rt = "mean_rt", var_rt = "var_rt", n_trials = "n_trials"
  ),
  # on the link scale (all three are logs). Default priors: drift
  # normal(0.5, 1), bound normal(0.5, 1), ndt normal(-1.5, 0.5). Drift 3.5
  # (log 1.25) is 0.75 SD above the prior centre, bound 1 (log 0) is 0.5 SD
  # below it, ndt 0.35 s (log -1.05) is 0.9 SD above it; the concentration
  # a * v = 3.5 gives a mean decision time of 0.24 s, so the RT statistics
  # carry real information on the three parameters
  pars = c(drift = log(3.5), bound = log(1), ndt = log(0.35)),
  # the SDs of the ezdm entry for drift and bound, a little less on ndt
  # because the RT variance gives it the least information
  sds = c(drift = 0.3, bound = 0.2, ndt = 0.15),
  # one row per subject (the statistics are per subject), so the grid's
  # n_trials is the number of trials summarised in each row; the default grid
  # (30, 100, 300) covers the range where the Gamma/normal approximation of
  # the RT statistics is worst (30) and good (300)
  generator = ezcdm_generator,
  sbc = list(
    generator = ezcdm_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
