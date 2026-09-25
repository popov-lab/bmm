# =============================================================================
# Stage 2 Stan/R kernel parity: the four LBA single-accumulator kernels
# (inst/stan_chunks/lba_*_functions.stan) compiled through a fixed_param
# generated-quantities block and checked against the R twins
# (.lba_lpdf_single()/.lba_lsurv_single()), plus an independent log-space
# quadrature reference for the normal kernel. Pattern follows
# local/lba_grid_check.R (the audit's own verification harness).
#
# Tolerances below are measured on THIS grid (see the worktree's
# local/report/lba_349_fixes_2026-09-25.md audit table) rather than assumed:
# normal/gamma/lognormal and the low-ratio Frechet rows agree with the R
# kernels to <1e-6 nats; a branch-boundary corner of the normal kernel
# (t = 5, A = gap = 1e-3, where delta = A/(t s) sits just above the 1e-4
# midpoint-rule cutoff) and the Frechet 16-point rule at A/gap > 6 are both
# documented as genuinely less accurate, so:
#  * every comparison is restricted to rows where the R reference is above
#    -40 nats (below that a density is representationally zero and comparing
#    two independent floating-point paths to high precision is meaningless
#    for either implementation);
#  * the Frechet kernel is compared at 1e-6 nats for A/gap <= 6 (matching the
#    chunk's own docstring) and at 5 nats for A/gap > 6 (the chunk documents
#    up to ~0.84-3 nats of drift there; 5 nats keeps slack while still
#    catching a regression an order of magnitude worse than what is shipped).
# =============================================================================

.lba_kernel_stan_model <- local({
  mod <- NULL
  function() {
    if (is.null(mod)) {
      sc_dir <- system.file("stan_chunks", package = "bmm")
      rl <- function(f) paste(readLines(file.path(sc_dir, paste0(f, ".stan"))), collapse = "\n")
      funs <- paste(
        .lba_stan_code("lba_normal_simple", c("driftc", "drifte"), "normal"),
        rl("lba_lognormal_functions"), rl("lba_gamma_functions"), rl("lba_frechet_functions"),
        sep = "\n"
      )
      prog <- paste0(
        "functions {\n", funs, "\n}\n",
        "data {\n",
        "  int N;\n  vector[N] t;\n  vector[N] v;\n  vector[N] b;\n  vector[N] A;\n",
        "  vector[N] s;\n  array[N] int dist;\n",
        "}\n",
        "generated quantities {\n",
        "  vector[N] bmm_lpdf;\n  vector[N] bmm_lsurv;\n",
        "  for (i in 1:N) {\n",
        "    if (dist[i] == 1) {\n",
        "      bmm_lpdf[i] = lba_normal_single_lpdf(t[i] | v[i], b[i], A[i], s[i]);\n",
        "      bmm_lsurv[i] = lba_normal_single_lccdf(t[i] | v[i], b[i], A[i], s[i]);\n",
        "    } else if (dist[i] == 2) {\n",
        "      bmm_lpdf[i] = lba_gamma_single_lpdf(t[i] | v[i], b[i], A[i], s[i]);\n",
        "      bmm_lsurv[i] = lba_gamma_single_lccdf(t[i] | v[i], b[i], A[i], s[i]);\n",
        "    } else if (dist[i] == 3) {\n",
        "      bmm_lpdf[i] = lba_lognormal_single_lpdf(t[i] | v[i], b[i], A[i], s[i]);\n",
        "      bmm_lsurv[i] = lba_lognormal_single_lccdf(t[i] | v[i], b[i], A[i], s[i]);\n",
        "    } else {\n",
        "      bmm_lpdf[i] = lba_frechet_single_lpdf(t[i] | v[i], b[i], A[i], s[i]);\n",
        "      bmm_lsurv[i] = lba_frechet_single_lccdf(t[i] | v[i], b[i], A[i], s[i]);\n",
        "    }\n",
        "  }\n",
        "}\n"
      )
      mod <<- cmdstanr::cmdstan_model(
        cmdstanr::write_stan_file(prog, dir = tempdir()),
        compile = TRUE, quiet = TRUE
      )
    }
    mod
  }
})

.lba_kernel_grid <- function() {
  ts <- c(0.001, 0.02, 0.5, 5)
  As <- c(1e-3, 0.5, 2)
  gaps <- c(1e-3, 0.5)
  vs <- list(
    normal = c(-3, -0.3, 1, 3), gamma = c(0.5, 2, 4),
    lognormal = c(-1, 0.5), frechet = c(0.8, 3)
  )
  do.call(rbind, lapply(names(vs), function(dist) {
    g <- expand.grid(t = ts, v = vs[[dist]], A = As, gap = gaps,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    g$b <- g$gap + g$A
    g$s <- 1
    g$dist <- dist
    g
  }))
}

# One fixed_param sample of the whole grid, cached: the grid and R twins never
# change within a test run, so every test_that below reuses this single call.
.lba_kernel_results <- local({
  res <- NULL
  function() {
    if (is.null(res)) {
      G <- .lba_kernel_grid()
      dist_code <- c(normal = 1L, gamma = 2L, lognormal = 3L, frechet = 4L)
      mod <- .lba_kernel_stan_model()
      dat <- list(N = nrow(G), t = G$t, v = G$v, b = G$b, A = G$A, s = G$s,
                 dist = unname(dist_code[G$dist]))
      fit <- mod$sample(data = dat, fixed_param = TRUE, iter_sampling = 1,
                        chains = 1, refresh = 0, sig_figs = 18,
                        show_messages = FALSE, show_exceptions = FALSE)
      dr <- as.data.frame(suppressWarnings(fit$draws(format = "df")))
      get_col <- function(nm) as.numeric(dr[1, grep(paste0("^", nm, "\\["), names(dr))])
      G$stan_lpdf <- get_col("bmm_lpdf")
      G$stan_lsurv <- get_col("bmm_lsurv")

      G$r_lpdf <- NA_real_
      G$r_lsurv <- NA_real_
      for (dist in unique(G$dist)) {
        i <- G$dist == dist
        G$r_lpdf[i] <- .lba_lpdf_single(G$t[i], G$v[i], G$b[i], G$A[i], G$s[i], dist)
        G$r_lsurv[i] <- .lba_lsurv_single(G$t[i], G$v[i], G$b[i], G$A[i], G$s[i], dist)
      }
      res <<- G
    }
    res
  }
})

# ref > -40: below this the reference density is representationally zero
# (exp(-40) ~ 4e-18) and a precise match between two independent
# floating-point paths carries no information either implementation is wrong.
.lba_kernel_meaningful <- function(x) is.finite(x) & x > -40

test_that("Stan normal LBA kernel matches the R kernel on the grid", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  G <- .lba_kernel_results()
  d <- G[G$dist == "normal", ]
  ok_lpdf <- .lba_kernel_meaningful(d$r_lpdf)
  ok_lsurv <- .lba_kernel_meaningful(d$r_lsurv)

  expect_gt(sum(ok_lpdf), 0)
  expect_lt(max(abs(d$stan_lpdf[ok_lpdf] - d$r_lpdf[ok_lpdf])), 1e-6)
  expect_lt(max(abs(d$stan_lsurv[ok_lsurv] - d$r_lsurv[ok_lsurv])), 1e-6)
})

test_that("Stan gamma LBA kernel matches the R kernel on the grid", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  G <- .lba_kernel_results()
  d <- G[G$dist == "gamma", ]
  ok_lpdf <- .lba_kernel_meaningful(d$r_lpdf)
  ok_lsurv <- .lba_kernel_meaningful(d$r_lsurv)

  expect_gt(sum(ok_lpdf), 0)
  expect_lt(max(abs(d$stan_lpdf[ok_lpdf] - d$r_lpdf[ok_lpdf])), 1e-6)
  expect_lt(max(abs(d$stan_lsurv[ok_lsurv] - d$r_lsurv[ok_lsurv])), 1e-6)
})

test_that("Stan lognormal LBA kernel matches the R kernel on the grid", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  G <- .lba_kernel_results()
  d <- G[G$dist == "lognormal", ]
  ok_lpdf <- .lba_kernel_meaningful(d$r_lpdf)
  ok_lsurv <- .lba_kernel_meaningful(d$r_lsurv)

  expect_gt(sum(ok_lpdf), 0)
  expect_lt(max(abs(d$stan_lpdf[ok_lpdf] - d$r_lpdf[ok_lpdf])), 1e-6)
  expect_lt(max(abs(d$stan_lsurv[ok_lsurv] - d$r_lsurv[ok_lsurv])), 1e-6)
})

test_that("Stan frechet LBA kernel matches the R kernel where A/gap <= 6", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  G <- .lba_kernel_results()
  d <- G[G$dist == "frechet" & (G$A / G$gap) <= 6, ]
  ok_lpdf <- .lba_kernel_meaningful(d$r_lpdf)
  ok_lsurv <- .lba_kernel_meaningful(d$r_lsurv)

  expect_gt(sum(ok_lpdf), 0)
  expect_lt(max(abs(d$stan_lpdf[ok_lpdf] - d$r_lpdf[ok_lpdf])), 1e-6)
  expect_lt(max(abs(d$stan_lsurv[ok_lsurv] - d$r_lsurv[ok_lsurv])), 1e-6)
})

test_that("Stan frechet LBA kernel stays within its documented envelope where A/gap > 6", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  G <- .lba_kernel_results()
  d <- G[G$dist == "frechet" & (G$A / G$gap) > 6, ]
  ok_lpdf <- .lba_kernel_meaningful(d$r_lpdf)
  ok_lsurv <- .lba_kernel_meaningful(d$r_lsurv)

  expect_gt(sum(ok_lpdf), 0)
  expect_lt(max(abs(d$stan_lpdf[ok_lpdf] - d$r_lpdf[ok_lpdf])), 5)
  expect_lt(max(abs(d$stan_lsurv[ok_lsurv] - d$r_lsurv[ok_lsurv])), 5)
})

# Independent of both the closed-form M assembly (R) and its Stan twin: the
# LBA density written directly as an integral over the start point k in
# (0, A) rather than over the implied drift u = (b - k) / t. Restricted to
# t >= 0.02 because the 401-point pre-scan cannot resolve the sliver the
# integrand occupies at faster responses (see local/lba_grid_check.R).
test_that("Stan normal LBA lpdf matches an independent start-point quadrature", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  lquad <- function(lg, A) {
    lo <- 0
    hi <- A
    for (round in 1:3) {
      ks <- seq(lo, hi, length.out = 401)
      lv <- vapply(ks, lg, numeric(1))
      fin <- is.finite(lv)
      if (!any(fin)) return(-Inf)
      m <- max(lv[fin])
      keep <- which(fin & lv > m - 80)
      lo <- ks[max(1, min(keep) - 1)]
      hi <- ks[min(401, max(keep) + 1)]
    }
    val <- stats::integrate(function(k) exp(vapply(k, lg, numeric(1)) - m), lo, hi,
                            rel.tol = 1e-11, subdivisions = 2000L,
                            stop.on.error = FALSE)$value
    m + log(val)
  }
  quad_lpdf <- function(t, v, b, A, s) {
    lg <- function(k) {
      u <- (b - k) / t
      log((b - k) / t^2) + stats::dnorm(u, mean = v, sd = s, log = TRUE)
    }
    -log(A) + lquad(lg, A) - stats::pnorm(v / s, log.p = TRUE)
  }

  G <- .lba_kernel_results()
  d <- G[G$dist == "normal" & G$t >= 0.02, ]
  quad <- mapply(quad_lpdf, d$t, d$v, d$b, d$A, d$s)

  expect_gt(nrow(d), 0)
  expect_lt(max(abs(d$stan_lpdf - quad)), 1e-6)
})
