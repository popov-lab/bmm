 .rdm_count_loglik_r <- function(rt, response, drift, n, gap, ndt, s, sp) {
   t <- rt - ndt
   b <- gap + sp
   A <- sp
   K <- length(drift)

   if (A == 0) {
     log_lik <- log(n[response]) +
       .dwald(t, drift = drift[response], bound = b, s = s, log = TRUE)
   } else {
     log_lik <- log(n[response]) +
       .dwald_full(t, drift = drift[response], bound = b, A = A,
                   s = s, log = TRUE)
   }

   for (j in seq_len(K)) {
     log_surv <- if (A == 0) {
       .pwald(t, drift = drift[j], bound = b, s = s,
              lower.tail = FALSE, log.p = TRUE)
     } else {
       .pwald_full(t, drift = drift[j], bound = b, A = A, s = s,
                   lower.tail = FALSE, log.p = TRUE)
     }

     if (j == response) {
       if (n[j] > 1) {
         log_lik <- log_lik + (n[j] - 1) * log_surv
       }
     } else {
       log_lik <- log_lik + n[j] * log_surv
     }
   }

   log_lik
 }

 .rdm_expose_stan_functions <- local({
   cache <- new.env(parent = emptyenv())

   function(kind = c("simple_fixed", "simple_free", "custom_free")) {
     kind <- match.arg(kind)
     if (exists(kind, envir = cache, inherits = FALSE)) {
       return(invisible(TRUE))
     }

     spec <- switch(
       kind,
       simple_fixed = list(
         data = data.frame(rt = c(0.5, 0.6), response = c(1L, 2L)),
         model = rdm(rt = "rt", response = "response", n_alternatives = 2),
         formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
       ),
       simple_free = list(
         data = data.frame(rt = c(0.5, 0.6), response = c(1L, 2L)),
         model = rdm(rt = "rt", response = "response", n_alternatives = 2),
         formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, s ~ 1, sp ~ 1)
       ),
       custom_free = list(
         data = data.frame(
           rt = c(0.5, 0.6, 0.7),
           resp = c("corr", "lure", "npl")
         ),
         model = rdm(rt = "rt", response = "resp", version = "custom"),
         formula = bmf(corr ~ 1, lure ~ 1, npl ~ 1, gap ~ 1, ndt ~ 1, s ~ 1, sp ~ 1)
       )
     )

     code <- suppressWarnings(
       stancode(spec$formula, spec$data, spec$model, backend = "rstan")
     )
     sm <- rstan::stan_model(model_code = code, verbose = FALSE)
     rstan::expose_stan_functions(sm)
     cache[[kind]] <- TRUE
     invisible(TRUE)
   }
 })

 test_that("RDM Stan fixed-sp wrapper matches R likelihood for simple grids", {
   skip_on_cran()
   skip_if_not_installed("rstan")

   .rdm_expose_stan_functions("simple_fixed")

   scenarios <- data.frame(
     rt = c(0.42, 0.58, 0.73, 0.91),
     response = c(1L, 2L, 1L, 2L),
     driftc = c(3.2, 2.8, 4.1, 3.4),
     drifte = c(1.2, 1.1, 1.8, 1.5),
     gap = c(0.9, 1.1, 0.8, 1.3),
     ndt = c(0.19, 0.21, 0.24, 0.28),
     n2 = c(1L, 3L, 1L, 3L)
   )

   for (i in seq_len(nrow(scenarios))) {
     row <- scenarios[i, ]
     stan_val <- rdm_simple_lpdf(
       rt = row$rt,
       mu = 0,
       driftc = row$driftc,
       drifte = row$drifte,
       gap = row$gap,
       ndt = row$ndt,
       s = 1,
       sp = exp(-100),
       response = row$response,
       n1 = 1L,
       n2 = as.integer(row$n2)
     )
     r_val <- .rdm_count_loglik_r(
       rt = row$rt,
       response = row$response,
       drift = c(row$driftc, row$drifte),
       n = c(1L, as.integer(row$n2)),
       gap = row$gap,
       ndt = row$ndt,
       s = 1,
       sp = 0
     )

     expect_equal(stan_val, r_val, tolerance = 1e-8)
   }
 })

 test_that("RDM Stan free-sp wrapper matches R likelihood for simple grids", {
   skip_on_cran()
   skip_if_not_installed("rstan")

   .rdm_expose_stan_functions("simple_free")

   scenarios <- data.frame(
     rt = c(0.405, 0.512, 0.688, 0.845),
     response = c(1L, 2L, 1L, 2L),
     driftc = c(3.0, 2.6, 3.8, 3.1),
     drifte = c(1.0, 1.2, 1.5, 1.4),
     gap = c(0.9, 1.1, 1.0, 1.2),
     ndt = c(0.2, 0.22, 0.25, 0.29),
     s = c(0.9, 1.0, 1.1, 0.95),
     sp = c(0.02, 0.03, 0.05, 0.04),
     n2 = c(1L, 3L, 1L, 3L)
   )

   for (i in seq_len(nrow(scenarios))) {
     row <- scenarios[i, ]
     stan_val <- rdm_simple_lpdf(
       rt = row$rt,
       mu = 0,
       driftc = row$driftc,
       drifte = row$drifte,
       gap = row$gap,
       ndt = row$ndt,
       s = row$s,
       sp = row$sp,
       response = row$response,
       n1 = 1L,
       n2 = as.integer(row$n2)
     )
     r_val <- .rdm_count_loglik_r(
       rt = row$rt,
       response = row$response,
       drift = c(row$driftc, row$drifte),
       n = c(1L, as.integer(row$n2)),
       gap = row$gap,
       ndt = row$ndt,
       s = row$s,
       sp = row$sp
     )

     expect_equal(stan_val, r_val, tolerance = 1e-6)
   }
 })

 test_that("RDM Stan free-sp wrapper matches R likelihood for custom grids", {
   skip_on_cran()
   skip_if_not_installed("rstan")

   .rdm_expose_stan_functions("custom_free")

   scenarios <- list(
     list(
       rt = 0.49,
       response = 1L,
       drift = c(3.4, 1.5, 0.9),
       n = c(1L, 2L, 1L),
       gap = 1.0,
       ndt = 0.21,
       s = 1.0,
       sp = 0.03
     ),
     list(
       rt = 0.63,
       response = 2L,
       drift = c(2.8, 1.8, 1.2),
       n = c(2L, 1L, 3L),
       gap = 1.1,
       ndt = 0.24,
       s = 0.95,
       sp = 0.04
     ),
     list(
       rt = 0.77,
       response = 3L,
       drift = c(2.5, 1.6, 1.4),
       n = c(1L, 1L, 2L),
       gap = 0.9,
       ndt = 0.26,
       s = 1.05,
       sp = 0.02
     )
   )

   for (scenario in scenarios) {
     stan_val <- rdm_custom_lpdf(
       rt = scenario$rt,
       mu = 0,
       corr = scenario$drift[1],
       lure = scenario$drift[2],
       npl = scenario$drift[3],
       gap = scenario$gap,
       ndt = scenario$ndt,
       s = scenario$s,
       sp = scenario$sp,
       response = scenario$response,
       n1 = scenario$n[1],
       n2 = scenario$n[2],
       n3 = scenario$n[3]
     )
     r_val <- .rdm_count_loglik_r(
       rt = scenario$rt,
       response = scenario$response,
       drift = scenario$drift,
       n = scenario$n,
       gap = scenario$gap,
       ndt = scenario$ndt,
       s = scenario$s,
       sp = scenario$sp
     )

     expect_equal(stan_val, r_val, tolerance = 1e-6)
   }
 })

 test_that("RDM count likelihood stays finite near ndt across the support boundary", {
   rt <- seq(0.2002, 0.205, length.out = 7)
   sp_grid <- c(1e-8, 1e-5, 1e-3, 0.01, 0.05)

   vals <- outer(rt, sp_grid, Vectorize(function(r, sp) {
     .rdm_count_loglik_r(
       rt = r,
       response = 1L,
       drift = c(3, 1.5),
       n = c(1L, 3L),
       gap = 1,
       ndt = 0.2,
       s = 1,
       sp = sp
     )
   }))

   expect_true(all(is.finite(vals)))
 })

 test_that("RDM log likelihood remains finite under small parameter perturbations", {
   ndt_max <- 0.45
   base <- c(
     driftc = log(3.2),
     drifte = log(1.4),
     gap = log(0.95),
     ndt_raw = stats::qlogis(0.19 / ndt_max),
     s = log(0.9),
     sp = log(0.03)
   )
   eps <- 1e-5

   ll_fun <- function(theta) {
     .rdm_count_loglik_r(
       rt = 0.54,
       response = 1L,
       drift = exp(theta[c("driftc", "drifte")]),
       n = c(1L, 3L),
       gap = exp(theta[["gap"]]),
       ndt = ndt_max * stats::plogis(theta[["ndt_raw"]]),
       s = exp(theta[["s"]]),
       sp = exp(theta[["sp"]])
     )
   }

   for (par in names(base)) {
     theta_minus <- theta_plus <- base
     theta_minus[[par]] <- theta_minus[[par]] - eps
     theta_plus[[par]] <- theta_plus[[par]] + eps

     ll_minus <- ll_fun(theta_minus)
     ll_mid <- ll_fun(base)
     ll_plus <- ll_fun(theta_plus)

     expect_true(all(is.finite(c(ll_minus, ll_mid, ll_plus))))

     slope_left <- (ll_mid - ll_minus) / eps
     slope_right <- (ll_plus - ll_mid) / eps

     expect_true(is.finite(slope_left))
     expect_true(is.finite(slope_right))
     expect_lt(abs(slope_right - slope_left), 1e4)
   }
 })
