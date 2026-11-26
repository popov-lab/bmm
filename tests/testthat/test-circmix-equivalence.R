# The asset holds develop's brms::mixture() log-likelihood (commit 53fd780d) on
# a subset of oberauer_lin_2017 for five posterior draws per default model
# version, together with the inputs it was computed from. It is build-ignored
# (#438), so the test is skipped when it is absent.

circmix_reference_models <- list(
  mixture2p_simple = function() mixture2p(resp_error = "dev_rad"),
  mixture3p_simple = function() {
    mixture3p(
      resp_error = "dev_rad", nt_features = "col_nt", set_size = "set_size",
      regex = TRUE
    )
  },
  imm_full = function() {
    imm(
      resp_error = "dev_rad", nt_features = "col_nt", nt_distances = "dist_nt",
      set_size = "set_size", regex = TRUE, version = "full"
    )
  },
  imm_bsc = function() {
    imm(
      resp_error = "dev_rad", nt_features = "col_nt", nt_distances = "dist_nt",
      set_size = "set_size", regex = TRUE, version = "bsc"
    )
  },
  imm_abc = function() {
    imm(
      resp_error = "dev_rad", nt_features = "col_nt", set_size = "set_size",
      regex = TRUE, version = "abc"
    )
  }
)

# a brmsprep with exactly what log_lik_<family>() reads: the standata the
# stack builds from the frozen rows and the reference's parameters on the
# natural scale, as brms hands them over after the inverse link
circmix_reference_prep <- function(frozen, rows, model) {
  dpars <- lapply(frozen$dpars, function(cell_values) {
    cell_values[, rows$cell, drop = FALSE]
  })
  dpars$b <- matrix(1, nrow = nrow(frozen$log_lik), ncol = nrow(rows))
  free <- setdiff(names(model$parameters), names(model$fixed_parameters))
  formula <- do.call(bmf, lapply(free, function(p) {
    stats::as.formula(paste(p, "~ 1"))
  }))
  structure(
    list(
      data = standata(formula, rows, model),
      dpars = dpars[intersect(names(dpars), names(model$links))],
      ndraws = nrow(frozen$log_lik),
      family = list(vp_nodes = model$vp_nodes)
    ),
    class = "brmsprep"
  )
}

test_that("the circmix families reproduce develop's mixture log-likelihood", {
  asset <- test_path("assets", "circmix_reference_loglik.rds")
  skip_if_not(file.exists(asset), "frozen develop reference is not in this build")
  ref <- readRDS(asset)

  for (name in names(circmix_reference_models)) {
    frozen <- ref$models[[name]]
    model <- circmix_reference_models[[name]]()
    prep <- circmix_reference_prep(frozen, ref$rows, model)
    log_lik <- get(paste0("log_lik_", name))
    candidate <- vapply(
      seq_len(nrow(ref$rows)),
      function(i) log_lik(i, prep),
      numeric(prep$ndraws)
    )
    expect_lt(max(abs(candidate - frozen$log_lik)), 1e-8, label = name)
  }
})
