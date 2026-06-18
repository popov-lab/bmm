# Tests for plot.bmm_sdt_roc() and plot.bmm_sdt_auc(). The ROC/AUC objects are
# built with mocked posterior draws (see helper-sdt-analysis.R); the plot
# methods themselves call no brms functions.

test_that("plot.bmm_sdt_roc() returns a ggplot for binary single-criterion", {
  skip_if_not_installed("ggplot2")
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  p <- plot(roc_sdt(fit, n_points = 15))
  expect_s3_class(p, "ggplot")
})

test_that("plot.bmm_sdt_roc() overlays operating points for binary multi-criteria", {
  skip_if_not_installed("ggplot2")
  fit <- fake_binary_fit(uv = TRUE, multi = TRUE)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(
      dprime = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8), sdratio = log(1.3))),
    ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
    variables = function(...) c("bsp_sdratio"),
    .package = "brms"
  )
  roc <- roc_sdt(fit)
  obs <- roc_observed(fit, conditions = "condition")
  p <- plot(roc, observed = obs)
  expect_s3_class(p, "ggplot")
  # ribbon + line (smooth curve) + 2 error bars + point (model) + point (observed)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_true("GeomRibbon" %in% geoms)
  expect_true("GeomPoint" %in% geoms)
})

test_that("plot.bmm_sdt_roc() returns a ggplot for a rating fit", {
  skip_if_not_installed("ggplot2")
  fit <- fake_rating_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  p <- plot(roc_sdt(fit))
  expect_s3_class(p, "ggplot")
})

test_that("plot.bmm_sdt_auc() returns a ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  p <- plot(auc_sdt(fit))
  expect_s3_class(p, "ggplot")
})


############################################################################# !
# Z-TRANSFORMED ROC                                                      ####
############################################################################# !

test_that("plot.bmm_sdt_roc(scale = 'z') drops endpoints and stays finite", {
  skip_if_not_installed("ggplot2")
  fit <- fake_binary_fit(uv = TRUE, multi = TRUE)
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(
      dprime = 1.2, criterion = c(-0.8, -0.3, 0, 0.3, 0.8), sdratio = log(1.3))),
    ranef = function(...) list(id = array(0, dim = c(1, 1, 1))),
    variables = function(...) c("bsp_sdratio"),
    .package = "brms"
  )
  roc <- roc_sdt(fit)
  obs <- roc_observed(fit, conditions = "condition")
  p   <- plot(roc, observed = obs, scale = "z")
  expect_s3_class(p, "ggplot")

  built <- ggplot2::ggplot_build(p)$data
  coord_cols <- c("x", "y", "xmin", "xmax", "ymin", "ymax")
  finite_ok <- vapply(built, function(d) {
    vals <- unlist(d[intersect(coord_cols, names(d))])
    all(is.finite(vals))
  }, logical(1))
  expect_true(all(finite_ok))
})

test_that("plot.bmm_sdt_roc(scale = 'z') returns a ggplot for a rating fit", {
  skip_if_not_installed("ggplot2")
  fit <- fake_rating_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  p <- plot(roc_sdt(fit), scale = "z")
  expect_s3_class(p, "ggplot")
})


############################################################################# !
# LATENT DISTRIBUTIONS                                                   ####
############################################################################# !

test_that("plot.bmm_sdt_latent() draws densities and boundary lines (binary)", {
  skip_if_not_installed("ggplot2")
  fit <- fake_binary_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0.3)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  p <- plot(latent_sdt(fit, n_grid = 100L))
  expect_s3_class(p, "ggplot")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_true(all(c("GeomArea", "GeomLine", "GeomVline") %in% geoms))
})

test_that("plot.bmm_sdt_latent() returns a ggplot for a rating fit", {
  skip_if_not_installed("ggplot2")
  fit <- fake_rating_fit()
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(dprime = 1.5, criterion = 0, spacing = 0)),
    ranef = function(...) list(),
    variables = function(...) character(0),
    .package = "brms"
  )
  p <- plot(latent_sdt(fit, n_grid = 80L))
  expect_s3_class(p, "ggplot")
})
