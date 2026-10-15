# A throwaway model registered as if another package defined it: a normal
# response with two parameters on brms' gaussian family, enough to run the
# pipeline on the mock backend. It has no model-case provider on purpose: the
# shared tests enumerate bmm_builtin_models(), so a registered model must not
# break them.
.model_regtest <- function(y = NULL, links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = list(y = y),
      other_vars = list(),
      domain = "Registry tests",
      task = "Testing the registry",
      name = "Registry test model",
      citation = "Test, A. (2026). A model for the registry tests.",
      version = "NA",
      requirements = "",
      parameters = list(mu = "Mean of the response", sigma = "Residual standard deviation"),
      links = list(mu = "identity", sigma = "log"),
      fixed_parameters = list(),
      default_priors = list(
        mu = list(main = "normal(0, 1)", effects = "normal(0, 0.5)"),
        sigma = list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
      )
    ),
    class = c("bmmodel", "regtest"),
    call = call
  )
  set_links(out, links)
}

# no stop_missing_args(): it looks the alias up by name, and a function
# defined in a test file is not on the search path
regtest <- function(y, links = NULL, ...) {
  .model_regtest(y = y, links = links, call = match.call(), ...)
}

configure_model_regtest <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)
  formula$family <- brms::brmsfamily("gaussian", link = model$links$mu, link_sigma = model$links$sigma)
  nlist(formula, data)
}

# registers the model and its configure_model() method the way an extension's
# .onLoad() and NAMESPACE would, and removes both when the calling test ends
local_regtest <- function(..., env = parent.frame()) {
  register_bmm_model("regtest", .model_regtest, regtest,
                     package = "bmmregtest", version = "0.1.0", ...)
  registerS3method("configure_model", "regtest", configure_model_regtest,
                   envir = asNamespace("bmm"))
  withr::defer({
    rm(list = "regtest", envir = .bmm_registry)
    rm(list = "configure_model.regtest", envir = asNamespace("bmm")[[".__S3MethodsTable__."]])
  }, envir = env)
}

regtest_data <- function() {
  withr::local_seed(442)
  data.frame(y = stats::rnorm(40), x = rep(c("a", "b"), 20), id = rep(1:2, each = 20))
}

test_that("a registered model is listed with its package, after bmm's own models", {
  local_regtest()
  expect_true("regtest" %in% bmm_models())
  expect_false("regtest" %in% bmm_builtin_models())
  expect_identical(head(model_names(), length(bmm_builtin_models())), bmm_builtin_models())
  expect_true("regtest" %in% model_names())
  expect_identical(get_model("regtest"), .model_regtest)
  expect_identical(get_model2("regtest"), regtest)
  expect_identical(model_package("regtest"), "bmmregtest")
  expect_identical(model_package("imm"), "bmm")
  expect_null(registry_entry(NA_character_))

  registry <- model_registry()
  row <- registry[registry$model == "regtest", ]
  expect_identical(row$package, "bmmregtest")
  expect_identical(row$url, NA_character_)
  expect_identical(registry$url[registry$model == "imm"], "https://venpopov.com/bmm/reference/imm.html")

  printed <- capture.output(print(bmm_models()))
  expect_true("- regtest() (bmmregtest): Registry test model" %in% printed)
  expect_true("- mixture2p(): Two-parameter mixture model by Zhang and Luck (2008)" %in% printed)
  md <- capture.output(print_pretty_models_md())
  expect_true("- `regtest()` (bmmregtest): Registry test model" %in% md)
})

test_that("a registered model with a url links to it in the markdown list", {
  local_regtest(url = "https://example.org/regtest.html")
  md <- capture.output(print_pretty_models_md(group = "Registry tests"))
  expect_identical(md[1], "- [`regtest()`](https://example.org/regtest.html) (bmmregtest): Registry test model")
  overview <- model_overview(group = "Registry tests")
  expect_match(overview$Model, "[`regtest()`](https://example.org/regtest.html)", fixed = TRUE)
})

test_that("a registered model fits on the mock backend and the fit records its package", {
  local_regtest()
  fit <- bmm(bmf(mu ~ x, sigma ~ 1), regtest_data(), regtest("y"),
             backend = "mock", mock_fit = 1, rename = FALSE)
  expect_s3_class(fit, "bmmfit")
  expect_identical(fit$bmm$package, "bmmregtest")
  expect_identical(fit$version$bmmregtest, package_version("0.1.0"))
  expect_identical(fit$version$bmm, utils::packageVersion("bmm"))
  expect_match(
    fit$model,
    paste0("^// generated with brms [^\n]* and bmm ", utils::packageVersion("bmm"), " and bmmregtest 0\\.1\\.0\n")
  )
})

test_that("a fit of a built-in model records bmm as its package and no other version", {
  fit <- stored_frame_fit(model_cases_mixture2p()$stored_frame$mixture2p)
  expect_identical(fit$bmm$package, "bmm")
  expect_false(any(c("bmmregtest", "bmmspike") %in% names(fit$version)))
  expect_match(fit$model, paste0(" and bmm ", utils::packageVersion("bmm"), "\n"), fixed = TRUE)
})

test_that("the extractors and bmm_data_check() accept a registered model", {
  local_regtest()
  formula <- bmf(mu ~ x + (1 | id), sigma ~ 1)
  data <- regtest_data()

  prior <- default_prior(formula, data, regtest("y"))
  expect_s3_class(prior, "brmsprior")
  expect_true(any(prior$prior == "normal(0, 1)" & prior$class == "Intercept"))

  code <- stancode(formula, data, regtest("y"))
  expect_match(code, "and bmmregtest 0.1.0", fixed = TRUE)
  expect_match(code, "normal_lpdf(Y | mu, sigma)", fixed = TRUE)

  sdata <- standata(formula, data, regtest("y"))
  expect_equal(sdata$N, nrow(data))
  expect_equal(as.numeric(sdata$Y), data$y)

  expect_no_error(suppressMessages(bmm_data_check(formula, data, regtest("y")))) # findings print as messages
})

test_that("register_bmm_model() refuses bmm's names and another package's model", {
  local_regtest()
  expect_error(
    register_bmm_model("regtest", .model_regtest, regtest, package = "othpkg", version = "1.0"),
    "already registered by package 'bmmregtest'"
  )
  expect_identical(model_package("regtest"), "bmmregtest")
  # the same package registering again replaces its entry, as a reload does
  expect_no_error(
    register_bmm_model("regtest", .model_regtest, regtest, package = "bmmregtest", version = "0.2.0")
  )
  expect_identical(registry_entry("regtest")$version, package_version("0.2.0"))
  register_bmm_model("regtest", .model_regtest, regtest,
                     package = "othpkg", version = "1.0", overwrite = TRUE)
  expect_identical(model_package("regtest"), "othpkg")

  expect_error(
    register_bmm_model("imm", .model_regtest, regtest, package = "othpkg", version = "1.0"),
    "model bmm defines itself"
  )
  expect_error(
    register_bmm_model("regtest", .model_regtest, regtest, package = "bmm"),
    "`package` must be the name"
  )
  expect_error(register_bmm_model("regtest", .model_regtest, package = "othpkg"), "alias")
})

test_that("register_bmm_model() checks that the constructor builds the named model", {
  needs_args <- function(y) y
  expect_error(
    register_bmm_model("regtest2", needs_args, regtest, package = "othpkg", version = "1.0"),
    "must run without arguments.*\"y\" is missing"
  )
  expect_error(
    register_bmm_model("regtest2", .model_regtest, regtest, package = "othpkg", version = "1.0"),
    "class vector ends in 'regtest2'.*ends in 'regtest'"
  )
  expect_error(
    register_bmm_model("regtest2", function() list(), regtest, package = "othpkg", version = "1.0"),
    "must return a `bmmodel`"
  )
  expect_error(
    register_bmm_model("regtest", .model_regtest, regtest, package = "othpkg", version = "x.y"),
    "invalid version"
  )
  expect_false("regtest2" %in% model_names())
})

test_that("check_model() says how to register a bmmodel bmm does not know, once", {
  err <- expect_error(check_model(regtest("y")), "register_bmm_model")
  expect_match(conditionMessage(err), "'bmmodel', 'regtest'", fixed = TRUE)
  expect_length(gregexpr("not known to bmm", conditionMessage(err))[[1]], 1)
  expect_error(check_model(structure(list(), class = "foo")), "object of class 'foo'")
})

test_that("add_bmm_version_to_stancode() names a registered model's package once", {
  local_regtest()
  body <- "\nfunctions {\n}\n"
  expected <- paste0(
    "// generated with brms 2.23.0 and bmm ", utils::packageVersion("bmm"), " and bmmregtest 0.1.0", body
  )
  once <- add_bmm_version_to_stancode(paste0("// generated with brms 2.23.0", body), "regtest")
  expect_identical(once, expected)
  expect_identical(add_bmm_version_to_stancode(once, "regtest"), expected)
  expect_identical(
    add_bmm_version_to_stancode(
      paste0("// generated with brms 2.23.0 and bmm 0.9.0 and bmmregtest 0.0.9", body), "regtest"
    ),
    expected
  )
  expect_identical(
    add_bmm_version_to_stancode(paste0("// generated with brms 2.23.0", body), "imm"),
    add_bmm_version_to_stancode(paste0("// generated with brms 2.23.0", body))
  )
})
