# model_spec() is the contract companion packages read (ADR 0005). The names
# and column sets below are written out on purpose: reading them from a
# constant in R/ would make these tests agree with whatever that constant says.

spec_fields <- c(
  "name", "package", "version", "domain", "task", "citation", "parameters",
  "default_priors", "init_ranges", "resp_vars", "other_vars", "data_columns",
  "functions", "contract_version"
)
parameter_columns <- c(
  "parameter", "description", "link", "fixed", "value", "alternative_group", "softmax"
)

# every constructor, once per version, plus an mpt() with trees, which is the
# only way an mpt has parameters
spec_test_models <- function() {
  ns <- asNamespace("bmm")
  models <- list()
  for (constructor in ls(ns, pattern = "^\\.model_", all.names = TRUE)) {
    for (version in model_versions(sub("^\\.model_", "", constructor))) {
      label <- if (is.na(version)) constructor else paste0(constructor, "(version = '", version, "')")
      models[[label]] <- if (is.na(version)) {
        get(constructor, envir = ns)()
      } else {
        get(constructor, envir = ns)(version = version)
      }
    }
  }
  models[["mpt() with trees"]] <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  models
}

# jsonlite writes a function as its source and an environment as {}, so the
# round trip alone would not catch them
is_plain_data <- function(x) {
  if (is.list(x)) all(vapply(x, is_plain_data, TRUE)) else is.atomic(x) || is.null(x)
}

test_that("model_spec() returns the contract fields for every model and version", {
  # 0.7 s over all models and versions, past the CRAN-tier limit of R10
  skip_on_cran()
  for (label in names(models <- spec_test_models())) {
    model <- models[[label]]
    spec <- model_spec(model)

    expect_identical(names(spec), spec_fields, label = label)
    expect_identical(spec$package, "bmm", label = label)
    expect_identical(spec$contract_version, "0.9", label = label)
    expect_identical(names(spec$parameters), parameter_columns, label = label)
    expect_identical(names(spec$default_priors), c("parameter", "main", "effects", "sd"), label = label)
    expect_identical(names(spec$init_ranges), c("parameter", "lower", "upper"), label = label)
    expect_identical(names(spec$data_columns), c("argument", "role"), label = label)
    expect_identical(names(spec$functions), c("density", "generator"), label = label)

    pars <- spec$parameters
    # `link` is what bmmtools hands to link_transform(); a name-only check let
    # a wrong link pass
    for (par in names(model$links)) {
      expect_identical(sum(pars$parameter == par), 1L, label = paste(label, par))
      expect_identical(pars$link[pars$parameter == par], model$links[[par]],
                       label = paste(label, par, "link"))
    }
    expect_identical(pars$parameter, names(model$parameters) %||% character(0), label = label)
    for (par in intersect(names(model$fixed_parameters), pars$parameter)) {
      expect_true(pars$fixed[pars$parameter == par], label = paste(label, par, "fixed"))
      expect_identical(pars$value[pars$parameter == par], model$fixed_parameters[[par]],
                       label = paste(label, par, "value"))
    }
    expect_true(all(is.na(pars$value[!pars$fixed])), label = label)
    expect_identical(spec$default_priors$parameter, names(model$default_priors) %||% character(0),
                     label = label)
    # a reader joins init_ranges to parameters, so brms plumbing such as `mu`
    # must not appear in one and not the other
    expect_true(all(spec$init_ranges$parameter %in% pars$parameter),
                label = paste(label, "init_ranges within parameters"))
    # other_vars is "the arguments the model was created with"; constructor
    # bookkeeping would become a contract item under ADR 0005
    alias <- sub("^\\.model_", "", sub("\\(.*$", "", label))
    expect_true(all(names(spec$other_vars) %in% names(formals(get(alias, envir = asNamespace("bmm"))))),
                label = paste(label, "other_vars within the arguments of", alias))

    exports <- getNamespaceExports("bmm")
    for (fun in Filter(Negate(is.na), spec$functions)) {
      expect_in(fun, exports)
    }
    expect_true(is_plain_data(spec), label = paste(label, "holds only atomic vectors"))
  }
})

test_that("model_spec() states each model's density and generator", {
  functions <- function(model) unlist(model_spec(model)$functions)
  ratings <- paste0("r", 1:4)
  expect_identical(functions(sdt_rating(ratings, "s", version = "dpsdt")),
                   c(density = "dsdt_dpsdt", generator = "rsdt_dpsdt"))
  expect_identical(functions(sdt_rating(ratings, "s", version = "metad")),
                   c(density = "dsdt_metad", generator = "rsdt_metad"))
  expect_identical(functions(sdt_rating(ratings, "s")),
                   c(density = "dsdt_rating", generator = "rsdt_rating"))
  expect_identical(functions(ezdm("m", "v", "u", "n", version = "4par")),
                   c(density = "dezdm", generator = "rezdm"))
  expect_identical(functions(imm("y", "nt", "d", 4, version = "abc")),
                   c(density = "dimm", generator = "rimm"))
})

test_that("model_spec() gives the data columns of the version, not of the model", {
  columns <- function(model) model_spec(model)$data_columns
  expect_identical(columns(imm("y", "nt", "d", 4, version = "abc"))$argument,
                   c("resp_error", "nt_features", "set_size"))
  expect_identical(columns(imm("y", "nt", "d", 4, version = "full"))$argument,
                   c("resp_error", "nt_features", "nt_distances", "set_size"))
  ezdm_columns <- columns(ezdm("m", "v", "u", "n", version = "4par"))
  expect_identical(ezdm_columns$argument, c("mean_rt", "var_rt", "n_upper", "n_trials"))
  expect_identical(ezdm_columns$role[ezdm_columns$argument == "mean_rt"],
                   "mean response time in seconds, one column per boundary (upper, lower)")
})

test_that("model_spec() round-trips through JSON", {
  skip_if_not_installed("jsonlite")
  for (label in names(models <- spec_test_models())) {
    spec <- model_spec(models[[label]])
    json <- expect_no_error(jsonlite::toJSON(spec, auto_unbox = TRUE))
    back <- jsonlite::fromJSON(json)
    expect_identical(names(back), spec_fields, label = label)
    if (nrow(spec$parameters) > 0) {
      expect_identical(back$parameters$parameter, spec$parameters$parameter, label = label)
      expect_identical(back$parameters$link, spec$parameters$link, label = label)
    }
  }
})

test_that("model_spec() of mixture2p() is the documented specification", {
  spec <- model_spec(mixture2p(resp_error = "y"))
  expect_s3_class(spec, "bmm_model_spec")
  expect_identical(unclass(spec), list(
    name = "Two-Parameter Mixture Model",
    package = "bmm",
    version = NA_character_,
    domain = "Working memory",
    task = "Continuous reproduction",
    citation = paste(
      "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution representations",
      "in visual working memory. Nature, 453(7192), 233-235.",
      "https://doi.org/10.1038/nature06860"
    ),
    parameters = data.frame(
      parameter = c("mu1", "kappa", "thetat"),
      description = c(
        paste(
          "Location parameter of the von Mises distribution for memory responses",
          "(in radians). Fixed internally to 0 by default."
        ),
        "Concentration parameter of the von Mises distribution",
        "Mixture weight for target responses"
      ),
      link = c("tan_half", "log", "logit"),
      fixed = c(TRUE, FALSE, FALSE),
      value = c(0, NA, NA),
      alternative_group = NA_character_,
      softmax = FALSE
    ),
    default_priors = data.frame(
      parameter = c("mu1", "kappa", "thetat"),
      main = c("normal(0, 0.5)", "normal(2, 1)", "logistic(0, 1)"),
      effects = c("normal(0, 0.25)", "normal(0, 1)", "normal(0, 2.5)"),
      sd = c("exponential(4)", "exponential(1)", "exponential(1)")
    ),
    init_ranges = data.frame(
      parameter = c("mu1", "kappa", "thetat"),
      lower = c(-0.65, 3.8, 0.25),
      upper = c(0.65, 15, 0.75)
    ),
    resp_vars = list(resp_error = "y"),
    other_vars = setNames(list(), character(0)),
    data_columns = data.frame(
      argument = "resp_error",
      role = "response error relative to the target, in radians"
    ),
    functions = list(density = "dmixture2p", generator = "rmixture2p"),
    contract_version = "0.9"
  ))
})

test_that("model_spec() marks alternative groups and softmax parameters", {
  pars <- model_spec(sdt_yn("r", "s", "n"))$parameters
  expect_identical(pars$parameter, c("d", "criterion", "bias", "sdratio"))
  expect_identical(pars$alternative_group, c(NA, "criterion|bias", "criterion|bias", NA))

  pars <- model_spec(mixture3p("y", nt_features = "nt", set_size = 4))$parameters
  expect_identical(pars$softmax, c(FALSE, FALSE, TRUE, TRUE))
})

test_that("model_spec() returns empty tables where a model declares nothing", {
  # m3 builds its inits in create_initfun.m3(), and an mpt() without trees has
  # no parameters
  m3_spec <- model_spec(m3(c("corr", "other", "npl"), num_options = c(1, 4, 5), version = "ss"))
  expect_identical(m3_spec$init_ranges,
                   data.frame(parameter = character(0), lower = numeric(0), upper = numeric(0)))
  mpt_spec <- model_spec(.model_mpt())
  expect_identical(nrow(mpt_spec$parameters), 0L)
  expect_identical(nrow(mpt_spec$init_ranges), 0L)
})

test_that("model_spec() of a fit adds the bmm version the fit was made with", {
  # the fixture predates the `functions` field and the #554 vocabulary; both
  # are bmm's, not the fit's, so they come from the installed constructor,
  # while the parameters stay what the fit was made with
  fit <- load_fixture_fit("bmmfit_mixture3p.rds")
  spec <- model_spec(fit)
  expect_identical(names(spec), c(spec_fields, "bmm_version"))
  expect_identical(spec$bmm_version, "1.4.1.9000")
  expect_identical(spec$name, "Three-Parameter Mixture Model")
  expect_identical(spec$domain, "Working memory")
  expect_identical(spec$task, "Continuous reproduction")
  expect_identical(spec$functions, list(density = "dmixture3p", generator = "rmixture3p"))
  expect_identical(spec$parameters$parameter, c("mu1", "kappa", "thetat", "thetant"))
})

test_that("model_spec() lists only the constructor's arguments in other_vars", {
  # mpt() stores its indicator matrices and the raw simplex next to its
  # arguments; sdt_yn() a trial_level flag; sdt_cdp() integer codes of dist
  # and threshold_type
  expect_identical(names(model_spec(.model_mpt())$other_vars),
                   c("tree_id", "covariates", "trees", "simplex", "restrictions"))
  expect_identical(names(model_spec(sdt_yn("r", "s"))$other_vars), c("stimulus", "n_trials", "dist"))
  expect_false(any(c("dist_int", "thresh_type_int") %in% names(model_spec(.model_sdt_cdp())$other_vars)))
})

test_that("model_spec() lists init ranges only for the parameters it declares", {
  # ddm() draws an initial value for `mu` too, which parameters leaves out
  # as brms plumbing
  expect_identical(model_spec(ddm("rt", "resp"))$init_ranges$parameter, c("drift", "bound", "ndt", "zr"))
})

test_that("model_spec() gives the expressions of mpt() trees as strings", {
  spec <- model_spec(mpt(mpt_2htm_trees(), tree_id = "item_type"))
  expect_identical(spec$other_vars$trees$old$branches,
                   list(old = "D + (1 - D) * g", new = "(1 - D) * (1 - g)"))
})

test_that("printing a model specification shows the model and its parameters", {
  out <- capture.output(print(model_spec(sdt_yn("r", "s", "n"))))
  expect_match(out[1], "Signal Detection", fixed = TRUE)
  expect_true(any(grepl("^Package: +bmm$", out)))
  expect_true(any(grepl("^Version: +none$", out)))
  expect_true(any(grepl("criterion\\|bias", out)))
})
