test_that("an unrecognized link target is refused, not appended", {
  expect_error(
    ddm(rt = "rt", response = "resp", links = list(zzz = "log")),
    "Unrecognized link target"
  )
  expect_error(
    m3(resp_cats = c("corr", "other"), num_options = c(1, 4), version = "ss",
       links = list(zzz = "log")),
    "Unrecognized link target"
  )
})

test_that("a link target one edit from a parameter is read as that parameter", {
  expect_warning(
    model <- ddm(rt = "rt", response = "resp", links = list(drif = "log")),
    "'drif' read as 'drift'"
  )
  expect_equal(model$links$drift, "log")
  expect_false("drif" %in% names(model$links))

  expect_warning(
    model <- ezdm("m", "v", "n", "t", links = list(bond = "softplus")),
    "'bond' read as 'bound'"
  )
  expect_equal(model$links$bound, "softplus")
})

test_that("a link target that names no single parameter is refused", {
  # 'b' is one edit from each of imm's a, c and s
  expect_error(
    imm(resp_error = "y", nt_features = "nt", nt_distances = "d",
        set_size = "ss", links = list(b = "log")),
    "Unrecognized link target"
  )
  # two entries resolving to the same parameter would silently drop one
  expect_error(
    ezdm("m", "v", "n", "t", links = list(bound = "log", boud = "softplus")),
    "name the same parameter"
  )
})

test_that("a link bmm does not implement is refused", {
  expect_error(
    ezdm("m", "v", "n", "t", links = list(bound = "logg")),
    "Unknown link function"
  )
  expect_error(
    m3(resp_cats = c("corr", "other"), num_options = c(1, 4),
       links = list(custom_par = "logg")),
    "Unknown link function"
  )
  expect_error(
    ddm(rt = "rt", response = "resp", links = stats::setNames(list("log"), "")),
    "must be a named list"
  )
})

test_that("a link the model cannot honour is refused at the constructor", {
  # brms writes no inverse-link code for loglog or softmax, so a custom family
  # built with either dies in stancode() with "argument is of length zero"
  expect_error(
    ddm(rt = "rt", response = "resp", links = list(bound = "loglog")),
    "Unknown link function\\(s\\): 'loglog'\\. ddm\\(\\) takes"
  )
  expect_error(
    ddm(rt = "rt", response = "resp", links = list(bound = "softmax")),
    "Unknown link function\\(s\\): 'softmax'\\. ddm\\(\\) takes"
  )
  expect_false(any(c("loglog", "softmax") %in% settable_link_functions(ddm("rt", "resp"))))

  # m3 substitutes the inverse link into its activation formulas, so it honours
  # inv_link()'s five rather than the ones a brms family can emit
  m3_ss <- function(links) {
    m3(resp_cats = c("corr", "other", "npl"),
       num_options = c("n_corr", "n_other", "n_npl"),
       choice_rule = "simple", version = "ss", links = links)
  }
  expect_error(m3_ss(list(c = "sqrt")), "'sqrt'\\. m3\\(\\) takes")
  expect_silent(m3_ss(list(c = "log")))
  expect_equal(
    settable_link_functions(m3_ss(NULL)),
    c("log", "softplus", "logit", "probit", "identity")
  )
})

test_that("a parameter fixed for scaling is refused with its own message", {
  expect_error(
    m3(resp_cats = c("corr", "other"), num_options = c(1, 4), version = "ss",
       links = list(b = "log")),
    "'b' has no link in m3\\(\\): the parameter is fixed for scaling"
  )
})

test_that("a model outside supported_models() is named by its class", {
  # what use_model_template() produces before the model is registered
  model <- structure(
    list(parameters = list(par1 = "", par2 = ""), links = list(par1 = "log")),
    class = c("bmmodel", "tmpl_demo")
  )
  expect_error(
    set_links(model, list(zzz = "log")),
    "tmpl_demo\\(\\) takes links for"
  )
})

test_that("a link wider than the model's default warns", {
  expect_warning(
    ezdm("m", "v", "n", "t", version = "4par", links = list(zr = "identity")),
    "allow values that the model's default"
  )
  # 1/eta is negative for every eta < 0, so inverse widens a positive parameter
  expect_warning(
    ddm(rt = "rt", response = "resp", links = list(bound = "inverse")),
    "allow values that the model's default"
  )
  # narrowing the range is what a link is for, and does not warn
  expect_silent(ezdm("m", "v", "n", "t", links = list(drift = "log")))
  expect_silent(ezdm("m", "v", "n", "t", links = list(bound = "softplus")))
})

test_that("a widening link warns once, not once per pipeline stage", {
  dat <- data.frame(rt = c(0.5, 0.6, 0.7), response = c(1, 0, 1))
  ff <- bmmformula(drift ~ 1, bound ~ 1, ndt ~ 1)
  model <- suppressWarnings(
    ddm(rt = "rt", response = "response", links = list(bound = "identity"))
  )
  expect_silent(check_links(model))
  expect_equal(check_model(model, dat, ff)$links$bound, "identity")
})

test_that("a link the pipeline fixed itself is not re-checked", {
  # resolve_fixed_links() swaps the link of a parameter that is fixed to a
  # constant, after set_links() recorded the estimation link, so the swap is
  # the pipeline's own doing and must not read as a user change
  model <- cswald(rt = "rt", response = "resp")
  model$parameters$p <- "a parameter the model can fix to a constant"
  model$links$p <- "log"
  attr(model, "links_default") <- model$links
  attr(model, "links_checked") <- model$links
  model$links_fixed <- list(p = "identity")
  model$links$p <- "identity"
  expect_silent(check_links(model))

  model$links_fixed <- NULL
  expect_warning(check_links(model), "allow values that the model's default")
})

test_that("a link the model does not pass on to the fit is refused", {
  expect_error(
    sdm(resp_error = "y", links = list(kappa = "softplus")),
    "cannot be changed in sdm\\(\\)"
  )
  expect_error(
    imm(resp_error = "y", nt_features = "nt", nt_distances = "d",
        set_size = "ss", links = list(c = "softplus")),
    "cannot be changed in imm\\(\\)"
  )
  expect_error(
    mixture2p(resp_error = "y", links = list(kappa = "softplus")),
    "cannot be changed in mixture2p\\(\\)"
  )
  expect_error(
    mixture3p(resp_error = "y", nt_features = "nt", set_size = "ss",
              links = list(kappa = "softplus")),
    "cannot be changed in mixture3p\\(\\)"
  )
  # naming the link the model already uses asks for no change
  expect_silent(sdm(resp_error = "y", links = list(kappa = "log")))
})

test_that("the links a model applies are exactly the settable ones", {
  # measured against the family each configure_model() builds: only these
  # models read model$links at fit time
  expect_equal(
    settable_links(ddm(rt = "rt", response = "resp")),
    c("drift", "bound", "ndt", "zr")
  )
  expect_equal(
    settable_links(ezdm("m", "v", "n", "t")),
    c("drift", "bound", "ndt", "s")
  )
  expect_equal(
    settable_links(cswald(rt = "rt", response = "resp")),
    names(cswald(rt = "rt", response = "resp")$links)
  )
  expect_null(settable_links(m3(resp_cats = c("a", "b"), num_options = c(1, 4))))
  expect_equal(settable_links(sdm(resp_error = "y")), character(0))
})

test_that("a custom link set on a model reaches the brms family", {
  dat <- data.frame(mean_rt = 0.6, var_rt = 0.03, n_upper = 80, n_trials = 100)
  ff <- bmmformula(drift ~ 1, bound ~ 1, ndt ~ 1)
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials",
    links = list(bound = "softplus")
  )
  family <- configure_model(model, check_data(model, dat, ff), ff)$formula$family
  expect_equal(family$link_bound, "softplus")
})

test_that("links set after construction are checked by the pipeline", {
  ff <- bmmformula(c ~ 1, kappa ~ 1)
  dat <- data.frame(y = c(0.1, -0.2, 0.3))

  model <- sdm(resp_error = "y")
  model$links$kappa <- "softplus"
  expect_error(check_model(model, dat, ff), "cannot be changed in sdm\\(\\)")

  model <- sdm(resp_error = "y")
  model$links$kapa <- "log"
  expect_warning(
    model <- check_model(model, dat, ff),
    "'kapa' read as 'kappa'"
  )
  expect_false("kapa" %in% names(model$links))

  # a typo that resolves to a parameter the user cannot set is reported as a
  # typo first, and is not offered back as a settable target
  model <- m3(resp_cats = c("corr", "other"), num_options = c(1, 4), version = "ss")
  model$links$bb <- "log"
  expect_warning(
    expect_error(
      check_model(model, NULL, bmf(c ~ 1, a ~ 1)),
      "'b' has no link in m3\\(\\).*Links can be set for 'c', 'a'"
    ),
    "'bb' read as 'b'"
  )

  # the documented m3 idiom of replacing the whole list keeps working
  model <- m3(resp_cats = c("corr", "other"), num_options = c(1, 4))
  model$links <- list(c = "log", a = "log")
  model <- suppressWarnings(
    check_model(model, NULL, bmf(corr ~ b + c, other ~ b, c ~ 1, a ~ 1))
  )
  expect_equal(model$links, list(c = "log", a = "log"))
})

test_that("an unchanged model passes the pipeline check untouched", {
  ff <- bmmformula(c ~ 1, kappa ~ 1)
  dat <- data.frame(y = c(0.1, -0.2, 0.3))
  model <- sdm(resp_error = "y")
  expect_equal(check_model(model, dat, ff)$links, model$links)
})
