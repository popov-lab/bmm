# the cases the shared tests run for sdm; see helper-model-cases.R. The
# default cell is "sdm"; every other activation x scaling cell (#558) is
# "sdm_<activation>_<scaling>" and runs the sdm_cont family
model_cases_sdm <- function() {
  re_formula <- bmf(kappa ~ 1 + (1 | ID), c ~ 1 + (1 | ID))
  cells <- list(sdm = sdm("dev_rad"))
  for (activation in c("vonmises", "laplace")) {
    for (scaling in c("density", "peak", "range")) {
      if (activation == "vonmises" && scaling == "density") next
      cells[[paste("sdm", activation, scaling, sep = "_")]] <- suppressMessages(
        sdm("dev_rad", activation = activation, scaling = scaling)
      )
    }
  }
  # the response options of K2: the full wheel in both activations, option
  # columns of two options (the target and the response, or the opposite
  # option where the response is the target), and four equally spaced options
  # with the response rounded onto them
  cells$sdm_grid <- sdm("dev_rad", n_options = 360)
  cells$sdm_grid_laplace <- suppressMessages(sdm("dev_rad", n_options = 360, activation = "laplace"))
  cells$sdm_options <- sdm("dev_rad", n_options = 360, options = c("opt1", "opt2"))
  cells$sdm_alternatives <- sdm("dev_rad", n_options = 360, n_alternatives = 4)
  data_of <- function(name) {
    dat <- small_data()
    if (name == "sdm_options") {
      dat$opt1 <- 0
      dat$opt2 <- deg2rad(round(rad2deg(dat$dev_rad)))
      dat$opt2[dat$opt2 == 0] <- pi
    }
    if (name == "sdm_alternatives") {
      dat$dev_rad <- .bin_options(dat$dev_rad, 4, "radians") * pi / 2
    }
    dat
  }
  cases <- function(build) Map(build, cells, names(cells))
  list(
    stored_frame = cases(function(model, name) {
      list(model = model, formula = bmf(c ~ 1, kappa ~ 1), data = data_of(name))
    }),
    sd_default = cases(function(model, name) {
      list(
        model = model, formula = re_formula, data = data_of(name),
        sd = c(kappa = "exponential(1)", c = "exponential(1)")
      )
    }),
    init = cases(function(model, name) list(model = model, formula = re_formula, data = data_of(name))),
    links = cases(function(model, name) list(model = model, settable = character(0))),
    fixture = list(sdm = list(file = "bmmfit_example1.rds"))
  )
}
