# the init function bmm() builds: the prior decides which parameters exist.
# A helper so that a model's own test file can check its initial values too
configured_initfun <- function(model, formula, data) {
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  formula <- check_formula(model, data, formula)
  config_args <- configure_model(model, data, formula)
  prior <- configure_prior(model, data, config_args$formula, NULL)
  create_initfun(model, data, config_args$formula, prior)
}
