# Test data for ezdm model
set.seed(123)
ezdm_test_data <- data.frame(
 subject = rep(1:5, each = 2),
 condition = rep(c("easy", "hard"), 5),
 mean_rt = runif(10, 0.4, 0.6),
 var_rt = runif(10, 0.01, 0.03),
 n_upper = sample(70:90, 10, replace = TRUE),
 n_trials = rep(100, 10)
)

# Test data for 4par version (separate upper/lower RTs)
ezdm_test_data_4par <- data.frame(
 subject = rep(1:5, each = 2),
 condition = rep(c("easy", "hard"), 5),
 mean_rt_upper = runif(10, 0.4, 0.55),
 mean_rt_lower = runif(10, 0.5, 0.7),
 var_rt_upper = runif(10, 0.01, 0.025),
 var_rt_lower = runif(10, 0.02, 0.04),
 n_upper = sample(70:90, 10, replace = TRUE),
 n_trials = rep(100, 10)
)

#############################################################################!
# MODEL SPECIFICATION TESTS                                               ####
#############################################################################!

test_that("ezdm() returns correct class for 3par version", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 expect_s3_class(model, "bmmodel")
 expect_s3_class(model, "ezdm")
 expect_s3_class(model, "ezdm_3par")
})

test_that("ezdm() returns correct class for 4par version", {
 model <- ezdm(
   mean_rt = c("mean_rt_upper", "mean_rt_lower"),
   var_rt = c("var_rt_upper", "var_rt_lower"),
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "4par"
 )
 expect_s3_class(model, "bmmodel")
 expect_s3_class(model, "ezdm")
 expect_s3_class(model, "ezdm_4par")
})

test_that("ezdm() stores correct response and other variables", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 expect_equal(model$resp_vars$mean_rt, "mean_rt")
 expect_equal(model$resp_vars$var_rt, "var_rt")
 expect_equal(model$resp_vars$n_upper, "n_upper")
 expect_equal(model$other_vars$n_trials, "n_trials")
})

test_that("ezdm() has correct parameters for 3par version", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 expect_true("drift" %in% names(model$parameters))
 expect_true("bound" %in% names(model$parameters))
 expect_true("ndt" %in% names(model$parameters))
 expect_true("s" %in% names(model$parameters))
 expect_false("zr" %in% names(model$parameters))
})

test_that("ezdm() has correct parameters for 4par version", {
 model <- ezdm(
   mean_rt = c("mean_rt_upper", "mean_rt_lower"),
   var_rt = c("var_rt_upper", "var_rt_lower"),
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "4par"
 )
 expect_true("drift" %in% names(model$parameters))
 expect_true("bound" %in% names(model$parameters))
 expect_true("ndt" %in% names(model$parameters))
 expect_true("zr" %in% names(model$parameters))
 expect_true("s" %in% names(model$parameters))
})

test_that("ezdm() errors on missing required arguments", {
 expect_error(ezdm(var_rt = "var_rt", n_upper = "n_upper", n_trials = "n_trials"))
 expect_error(ezdm(mean_rt = "mean_rt", n_upper = "n_upper", n_trials = "n_trials"))
 expect_error(ezdm(mean_rt = "mean_rt", var_rt = "var_rt", n_trials = "n_trials"))
 expect_error(ezdm(mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper"))
})

#############################################################################!
# CHECK_DATA TESTS                                                        ####
#############################################################################!

test_that("check_data.ezdm() passes with valid data", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
 expect_silent(check_data(model, ezdm_test_data, formula))
})

test_that("check_data.ezdm() errors on missing variables", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

 # Missing mean_rt
 bad_data <- ezdm_test_data[, c("var_rt", "n_upper", "n_trials")]
 expect_error(
   check_data(model, bad_data, formula),
   "missing from the data"
 )

 # Missing n_trials
 bad_data <- ezdm_test_data[, c("mean_rt", "var_rt", "n_upper")]
 expect_error(
   check_data(model, bad_data, formula),
   "missing from the data"
 )
})

test_that("check_data.ezdm() errors on invalid mean_rt values", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

 # Negative mean_rt
 bad_data <- ezdm_test_data
 bad_data$mean_rt[1] <- -0.5
 expect_error(
   check_data(model, bad_data, formula),
   "must be positive"
 )

 # Zero mean_rt
 bad_data <- ezdm_test_data
 bad_data$mean_rt[1] <- 0
 expect_error(
   check_data(model, bad_data, formula),
   "must be positive"
 )
})

test_that("check_data.ezdm() errors on invalid var_rt values", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

 # Negative variance
 bad_data <- ezdm_test_data
 bad_data$var_rt[1] <- -0.01
 expect_error(
   check_data(model, bad_data, formula),
   "must be positive"
 )
})

test_that("check_data.ezdm() errors on n_upper > n_trials", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

 bad_data <- ezdm_test_data
 bad_data$n_upper[1] <- 150
 expect_error(
   check_data(model, bad_data, formula),
   "cannot exceed total trials"
 )
})

test_that("check_data.ezdm() warns on likely millisecond RTs", {
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

 bad_data <- ezdm_test_data
 bad_data$mean_rt <- bad_data$mean_rt * 1000  # milliseconds
 expect_warning(
   check_data(model, bad_data, formula),
   "greater than 10"
 )
})

#############################################################################!
# MODEL FITTING TESTS (MOCK BACKEND)                                      ####
#############################################################################!

test_that("ezdm 3par model runs with mock backend", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
 expect_silent(
   bmm(formula, ezdm_test_data, model, backend = "mock", mock = 1, rename = FALSE)
 )
})

test_that("ezdm 4par model runs with mock backend", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = c("mean_rt_upper", "mean_rt_lower"),
   var_rt = c("var_rt_upper", "var_rt_lower"),
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "4par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)
 expect_silent(
   bmm(formula, ezdm_test_data_4par, model, backend = "mock", mock = 1, rename = FALSE)
 )
})

test_that("ezdm 3par with condition predictor runs with mock backend", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 0 + condition, bound ~ 1, ndt ~ 1)
 expect_silent(
   bmm(formula, ezdm_test_data, model, backend = "mock", mock = 1, rename = FALSE)
 )
})

test_that("ezdm 3par with random effects runs with mock backend", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(
   drift ~ 0 + condition + (1 | subject),
   bound ~ 1 + (1 | subject),
   ndt ~ 1
 )
 expect_silent(
   bmm(formula, ezdm_test_data, model, backend = "mock", mock = 1, rename = FALSE)
 )
})

test_that("ezdm 4par with zr predictor runs with mock backend", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = c("mean_rt_upper", "mean_rt_lower"),
   var_rt = c("var_rt_upper", "var_rt_lower"),
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "4par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 0 + condition)
 expect_silent(
   bmm(formula, ezdm_test_data_4par, model, backend = "mock", mock = 1, rename = FALSE)
 )
})

#############################################################################!
# CONFIGURE_MODEL TESTS                                                   ####
#############################################################################!

test_that("configure_model.ezdm_3par returns correct family", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
 fit <- bmm(formula, ezdm_test_data, model, backend = "mock", mock = 1, rename = FALSE)
 expect_equal(fit$family$name, "ezdm_3par")
})

test_that("configure_model.ezdm_4par returns correct family", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = c("mean_rt_upper", "mean_rt_lower"),
   var_rt = c("var_rt_upper", "var_rt_lower"),
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "4par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)
 fit <- bmm(formula, ezdm_test_data_4par, model, backend = "mock", mock = 1, rename = FALSE)
 expect_equal(fit$family$name, "ezdm_4par")
})

test_that("ezdm 3par has correct distributional parameters", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
 fit <- bmm(formula, ezdm_test_data, model, backend = "mock", mock = 1, rename = FALSE)
 expect_true(all(c("mu", "drift", "bound", "ndt", "s") %in% fit$family$dpars))
})

test_that("ezdm 4par has correct distributional parameters", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = c("mean_rt_upper", "mean_rt_lower"),
   var_rt = c("var_rt_upper", "var_rt_lower"),
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "4par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)
 fit <- bmm(formula, ezdm_test_data_4par, model, backend = "mock", mock = 1, rename = FALSE)
 expect_true(all(c("mu", "drift", "bound", "ndt", "zr", "s") %in% fit$family$dpars))
})

#############################################################################!
# DEFAULT PRIORS TESTS                                                    ####
#############################################################################!

test_that("ezdm 3par has default priors for all parameters", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = "mean_rt",
   var_rt = "var_rt",
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "3par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
 fit <- bmm(formula, ezdm_test_data, model, backend = "mock", mock = 1, rename = FALSE)
 prior <- brms::prior_summary(fit)

 # Check that priors exist for intercepts (uses dpar column for parameter names)
 expect_true(any(prior$dpar == "drift" & prior$class == "Intercept"))
 expect_true(any(prior$dpar == "bound" & prior$class == "Intercept"))
 expect_true(any(prior$dpar == "ndt" & prior$class == "Intercept"))
})

test_that("ezdm 4par has default priors for zr parameter", {
 skip_on_cran()
 model <- ezdm(
   mean_rt = c("mean_rt_upper", "mean_rt_lower"),
   var_rt = c("var_rt_upper", "var_rt_lower"),
   n_upper = "n_upper",
   n_trials = "n_trials",
   version = "4par"
 )
 formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)
 fit <- bmm(formula, ezdm_test_data_4par, model, backend = "mock", mock = 1, rename = FALSE)
 prior <- brms::prior_summary(fit)

 # Check that prior exists for zr
 expect_true(any(prior$dpar == "zr" & prior$class == "Intercept"))
})