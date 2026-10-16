# conditional_effects() and its plots on a real fit of every model, one
# parameter at a time; the fits are the reference cases of
# helper-reference-fits.R

skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")

# the native range a parameter's link maps onto, when it has one
native_range <- function(par_info) {
  if (par_info$softmax || par_info$link %in% c("logit", "probit", "cloglog")) {
    c(0, 1)
  } else if (par_info$link == "log") {
    c(0, Inf)
  }
}

for (name in names(reference_cases())) {
  test_that(paste("conditional_effects() works on a", name, "fit"), {
    fit <- reference_fit(name)
    for (par in reference_ce_pars(reference_cases()[[name]])) {
      range <- native_range(.get_parameter_info(fit, par))
      for (scale in c("native", "sampling")) {
        label <- paste(name, par, scale)
        ce <- conditional_effects(fit, par = par, scale = scale)
        expect_s3_class(ce, "brms_conditional_effects")
        expect_gt(length(ce), 0)
        for (effect in ce) {
          expect_true(all(is.finite(effect$estimate__)), info = label)
          expect_true(all(effect$lower__ <= effect$estimate__ &
                            effect$estimate__ <= effect$upper__), info = label)
          if (scale == "native" && !is.null(range)) {
            expect_range(c(effect$lower__, effect$upper__), range[1], range[2], info = label)
          }
        }
        # brms labels the colour of a one-variable effect NA, which ggplot2
        # >= 4.0 reports as an unknown label on every build, plot() included
        for (p in SM(plot(ce, plot = FALSE))) {
          expect_s3_class(p, "ggplot")
          expect_no_error(SM(ggplot2::ggplot_build(p)))
        }
      }
    }
  })
}
