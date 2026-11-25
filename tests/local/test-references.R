# Refits of the reference cases against tests/local/references.rds (ADR 0004,
# issue #123); regenerate with tests/local/make-references.R

skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")

references <- if (file.exists("references.rds")) readRDS("references.rds") else list()

# The circular mixtures were referenced on 1.4, whose brms::mixture() family
# made mu1 the location and carried the constant components mu2 and kappa2. Since
# 1.5.0 they have their own families (#434): mu is the main parameter, so its
# intercept is brms's `Intercept`, and imm's fixed b = 0 reaches brms as a
# constant. The 1.4 values are kept, so that a refit is compared with what 1.4
# estimated; only the names change.
circmix_names_1_5 <- function(ref, name) {
  if (name %in% c("mixture2p", "mixture3p", "imm")) {
    ref <- ref[!ref$parameter %in% c("mu2_Intercept", "kappa2_Intercept"), ]
    ref$parameter[ref$parameter == "mu1_Intercept"] <- "Intercept"
  }
  if (name == "imm" && !"b_Intercept" %in% ref$parameter) {
    ref <- rbind(ref, data.frame(parameter = "b_Intercept", Estimate = 0, Est.Error = 0))
  }
  ref
}

for (name in names(reference_cases())) {
  test_that(paste(name, "reproduces its reference fixed effects"), {
    skip_if(
      is.null(references[[name]]),
      paste0("references.rds has no '", name, "' entry: run ",
             "Rscript tests/local/make-references.R ", name)
    )
    ref <- circmix_names_1_5(references[[name]], name)
    fit <- reference_fit(name)
    fe <- brms::fixef(fit)
    expect_setequal(rownames(fe), ref$parameter)

    # At 1000 draws the Monte Carlo error of a posterior mean is a few
    # hundredths of a posterior SD, so 0.5 SD leaves ample room for platform
    # differences, yet is far below the shift a changed likelihood or default
    # prior gives these small data sets.
    for (i in seq_len(nrow(ref))) {
      par <- ref$parameter[i]
      expect_range(fe[par, "Estimate"],
                   ref$Estimate[i] - 0.5 * ref$Est.Error[i],
                   ref$Estimate[i] + 0.5 * ref$Est.Error[i],
                   info = paste(name, par, "estimate"))
      expect_range(fe[par, "Est.Error"], 0.7 * ref$Est.Error[i], 1.4 * ref$Est.Error[i],
                   info = paste(name, par, "posterior SD"))
    }

    np <- brms::nuts_params(fit)
    expect_lte(sum(np$Value[np$Parameter == "divergent__"]), 5)
    # fixed parameters are constants in fixef(), with SD 0 and no Rhat
    sampled <- ref$parameter[ref$Est.Error > 0]
    expect_range(brms::rhat(fit)[paste0("b_", sampled)], 0, 1.05,
                 info = paste(name, "Rhat"))
  })
}
