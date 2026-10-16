# Refits of the reference cases against tests/local/references.rds (ADR 0004,
# issue #123); regenerate with tests/local/make-references.R

skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")

references <- if (file.exists("references.rds")) readRDS("references.rds") else list()

for (name in names(reference_cases())) {
  test_that(paste(name, "reproduces its reference fixed effects"), {
    skip_if(
      is.null(references[[name]]),
      paste0("references.rds has no '", name, "' entry: run ",
             "Rscript tests/local/make-references.R ", name)
    )
    ref <- references[[name]]
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
