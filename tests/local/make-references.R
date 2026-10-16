# Writes tests/local/references.rds, the reference summaries that
# test-references.R compares refits with (ADR 0004, issue #123): per model the
# posterior mean and SD of every fixed effect of the case in
# helper-reference-fits.R. From the package root:
#
#   Rscript tests/local/make-references.R               # every model
#   Rscript tests/local/make-references.R sdm sdt_yn    # only these
#
# Regenerate a model's reference only when a change to its likelihood, priors
# or case is intended; the point of the file is to catch the unintended ones.
devtools::load_all(quiet = TRUE)
source("tests/local/setup.R")
source("tests/local/helper-reference-fits.R")

path <- "tests/local/references.rds"
models <- commandArgs(trailingOnly = TRUE)
if (length(models) == 0) models <- names(reference_cases())
unknown <- setdiff(models, names(reference_cases()))
if (length(unknown) > 0) stop("No reference case for: ", paste(unknown, collapse = ", "))

refs <- if (file.exists(path)) readRDS(path) else list()
versions <- c(
  bmm = read.dcf("DESCRIPTION", "Version")[[1]],
  brms = as.character(utils::packageVersion("brms")),
  cmdstan = cmdstanr::cmdstan_version()
)

for (name in models) {
  message("Fitting ", name)
  start <- Sys.time()
  fe <- brms::fixef(reference_fit(name))
  ref <- data.frame(parameter = rownames(fe), Estimate = fe[, "Estimate"],
                    Est.Error = fe[, "Est.Error"], row.names = NULL)
  # per model, because a partial rerun leaves the other models' provenance as it was
  attr(ref, "versions") <- versions
  attr(ref, "date") <- Sys.Date()
  refs[[name]] <- ref
  message(sprintf("  %s: %d fixed effects, %.0f s", name, nrow(ref),
                  as.numeric(difftime(Sys.time(), start, units = "secs"))))
}

saveRDS(refs[sort(names(refs))], path)
message("Wrote ", path, " (", file.size(path), " bytes)")
