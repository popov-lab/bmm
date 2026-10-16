# the blanket sd priors of a brmsprior, and the one of a parameter; shared so a
# model's own test file can check its sd defaults
sd_rows <- function(pr) {
  pr[pr$class == "sd" & pr$coef == "" & pr$group == "" & pr$prior != "", ]
}

sd_default <- function(pr, par) {
  rows <- sd_rows(pr)
  rows[rows$nlpar == par | rows$dpar == par, ]$prior
}
