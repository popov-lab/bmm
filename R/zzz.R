.onLoad <- function(libname, pkgname) {
  suppressMessages(bmm_options(reset_options = TRUE))
}

.onAttach <- function(libname, pkgname) {
  # add banner of package
  banner <- " _
| |_ _____ _____
| . |     |     |
|___|_|_|_|_|_|_|
"
  local_version <- utils::packageVersion("bmm")
  versionMsg <- paste0("Loading bmm (version: ", local_version, ").\n")

  startUpMsg <- c(
    paste0(
      "A short introduction to package is available by calling help(\"bmm\"). \n",
      "More detailed articles on how to fit different models are available online at https://venpopov.com/bmm/articles \n",
      "You can view the list of currently available models by calling supported_models().\n"
    )
  )

  optionsMsg <- tryCatch2(bmm_options())$message

  if (interactive()) {
    packageStartupMessage(banner, versionMsg, startUpMsg, optionsMsg)
  }
}
