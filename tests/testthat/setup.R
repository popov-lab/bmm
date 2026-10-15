# pp_check(resp_var = "all") builds a grob grid, and grob construction needs a
# device for font metrics. Without one open, R starts its default device, which
# is pdf() in a non-interactive session and leaves an Rplots.pdf behind.
grDevices::pdf(NULL)
withr::defer(grDevices::dev.off(), teardown_env())
