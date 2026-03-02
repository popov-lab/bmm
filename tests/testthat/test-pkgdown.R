test_that("_pkgdown.yml reference section covers all non-internal topics", {
  skip_if_not_installed("yaml")

  # Find package root: try test_path() first, then fall back to getwd()-relative paths.
  # _pkgdown.yml is excluded from the R CMD CHECK tarball via .Rbuildignore, so the
  # test skips gracefully when running inside R CMD CHECK.
  pkg_root <- local({
    candidates <- c(
      tryCatch(normalizePath(test_path("..", "..")), error = function(e) ""),
      normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
      normalizePath(file.path(getwd(), ".."), mustWork = FALSE),
      normalizePath(getwd(), mustWork = FALSE)
    )
    found <- candidates[file.exists(file.path(candidates, "_pkgdown.yml"))]
    if (length(found) > 0) found[1] else NULL
  })
  skip_if(is.null(pkg_root), "_pkgdown.yml not found (skipping outside source tree)")

  pkgdown_file <- file.path(pkg_root, "_pkgdown.yml")

  # Find man directory
  man_dir <- file.path(pkg_root, "man")
  skip_if(!dir.exists(man_dir), "man directory not found")

  has_keyword_pattern <- '^has_keyword\\("(.+)"\\)$'

  # Parse each .Rd file to get topic name and keywords
  rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)

  rd_info <- lapply(rd_files, function(rd_file) {
    content <- readLines(rd_file, warn = FALSE)
    name_lines <- grep("^\\\\name\\{", content, value = TRUE)
    if (length(name_lines) == 0L) return(NULL)
    name <- sub("\\\\name\\{(.+)\\}", "\\1", name_lines[1])
    keywords <- gsub("\\\\keyword\\{(.+)\\}", "\\1",
                     grep("^\\\\keyword\\{", content, value = TRUE))
    list(name = name, keywords = keywords)
  })
  rd_info <- Filter(Negate(is.null), rd_info)

  # Keep only non-internal public topics
  public_topics <- Filter(function(x) !("internal" %in% x$keywords), rd_info)
  public_topic_names <- vapply(public_topics, `[[`, character(1), "name")

  # Parse _pkgdown.yml reference section
  pkgdown_config <- yaml::read_yaml(pkgdown_file)

  # Build the set of topics covered by the reference section
  covered_topics <- character(0)
  for (section in pkgdown_config$reference) {
    for (item in section$contents) {
      if (!is.character(item)) next
      if (grepl(has_keyword_pattern, item)) {
        # Expand has_keyword("xxx") to all topics carrying that keyword
        keyword <- sub(has_keyword_pattern, "\\1", item)
        matching_names <- vapply(
          Filter(function(x) keyword %in% x$keywords, public_topics),
          `[[`, character(1), "name"
        )
        covered_topics <- c(covered_topics, matching_names)
      } else {
        covered_topics <- c(covered_topics, item)
      }
    }
  }

  missing_topics <- setdiff(public_topic_names, covered_topics)

  expect_equal(
    missing_topics, character(0),
    info = paste(
      "Topics missing from _pkgdown.yml reference section:",
      paste(missing_topics, collapse = ", ")
    )
  )
})
