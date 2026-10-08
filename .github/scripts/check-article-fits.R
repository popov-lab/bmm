open_fence <- "^[\t >]*(`{3,})\\s*\\{\\s*[rR]([ ,}]|$)"
any_fence <- "^[\t >]*(`{3,})"
close_fence <- "^[\t >]*(`{3,})\\s*$"
fits_path <- "(^|/)fits(/|$)"

fence_length <- function(line, pattern) {
  found <- regmatches(line, regexec(pattern, line))[[1]]
  if (length(found) == 0L) 0L else nchar(found[2])
}

# Follows knitr: a chunk opens with a fence and a header such as
# {r setup}, and closes with a fence at least as long; any other
# fence, such as a plain r fence or a four-backtick block that
# shows chunk syntax, is skipped as a whole.
scan_lines <- function(lines) {
  chunks <- list()
  inline <- integer()
  fence <- 0L
  first <- NULL
  header <- ""
  for (i in seq_along(lines)) {
    if (fence == 0L) {
      fence <- fence_length(lines[i], any_fence)
      if (fence > 0L) {
        header <- lines[i]
        first <- if (grepl(open_fence, lines[i])) i + 1L
      } else if (grepl("`r [^`]*(fits/|[\"']fits[\"'])", lines[i])) {
        inline <- c(inline, i)
      }
    } else if (fence_length(lines[i], close_fence) >= fence) {
      if (!is.null(first)) {
        code <- if (i - 1L >= first) sub("^[\t >]+", "", lines[first:(i - 1L)]) else character()
        chunks[[length(chunks) + 1L]] <- list(first = first, header = header, code = code)
      }
      fence <- 0L
      first <- NULL
    }
  }
  list(chunks = chunks, inline = inline)
}

not_evaluated <- function(chunk) {
  options <- paste(c(chunk$header, grep("^#\\|", chunk$code, value = TRUE)), collapse = " ")
  grepl("eval\\s*(=|:)\\s*(FALSE|F|false)\\b", options)
}

after_file_equals <- function(pd, id) {
  parent <- pd$parent[pd$id == id]
  siblings <- pd[pd$parent == parent, ]
  siblings <- siblings[order(siblings$line1, siblings$col1), ]
  at <- which(siblings$id == id)
  parent > 0 && at > 2L &&
    siblings$token[at - 1L] == "EQ_SUB" &&
    siblings$token[at - 2L] == "SYMBOL_SUB" && gsub("`", "", siblings$text[at - 2L]) == "file"
}

inside_file_argument <- function(pd, id) {
  while (id > 0L) {
    if (after_file_equals(pd, id)) return(TRUE)
    id <- pd$parent[pd$id == id]
  }
  FALSE
}

references <- data.frame(article = character(), line = integer(), fit = character())
problems <- data.frame(article = character(), line = integer(), message = character())
note <- function(article, line, message) {
  problems[nrow(problems) + 1L, ] <<- list(article, line, gsub("[\r\n]+", " ", message))
}

articles <- list.files("vignettes", pattern = "\\.[Rrq]md$", recursive = TRUE, full.names = TRUE)
for (article in articles) {
  lines <- iconv(readLines(article, warn = FALSE, encoding = "UTF-8"), "UTF-8", "UTF-8", sub = "?")
  scanned <- scan_lines(lines)
  for (i in scanned$inline) {
    note(article, i, "Inline code names a fits/ path, which cannot be checked. Load the fit in a chunk with file = \"fits/<name>\".")
  }
  for (chunk in scanned$chunks) {
    if (not_evaluated(chunk)) next
    code <- paste(chunk$code, collapse = "\n")
    parsed <- tryCatch(parse(text = code, keep.source = TRUE, encoding = "UTF-8"), error = function(e) NULL)
    if (is.null(parsed)) {
      if (grepl("[\"']fits[/\"']", code)) {
        note(article, chunk$first, "A code chunk that mentions fits/ cannot be parsed, so its fits cannot be checked.")
      }
      next
    }
    pd <- utils::getParseData(parsed)
    if (is.null(pd)) next
    for (k in which(pd$token == "STR_CONST")) {
      value <- tryCatch(str2lang(pd$text[k]), error = function(e) "")
      if (!is.character(value) || !grepl(fits_path, value)) next
      if (!grepl("fits/", value, fixed = TRUE) && !inside_file_argument(pd, pd$id[k])) next
      line <- chunk$first + pd$line1[k] - 1L
      literal <- pd$parent[k]
      if (after_file_equals(pd, literal) && sum(pd$parent == literal) == 1L && grepl("^fits/.+", value)) {
        references[nrow(references) + 1L, ] <- list(article, line, value)
      } else {
        note(article, line, paste0("The string \"", value, "\" cannot be checked. Pass the fit as a literal, file = \"fits/<name>\", so that the build can check it."))
      }
    }
  }
}

fit_file <- function(article, fit) {
  file.path(dirname(article), if (identical(tools::file_ext(fit), "rds")) fit else paste0(fit, ".rds"))
}
readable <- list()
for (k in seq_len(nrow(references))) {
  path <- fit_file(references$article[k], references$fit[k])
  if (is.null(readable[[path]])) {
    readable[[path]] <- if (!file.exists(path)) {
      "is not in the article-fits release"
    } else {
      tryCatch({ readRDS(path); "ok" }, error = function(e) paste0("cannot be read (", conditionMessage(e), ")"), warning = function(w) paste0("cannot be read (", conditionMessage(w), ")"))
    }
  }
  if (readable[[path]] != "ok") {
    note(references$article[k], references$line[k],
         paste0(basename(references$article[k]), " loads ", references$fit[k], ", which ", readable[[path]], ". It is expected at ", path, ". Upload it to the article-fits release before the site is built."))
  }
}

for (k in seq_len(nrow(problems))) {
  cat(sprintf("::error file=%s,line=%d::%s\n", problems$article[k], problems$line[k], problems$message[k]))
}
cat(sprintf("Checked %d fit references in %d articles; %d problems\n", nrow(references), length(articles), nrow(problems)))
if (nrow(problems) > 0L) quit(save = "no", status = 1L)
