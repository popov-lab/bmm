# Run from the repository root. Lists every string literal that an evaluated R
# chunk of a pkgdown article passes as file = "<relative path>", and fails when
# that file is missing next to the article or cannot be read, because bmm() and
# brm() refit a missing fit during the build instead of failing it. Code that an
# article sources, includes as a child document or computes inline is not seen;
# a fits/ path that is not a literal file argument fails, inline or in a chunk.

chunk_begin <- "^[\t >]*```+\\s*\\{([a-zA-Z0-9_]+( *[ ,].*)?)\\}\\s*$"
chunk_end <- "^[\t >]*```+\\s*$"
inline_fits <- "`r[ #][^`]*(fits/|[\"']fits[\"'])"
fits_path <- "(^|/)fits(/|$)"
not_evaluated <- "eval\\s*(=|:)\\s*(FALSE|F|false)\\b"
writers <- c("saveRDS", "save", "write_rds")

finding <- function(article, line, kind, text) {
  n <- length(line)
  data.frame(article = rep_len(article, n), line = as.integer(line), kind = rep_len(kind, n), text = rep_len(text, n))
}

# knitr 1.52 (group_indices(), match_chunk_end()): plain fences are text; a
# header such as ```{r} opens a chunk wherever it stands, and only a fence with
# the header's own prefix and backtick count closes it, or the end of the file.
scan_chunks <- function(lines) {
  header <- last <- integer()
  in_chunk <- logical(length(lines))
  fence <- NA_character_
  for (i in seq_along(lines)) {
    opens <- grepl(chunk_begin, lines[i])
    if (is.na(fence)) {
      if (opens) {
        header <- c(header, i)
        # knitr sets no closing pattern for a chunk that opens on line 1
        fence <- if (i == 1L) chunk_end else sub("(^[\t >]*```+).*", "^\\1\\\\s*$", lines[i])
      }
    } else if (opens && grepl(sub("^([^`]*`+).*", "^\\1\\\\{", fence), lines[i])) {
      last <- c(last, i - 1L)
      header <- c(header, i)
    } else if (grepl(fence, lines[i])) {
      last <- c(last, i - 1L)
      fence <- NA_character_
      in_chunk[i] <- TRUE
    }
    in_chunk[i] <- in_chunk[i] || !is.na(fence)
  }
  if (!is.na(fence)) last <- c(last, length(lines))
  list(chunks = data.frame(header = header, last = last), text = !in_chunk)
}

chunk_engine <- function(header) sub("^([a-zA-Z0-9_]+).*$", "\\1", trimws(sub(chunk_begin, "\\1", header)))

strip_indent <- function(header, code) {
  indent <- sub("^([\t >]*).*", "\\1", header)
  if (nzchar(indent)) sub(paste0("^", sub("\\s+$", "", indent)), "", sub(paste0("^", indent), "", code)) else code
}

evaluated <- function(header, code) {
  options <- code[seq_len(match(FALSE, startsWith(code, "#| "), nomatch = length(code) + 1L) - 1L)]
  !grepl(not_evaluated, paste(c(header, options), collapse = " "))
}

parent_of <- function(pd, id) pd$parent[pd$id == id]

after_file_equals <- function(pd, id) {
  parent <- parent_of(pd, id)
  siblings <- pd[pd$parent == parent, ]
  siblings <- siblings[order(siblings$line1, siblings$col1), ]
  at <- which(siblings$id == id)
  parent > 0 && at > 2L &&
    siblings$token[at - 1L] == "EQ_SUB" &&
    siblings$token[at - 2L] == "SYMBOL_SUB" && gsub("`", "", siblings$text[at - 2L]) == "file"
}

file_argument_call <- function(pd, id) {
  while (id > 0L && !after_file_equals(pd, id)) id <- parent_of(pd, id)
  if (id > 0L) parent_of(pd, id) else NA_integer_
}

call_name <- function(pd, id) {
  heads <- pd$id[pd$parent == id & pd$token == "expr"]
  c(pd$text[pd$token == "SYMBOL_FUNCTION_CALL" & pd$parent %in% heads], "")[1]
}

enclosing_call_name <- function(pd, id) {
  while (id > 0L && !nzchar(call_name(pd, id))) id <- parent_of(pd, id)
  if (id > 0L) call_name(pd, id) else ""
}

normalize_path <- function(path) sub("^(\\./)+", "", path)

relative_path <- function(path) {
  path <- normalize_path(path)
  nzchar(path) && !endsWith(path, "/") && !grepl("^(/|~|\\\\|[A-Za-z]:)", path) &&
    !grepl("://", path, fixed = TRUE) && !grepl("(^|/)\\.{1,2}(/|$)", path)
}

check_string <- function(article, line, pd, k) {
  value <- tryCatch(str2lang(pd$text[k]), error = function(e) "")
  literal <- pd$parent[k]
  holder <- file_argument_call(pd, literal)
  caller <- if (is.na(holder)) enclosing_call_name(pd, literal) else call_name(pd, holder)
  if (caller %in% writers) {
    NULL
  } else if (after_file_equals(pd, literal) && sum(pd$parent == literal) == 1L && relative_path(value)) {
    finding(article, line, "reference", value)
  } else if (grepl(fits_path, value) && (grepl("fits/", value, fixed = TRUE) || !is.na(holder))) {
    finding(article, line, "problem", paste0(
      "The string \"", value, "\" cannot be checked. Pass the fit as a literal file argument, ",
      "readRDS(file = \"fits/<name>\") or bmm(..., file = \"fits/<name>\"), so that the build can check it."
    ))
  }
}

check_chunk <- function(article, header, code) {
  parsed <- tryCatch(parse(text = paste(code, collapse = "\n"), keep.source = TRUE, encoding = "UTF-8"), error = function(e) NULL)
  pd <- if (!is.null(parsed)) utils::getParseData(parsed)
  if (is.null(parsed) && any(grepl("[\"']fits[/\"']", code))) {
    finding(article, header + 1L, "problem", "A code chunk that mentions fits/ cannot be parsed, so its fits cannot be checked.")
  } else if (!is.null(pd)) {
    do.call(rbind, lapply(which(pd$token == "STR_CONST"), function(k) check_string(article, header + pd$line1[k], pd, k)))
  }
}

check_article <- function(article) {
  lines <- iconv(readLines(article, warn = FALSE, encoding = "UTF-8"), "UTF-8", "UTF-8", sub = "?")
  scanned <- scan_chunks(lines)
  inline <- which(scanned$text & grepl(inline_fits, lines))
  chunks <- scanned$chunks
  rbind(
    finding(article, inline, "problem", "Inline code names a fits/ path, which cannot be checked. Load the fit in a chunk with file = \"fits/<name>\"."),
    do.call(rbind, lapply(seq_len(nrow(chunks)), function(k) {
      header <- lines[chunks$header[k]]
      code <- strip_indent(header, lines[seq_len(chunks$last[k] - chunks$header[k]) + chunks$header[k]])
      if (tolower(chunk_engine(header)) == "r" && evaluated(header, code)) check_chunk(article, chunks$header[k], code)
    }))
  )
}

resolve <- function(article, path) {
  path <- normalize_path(path)
  file.path(dirname(article), ifelse(tools::file_ext(path) == "rds", path, paste0(path, ".rds")))
}

file_status <- function(file) {
  if (!file.exists(file)) {
    "missing"
  } else {
    unreadable <- function(e) paste0("cannot be read (", conditionMessage(e), ")")
    tryCatch({ readRDS(file); "ok" }, error = unreadable, warning = unreadable)
  }
}

reference_problem <- function(article, path, file, status) {
  release <- startsWith(normalize_path(path), "fits/")
  missing <- ifelse(release, "is not in the article-fits release", "is not in the tree")
  advice <- ifelse(release, " Upload it to the article-fits release before the site is built.",
                   ifelse(status == "missing", " Commit it, or load it from fits/ and upload it to the article-fits release.", ""))
  paste0(basename(article), " loads ", path, ", which ", ifelse(status == "missing", missing, status),
         ". It is expected at ", file, ".", advice)
}

escape_data <- function(x) gsub("\n", "%0A", gsub("\r", "%0D", gsub("%", "%25", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)
escape_property <- function(x) gsub(",", "%2C", gsub(":", "%3A", escape_data(x), fixed = TRUE), fixed = TRUE)

articles <- list.files("vignettes", pattern = "\\.[Rrq]md$", recursive = TRUE)
articles <- file.path("vignettes", articles[!startsWith(basename(articles), "_") & !grepl("^tutorials", dirname(articles))])
findings <- do.call(rbind, c(list(finding(character(), integer(), character(), character())), lapply(articles, check_article)))

references <- findings[findings$kind == "reference", ]
files <- resolve(references$article, references$text)
checked <- vapply(unique(files), file_status, character(1))
status <- unname(checked[match(files, unique(files))])
bad <- status != "ok"
problems <- rbind(
  findings[findings$kind == "problem", ],
  finding(references$article[bad], references$line[bad], "problem",
          reference_problem(references$article[bad], references$text[bad], files[bad], status[bad]))
)

cat(sprintf("::error file=%s,line=%d::%s\n", escape_property(problems$article), problems$line, escape_data(problems$text)), sep = "")
cat(sprintf("Checked %d fit references in %d articles; %d problems\n", nrow(references), length(articles), nrow(problems)))
if (nrow(problems) > 0L) quit(save = "no", status = 1L)
