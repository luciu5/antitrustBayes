.ab_test_path <- function(...) {
  parts <- as.character(c(...))
  source_root <- testthat::test_path("..", "..")
  source_path <- do.call(file.path, c(list(source_root), as.list(parts)))
  if (file.exists(source_path) || dir.exists(source_path)) {
    return(normalizePath(source_path, mustWork = TRUE))
  }

  if (length(parts) > 1L && identical(parts[[1L]], "inst")) {
    installed_path <- do.call(system.file, c(
      as.list(parts[-1L]),
      list(package = "antitrustBayes", mustWork = FALSE)
    ))
    if (nzchar(installed_path) &&
        (file.exists(installed_path) || dir.exists(installed_path))) {
      return(installed_path)
    }
  }
  source_path
}
