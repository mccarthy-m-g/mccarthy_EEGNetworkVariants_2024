#' Create a table for environment session info
#'
#' @return A tibble.
session_info_environment <- function() {

  session_info_environment <- sessioninfo::session_info()$platform

  session_info_environment$python <- sessioninfo::python_info()$version

  session_info_environment <- session_info_environment |>
    unlist() |>
    as.data.frame() |>
    tibble::rownames_to_column() |>
    dplyr::rename_with(\(.x) c("Setting", "Value")) |>
    tibble::add_row(
      Setting = c("zotero", "better bibtex"),
      Value   = c("6.0.30", "6.7.140")
    )

  session_info_environment

}

#' Create a table for R package session info
#'
#' @return A tibble.
session_info_packages <- function() {

  project_dependencies <- renv::dependencies(here::here()) |>
    tibble::as_tibble() |>
    dplyr::pull(Package) |>
    unique()

  session_info_packages <- sessioninfo::session_info(project_dependencies) |>
    purrr::pluck("packages") |>
    tibble::as_tibble() |>
    dplyr::filter(package %in% project_dependencies) |>
    dplyr::select(package, ondiskversion, source)

  colnames(session_info_packages) <- c("Package", "Version", "Source")

  session_info_packages

}
