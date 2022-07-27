#' Check if the Zotero library with citation keys is available
#'
#' @param public_library_id
#'
#' @return
bbt_zotero_library_unavailable <- function(public_library_id) {

  # Check if Zotero (with Better BibTeX) is unavailable. If Zotero (with Better
  # BibTeX) is available then check if the public Zotero library for the project
  # exists in the user's Zotero library.
  if (!rbbt::has_bbt()) {
    library_unavailable <- TRUE
  } else {
    library_unavailable <- ifelse(
      TRUE == tryCatch(
        rbbt::bbt_library_id(public_library_id),
        error = function(e) TRUE
      ),
      TRUE,
      FALSE
    )
  }

  library_unavailable

}

#' Write `references.json` file using Zotero library
#'
#' @param path
#' @param keys
#' @param ignore
#' @param public_library_id
#'
#' @return
write_bib <- function(path, keys, ignore, public_library_id) {

  if (bbt_zotero_library_unavailable(public_library_id)) {
    path
  } else {
    rbbt::bbt_write_bib(
      path,
      keys = rbbt::bbt_detect_citations(
        path = keys
      ),
      ignore = paste0("R-", ignore),
      overwrite = TRUE
    )
    path
  }

}

#' Write `packages.bib` file
#'
#' @param path
#' @param packages
#'
#' @return
write_packages_bib <- function(path, packages) {

  knitr::write_bib(
    x = packages,
    file = path
  )

  path

}
