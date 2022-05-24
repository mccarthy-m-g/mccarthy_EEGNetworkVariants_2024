#' Title
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

#' Title
#'
#' @param path
#' @param keys
#' @param public_library_id
#'
#' @return
write_bib <- function(path, keys, public_library_id) {

  if (bbt_zotero_library_unavailable(public_library_id)) {
    path
  } else {
    rbbt::bbt_write_bib(
      path,
      keys = rbbt::bbt_detect_citations(
        path = keys
      ),
      overwrite = TRUE
    )
    path
  }

}
