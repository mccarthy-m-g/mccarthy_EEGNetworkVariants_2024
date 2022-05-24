# TODO: test this function
#' Title
#'
#' @param input_file
#' @param participants
#'
#' @return
summarize_participant_descriptives <- function(input_file, participants = NULL) {

  descriptives <- readr::read_csv(input_file, col_types = readr::cols())

  descriptives |>
    # Only filter when a value is supplied to participants, otherwise return
    # the full data set
    dplyr::filter(
      if (is.null(participants)) TRUE else participant %in% participants
    ) |>
    dplyr::summarise(
      n = dplyr::n(),
      # TODO: maybe turn the set of four functions, min, max, mean, sd, into
      # their own little function to reduce repetition
      years_age_min              = min(years_age, na.rm = TRUE),
      years_age_max              = max(years_age, na.rm = TRUE),
      years_age_mean             = mean(years_age, na.rm = TRUE),
      years_age_sd               = sd(years_age, na.rm = TRUE),
      years_education_min        = min(years_education, na.rm = TRUE),
      years_education_max        = max(years_education, na.rm = TRUE),
      years_education_mean       = mean(years_education, na.rm = TRUE),
      years_education_sd         = sd(years_education, na.rm = TRUE),
      days_pre_to_post_min       = min(days_pre_to_post, na.rm = TRUE),
      days_pre_to_post_max       = max(days_pre_to_post, na.rm = TRUE),
      days_pre_to_post_mean      = mean(days_pre_to_post, na.rm = TRUE),
      days_pre_to_post_sd        = sd(days_pre_to_post, na.rm = TRUE),
      days_post_to_followup_min  = min(days_post_to_followup, na.rm = TRUE),
      days_post_to_followup_max  = max(days_post_to_followup, na.rm = TRUE),
      days_post_to_followup_mean = mean(days_post_to_followup, na.rm = TRUE),
      days_post_to_followup_sd   = sd(days_post_to_followup, na.rm = TRUE),
    ) |>
    cbind(
      descriptives |>
        dplyr::filter(
          if (is.null(participants)) TRUE else participant %in% participants
        ) |>
        dplyr::count(gender) |>
        tidyr::pivot_wider(names_from = gender, values_from = n)
    )

}

# TODO: Write these
get_preprocessing_descriptives <- function(input_files) {

  # TODO: Map over recordings and extract metadata, then save to a csv
  # Or just use the file handler, in which case probably don't need this function
  # and can just use the next one

}

summarize_preprocessing_descriptives <- function(input_file) {

}
