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

#' Title
#'
#' @param x
#'
#' @return
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Title
#'
#' @param input_file
#'
#' @return A list of tibbles.
summarize_bad_channels <- function(input_file) {

  file_handler <- input_file |>
    dplyr::mutate(
      bad_channels = stringr::str_split(bad_channels, "[[:space:]]"),
      bad_channels_n = ifelse(is.na(bad_channels), 0, lengths(bad_channels))
    ) |>
    dplyr::select(participant:task, bad_channels, bad_channels_n)

  # TODO: Considering filtering first to remove 0 bad channels
  descriptives <- file_handler |>
    dplyr::filter(bad_channels_n != 0) |>
    dplyr::summarise(
      bad_channels_n_min = min(bad_channels_n),
      bad_channels_n_max = max(bad_channels_n),
      bad_channels_n_mode = mode(bad_channels_n)
    )

  count_per_participant <- file_handler |>
    dplyr::count(bad_channels_n) |>
    # Reverse the arrange order so cumsum() goes from max to min.
    dplyr::arrange(dplyr::desc(bad_channels_n)) |>
    dplyr::mutate(
      n_cumulative = cumsum(n),
      percent = n / sum(n) * 100,
      percent_cumulative = cumsum(percent)
    )

  count_per_channel <- file_handler |>
    dplyr::select(participant:task, bad_channels) |>
    tidyr::unnest(cols = bad_channels) |>
    na.omit() |>
    dplyr::count(bad_channels) |>
    # Reverse the arrange order so cumsum() goes from the most to least marked
    # channels.
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::mutate(
      n_cumulative = cumsum(n),
      percent = n / sum(n) * 100,
      percent_cumulative = cumsum(percent)
    )

  list(
    descriptives = descriptives,
    count_per_participant = count_per_participant,
    count_per_channel = count_per_channel
  )

}

summarize_bad_segments <- function(input_file) {

  file_handler <- input_file |>
    dplyr::mutate(
      annotation_durations = stringr::str_split(annotation_durations, "[[:space:]]"),
      bad_segments_n = ifelse(is.na(annotation_durations), 0, lengths(annotation_durations)),
    ) |>
    dplyr::select(participant:task, annotation_durations, bad_segments_n)

  descriptives <- file_handler |>
    dplyr::summarise(
      bad_segments_n_min = min(bad_segments_n),
      bad_segments_n_max = max(bad_segments_n),
      bad_segments_n_mode = mode(bad_segments_n)
    )

  count_per_participant <- file_handler |>
    dplyr::count(bad_segments_n) |>
    # Reverse the arrange order so cumsum() goes from max to min.
    dplyr::arrange(dplyr::desc(bad_segments_n)) |>
    dplyr::mutate(
      n_cumulative = cumsum(n),
      percent = n / sum(n) * 100,
      percent_cumulative = cumsum(percent)
    )

  bad_segment_durations <- file_handler |>
    dplyr::select(participant:task, annotation_durations) |>
    tidyr::unnest(cols = annotation_durations) |>
    dplyr::mutate(annotation_durations = as.numeric(annotation_durations))
  na.omit()

  bad_segment_durations_descriptives <- bad_segment_durations |>
    dplyr::summarise(
      min = min(annotation_durations),
      max = max(annotation_durations),
      mean = mean(annotation_durations),
      sd = sd(annotation_durations),
      median = median(annotation_durations)
    )

  bad_segment_duration_totals <- bad_segment_durations |>
    dplyr::group_by(label = paste(participant, session, task, sep = "_")) |>
    dplyr::summarise(
      sum = sum(annotation_durations)
    )

  list(
    descriptives = descriptives,
    count_per_participant = count_per_participant,
    bad_segment_durations_descriptives = bad_segment_durations_descriptives,
    bad_segment_duration_totals = bad_segment_duration_totals
  )

}

