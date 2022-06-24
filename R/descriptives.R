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

#' Calculate the mode of a vector
#'
#' @param x
#'
#' @return
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Get summary statistics about bad EEG channels
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

  count_per_recording <- file_handler |>
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
    count_per_recording = count_per_recording,
    count_per_channel = count_per_channel
  )

}

#' Plot summary statistics about bad EEG channels
#'
#' @param input
#'
#' @return
plot_bad_channel_counts <- function(input) {

  plot_count_per_recording <- input |>
    purrr::pluck("count_per_recording") |>
    ggplot2::ggplot(ggplot2::aes(x = n, y = factor(bad_channels_n))) +
      ggplot2::geom_col() +
      ggplot2::labs(
        x = "Number of recordings with n bad channels",
        y = "Number of bad channels"
      )

  plot_count_per_channel <- input |>
    purrr::pluck("count_per_channel") |>
    ggplot2::ggplot(
      ggplot2::aes(x = n, y = forcats::fct_reorder(bad_channels, n))
    ) +
      ggplot2::geom_col() +
      ggplot2::labs(
        x = "Number of recordings marked bad in",
        y = "Bad channel"
      )

  list(
    plot_count_per_recording = plot_count_per_recording,
    plot_count_per_channel   = plot_count_per_channel
  )

}

#' Get summary statistics about bad segments flagged in EEG data
#'
#' @param input_file
#'
#' @return
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

  count_per_recording <- file_handler |>
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
    dplyr::mutate(annotation_durations = as.numeric(annotation_durations)) |>
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
    count_per_recording = count_per_recording,
    bad_segment_durations = bad_segment_durations,
    bad_segment_durations_descriptives = bad_segment_durations_descriptives,
    bad_segment_duration_totals = bad_segment_duration_totals
  )

}

#' Plot summary statistics about bad segments
#'
#' @param input
#'
#' @return List.
plot_bad_segment_descriptives <- function(input) {

  plot_count_per_recording <- input |>
    purrr::pluck("count_per_recording") |>
    ggplot2::ggplot(ggplot2::aes(x = n, y = factor(bad_segments_n))) +
    ggplot2::geom_col() +
    ggplot2::labs(
      x = "Number of recordings with n bad segments",
      y = "Number of bad segments"
    )

  plot_duration_totals <- input |>
    purrr::pluck("bad_segment_duration_totals") |>
    ggplot2::ggplot(ggplot2::aes(x = sum)) +
      ggplot2::geom_histogram(bins = 60) +
      ggplot2::coord_flip() +
      ggplot2::xlab("Bad segment duration total per recording (seconds)")

  plot_durations <- input |>
    purrr::pluck("bad_segment_durations") |>
    ggplot2::ggplot(ggplot2::aes(x = annotation_durations)) +
      ggplot2::geom_histogram(bins = 60) +
      ggplot2::coord_flip() +
      ggplot2::xlab("Bad segment duration (seconds)")

  list(
    plot_count_per_recording = plot_count_per_recording,
    plot_duration_totals = plot_duration_totals,
    plot_durations = plot_durations
  )

}

# TODO: Summarize ICA bad indices
