#' Create a vector of channels names ordered from left to centre to right, front to back
#'
#' @param input A connectivity matrix.
#'
#' @return Character vector of ordered EEG channels names.
channel_order <- function(input) {

  # Define helper functions for sorting
  get_channel_lobe <- function(channel) {
    stringr::str_extract(channel, "[A-Y|a-y]+")
  }

  get_channel_lobe("Fpz")

  get_channel_number <- function(channel) {
    channel_number <- channel |>
      stringr::str_extract("[[:digit:]]+") |>
      as.numeric()
    channel_number <- ifelse(
      is.na(channel_number),
      0,
      channel_number
    )
    channel_number
  }

  # Set the channel order
  channels <- tibble::tibble(channel = colnames(input))

  channel_order <- channels |>
    dplyr::mutate(
      channel_lobe = get_channel_lobe(channel),
      # Order the channels from anterior to posterior
      channel_lobe = factor(
        channel_lobe,
        levels = c(
          "Fp", "AF", "F", "FT", "FC", "T", "C", "TP", "CP", "P", "PO", "O"
        )
      ),
      channel_number = get_channel_number(channel),
      # In the standard 10-10 system odd channels are on the left and even
      # channels are on the right. "z" or 0 channels are in the centre.
      channel_hemisphere = dplyr::case_when(
        # A number is odd if division by 2 gives a remainder of 1.
        channel_number %% 2 == 1 ~ "left",
        channel_number == 0      ~ "centre",
        # A number is even if division by 2 gives a remainder of 0.
        channel_number %% 2 == 0 ~ "right"
      ),
      channel_hemisphere = factor(
        channel_hemisphere,
        levels = c("left", "centre", "right")
      ),
      channel_number = factor(
        channel_number,
        levels = c(9, 7, 5, 3, 1, 0, 2, 4, 6, 8, 10)
      )
    ) %>%
    dplyr::group_by(channel_hemisphere) %>%
    dplyr::arrange(channel_lobe, channel_number, .by_group = TRUE)

  dplyr::pull(channel_order, channel)

}

#' Estimate phase coupling between pairs of channels using PLI
#'
#' @param input_file
#'
#' @details
#' This function is a wrapper for the `mne.connectivity.spectral_connectivity()`
#' function.
#'
#' @references
#'
#' - <https://mne.tools/0.22/generated/mne.read_epochs.html#mne.read_epochs>
#' - <https://mne.tools/0.22/generated/mne.connectivity.spectral_connectivity.html#mne.connectivity.spectral_connectivity>
#'
#' @return A character vector of file paths.
estimate_phase_coupling <- function(input_file) {

  # Read
  epochs_data <-  mne$read_epochs(input_file, preload = TRUE, verbose = FALSE)

  # Process

  ## Estimate phase coupling using the PLI
  phase_coupling <- mne$connectivity$spectral_connectivity(
    epochs_data,
    method = 'pli',
    mode = 'multitaper',
    sfreq = epochs_data$info$sfreq,
    fmin = epochs_data$info$highpass,
    fmax = epochs_data$info$lowpass,
    faverage = TRUE,
    mt_adaptive = TRUE,
    n_jobs = 1L,
    verbose = FALSE
  )
  names(phase_coupling) <- c("con", "freqs", "times", "n_epochs", "n_tapers")

  ## Tidy the connectivity matrix
  connectivity_matrix <- matrix(
    phase_coupling$con,
    nrow = 64,
    ncol = 64,
    dimnames = list(
      epochs_data$info$ch_names, # rownames
      epochs_data$info$ch_names  # colnames
    )
  )
  ### The diagonal should be 1.0 to reflect perfect similarity of a label with
  # itself. MNE zeroes out the diagonals so it needs to be fixed here.
  diag(connectivity_matrix) <- 1.0
  ### The matrix should be symmetric for visualizations. MNE only fills the lower
  # triangle of the matrix so it needs to be fixed here.
  connectivity_matrix <- as.matrix(
    Matrix::forceSymmetric(connectivity_matrix, uplo = "L")
  )
  ### The rows and columns should be in an order that makes sense.
  reorder_channels <- channel_order(connectivity_matrix)
  connectivity_matrix <- connectivity_matrix[reorder_channels, reorder_channels]

  # Write
  output_file <- input_file |>
    # Directory
    stringr::str_replace("recordings", "sensor-connectivity") |>
    stringr::str_replace("filt-ds_reref_ica-interp_epoch-filt", "phase") |>
    # File name
    stringr::str_replace("filt-ds_reref_ica-interp_epoch-filt", "phase") |>
    stringr::str_replace("-epo.fif.gz", ".csv")

  fs::dir_create(fs::path_dir(output_file))
  connectivity_matrix |>
    tibble::as_tibble() |>
    readr::write_csv(output_file)

  output_file

}


#' Estimate PLI with Hilbert transform
#'
#' @param input_file
#' @param absolute_last Whether to take the absolute value after averaging
#' across epochs (TRUE; this is what MNE does), or to take the absolute value
#' at each epoch before averaging across epochs (FALSE).
#'
#' @details
#' This function calculates the Phase Lag Index following Stam's original
#' definition using the instantaneous phase angle differences obtained using
#' the Hilbert Transform.
#'
#' @references
#'
#' @return A character vector of file paths.
estimate_phase_coupling_hilbert <- function(input_file, absolute_last) {

  # Read
  epochs_data <-  mne$read_epochs(input_file, preload = TRUE, verbose = FALSE)

  # Process

  ## The absolute value of the signed PLI can either be returned for the
  ## connectivity matrix of each epoch, or for the averaged across epochs
  ## connectivity matrix. In general, taking the absolute value of each
  ## epoch will result in higher mean connectivity. Both methods make
  ## different assumptions about the nature of the phase coupling.
  absolute_by_epoch <- ifelse(absolute_last, FALSE, TRUE)

  ## Estimate phase coupling using the PLI
  phase_coupling <- pli$phase_lag_index(
    epochs_data, average = TRUE, absolute = absolute_by_epoch
  )
  ### Take the absolute value after getting the mean over epochs
  if (absolute_last) phase_coupling <- abs(phase_coupling)

  ## Tidy the connectivity matrix
  connectivity_matrix <- matrix(
    phase_coupling,
    nrow = 64,
    ncol = 64,
    dimnames = list(
      epochs_data$info$ch_names, # rownames
      epochs_data$info$ch_names  # colnames
    )
  )
  ### The diagonal should be 1.0 to reflect perfect similarity of a label with
  # itself. MNE zeroes out the diagonals so it needs to be fixed here.
  diag(connectivity_matrix) <- 1.0
  ### The rows and columns should be in an order that makes sense.
  reorder_channels <- channel_order(connectivity_matrix)
  connectivity_matrix <- connectivity_matrix[reorder_channels, reorder_channels]

  # Write
  output_file <- input_file |>
    # Directory
    stringr::str_replace("recordings", "sensor-connectivity") |>
    stringr::str_replace("filt-ds_reref_ica-interp_epoch-filt", "phase-hilbert") |>
    # File name
    stringr::str_replace("filt-ds_reref_ica-interp_epoch-filt", "phase-hilbert") |>
    stringr::str_replace("-epo.fif.gz", ".csv")

  fs::dir_create(fs::path_dir(output_file))
  connectivity_matrix |>
    tibble::as_tibble() |>
    readr::write_csv(output_file)

  output_file

}

#' Estimate amplitude coupling between pairs of channels using AEC
#'
#' @param input_file
#' @param absolute_last Whether to take the absolute value after averaging
#' across epochs (TRUE; this is what the `mne.connectivity.spectral_connectivity()`
#' function does), or to take the absolute value at each epoch before averaging
#' across epochs (FALSE).
#'
#' @references
#'
#' - <https://mne.tools/0.22/generated/mne.read_epochs.html#mne.read_epochs>
#' - <https://mne.tools/0.22/generated/mne.connectivity.envelope_correlation.html#mne.connectivity.envelope_correlation>
#'
#' @return A character vector of file paths.
estimate_amplitude_coupling <- function(input_file, absolute_last) {

  # Read
  epochs_data <- mne$read_epochs(input_file, preload = TRUE, verbose = FALSE)

  # Process

  ## The absolute value of the signed corr can either be returned for the
  ## connectivity matrix of each epoch, or for the averaged across epochs
  ## connectivity matrix. In general, taking the absolute value of each
  ## epoch will result in higher mean connectivity. Both methods make
  ## different assumptions about the nature of the phase coupling.
  absolute_by_epoch <- ifelse(absolute_last, FALSE, TRUE)

  ## Estimate amplitude coupling using the AEC
  amplitude_coupling <- mne$connectivity$envelope_correlation(
    epochs_data,
    combine = 'mean',
    orthogonalize = 'pairwise',
    log = FALSE,
    absolute = absolute_by_epoch
  )
  ### Take the absolute value after getting the mean over epochs
  if (absolute_last) amplitude_coupling <- abs(amplitude_coupling)

  ## Tidy the connectivity matrix
  connectivity_matrix <- matrix(
    amplitude_coupling,
    nrow = 64,
    ncol = 64,
    dimnames = list(
      epochs_data$info$ch_names, # rownames
      epochs_data$info$ch_names  # colnames
    )
  )
  ### The diagonal should be 1.0 to reflect perfect similarity of a label with
  # itself. MNE zeroes out the diagonals so it needs to be fixed here.
  diag(connectivity_matrix) <- 1.0
  ### The rows and columns should be in an order that makes sense.
  reorder_channels <- channel_order(connectivity_matrix)
  connectivity_matrix <- connectivity_matrix[reorder_channels, reorder_channels]

  # Write
  output_file <- input_file |>
    # Directory
    stringr::str_replace("recordings", "sensor-connectivity") |>
    stringr::str_replace("filt-ds_reref_ica-interp_epoch-filt", "amplitude") |>
    # File name
    stringr::str_replace("filt-ds_reref_ica-interp_epoch-filt", "amplitude") |>
    stringr::str_replace("-epo.fif.gz", ".csv")

  fs::dir_create(fs::path_dir(output_file))
  connectivity_matrix |>
    tibble::as_tibble() |>
    readr::write_csv(output_file)

  output_file
}

#' Collect all connectivity matrices into a list
#'
#' Collects all the connectivity matrices into a named list, including
#' metadata for each matrix.
#'
#' @param input_file
#'
#' @return A list of lists. Each list contains two elements: The connectivity
#' matrix, and the metadata for that matrix.
get_connectivity_matrix <- function(input_file) {

  # Connectivity matrix
  connectivity_matrix <- input_file |>
    readr::read_csv(col_types = readr::cols()) |>
    as.matrix()
  rownames(connectivity_matrix) <- colnames(connectivity_matrix)

  # Metadata
  case <- input_file |>
    fs::path_file() |>
    stringr::str_extract("[^_]*_[^_]*_[^_]*")
  participant <- case |>
    stringr::str_extract("[^_]*")

  # The branch name does not contain metadata about the case so it needs to
  # be stored for retrieval by other functions
  list(
    connectivity_matrix = connectivity_matrix,
    metadata = list(
      case = case,
      participant = participant,
      input_file = input_file
    )
  )

}

#' Get participant IDs from connectivity matrix metadata
#'
#' @param input List of connectivity matrices.
#'
#' @return A character vector of participant IDs in ascending order.
get_participant_levels <- function(input) {

  input |>
    purrr::map(~ .x$metadata$participant) |>
    unlist() |>
    unique() |>
    sort()

}

#' Get connectivity profiles from connectivity matrices
#'
#' @param input List of connectivity matrices.
#'
#' @return A tibble containing the connectivity profiles from each recording.
get_connectivity_profiles <- function(input) {

  # Get levels for factors
  participant_levels <- get_participant_levels(input)
  label_levels <- case_order(participant_levels)

  # Get lower triangle of connectivity matrices and store in a tibble to prepare
  # for plotting
  input |>
    purrr::map_dfr(
      ~{

        # Get the lower triangle of the connectivity matrix as a named vector, which
        # can be used as a connectivity profile
        connectivity_matrix <- .x$connectivity_matrix
        connectivity_matrix_cellnames <- outer(
          colnames(connectivity_matrix),
          rownames(connectivity_matrix),
          "paste",
          sep = "_x_"
        )

        connectivity_profile <- connectivity_matrix[lower.tri(connectivity_matrix)]
        names(connectivity_profile) <- connectivity_matrix_cellnames[
          lower.tri(connectivity_matrix_cellnames)
        ]

        tibble::tibble(
          participant = .x$metadata$participant,
          case = .x$metadata$case,
          pair = names(connectivity_profile),
          value = connectivity_profile
        )
      },
      .id = "branch"
    ) |>
    dplyr::mutate(case = factor(case, levels = label_levels))

}

#' Calculate group-wide connectivity summary statistics
#'
#' @param input List of connectivity matrices.
#'
#' @return A tibble containing summary statistics across all connectomes.
summarize_connectivity <- function(input) {

  connectivity_profiles <- get_connectivity_profiles(input)

  connectivity_profiles |>
    dplyr::summarise(
      min = min(value),
      max = max(value),
      mean = mean(value),
      sd = sd(value),
      mode = mode(value)
    )

}

#' Plot a connectivity matrix
#'
#' @param input List of connectivity matrices.
#' @param method A character string for the legend title.
#'
#' @return A list of lists. Each list contains two elements: The connectivity
#' matrix plot, and the metadata for that plot.
plot_connectivity <- function(input, method) {

  connectivity_matrix_plot <- input$connectivity_matrix |>
    as.data.frame.table() |>
    ggplot2::ggplot(ggplot2::aes(x = Var1, y = Var2, fill = Freq)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, 1),
      option = "cividis"
    ) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
    # Reversing the y-axis is needed for the diagonal to go from the upper-left
    # to the bottom-right.
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(
      title = input$metadata$case,
      x = "EEG Channels",
      y = "EEG Channels",
      fill = method
    )

  list(
    plot = connectivity_matrix_plot,
    metadata = append(input$metadata, list(method = method))
  )

}

#' Plot connectivity profiles
#'
#' @param input List of connectivity matrices.
#' @param method A character string for the legend title.
#'
#' @return A ggplot of connectivity profiles.
plot_connectivity_profiles <- function(input, method) {

  connectivity_profiles <- get_connectivity_profiles(input)

  # Get levels for factors
  participant_levels <- get_participant_levels(input)
  session_levels <- c("pre", "post", "fu")
  state_levels <- c("rc1", "rc2", "ro1", "ro2")

  # Plot

  ## The axis-tick colours are being used as annotations for the different
  ## sessions and states. This is a bit hackish, but works for our purposes.
  axis_tick_colours <- colorspace::diverging_hcl(
    12, h = c(299, 135), c = 60, l = c(20, 80), power = c(0.7, 1.3)
  )

  case_colours <- rep(
    c(rev(axis_tick_colours[7:12]), axis_tick_colours[1:6]),
    times = length(participant_levels)
  )

  state_legend <- tibble::tibble(
    State = factor(c("Eyes closed", "Eyes open"))
  )

  session_legend_colours <- colorspace::diverging_hcl(
    12, h = c(299, 135), c = 0, l = c(20, 80), power = c(0.7, 1.3)
  )

  session_levels <- c("1-1", "1-2", "2-1", "2-2", "3-1", "3-2")

  session_legend <- tibble::tibble(
    Session = factor(session_levels, levels = session_levels)
  )

  ggplot2::ggplot(connectivity_profiles, ggplot2::aes(x = pair, y = case)) +
    # These rectangles are purely here to make the legend show up; they aren't
    # actually visible in the plot.
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = Session),
      data = session_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = session_legend_colours[6:1],
      guide = ggplot2::guide_legend(
        title = "Session-Recording", order = 3
      )
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = State),
      data = state_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(case_colours[9], case_colours[3]),
      guide = ggplot2::guide_legend(order = 2)
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, 1),
      option = "cividis",
      guide  = ggplot2::guide_colourbar(order = 1)
    ) +
    ggh4x::facet_nested(
      rows = ggplot2::vars(participant),
      independent = TRUE,
      scales = "free",
      switch = "y"
    ) +
    # Expand reduces the spacing between facets
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev) +
    ggplot2::labs(
      x = "EEG Channel Pairs",
      y = NULL,
      fill = method
    ) +
    ggplot2::theme(
      panel.spacing = grid::unit(0, "lines"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      strip.text.y.left = ggplot2::element_text(angle = 0),
      strip.placement = "outside",
      strip.switch.pad.grid = grid::unit(0, "lines"),
      axis.text = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_line(colour = case_colours, size = 1),
      axis.ticks.length.y = grid::unit(0.5, "lines")
    )

}

# TODO: This is the patchwork plot of all the functional connectivity matrices for a subject
plot_connectivity_patchwork <- function() {

}
