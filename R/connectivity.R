# TODO: Decide on the order I want.
#' Title
#'
#' @param input
#'
#' @return Character vector.
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

#' Title
#'
#' @param input_file
#'
#' @return
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

#' Title
#'
#' @param input_file
#'
#' @return A List.
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

# TODO: Finish writing then add target
plot_connectivity <- function(input, method) {

  connectivity_matrix_plot <- input$connectivity_matrix |>
    as.data.frame.table() |>
    ggplot2::ggplot(ggplot2::aes(x=Var1, y=Var2, fill=Freq)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_continuous(limits = c(0,1), type = "viridis") +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
    # Reversing the y-axis is needed for the diagonal to go from the upper-left
    # to the bottom-right.
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(
      x = "EEG Channels",
      y = "EEG Channels",
      fill = method
    )

  list(
    plot = connectivity_matrix_plot,
    metadata = append(input$metadata, list(method = method))
  )
}

# TODO: This is the patchwork plot of all the functional connectivity matrices for a subject
plot_connectivity_patchwork <- function() {

}
