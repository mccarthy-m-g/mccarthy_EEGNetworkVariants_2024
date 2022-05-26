#' Title
#'
#' @param input_file
#'
#' @return
#' @export
#'
#' @examples
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
#' @return
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
plot_connectivity <- function(input_file) {

  # Read
  connectivity_matrix <- get_connectivity_matrix(input_file)

  # Process
  # TODO: Modify/write plotting code
  con_lower_tri <- con_mat[lower.tri(con_mat)]

  # Some descriptive stats
  con_mean <- sprintf("%.2f", mean(con_lower_tri))
  con_min  <- sprintf("%.2f", min(con_lower_tri))
  con_max  <- sprintf("%.2f", max(con_lower_tri))

  # Plot
  con_mat %>%
    as.data.frame.table() %>%
    ggplot(aes(x=Var1, y=Var2, fill=Freq)) +
    ggplot2::geom_raster() +
    scale_fill_continuous(limits = c(0,1), type = "viridis") +
    scale_x_discrete(guide = guide_axis(angle = 90), limits = rev) +
    labs(
      title = case_name,
      subtitle = paste0(
        "Mean: ", con_mean, "\n",
        "Min: ", con_min, "\n",
        "Max: ", con_max
      ),
      x = "",
      y = "",
      fill = "Con Strength"
    )

  # Write
  output_file <- input_file |>
    # Directory
    stringr::str_replace("data", "figures") |>
    # File name
    stringr::str_replace(".csv", ".png")

  fs::dir_create(fs::path_dir(output_file))
  # TODO: Add code for saving

  output_file
}

# TODO: This is the patchwork plot of all the functional connectivity matrices for a subject
plot_connectivity_patchwork <- function() {

}
