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

#' Plot connectivity histograms
#'
#' @param input List of connectivity matrices.
#' @param method A character string for the x-axis title.
#'
#' @return A ggplot of connectivity histograms.
plot_connectivity_histogram <- function(input, method) {

  connectivity_profiles <- get_connectivity_profiles(input)

  ggplot2::ggplot(
    connectivity_profiles, ggplot2::aes(x = value, fill = participant)
  ) +
    ggplot2::geom_histogram(binwidth = .01) +
    ggplot2::scale_fill_viridis_d(option = "D") +
    ggplot2::labs(
      x = "AEC",
      y = "Count",
      fill = "Participant"
    ) +
    ggplot2::theme_classic()

}

#' Plot all connectivity matrices per participant
#'
#' @param input List of connectivity matrices.
#' @param method A character string for the legend title.
#'
#' @return A named list of patchwork plots, where each index is for all matrices
#'   of a given participant.
plot_connectivity_patchwork <- function(input, method) {

  # Get levels for factors
  participant_levels <- get_participant_levels(input)
  label_levels <- tidyr::expand_grid(
    participant = participant_levels,
    session_state = factor(
      c(
        "pre_rc1",
        "pre_rc2",
        "pre_ro1",
        "pre_ro2",
        "post_rc1",
        "post_rc2",
        "post_ro1",
        "post_ro2",
        "fu_rc1",
        "fu_rc2",
        "fu_ro1",
        "fu_ro2"
      )
    )
  ) |>
    dplyr::arrange(participant) |>
    dplyr::mutate(label = paste0(participant, "_", session_state)) |>
    purrr::pluck("label")

  session_levels <- c("pre", "post", "fu")
  state_levels <- c("rc1", "rc2", "ro1", "ro2")

  sensor_levels <- c(
    "Fp1",
    "Fp2",
    "AF7",
    "AF3",
    "AFz",
    "AF4",
    "AF8",
    "F7",
    "F5",
    "F3",
    "F1",
    "Fz",
    "F2",
    "F4",
    "F6",
    "F8",
    "FT9",
    "FT7",
    "FC5",
    "FC3",
    "FC1",
    "FCz",
    "FC2",
    "FC4",
    "FC6",
    "FT8",
    "FT10",
    "T7",
    "C5",
    "C3",
    "C1",
    "Cz",
    "C2",
    "C4",
    "C6",
    "T8",
    "TP9",
    "TP7",
    "CP5",
    "CP3",
    "CP1",
    "CPz",
    "CP2",
    "CP4",
    "CP6",
    "TP8",
    "TP10",
    "P7",
    "P5",
    "P3",
    "P1",
    "Pz",
    "P2",
    "P4",
    "P6",
    "P8",
    "PO7",
    "PO3",
    "POz",
    "PO4",
    "PO8",
    "O1",
    "Oz",
    "O2"
  )

  # Get the connectivity matrices as a list of tibbles
  connectivity_matrices <- input |>
    purrr::map_dfr(
      ~{
        connectivity_matrix <- .x$connectivity_matrix
        diag(connectivity_matrix) <- NA
        connectivity_matrix |>
          as.table() |>
          as.data.frame(stringsAsFactors = FALSE) |>
          tibble::as_tibble() |>
          ####
          dplyr::mutate(
            Var1 = factor(Var1, levels = sensor_levels),
            Var2 = factor(Var2, levels = sensor_levels)
          ) |>
          ####
          tibble::add_column(
            participant = .x$metadata$participant,
            case = .x$metadata$case,
            .before = 1
          )
      },
      .id = "branch"
    ) |>
    dplyr::mutate(
      session = dplyr::case_when(
        stringr::str_detect(case, "pre")  ~ "Session 1",
        stringr::str_detect(case, "post") ~ "Session 2",
        stringr::str_detect(case, "fu")   ~ "Session 3",
      ),
      state = dplyr::case_when(
        stringr::str_detect(case, "rc1") ~ "Eyes closed 1",
        stringr::str_detect(case, "rc2") ~ "Eyes closed 2",
        stringr::str_detect(case, "ro1") ~ "Eyes open 1",
        stringr::str_detect(case, "ro2") ~ "Eyes open 2"
      ),
      title = paste0(participant, ", ", session, ", ", state),
      case = factor(case, levels = label_levels)
    ) |>
    dplyr::group_by(participant)

  connectivity_matrices |>
    dplyr::group_split() |>
    rlang::set_names(unlist(dplyr::group_keys(connectivity_matrices))) |>
    purrr::map(
      ~ {
        connectivity_matrix <- dplyr::group_by(.x, case)

        connectivity_plots <- connectivity_matrix |>
          dplyr::group_split() |>
          rlang::set_names(unlist(dplyr::group_keys(connectivity_matrix))) |>
          purrr::map(
            ~ {
              ggplot2::ggplot(.x, ggplot2::aes(x = Var1, y = Var2, fill = Freq)) +
                ggplot2::geom_raster() +
                ggplot2::scale_fill_viridis_c(
                  limits = c(0, NA)
                ) +
                ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
                # Reversing the y-axis is needed for the diagonal to go from the upper-left
                # to the bottom-right.
                ggplot2::scale_y_discrete(limits = rev) +
                ggplot2::labs(
                  title = .x$title,
                  x = "EEG Channels",
                  y = "EEG Channels",
                  fill = method
                ) +
                ggplot2::theme(
                  axis.text = ggplot2::element_blank(),
                  axis.ticks = ggplot2::element_blank()
                )
            }
          )

        patchwork::wrap_plots(connectivity_plots, ncol = 3, byrow = FALSE)
      }
    )
}

#' Save connectivity patchwork figures
#'
#' @param connectivity_patchwork The connectivity patchwork plots
#' @param freq_band The frequency band
#'
#' @return A character vector of the file paths.
save_connectivity_patchwork <- function(connectivity_patchwork, freq_band, mode) {

  purrr::map(
    names(connectivity_patchwork),
    ~ {

      filename <- paste0(
        "figures/", freq_band, "/", mode,"-connectivity-patchwork/", .x, ".png"
      )
      names(filename) <- .x

      ggplot2::ggsave(
        filename = filename,
        plot = connectivity_patchwork[[.x]],
        device = ragg::agg_png,
        width = 21.59,
        height = 18,
        units = "cm",
        dpi = "retina",
        scaling = 0.65
      )

      filename

    }
  ) |>
    unlist()

}

#' Plot connectivity histograms
#'
#' @param delta,theta,alpha,beta,gamma The connectivity matrices for a given
#'   frequency band.
#' @param method String. "PLI" or "AEC".
#'
#' @return A ggplot2 object.
plot_connectivity_histograms <- function(
  delta,
  theta,
  alpha,
  beta,
  gamma,
  method
) {

  connectivity_profiles <- purrr::map_dfr(
    list(
      Delta = delta,
      Theta = theta,
      Alpha = alpha,
      Beta  = beta,
      Gamma = gamma
    ),
    get_connectivity_profiles,
    .id = "frequency_band"
  ) |>
    mutate(
      frequency_band = factor(
        frequency_band,
        levels = c("Delta", "Theta", "Alpha", "Beta", "Gamma")
      )
    )

  ggplot2::ggplot(
    connectivity_profiles, ggplot2::aes(x = value, fill = participant)
  ) +
    ggplot2::geom_histogram(binwidth = .01) +
    ggplot2::geom_rug(alpha = 0.1) +
    ggplot2::scale_fill_viridis_d(option = "D") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::facet_wrap(
      ggplot2::vars(frequency_band), ncol = 1, scales = "free_y"
    ) +
    ggplot2::labs(
      x = method,
      y = "Count",
      fill = "Participant"
    ) +
    ggplot2::theme_classic()

}

#' Save connectivity histogram patchwork figure
#'
#' @param phase The phase connectivity histograms plot
#' @param amplitude The amplitude connectivity histograms plot
#'
#' @return A character vector of the file path.
save_connectivity_histogram_patchwork <- function(
  phase, amplitude
) {

  filename <- "figures/connectivity-histogram-patchwork.png"

  connectivity_histogram_patchwork <- phase + amplitude +
    patchwork::plot_layout(ncol = 2, guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")

  ggplot2::ggsave(
    filename = filename,
    plot = connectivity_histogram_patchwork,
    device = ragg::agg_png,
    width = 21.59,
    height = 18,
    units = "cm",
    dpi = "retina",
    scaling = 0.65
  )

  filename

}

#' Save illustrative connectomes figure
#'
#' Save a patchwork figure with the PLI and AEC functional connectome
#' corresponding to a single participant/recording.
#'
#' @param phase The phase connectivity histograms plot
#' @param amplitude The amplitude connectivity histograms plot
#' @param participant The participant ID.
#' @param recording The recording ID. One of "pre_rc1", "pre_rc2", "pre_ro1",
#'   "pre_ro2", "post_rc1", "post_rc2", "post_ro1", "post_ro2", "fu_rc1",
#'   "fu_rc2", "fu_ro1", "fu_ro2".
#'
#' @return A character vector of the file path.
save_illustrative_connectomes_figure <- function(
  phase_matrices,
  amplitude_matrices,
  participant,
  recording
) {

  filename <- paste0(
    "figures/illustrative-connectomes-", participant, "_", recording, ".png"
  )

  # Get levels for factors
  sensor_levels <- c(
    "Fp1",
    "Fp2",
    "AF7",
    "AF3",
    "AFz",
    "AF4",
    "AF8",
    "F7",
    "F5",
    "F3",
    "F1",
    "Fz",
    "F2",
    "F4",
    "F6",
    "F8",
    "FT9",
    "FT7",
    "FC5",
    "FC3",
    "FC1",
    "FCz",
    "FC2",
    "FC4",
    "FC6",
    "FT8",
    "FT10",
    "T7",
    "C5",
    "C3",
    "C1",
    "Cz",
    "C2",
    "C4",
    "C6",
    "T8",
    "TP9",
    "TP7",
    "CP5",
    "CP3",
    "CP1",
    "CPz",
    "CP2",
    "CP4",
    "CP6",
    "TP8",
    "TP10",
    "P7",
    "P5",
    "P3",
    "P1",
    "Pz",
    "P2",
    "P4",
    "P6",
    "P8",
    "PO7",
    "PO3",
    "POz",
    "PO4",
    "PO8",
    "O1",
    "Oz",
    "O2"
  )

  connectomes <- purrr::map(
    list(phase_matrices, amplitude_matrices),
    function(.x) {
      .x |>
        purrr::keep(
          function(.y) {
            .y$metadata$participant == participant &
              stringr::str_detect(.y$metadata$case, recording)
          }
        ) |>
        purrr::pluck(1, "connectivity_matrix")
    }
  )

  connectivity_plots <- purrr::map2(
    connectomes, c("PLI", "AEC"),
    function(.x, .y) {

      diag(.x) <- NA
      .x <- .x |>
        as.table() |>
        as.data.frame(stringsAsFactors = FALSE) |>
        tibble::as_tibble() |>
        dplyr::mutate(
          Var1 = factor(Var1, levels = sensor_levels),
          Var2 = factor(Var2, levels = sensor_levels)
        )

      ggplot2::ggplot(.x, ggplot2::aes(x = Var1, y = Var2, fill = Freq)) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_viridis_c(
          limits = c(0, NA)
        ) +
        ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
        # Reversing the y-axis is needed for the diagonal to go from the upper-left
        # to the bottom-right.
        ggplot2::scale_y_discrete(limits = rev) +
        ggplot2::labs(
          x = "EEG Channels",
          y = "EEG Channels",
          fill = .y
        ) +
        ggplot2::theme(
          axis.text  = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank()
        )
    }
  )

  connectivity_patch <- connectivity_plots |>
    patchwork::wrap_plots() +
    patchwork::plot_annotation(tag_levels = "A")

  ggplot2::ggsave(
    filename = filename,
    plot = connectivity_patch,
    device = ragg::agg_png,
    width = 21.59,
    height = 9,
    units = "cm",
    dpi = "retina",
    scaling = 1
  )

  filename

}
