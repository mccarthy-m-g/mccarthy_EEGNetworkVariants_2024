#' Check if EEG recordings are unavailable
#'
#' This is used to determine whether or not to skip preprocessing targets.
#'
#' @return Logical.
eeg_recordings_unavailable <- function() {
  !fs::dir_exists(
    here::here("data", "EEG-recordings", "raw")
  )
}

#' Filter and downsample Raw EEG data
#'
#' This function is a wrapper for the
#'
#' @param input_file
#' @param l_freq
#' @param h_freq
#' @param phase
#' @param resample_rate
#'
#' @references
#'
#' - <https://mne.tools/0.22/generated/mne.io.read_raw_brainvision.html#mne.io.read_raw_brainvision>
#' - <https://mne.tools/0.22/generated/mne.io.Raw.html#mne.io.Raw.filter>
#' - <https://mne.tools/0.22/generated/mne.io.Raw.html#mne.io.Raw.resample>
#' - <https://mne.tools/0.22/generated/mne.io.Raw.html#mne.io.Raw.save>
#'
#' @return
preprocess_filter_downsample <- function(
  input_file,
  l_freq,
  h_freq,
  phase,
  resample_rate
) {

  # Read
  raw_data <- mne$io$read_raw_brainvision(
    input_file,
    preload = TRUE,
    verbose = FALSE
  )

  # Preprocess
  raw_data$filter(
    l_freq = l_freq,
    h_freq = h_freq,
    picks = "data",
    phase = phase,
    verbose = FALSE
  )
  raw_data$resample(resample_rate, npad = "auto")

  # Write
  output_file <- input_file |>
    stringr::str_replace("raw", "filt-ds") |>
    stringr::str_replace(".vhdr", "_filt-ds_raw.fif.gz")

  fs::dir_create(fs::path_dir(output_file))
  raw_data$save(output_file, overwrite = TRUE)

  output_file

}

#' Get EEG channel dictionary for renaming
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
get_channel_dictionary <- function(file) {
  channel_dictionary <- readr::read_csv(file, col_types = readr::cols())
  channel_dictionary_new <- channel_dictionary |>
    dplyr::select(new_name) |>
    dplyr::group_split(dplyr::row_number(), .keep = FALSE) |>
    purrr::map(unlist)
  names(channel_dictionary_new) <- channel_dictionary$channel
  dict(channel_dictionary_new)
}

#' Re-reference Raw EEG data to average and set montage
#'
#' Performs the following...
#'
#' @param input_file A path to a file.
#' @param channel_dictionary_file A path to a file.
#'
#' @references
#'
#' - <https://mne.tools/0.22/generated/mne.io.read_raw_fif.html#mne.io.read_raw_fif>
#' - <https://mne.tools/0.22/generated/mne.io.Raw.html#mne.io.Raw.rename_channels>
#' - <https://mne.tools/0.22/generated/mne.add_reference_channels.html#mne.add_reference_channels>
#' - <https://mne.tools/0.22/generated/mne.io.Raw.html#mne.io.Raw.set_eeg_reference>
#' - <https://mne.tools/0.22/generated/mne.channels.make_standard_montage.html#mne.channels.make_standard_montage>
#' - <https://mne.tools/0.22/generated/mne.io.Raw.html#mne.io.Raw.set_montage>
#'
#' @return A character vector of file paths for files written by the function.
preprocess_rereference <- function( # maybe rename to rereference
  input_file,
  channel_dictionary_file
) {

  # Read
  raw_data <- mne$io$read_raw_fif(input_file, preload = TRUE, verbose = FALSE)

  # Preprocess
  channel_dictionary  <- get_channel_dictionary(channel_dictionary_file)
  raw_data <- raw_data$rename_channels(channel_dictionary)

  mne$add_reference_channels(
    raw_data,
    ref_channels = "Cz",
    copy = FALSE
  )
  raw_data$set_eeg_reference(ref_channels = "average", projection = FALSE)

  easycap_1010_montage <- mne$channels$make_standard_montage("easycap-M1")
  raw_data <- raw_data$set_montage(easycap_1010_montage, verbose = FALSE)

  # Write
  output_file <- input_file |>
    stringr::str_replace("filt-ds", "filt-ds_reref") |>
    stringr::str_replace("filt-ds_raw", "filt-ds_reref_raw")

  fs::dir_create(fs::path_dir(output_file))
  raw_data$save(output_file, overwrite = TRUE)

  output_file
}

#' Fit ICA to Raw EEG data then interpolate bad channels
#'
#'
#' @param input_file
#' @param bad_channels
#' @param annotation_onsets
#' @param annotation_durations
#' @param annotation_labels
#'
#' @references
#'
#' - <https://mne.tools/0.22/generated/mne.Annotations.html#mne.Annotations>
#' - <https://mne.tools/0.22/generated/mne.io.Raw.html#mne.io.Raw.set_annotations>
#' - <https://mne.tools/0.22/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA>
#' - <https://mne.tools/0.22/generated/mne.io.Raw.html#mne.io.Raw.interpolate_bads>
#'
#' @return
preprocess_ica <- function(
  input_file,
  bad_channels,
  annotation_onsets,
  annotation_durations,
  annotation_labels,
  ica_n_components,
  ica_method,
  ica_seed,
  bad_ica_indices
) {

  # Read
  raw_data <- mne$io$read_raw_fif(input_file, preload = TRUE, verbose = FALSE)

  # Preprocess

  ## Mark bad channels
  if (!is.na(bad_channels)) {
    raw_data$info["bads"] <- bad_channels |> # turn this into a helper function
      stringr::str_split("[[:space:]]") |>
      unlist() |>
      array()
  }

  ## Annotate bad segments
  if (!is.na(annotation_onsets)) {
    annotation_onsets <- annotation_onsets |>
      stringr::str_split("[[:space:]]") |>
      unlist() |>
      as.numeric() |>
      array()
    annotation_durations <- annotation_durations |>
      stringr::str_split("[[:space:]]") |>
      unlist() |>
      as.numeric() |>
      array()
    annotation_labels <- annotation_labels |>
      stringr::str_split("[[:space:]]") |>
      unlist() |>
      array()
    bad_segments <- mne$Annotations(
      onset = annotation_onsets,
      duration = annotation_durations,
      description = annotation_labels
    )
    # A warning is thrown if an annotation expands outside the data range. These
    # can be ignored as they just represent bad segments I wanted to reach the
    # very beginning or end.
    raw_data$set_annotations(bad_segments, emit_warning = FALSE)
  }

  ## ICA

  ### Fit ICA
  ica_decomposition <- mne$preprocessing$ICA(
    n_components = ica_n_components,
    method = ica_method,
    random_state = ica_seed,
    verbose = FALSE
  )
  ica_decomposition$fit(raw_data)

  ### Retrieve bad component indices
  bad_ica_indices <- bad_ica_indices |>
    stringr::str_split("[[:space:]]") |>
    unlist() |>
    as.integer() |>
    array()

  ### Apply ICA
  ica_decomposition$apply(raw_data, exclude = bad_ica_indices)

  ## Interpolate bad channels
  if (!is.na(bad_channels)) {
    raw_data$interpolate_bads(reset_bads = TRUE, verbose = FALSE)
  }

  # Write
  output_file <- input_file |>
    stringr::str_replace("filt-ds_reref", "filt-ds_reref_ica-interp") |> # Directory
    stringr::str_replace("filt-ds_reref_raw", "filt-ds_reref_ica-interp_raw") # File name end

  fs::dir_create(fs::path_dir(output_file))
  raw_data$save(output_file, overwrite = TRUE)

  output_file
}


#' Epoch Raw EEG data then filter into five frequency bands
#'
#' @param input_file
#' @param filter_freq_band
#' @param filter_l_freq
#' @param filter_h_freq
#'
#' @references
#'
#' - <https://mne.tools/0.22/generated/mne.make_fixed_length_events.html#mne.make_fixed_length_events>
#' - <https://mne.tools/0.22/generated/mne.Epochs.html#mne.Epochs>
#' - <https://mne.tools/0.22/generated/mne.Epochs.html#mne.Epochs.filter>
#' - <https://mne.tools/0.22/generated/mne.Epochs.html#mne.Epochs.crop>
#'
#' @return
preprocess_filter_epoch <- function(
  input_file,
  filter_freq_band,
  filter_l_freq,
  filter_h_freq
  ) {

  # Read
  raw_data <- mne$io$read_raw_fif(input_file, preload = TRUE, verbose = FALSE)

  # Preprocess

  ## This creates an event every five seconds, with times (s) at:
  # 0-8, 5-13, 10-18, ..., etc.
  # The purpose of overlapping events here is to then create a
  # non-overlapping but contiguous series of 5 second epochs after
  # filtering. The overlap and durations are to deal with edge
  # artifacts.
  events <- raw_data |>
    mne$make_fixed_length_events(
    duration = 8.0,
    overlap = 3.0
  ) |>
    reticulate::r_to_py()
  # r_to_py() turns this into float so need to force int
  events <- events$astype("int")


  ## Preloading is necessary for filtering and cropping the epochs data
  epochs_data <- mne$Epochs(
    raw_data,
    events,
    event_id = 1L,
    tmin = 0.0,
    tmax = 8.0,
    baseline = NULL,
    detrend = NULL, # TODO: See if detrending improves results: 1L # Linear detrend
    reject_by_annotation = TRUE, # Epochs with BAD segments should be dropped
    preload = TRUE,
    verbose = FALSE
  )

  ## Frequency bands:
  # - delta (1-4 Hz)
  # - theta (4-8 Hz)
  # - alpha (8-13 Hz)
  # - beta (13-30 Hz)
  # - gamma (30-50 Hz)
  epochs_data <- epochs_data$filter(
    l_freq = filter_l_freq,
    h_freq = filter_h_freq,
    picks = "data",
    phase = "zero-double",
    verbose = FALSE
  )

  ## This creates non-overlapping, contiguous, 5 second slices with times (s) at:
  # 1-6, 6-11, 11-16, ..., etc. The length of 5 seconds was selected based on
  # the minimum epoch length required to get reliable phase coupling estimates
  # in the Delta band. Originally, 2 second slices were chosen (in which case
  # the events object above should be changed to
  # `mne.make_fixed_length_events(..., overlap=6.0)`); but these were too short
  # and lead to a warning about unreliable estimates in Delta band connectivity.
  epochs_data <- epochs_data$crop(tmin=1.0, tmax=6.0)

  # Write
  output_file <- input_file |>
    # Directory
    stringr::str_replace(
      "filt-ds_reref_ica-interp",
      paste0("filt-ds_reref_ica-interp_epoch-filt/", filter_freq_band)
    ) |>
    # File name end
    stringr::str_replace(
      "filt-ds_reref_ica-interp_raw",
      paste0("filt-ds_reref_ica-interp_epoch-filt-", filter_freq_band, "-epo")
    )

  fs::dir_create(fs::path_dir(output_file))
  epochs_data$save(output_file, overwrite = TRUE)

  output_file
}
