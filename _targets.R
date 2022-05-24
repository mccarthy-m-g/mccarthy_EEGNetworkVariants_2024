# Follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load R packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(here)
library(tibble)
library(rbbt)
library(dplyr)

# Load Python packages required to define the pipeline:
library(reticulate)
use_miniconda("r-reticulate-mne")
mne <- import("mne")

# Set target options:
tar_option_set(
  packages = c("reticulate", "tidyverse"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# Load the R scripts with your custom functions:
#for (file in list.files("R", full.names = TRUE)) source(file) # TODO: Uncomment this once all done
source("R/preprocessing.R")
source("R/connectivity.R")
source("R/similarity.R")
source("R/descriptives.R")
source("R/references.R")

# Pipeline params:
eeg_recordings_unavailable <- FALSE # TODO: Figure out how to implement this, either with tar_skip() or an ifelse
slices <- 17:24 # currently inactive

participants_final <- c(
  # P01 Excluded for not having a follow-up session
  # P02 Excluded for not having a follow-up session
  "P03",
  "P04",
  # P05 Excluded for not having a follow-up session
  "P06",
  "P07",
  "P08",
  "P09",
  # P10 Excluded for having an unusable recording
  "P11",
  # P12 Did not participate in the study
  # P13 Excluded for having two unusable recordings
  "P14",
  # P15 Excluded for not having data from the follow-up session
  "P16",
  "P17",
  # P18 Excluded for having two unusable recordings
  "P19",
  "P20",
  "P21",
  "P22"
)

# FIXME: Some of the changes I made to preprocessing while testing the targets
# pipeline changed the results. Should revert back to the original steps to fix.

# Static branching
preprocessing_targets <- list(
  # Get file handler containing preprocessing parameters ----
  tar_target(
    file_handler_path,
    here("data", "file-handler.csv"),
    format = "file"
  ),
  tar_target(
    file_handler,
    read_csv(file_handler_path, col_types = cols())
  ),
  tar_target(
    file_handler_final, # Exclusions applied
    filter(file_handler, participant %in% participants_final)
  ),
  # Reference ----
  tar_target(
    channel_dictionary,
    here("data", "EEG-recordings", "channel-dictionary.csv"),
    format = "file"
  ),
  # FIXME This doesn't work with tar_files_input
  tar_files(
    raw_ref_input,
    file_handler_final[["raw"]] # [slices]
  ),
  tar_target(
    raw_ref,
    preprocess_reference(
      raw_ref_input,
      channel_dictionary
    ),
    pattern = map(raw_ref_input),
    format = "file"
  ),
  # Filter and downsample ----
  tar_target(
    raw_ref_filt_ds,
    preprocess_filter_downsample(
      raw_ref,
      l_freq = 0.1,
      h_freq = 50,
      phase = "zero-double",
      resample_rate = 200
    ),
    pattern = map(raw_ref),
    format = "file"
  ),
  # ICA ----
  tar_target(
    bad_channels,
    file_handler_final[["bad_channels"]] # [slices]
  ),
  tar_target(
    annotation_onsets,
    file_handler_final[["annotation_onsets"]] # [slices]
  ),
  tar_target(
    annotation_durations,
    file_handler_final[["annotation_durations"]] # [slices]
  ),
  tar_target(
    annotation_labels,
    file_handler_final[["annotation_labels"]] # [slices]
  ),
  tar_target(
    bad_ica_indices,
    file_handler_final[["bad_indices"]] # [slices]
  ),
  tar_target(
    raw_ref_filt_ds_ica,
    preprocess_ica(
      raw_ref_filt_ds,
      bad_channels,
      annotation_onsets,
      annotation_durations,
      annotation_labels,
      ica_n_components = 0.995,
      ica_method = "picard",
      ica_seed = 666L,
      bad_ica_indices
    ),
    pattern = map(
      raw_ref_filt_ds,
      bad_channels,
      annotation_onsets,
      annotation_durations,
      annotation_labels,
      bad_ica_indices
    ),
    format = "file"
  )
)

# TODO: Calculate amplitude connectivity.
# Might be able to use tar_map() and metaprogramming to do both phase and amplitude in a single target
# TODO: See if this can be split up into separate tar_map functions
connectivity_estimation_targets <- list(
  tar_map(
    values = tibble(
      filter_freq_band = c("delta", "theta", "alpha", "beta", "gamma"),
      filter_l_freq    = c(1, 4,  8, 13, 30),
      filter_h_freq    = c(4, 8, 13, 30, 50)
    ),
    names = filter_freq_band,
    # Epoch and filter ----
    tar_target(
      raw_ref_filt_ds_ica_epoch_filt,
      preprocess_filter_epoch(
        raw_ref_filt_ds_ica,
        filter_freq_band,
        filter_l_freq,
        filter_h_freq
      ),
      pattern = map(raw_ref_filt_ds_ica),
      format = "file"
    ),
    # Estimate phase coupling ----
    # FIXME: This is giving different values from my original estimates.
    # Might be caused by:
    # - changing the montage
    # - changing the order of some steps
    # - screwing up the code somewhere
    # The recordings are 99.99% identical but there's maybe a slight difference
    # in amplitude scaling of some channels, so it's probably the montage or referencing.
    tar_target(
      phase_coupling,
      estimate_phase_coupling(raw_ref_filt_ds_ica_epoch_filt),
      pattern = map(raw_ref_filt_ds_ica_epoch_filt),
      format = "file"
    ),
    tar_target(
      phase_connectivity_matrix,
      get_connectivity_matrix(phase_coupling),
      pattern = map(phase_coupling),
      iteration = "list"
    ),
    ########## TODO: Plotting connectivity!
    # Estimate similarity ----
    tar_target(
      phase_similarity,
      estimate_similarity(phase_connectivity_matrix)
    ),
    tar_target(
      phase_similarity_plot,
      plot_similarity(phase_similarity, rv)
    ),
    tar_map(
      values = tibble(
        participants = participants_final
      ),
      tar_target(
        phase_similarity_key_plot,
        plot_similarity_key(phase_similarity, participants)
      ),
      tar_target(
        phase_similarity_highlight_plot,
        plot_similarity_highlight(phase_similarity, participants)
      )
    )
  ),
  tar_target(
    similarity_archetype_plot,
    plot_similarity_archetype(phase_similarity_gamma)
  )
)

descriptives_targets <- list(
  # Participant descriptive statistics ----
  tar_target(
    participant_descriptives_path,
    here("data", "participant-descriptives.csv"),
    format = "file"
  ),
  ## Before exclusions
  tar_target(
    participant_descriptives,
    summarize_participant_descriptives(participant_descriptives_path)
  ),
  ## After exclusions
  tar_target(
    participant_descriptives_final,
    summarize_participant_descriptives(
      participant_descriptives_path,
      participants = participants_final
    )
  )
  # TODO: Preprocessing descriptive statistics
)

manuscripts_targets <- list(
  # Child documents for different sections of the manuscript ----
  tar_target(
    introduction,
    here("manuscripts", "child-documents", "introduction.Rmd"),
    format = "file"
  ),
  tar_target(
    methods,
    here("manuscripts", "child-documents", "methods.Rmd"),
    format = "file"
  ),
  # Bibliography file ----
  tar_target(
    references,
    write_bib(
      here("manuscripts", "references.json"),
      keys = c(
        introduction,
        methods
      ),
      public_library_id = "mccarthymg"
    ),
    format = "file"
  ),
  # Thesis (papaja) ----
  ## Note: The following error is thrown at the start of the pipeline because
  ## this target uses params and child documents: `Error in eval(x, envir =
  ## envir) : object 'params' not found`. The source of the error is explained
  ## in the following GitHub issue
  ## https://github.com/ropensci/targets/issues/256. Importantly, this error can
  ## be ignored as the report renders as intended in spite of it.
  tar_render(
    thesis,
    here("manuscripts", "thesis", "index.Rmd"),
    params = list(
      introduction_path = introduction,
      methods_path      = methods,
      references_path   = references
    )
  )
)

list(
  preprocessing_targets,
  connectivity_estimation_targets,
  descriptives_targets,
  manuscripts_targets
)
