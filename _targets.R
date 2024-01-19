# Follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load R packages required to define the pipeline:
library(targets)
library(tarchetypes)

library(here)
library(fs)
library(tidyverse)

library(FactoMineR)
library(glmmTMB)
library(emmeans)
library(performance)

library(ggdist)
library(ggh4x)
library(ggnewscale)
library(patchwork)
library(see)

library(flextable)
library(ftExtra)
library(gtsummary)

library(rmarkdown)
library(knitr)
library(officer)
library(officedown)
library(papaja)
library(rbbt)

# Load Python packages required to define the pipeline:
library(reticulate)
use_miniconda("r-reticulate-mne")
mne <- import("mne")
pli <- import_from_path("pli", "python")

# Set target options:
tar_option_set(
  packages = c("reticulate", "tidyverse"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# Load the R scripts with your custom functions:
source("R/preprocessing.R")
source("R/connectivity.R")
source("R/similarity.R")
source("R/descriptives.R")
source("R/references.R")
source("R/session_info.R")
source("R/utils.R")

# Pipeline params:
eeg_recordings_unavailable <- FALSE # TODO: Figure out how to implement this, either with tar_skip() or an ifelse
participant_descriptives_available <- TRUE
slices <- 17:24 # currently inactive

participants_final <- c(
  # P01 Excluded for not having a follow-up session
  # P02 Excluded for not having a follow-up session
  "P03",
  "P04", # Consider excluding because fu rc1, fu ro1 have bad channels next to each other
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
  "P20", # Consider excluding because fu ro2, post rc1, post ro1, post ro1, pre rc2, pre ro1, pre ro2 have bad channels next to each other
  "P21",
  "P22"
)

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
  tar_files(
    raw_input,
    file_handler_final[["raw"]] # [slices]
  ),
  # Filter and downsample ----
  tar_target(
    raw_filt_ds,
    preprocess_filter_downsample(
      raw_input,
      l_freq = 0.1,
      h_freq = 50,
      phase = "zero-double",
      resample_rate = 200
    ),
    pattern = map(raw_input),
    format = "file"
  ),
  # Re-reference and rename channels ----
  tar_target(
    channel_dictionary,
    here("data", "EEG-recordings", "channel-dictionary.csv"),
    format = "file"
  ),
  tar_target(
    raw_filt_ds_reref,
    preprocess_rereference(
      raw_filt_ds,
      channel_dictionary
    ),
    pattern = map(raw_filt_ds),
    format = "file"
  ),
  # ICA then interpolate bad channels ----
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
    raw_filt_ds_reref_ica_interp,
    preprocess_ica(
      raw_filt_ds_reref,
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
      raw_filt_ds_reref,
      bad_channels,
      annotation_onsets,
      annotation_durations,
      annotation_labels,
      bad_ica_indices
    ),
    format = "file"
  )
)

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
    # Epoch then filter into frequency bands ----
    tar_target(
      raw_filt_ds_reref_ica_interp_epoch_filt,
      preprocess_filter_epoch(
        raw_filt_ds_reref_ica_interp,
        filter_freq_band,
        filter_l_freq,
        filter_h_freq
      ),
      pattern = map(raw_filt_ds_reref_ica_interp),
      format = "file"
    ),
    # Estimate phase coupling ----
    tar_target(
      phase_coupling,
      estimate_phase_coupling(raw_filt_ds_reref_ica_interp_epoch_filt),
      pattern = map(raw_filt_ds_reref_ica_interp_epoch_filt),
      format = "file"
    ),
    tar_target(
      phase_connectivity_matrix,
      get_connectivity_matrix(phase_coupling),
      pattern = map(phase_coupling),
      iteration = "list"
    ),
    tar_target(
      phase_connectivity_summary,
      summarize_connectivity(phase_connectivity_matrix)
    ),
    tar_target(
      phase_connectivity_plot,
      plot_connectivity(phase_connectivity_matrix, "PLI"),
      pattern = map(phase_connectivity_matrix),
      iteration = "list"
    ),
    tar_target(
      phase_connectivity_profile_plot,
      plot_connectivity_profiles(phase_connectivity_matrix, "PLI")
    ),
    tar_target(
      phase_connectivity_histogram_plot,
      plot_connectivity_histogram(phase_connectivity_matrix, "PLI")
    ),
    tar_target(
      phase_connectivity_patchwork_plot,
      plot_connectivity_patchwork(phase_connectivity_matrix, "PLI")
    ),
    # Estimate phase connectivity similarity ----
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
    ),
    # Model phase connectivity similarity ----
    tar_target(
      phase_similarity_glmmTMB,
      glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | x_label) + (1 | y_label), # Random effects
        data = phase_similarity
      )
    ),
    tar_target(
      phase_similarity_emmeans,
      emmeans_similarity(phase_similarity_glmmTMB)
    ),
    tar_target(
      phase_similarity_emmeans_tidy,
      tidy_emmeans_similarity(phase_similarity_emmeans)
    ),
    tar_target(
      phase_similarity_contrasts,
      contrast_similarity(phase_similarity_emmeans)
    ),
    tar_target(
      phase_similarity_contrasts_tidy,
      tidy_contrast_similarity(phase_similarity_contrasts)
    ),
    tar_target(
      phase_similarity_contrasts_cohens_d,
      contrast_similarity_cohens_d(
        phase_similarity_glmmTMB,
        phase_similarity_contrasts
      )
    ),
    tar_target(
      phase_similarity_contrasts_cohens_d_tidy,
      tidy_contrast_similarity_cohens_d(
        phase_similarity_contrasts_cohens_d
      )
    ),
    tar_target(
      phase_similarity_contrasts_plot,
      plot_similarity_contrasts(phase_similarity_contrasts)
    ),
    tar_target(
      phase_similarity_contrasts_table_nhst,
      make_contrast_results_table_nhst(
        phase_similarity_emmeans_tidy,
        phase_similarity_contrasts_tidy
      )
    ),
    tar_target(
      phase_similarity_contrasts_table_smd,
      make_contrast_results_table_smd(
        phase_similarity_glmmTMB,
        phase_similarity_emmeans_tidy,
        phase_similarity_contrasts_tidy,
        phase_similarity_contrasts_cohens_d_tidy
      )
    ),
    tar_target(
      phase_similarity_model_fit_table,
      make_model_fit_table(phase_similarity_glmmTMB)
    ),
    ## Model diagnostics
    tar_target(
      phase_similarity_glmmTMB_ppreds,
      check_predictions(phase_similarity_glmmTMB)
    ),
    tar_target(
      phase_similarity_glmmTMB_uniformity,
      check_uniformity(phase_similarity_glmmTMB)
    ),
    tar_target(
      phase_similarity_glmmTMB_randnorm,
      check_model(phase_similarity_glmmTMB, check = c("reqq"), panel = FALSE)
    ),
    # TODO: Upgrade `performance` to 9.1 so I get CIs for these
    tar_target(
      phase_similarity_glmmTMB_collinearity,
      check_collinearity(phase_similarity_glmmTMB)
    ),
    # Fit phase connectivity similarity submodels ----
    tar_target(
      phase_similarity_subset,
      subset_similarity_results(phase_similarity, participants_final)
    ),
    tar_target(
      phase_similarity_subset_glmmTMB,
      subset_glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | x_label) + (1 | y_label), # Random effects
        data = phase_similarity_subset
      )
    ),
    tar_target(
      phase_similarity_subset_emmeans,
      subset_emmeans_similarity(phase_similarity_subset_glmmTMB)
    ),
    tar_target(
      phase_similarity_subset_contrasts,
      subset_contrast_similarity(phase_similarity_subset_emmeans)
    ),
    tar_target(
      phase_similarity_subset_contrasts_plot,
      plot_subset_similarity_contrasts(
        phase_similarity_contrasts,
        phase_similarity_subset_contrasts,
        phase_similarity_subset_glmmTMB
      )
    ),
    tar_target(
      phase_similarity_subset_contrasts_plot_2,
      plot_subset_similarity_contrasts_2(
        phase_similarity_contrasts,
        phase_similarity_subset_contrasts,
        phase_similarity_subset_glmmTMB
      )
    ),
    # Estimate amplitude coupling ----
    tar_target(
      amplitude_coupling,
      estimate_amplitude_coupling(
        raw_filt_ds_reref_ica_interp_epoch_filt,
        absolute_last = TRUE
      ),
      pattern = map(raw_filt_ds_reref_ica_interp_epoch_filt),
      format = "file"
    ),
    tar_target(
      amplitude_connectivity_matrix,
      get_connectivity_matrix(amplitude_coupling),
      pattern = map(amplitude_coupling),
      iteration = "list"
    ),
    tar_target(
      amplitude_connectivity_summary,
      summarize_connectivity(amplitude_connectivity_matrix)
    ),
    tar_target(
      amplitude_connectivity_plot,
      plot_connectivity(amplitude_connectivity_matrix, "AEC"),
      pattern = map(amplitude_connectivity_matrix),
      iteration = "list"
    ),
    tar_target(
      amplitude_connectivity_profile_plot,
      plot_connectivity_profiles(amplitude_connectivity_matrix, "AEC")
    ),
    tar_target(
      amplitude_connectivity_histogram_plot,
      plot_connectivity_histogram(amplitude_connectivity_matrix, "AEC")
    ),
    tar_target(
      amplitude_connectivity_patchwork_plot,
      plot_connectivity_patchwork(amplitude_connectivity_matrix, "AEC")
    ),
    # Estimate amplitude coupling similarity ----
    tar_target(
      amplitude_similarity,
      estimate_similarity(amplitude_connectivity_matrix)
    ),
    tar_target(
      amplitude_similarity_plot,
      plot_similarity(amplitude_similarity, rv)
    ),
    tar_map(
      values = tibble(
        participants = participants_final
      ),
      tar_target(
        amplitude_similarity_key_plot,
        plot_similarity_key(amplitude_similarity, participants)
      ),
      tar_target(
        amplitude_similarity_highlight_plot,
        plot_similarity_highlight(amplitude_similarity, participants)
      )
    ),
    # Model amplitude connectivity similarity ----
    tar_target(
      amplitude_similarity_glmmTMB,
      glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | x_label) + (1 | y_label), # Random effects
        data = amplitude_similarity
      )
    ),
    tar_target(
      amplitude_similarity_emmeans,
      emmeans_similarity(amplitude_similarity_glmmTMB)
    ),
    tar_target(
      amplitude_similarity_emmeans_tidy,
      tidy_emmeans_similarity(amplitude_similarity_emmeans)
    ),
    tar_target(
      amplitude_similarity_contrasts,
      contrast_similarity(amplitude_similarity_emmeans)
    ),
    tar_target(
      amplitude_similarity_contrasts_tidy,
      tidy_contrast_similarity(amplitude_similarity_contrasts)
    ),
    tar_target(
      amplitude_similarity_contrasts_cohens_d,
      contrast_similarity_cohens_d(
        amplitude_similarity_glmmTMB,
        amplitude_similarity_contrasts
      )
    ),
    tar_target(
      amplitude_similarity_contrasts_cohens_d_tidy,
      tidy_contrast_similarity_cohens_d(
        amplitude_similarity_contrasts_cohens_d
      )
    ),
    tar_target(
      amplitude_similarity_contrasts_plot,
      plot_similarity_contrasts(amplitude_similarity_contrasts)
    ),
    tar_target(
      amplitude_similarity_contrasts_table_nhst,
      make_contrast_results_table_nhst(
        amplitude_similarity_emmeans_tidy,
        amplitude_similarity_contrasts_tidy
      )
    ),
    tar_target(
      amplitude_similarity_contrasts_table_smd,
      make_contrast_results_table_smd(
        amplitude_similarity_glmmTMB,
        amplitude_similarity_emmeans_tidy,
        amplitude_similarity_contrasts_tidy,
        amplitude_similarity_contrasts_cohens_d_tidy
      )
    ),
    tar_target(
      amplitude_similarity_model_fit_table,
      make_model_fit_table(amplitude_similarity_glmmTMB)
    ),
    ## TODO: Model diagnostics


    # Fit amplitude connectivity similarity submodels ----
    tar_target(
      amplitude_similarity_subset,
      subset_similarity_results(amplitude_similarity, participants_final)
    ),
    tar_target(
      amplitude_similarity_subset_glmmTMB,
      subset_glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | x_label) + (1 | y_label), # Random effects
        data = amplitude_similarity_subset
      )
    ),
    tar_target(
      amplitude_similarity_subset_emmeans,
      subset_emmeans_similarity(amplitude_similarity_subset_glmmTMB)
    ),
    tar_target(
      amplitude_similarity_subset_contrasts,
      subset_contrast_similarity(amplitude_similarity_subset_emmeans)
    ),
    tar_target(
      amplitude_similarity_subset_contrasts_plot,
      plot_subset_similarity_contrasts(
        amplitude_similarity_contrasts,
        amplitude_similarity_subset_contrasts,
        amplitude_similarity_subset_glmmTMB
      )
    ),
    tar_target(
      amplitude_similarity_subset_contrasts_plot_2,
      plot_subset_similarity_contrasts_2(
        amplitude_similarity_contrasts,
        amplitude_similarity_subset_contrasts,
        amplitude_similarity_subset_glmmTMB
      )
    ),
    # Estimate phase coupling (Hilbert transform) ----
    tar_target(
      phase_coupling_hilbert,
      estimate_phase_coupling_hilbert(
        raw_filt_ds_reref_ica_interp_epoch_filt,
        absolute_last = TRUE
      ),
      pattern = map(raw_filt_ds_reref_ica_interp_epoch_filt),
      format = "file"
    ),
    tar_target(
      phase_connectivity_matrix_hilbert,
      get_connectivity_matrix(phase_coupling_hilbert),
      pattern = map(phase_coupling_hilbert),
      iteration = "list"
    ),
    tar_target(
      phase_connectivity_profile_plot_hilbert,
      plot_connectivity_profiles(phase_connectivity_matrix_hilbert, "PLI (Hilbert)")
    ),
    tar_target(
      phase_connectivity_patchwork_plot_hilbert,
      plot_connectivity_patchwork(phase_connectivity_matrix_hilbert, "PLI")
    ),
    # Estimate phase connectivity similarity (Hilbert transform) ----
    tar_target(
      phase_similarity_hilbert,
      estimate_similarity(phase_connectivity_matrix_hilbert)
    ),
    tar_target(
      phase_similarity_plot_hilbert,
      plot_similarity(phase_similarity_hilbert, rv)
    ),
    # Model phase connectivity similarity (Hilbert transform) ----
    tar_target(
      phase_similarity_glmmTMB_hilbert,
      glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | x_label) + (1 | y_label), # Random effects
        data = phase_similarity_hilbert
      )
    ),
    tar_target(
      phase_similarity_emmeans_hilbert,
      emmeans_similarity(phase_similarity_glmmTMB_hilbert)
    ),
    tar_target(
      phase_similarity_emmeans_tidy_hilbert,
      tidy_emmeans_similarity(phase_similarity_emmeans_hilbert)
    ),
    tar_target(
      phase_similarity_contrasts_hilbert,
      contrast_similarity(phase_similarity_emmeans_hilbert)
    ),
    tar_target(
      phase_similarity_contrasts_tidy_hilbert,
      tidy_contrast_similarity(phase_similarity_contrasts_hilbert)
    ),
    tar_target(
      phase_similarity_contrasts_cohens_d_hilbert,
      contrast_similarity_cohens_d(
        phase_similarity_glmmTMB_hilbert,
        phase_similarity_contrasts_hilbert
      )
    ),
    tar_target(
      phase_similarity_contrasts_cohens_d_tidy_hilbert,
      tidy_contrast_similarity_cohens_d(
        phase_similarity_contrasts_cohens_d_hilbert
      )
    ),
    tar_target(
      phase_similarity_contrasts_plot_hilbert,
      plot_similarity_contrasts(phase_similarity_contrasts_hilbert)
    ),
    tar_target(
      phase_similarity_contrasts_table_nhst_hilbert,
      make_contrast_results_table_nhst(
        phase_similarity_emmeans_tidy_hilbert,
        phase_similarity_contrasts_tidy_hilbert
      )
    ),
    tar_target(
      phase_similarity_contrasts_table_smd_hilbert,
      make_contrast_results_table_smd(
        phase_similarity_glmmTMB_hilbert,
        phase_similarity_emmeans_tidy_hilbert,
        phase_similarity_contrasts_tidy_hilbert,
        phase_similarity_contrasts_cohens_d_tidy_hilbert
      )
    ),
    tar_target(
      phase_similarity_model_fit_table_hilbert,
      make_model_fit_table(phase_similarity_glmmTMB_hilbert)
    ),
    # Fit phase connectivity similarity submodels (Hilbert transform) ----
    tar_target(
      phase_similarity_subset_hilbert,
      subset_similarity_results(phase_similarity_hilbert, participants_final)
    ),
    tar_target(
      phase_similarity_subset_glmmTMB_hilbert,
      subset_glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | x_label) + (1 | y_label), # Random effects
        data = phase_similarity_subset_hilbert
      )
    ),
    tar_target(
      phase_similarity_subset_emmeans_hilbert,
      subset_emmeans_similarity(phase_similarity_subset_glmmTMB_hilbert)
    ),
    tar_target(
      phase_similarity_subset_contrasts_hilbert,
      subset_contrast_similarity(phase_similarity_subset_emmeans_hilbert)
    ),
    tar_target(
      phase_similarity_subset_contrasts_plot_hilbert,
      plot_subset_similarity_contrasts(
        phase_similarity_contrasts_hilbert,
        phase_similarity_subset_contrasts_hilbert,
        phase_similarity_subset_glmmTMB_hilbert
      )
    ),
    # Model phase connectivity similarity (maximal model) ----
    tar_target(
      phase_similarity_glmmTMB_maximal,
      glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | pair_participant) + (1 | x_label) + (1 | y_label), # Random effects
        data = phase_similarity
      )
    ),
    tar_target(
      phase_similarity_emmeans_maximal,
      emmeans_similarity(phase_similarity_glmmTMB_maximal)
    ),
    tar_target(
      phase_similarity_emmeans_tidy_maximal,
      tidy_emmeans_similarity(phase_similarity_emmeans_maximal)
    ),
    tar_target(
      phase_similarity_contrasts_maximal,
      contrast_similarity(phase_similarity_emmeans_maximal)
    ),
    tar_target(
      phase_similarity_contrasts_tidy_maximal,
      tidy_contrast_similarity(phase_similarity_contrasts_maximal)
    ),
    tar_target(
      phase_similarity_contrasts_cohens_d_maximal,
      contrast_similarity_cohens_d(
        phase_similarity_glmmTMB_maximal,
        phase_similarity_contrasts_maximal
      )
    ),
    tar_target(
      phase_similarity_contrasts_cohens_d_tidy_maximal,
      tidy_contrast_similarity_cohens_d(
        phase_similarity_contrasts_cohens_d_maximal
      )
    ),
    tar_target(
      phase_similarity_contrasts_plot_maximal,
      plot_similarity_contrasts(phase_similarity_contrasts_maximal)
    ),
    tar_target(
      phase_similarity_contrasts_table_nhst_maximal,
      make_contrast_results_table_nhst(
        phase_similarity_emmeans_tidy_maximal,
        phase_similarity_contrasts_tidy_maximal
      )
    ),
    tar_target(
      phase_similarity_contrasts_table_smd_maximal,
      make_contrast_results_table_smd(
        phase_similarity_glmmTMB_maximal,
        phase_similarity_emmeans_tidy_maximal,
        phase_similarity_contrasts_tidy_maximal,
        phase_similarity_contrasts_cohens_d_tidy_maximal
      )
    ),
    tar_target(
      phase_similarity_model_fit_table_maximal,
      make_model_fit_table(
        phase_similarity_glmmTMB_maximal, maximal_model = TRUE
      )
    ),
    # Fit phase connectivity similarity submodels (maximal models) ----
    tar_target(
      phase_similarity_subset_glmmTMB_maximal,
      subset_glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | pair_participant) + (1 | x_label) + (1 | y_label), # Random effects
        data = phase_similarity_subset
      )
    ),
    tar_target(
      phase_similarity_subset_emmeans_maximal,
      subset_emmeans_similarity(phase_similarity_subset_glmmTMB_maximal)
    ),
    tar_target(
      phase_similarity_subset_contrasts_maximal,
      subset_contrast_similarity(phase_similarity_subset_emmeans_maximal)
    ),
    tar_target(
      phase_similarity_subset_contrasts_plot_maximal,
      plot_subset_similarity_contrasts(
        phase_similarity_contrasts_maximal,
        phase_similarity_subset_contrasts_maximal,
        phase_similarity_subset_glmmTMB_maximal
      )
    ),
    # Model amplitude connectivity similarity (maximal model) ----
    tar_target(
      amplitude_similarity_glmmTMB_maximal,
      glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | pair_participant) + (1 | x_label) + (1 | y_label), # Random effects
        data = amplitude_similarity
      )
    ),
    tar_target(
      amplitude_similarity_emmeans_maximal,
      emmeans_similarity(amplitude_similarity_glmmTMB_maximal)
    ),
    tar_target(
      amplitude_similarity_emmeans_tidy_maximal,
      tidy_emmeans_similarity(amplitude_similarity_emmeans_maximal)
    ),
    tar_target(
      amplitude_similarity_contrasts_maximal,
      contrast_similarity(amplitude_similarity_emmeans_maximal)
    ),
    tar_target(
      amplitude_similarity_contrasts_tidy_maximal,
      tidy_contrast_similarity(amplitude_similarity_contrasts_maximal)
    ),
    tar_target(
      amplitude_similarity_contrasts_cohens_d_maximal,
      contrast_similarity_cohens_d(
        amplitude_similarity_glmmTMB_maximal,
        amplitude_similarity_contrasts_maximal
      )
    ),
    tar_target(
      amplitude_similarity_contrasts_cohens_d_tidy_maximal,
      tidy_contrast_similarity_cohens_d(
        amplitude_similarity_contrasts_cohens_d_maximal
      )
    ),
    tar_target(
      amplitude_similarity_contrasts_plot_maximal,
      plot_similarity_contrasts(amplitude_similarity_contrasts_maximal)
    ),
    tar_target(
      amplitude_similarity_contrasts_table_nhst_maximal,
      make_contrast_results_table_nhst(
        amplitude_similarity_emmeans_tidy_maximal,
        amplitude_similarity_contrasts_tidy_maximal
      )
    ),
    tar_target(
      amplitude_similarity_contrasts_table_smd_maximal,
      make_contrast_results_table_smd(
        amplitude_similarity_glmmTMB_maximal,
        amplitude_similarity_emmeans_tidy_maximal,
        amplitude_similarity_contrasts_tidy_maximal,
        amplitude_similarity_contrasts_cohens_d_tidy_maximal
      )
    ),
    tar_target(
      amplitude_similarity_model_fit_table_maximal,
      make_model_fit_table(
        amplitude_similarity_glmmTMB_maximal, maximal_model = TRUE
      )
    ),
    # Fit amplitude connectivity similarity submodels (maximal models) ----
    tar_target(
      amplitude_similarity_subset_glmmTMB_maximal,
      subset_glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | pair_participant) + (1 | x_label) + (1 | y_label), # Random effects
        data = amplitude_similarity_subset
      )
    ),
    tar_target(
      amplitude_similarity_subset_emmeans_maximal,
      subset_emmeans_similarity(amplitude_similarity_subset_glmmTMB_maximal)
    ),
    tar_target(
      amplitude_similarity_subset_contrasts_maximal,
      subset_contrast_similarity(amplitude_similarity_subset_emmeans_maximal)
    ),
    tar_target(
      amplitude_similarity_subset_contrasts_plot_maximal,
      plot_subset_similarity_contrasts(
        amplitude_similarity_contrasts_maximal,
        amplitude_similarity_subset_contrasts_maximal,
        amplitude_similarity_subset_glmmTMB_maximal
      )
    ),
    # Model phase connectivity similarity (Hilbert transform, maximal model) ----
    tar_target(
      phase_similarity_glmmTMB_hilbert_maximal,
      glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | pair_participant) + (1 | x_label) + (1 | y_label), # Random effects
        data = phase_similarity_hilbert
      )
    ),
    tar_target(
      phase_similarity_emmeans_hilbert_maximal,
      emmeans_similarity(phase_similarity_glmmTMB_hilbert_maximal)
    ),
    tar_target(
      phase_similarity_emmeans_tidy_hilbert_maximal,
      tidy_emmeans_similarity(phase_similarity_emmeans_hilbert_maximal)
    ),
    tar_target(
      phase_similarity_contrasts_hilbert_maximal,
      contrast_similarity(phase_similarity_emmeans_hilbert_maximal)
    ),
    tar_target(
      phase_similarity_contrasts_tidy_hilbert_maximal,
      tidy_contrast_similarity(phase_similarity_contrasts_hilbert_maximal)
    ),
    tar_target(
      phase_similarity_contrasts_cohens_d_hilbert_maximal,
      contrast_similarity_cohens_d(
        phase_similarity_glmmTMB_hilbert_maximal,
        phase_similarity_contrasts_hilbert_maximal
      )
    ),
    tar_target(
      phase_similarity_contrasts_cohens_d_tidy_hilbert_maximal,
      tidy_contrast_similarity_cohens_d(
        phase_similarity_contrasts_cohens_d_hilbert_maximal
      )
    ),
    tar_target(
      phase_similarity_contrasts_plot_hilbert_maximal,
      plot_similarity_contrasts(phase_similarity_contrasts_hilbert_maximal)
    ),
    tar_target(
      phase_similarity_contrasts_table_nhst_hilbert_maximal,
      make_contrast_results_table_nhst(
        phase_similarity_emmeans_tidy_hilbert_maximal,
        phase_similarity_contrasts_tidy_hilbert_maximal
      )
    ),
    tar_target(
      phase_similarity_contrasts_table_smd_hilbert_maximal,
      make_contrast_results_table_smd(
        phase_similarity_glmmTMB_hilbert_maximal,
        phase_similarity_emmeans_tidy_hilbert_maximal,
        phase_similarity_contrasts_tidy_hilbert_maximal,
        phase_similarity_contrasts_cohens_d_tidy_hilbert_maximal
      )
    ),
    tar_target(
      phase_similarity_model_fit_table_hilbert_maximal,
      make_model_fit_table(
        phase_similarity_glmmTMB_hilbert_maximal, maximal_model = TRUE
      )
    ),
    # Fit phase connectivity similarity submodels (Hilbert transform, maximal models) ----
    tar_target(
      phase_similarity_subset_glmmTMB_hilbert_maximal,
      subset_glmmTMB_similarity(
        rv ~
          within_participant * within_session * within_state + # Fixed effects
          (1 | pair_participant) + (1 | x_label) + (1 | y_label), # Random effects
        data = phase_similarity_subset_hilbert
      )
    ),
    tar_target(
      phase_similarity_subset_emmeans_hilbert_maximal,
      subset_emmeans_similarity(phase_similarity_subset_glmmTMB_hilbert_maximal)
    ),
    tar_target(
      phase_similarity_subset_contrasts_hilbert_maximal,
      subset_contrast_similarity(phase_similarity_subset_emmeans_hilbert_maximal)
    ),
    tar_target(
      phase_similarity_subset_contrasts_plot_hilbert_maximal,
      plot_subset_similarity_contrasts(
        phase_similarity_contrasts_hilbert_maximal,
        phase_similarity_subset_contrasts_hilbert_maximal,
        phase_similarity_subset_glmmTMB_hilbert_maximal
      )
    ),
    # Save manuscript figures ----
    tar_target(
      phase_similarity_matrix_figure,
      save_similarity_matrix_figure(
        paste0("figures/", filter_freq_band, "/phase-similarity-matrix.png"),
        phase_similarity_plot
      ),
      format = "file"
    ),
    tar_target(
      phase_similarity_results_figure,
      save_results_figure(
        paste0("figures/", filter_freq_band, "/phase_similarity_results.png"),
        phase_connectivity_profile_plot,
        phase_similarity_plot,
        phase_similarity_contrasts_plot
      ),
      format = "file"
    ),
    tar_target(
      phase_similarity_subset_results_figure,
      save_figure(
        paste0("figures/", filter_freq_band, "/phase_similarity_subset_results.png"),
        phase_similarity_subset_contrasts_plot_2,
        height = 9
      ),
      format = "file"
    ),
    tar_target(
      phase_similarity_results_figure_maximal,
      save_results_figure(
        paste0("figures/", filter_freq_band, "/phase_similarity_results_maximal.png"),
        phase_connectivity_profile_plot,
        phase_similarity_plot,
        phase_similarity_contrasts_plot_maximal,
        phase_similarity_subset_contrasts_plot_maximal
      ),
      format = "file"
    ),
    tar_target(
      amplitude_similarity_matrix_figure,
      save_similarity_matrix_figure(
        paste0("figures/", filter_freq_band, "/amplitude-similarity-matrix.png"),
        amplitude_similarity_plot
      ),
      format = "file"
    ),
    tar_target(
      amplitude_similarity_results_figure,
      save_results_figure(
        paste0("figures/", filter_freq_band, "/amplitude_similarity_results.png"),
        amplitude_connectivity_profile_plot,
        amplitude_similarity_plot,
        amplitude_similarity_contrasts_plot
      ),
      format = "file"
    ),
    tar_target(
      amplitude_similarity_subset_results_figure,
      save_figure(
        paste0("figures/", filter_freq_band, "/amplitude_similarity_subset_results.png"),
        amplitude_similarity_subset_contrasts_plot_2,
        height = 9
      ),
      format = "file"
    ),
    tar_target(
      amplitude_similarity_results_figure_maximal,
      save_results_figure(
        paste0("figures/", filter_freq_band, "/amplitude_similarity_results_maximal.png"),
        amplitude_connectivity_profile_plot,
        amplitude_similarity_plot,
        amplitude_similarity_contrasts_plot_maximal,
        amplitude_similarity_subset_contrasts_plot_maximal
      ),
      format = "file"
    ),
    tar_target(
      phase_similarity_matrix_figure_hilbert,
      save_similarity_matrix_figure(
        paste0("figures/", filter_freq_band, "/phase-similarity-matrix-hilbert.png"),
        phase_similarity_plot_hilbert
      ),
      format = "file"
    ),
    tar_target(
      phase_similarity_results_figure_hilbert,
      save_results_figure(
        paste0("figures/", filter_freq_band, "/phase_similarity_results_hilbert.png"),
        phase_connectivity_profile_plot_hilbert,
        phase_similarity_plot_hilbert,
        phase_similarity_contrasts_plot_hilbert,
        phase_similarity_subset_contrasts_plot_hilbert
      ),
      format = "file"
    ),
    tar_target(
      phase_similarity_results_figure_hilbert_maximal,
      save_results_figure(
        paste0("figures/", filter_freq_band, "/phase_similarity_results_hilbert_maximal.png"),
        phase_connectivity_profile_plot_hilbert,
        phase_similarity_plot_hilbert,
        phase_similarity_contrasts_plot_hilbert_maximal,
        phase_similarity_subset_contrasts_plot_hilbert_maximal
      ),
      format = "file"
    ),
    tar_target(
      phase_connectivity_patchwork_figure,
      save_connectivity_patchwork(
        phase_connectivity_patchwork_plot, filter_freq_band, "phase"
      ),
      format = "file"
    ),
    tar_target(
      amplitude_connectivity_patchwork_figure,
      save_connectivity_patchwork(
        amplitude_connectivity_patchwork_plot, filter_freq_band, "amplitude"
      ),
      format = "file"
    )
    ,
    tar_target(
      phase_connectivity_patchwork_figure_hilbert,
      save_connectivity_patchwork(
        phase_connectivity_patchwork_plot_hilbert, filter_freq_band, "phase-hilbert"
      ),
      format = "file"
    )
  ),
  tar_target(
    similarity_archetype_plot,
    plot_similarity_archetype(phase_similarity_gamma)
  ),
  tar_target(
    similarity_archetypes_figure,
    save_similarity_archetypes_figure(
      "figures/similarity_archetypes.png",
      similarity_archetype_plot
    ),
    format = "file"
  ),
  tar_target(
    amplitude_similarity_patchwork_figure,
    save_amplitude_similarity_patchwork_figure(
      amplitude_similarity_plot_delta,
      amplitude_similarity_plot_theta,
      amplitude_similarity_plot_beta,
      amplitude_similarity_plot_gamma
    ),
    format = "file"
  ),
  tar_target(
    jet_similarity_patchwork_figure,
    save_jet_similarity_patchwork_figure(
      phase_similarity_plot_delta,
      phase_similarity_plot_theta,
      phase_similarity_plot_alpha,
      phase_similarity_plot_beta,
      phase_similarity_plot_gamma,
      amplitude_similarity_plot_alpha
    ),
    format = "file"
  ),
  tar_target(
    phase_connectivity_histograms_plot,
    plot_connectivity_histograms(
      phase_connectivity_matrix_delta,
      phase_connectivity_matrix_theta,
      phase_connectivity_matrix_alpha,
      phase_connectivity_matrix_beta,
      phase_connectivity_matrix_gamma,
      "PLI"
    )
  ),
  tar_target(
    amplitude_connectivity_histograms_plot,
    plot_connectivity_histograms(
      amplitude_connectivity_matrix_delta,
      amplitude_connectivity_matrix_theta,
      amplitude_connectivity_matrix_alpha,
      amplitude_connectivity_matrix_beta,
      amplitude_connectivity_matrix_gamma,
      "AEC"
    )
  ),
  tar_target(
    connectivity_histogram_patchwork_figure,
    save_connectivity_histogram_patchwork(
      phase_connectivity_histograms_plot,
      amplitude_connectivity_histograms_plot
    ),
    format = "file"
  ),
  tar_target(
    illustrative_connectomes_figure,
    save_illustrative_connectomes_figure(
      phase_connectivity_matrix_delta,
      amplitude_connectivity_matrix_delta,
      participant = "P16",
      recording = "fu_ro2"
    ),
    format = "file"
  ),
  tar_target(
    phase_subset_similarity_contrast_barplot,
    plot_subset_similarity_contrast_bars(
      phase_similarity_subset_contrasts_delta,
      phase_similarity_subset_contrasts_theta,
      phase_similarity_subset_contrasts_alpha,
      phase_similarity_subset_contrasts_beta,
      phase_similarity_subset_contrasts_gamma
    )
  ),
  tar_target(
    phase_subset_similarity_contrast_barplot_figure,
    save_figure(
      paste0("figures/phase-similarity-subset-contrast-barplot.png"),
      phase_subset_similarity_contrast_barplot,
      height = 9,
      scaling = 1
    ),
    format = "file"
  ),
  tar_target(
    amplitude_subset_similarity_contrast_barplot,
    plot_subset_similarity_contrast_bars(
      alpha = amplitude_similarity_subset_contrasts_alpha,
    )
  ),
  tar_target(
    amplitude_subset_similarity_contrast_barplot_figure,
    save_figure(
      paste0("figures/amplitude-similarity-subset-contrast-barplot.png"),
      amplitude_subset_similarity_contrast_barplot,
      width = 21.59/2,
      height = 9,
      scaling = 1
    ),
    format = "file"
  ),
  tar_target(
    subset_similarity_contrast_barplot_figure,
    save_subset_similarity_contrast_barplot_figure(
      paste0("figures/similarity-subset-contrast-barplot.png"),
      phase_subset_similarity_contrast_barplot,
      amplitude_subset_similarity_contrast_barplot
    ),
    format = "file"
  ),
  # Make tables ----
  tar_target(
    phase_connectivity_summaries,
    summarize_connectivity_distributions(
      phase_connectivity_matrix_delta,
      phase_connectivity_matrix_theta,
      phase_connectivity_matrix_alpha,
      phase_connectivity_matrix_beta,
      phase_connectivity_matrix_gamma,
      "Phase"
    )
  ),
  tar_target(
    amplitude_connectivity_summaries,
    summarize_connectivity_distributions(
      amplitude_connectivity_matrix_delta,
      amplitude_connectivity_matrix_theta,
      amplitude_connectivity_matrix_alpha,
      amplitude_connectivity_matrix_beta,
      amplitude_connectivity_matrix_gamma,
      "Amplitude"
    )
  ),
  tar_target(
    connectivity_summary_table,
    make_connectivity_summary_table(
      phase_connectivity_summaries,
      amplitude_connectivity_summaries
    )
  ),
  tar_target(
    similarity_glmmTMB_convergence_table,
    make_convergence_table(
      phase_similarity_glmmTMB_delta,
      phase_similarity_glmmTMB_maximal_delta,
      phase_similarity_glmmTMB_theta,
      phase_similarity_glmmTMB_maximal_theta,
      phase_similarity_glmmTMB_alpha,
      phase_similarity_glmmTMB_maximal_alpha,
      phase_similarity_glmmTMB_beta,
      phase_similarity_glmmTMB_maximal_beta,
      phase_similarity_glmmTMB_gamma,
      phase_similarity_glmmTMB_maximal_gamma,
      amplitude_similarity_glmmTMB_alpha,
      amplitude_similarity_glmmTMB_maximal_alpha,
      phase_similarity_glmmTMB_hilbert_delta,
      phase_similarity_glmmTMB_hilbert_maximal_delta,
      phase_similarity_glmmTMB_hilbert_theta,
      phase_similarity_glmmTMB_hilbert_maximal_theta,
      phase_similarity_glmmTMB_hilbert_alpha,
      phase_similarity_glmmTMB_hilbert_maximal_alpha,
      phase_similarity_glmmTMB_hilbert_beta,
      phase_similarity_glmmTMB_hilbert_maximal_beta,
      phase_similarity_glmmTMB_hilbert_gamma,
      phase_similarity_glmmTMB_hilbert_maximal_gamma
    )
  )
)

descriptives_targets <- list(
  # Participant descriptive statistics ----
  tar_target(
    participant_descriptives_path,
    ifelse(
      participant_descriptives_available,
      here("data", "participant-descriptives.csv"),
      here("data", "participant-descriptives-synthetic.csv")
    ),
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
  ),
  # Preprocessing descriptive statistics ----
  tar_target(
    bad_channels_descriptives_final,
    summarize_bad_channels(file_handler_final)
  ),
  tar_target(
    bad_channels_plots_final,
    plot_bad_channel_counts(bad_channels_descriptives_final)
  ),
  tar_target(
    bad_channels_figure,
    save_figure(
      "figures/eeg-diagnostics/bad-eeg-channels.png",
      bad_channels_plots_final,
      height = 9,
      scaling = 1
    ),
    format = "file"
  ),
  tar_target(
    bad_segments_descriptives_final,
    summarize_bad_segments(file_handler_final)
  ),
  tar_target(
    bad_segments_plots_final,
    plot_bad_segment_descriptives(bad_segments_descriptives_final)
  ),
  tar_target(
    bad_segments_figure,
    save_figure(
      "figures/eeg-diagnostics/bad-eeg-segments.png",
      bad_segments_plots_final,
      height = 18,
      scaling = 1
    ),
    format = "file"
  )
  # TODO: Add targets for bad ICA indices
)

manuscripts_targets <- list(
  # Child documents for different sections of the manuscript ----
  tar_target(
    frontmatter,
    here("manuscripts", "child-documents", "frontmatter.Rmd"),
    format = "file"
  ),
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
  tar_target(
    results,
    here("manuscripts", "child-documents", "results.Rmd"),
    format = "file"
  ),
  tar_target(
    discussion,
    here("manuscripts", "child-documents", "discussion.Rmd"),
    format = "file"
  ),
  tar_target(
    appendix,
    here("manuscripts", "child-documents", "appendix.Rmd"),
    format = "file"
  ),
  # Bibliography file ----
  tar_target(
    packages,
    c(
      "base",
      "targets",
      "renv",
      "reticulate",
      "FactoMineR",
      "glmmTMB",
      "emmeans",
      "DHARMa",
      "performance",
      "ggplot2",
      "ggdist",
      "ggh4x",
      "ggnewscale",
      "patchwork",
      "tidyverse",
      "rmarkdown",
      "papaja",
      "officer",
      "officedown",
      "flextable",
      "ftExtra",
      "gtsummary"
    )
  ),
  tar_target(
    package_references,
    write_packages_bib(
      here("manuscripts", "packages.bib"),
      packages = packages
    ),
    format = "file"
  ),
  tar_target(
    references,
    write_bib(
      here("manuscripts", "references.json"),
      keys = c(
        introduction,
        methods,
        results,
        discussion
      ),
      ignore = packages,
      public_library_id = "mccarthymg"
    ),
    format = "file"
  ),
  # Session info ----
  tar_target(
    session_information_environment,
    session_info_environment()
  ),
  tar_target(
    session_information_packages,
    session_info_packages(),
    cue = tar_cue("always")
  ),
  # Thesis (officedown) ----
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
      references_path = references,
      packages_path   = package_references
    )
  ),
  # Supplement ----
  tar_render(
    supplement_phase_coupling_functional_connectomes,
    here(
      "manuscripts", "supplement", "phase-coupling-functional-connectomes.Rmd"
    )
  ),
  tar_render(
    supplement_amplitude_coupling_functional_connectomes,
    here(
      "manuscripts", "supplement", "amplitude-coupling-functional-connectomes.Rmd"
    )
  ),
  tar_render(
    supplement_phase_coupling_functional_connectomes_hilbert,
    here(
      "manuscripts", "supplement", "phase-coupling-functional-connectomes-hilbert.Rmd"
    )
  ),
  tar_render(
    supplement_similarity_matrices,
    here(
      "manuscripts", "supplement", "similarity-matrices.Rmd"
    )
  ),
  tar_render(
    supplement_model_summaries,
    here(
      "manuscripts", "supplement", "model-summaries.Rmd"
    )
  ),
  tar_render(
    supplement_maximal_model_results,
    here(
      "manuscripts", "supplement", "maximal-model-results.Rmd"
    )
  ),
  tar_render(
    supplement_eeg_preprocessing_diagnostics,
    here(
      "manuscripts", "supplement", "eeg-preprocessing-diagnostics.Rmd"
    )
  )
)

list(
  preprocessing_targets,
  connectivity_estimation_targets,
  descriptives_targets,
  manuscripts_targets
)
