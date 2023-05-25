library(targets)
library(tarchetypes)

library(here)
library(tidyverse)

library(ggdist)
library(ggh4x)
library(ggnewscale)
library(patchwork)

source(here("R/connectivity.R"))
source(here("R/similarity.R"))

# Connectivity ----
targets::tar_load(phase_connectivity_matrix_alpha)
targets::tar_load(phase_connectivity_plot_alpha)

example_case <- phase_connectivity_matrix_alpha$phase_connectivity_matrix_alpha_b60434fe
which_recording <- example_case$metadata$case
which_participant <- example_case$metadata$participant

functional_connectome <- plot_connectivity(example_case, "PLI or AEC")$plot +
  ggplot2::labs(
    title = "Example Functional Connectome"
  ) +
  ggplot2::theme(
    legend.title = ggplot2::element_text(vjust = 3),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )

connectivity_profile <- get_connectivity_profiles(phase_connectivity_matrix_alpha) |>
  dplyr::filter(case == which_recording) |>
  ggplot2::ggplot(ggplot2::aes(x = pair, y = case, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "cividis") +
    ggplot2::labs(
      title = "Example Connectivity Profile",
      x = "EEG Channel Pairs",
      y = NULL,
      fill = "PLI or AEC"
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "none"
    )

connectivity_patchwork <- functional_connectome + connectivity_profile +
  patchwork::plot_layout(ncol = 1, heights = c(1, 0.1), guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A")

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/connectivity.png"),
  plot = connectivity_patchwork,
  device = ragg::agg_png,
  width = 26.4583,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 2
)

# Similarity ----
targets::tar_load(phase_similarity_alpha)

similarity_matrix <- phase_similarity_alpha |>
  dplyr::filter(
    x_participant == which_participant & y_participant == which_participant
  ) |>
  dplyr::mutate(
    x_participant = forcats::fct_recode(x_participant, P00 = which_participant),
    y_participant = forcats::fct_recode(y_participant, P00 = which_participant)
  ) |>
  plot_similarity(rv) +
    ggplot2::labs(
      title = "Example Similarity Matrix"
    ) +
    ggplot2::theme(
      axis.ticks.x.top = ggplot2::element_line(size = 7),
      axis.ticks.y.left = ggplot2::element_line(size = 7.6)
    )

connectivity_profiles <- plot_connectivity_profiles(
  phase_connectivity_matrix_alpha[1:12], "PLI or AEC"
  ) +
  ggplot2::labs(title = "Example Connectivity Profiles") +
  ggplot2::theme(
    axis.ticks.y.left = ggplot2::element_line(size = 0.9),
  )
connectivity_profiles$data$participant <- "P00"

similarity_patchwork <- connectivity_profiles + similarity_matrix +
  patchwork::plot_layout(ncol = 1, heights = c(0.1, 1), guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    legend.key.height = grid::unit(0.38, "cm"),
    legend.key.width  = grid::unit(0.4, "cm")
  )

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/similarity.png"),
  plot = similarity_patchwork,
  device = ragg::agg_png,
  width = 26.4583,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 2
)

# Archetypal outcomes ----
save_similarity_archetypes_figure <- function(filename, plots) {

  # Hide the legends so we can only show one for the entire patchwork
  group_effect <- plots$group_effect +
    ggplot2::ggtitle("Group effect") +
    ggplot2::theme(legend.position = "none")
  session_effect <- plots$session_effect +
    ggplot2::ggtitle("Session effect") +
    ggplot2::theme(legend.position = "none")
  state_effect <- plots$state_effect +
    ggplot2::ggtitle("State effect") +
    ggplot2::theme(legend.position = "none")
  individual_session_effect <- plots$individual_session_effect +
    ggplot2::ggtitle("Individual-session effect") +
    ggplot2::theme(legend.position = "none")
  individual_state_effect <- plots$individual_state_effect +
    ggplot2::ggtitle("Individual-state effect") +
    ggplot2::theme(legend.position = "none")

  # Modify the legend so that it's discrete instead of continuous. This legend
  # will be extracted to represent the entire patchwork, since collecting
  # legends does not work with these plots. Also add legends for session and
  # state.
  axis_tick_colours <- colorspace::diverging_hcl(
    12, h = c(299, 135), c = 60, l = c(20, 80), power = c(0.7, 1.3)
  )

  case_colours <- rep(
    c(rev(axis_tick_colours[7:12]), axis_tick_colours[1:6]),
    times = 14 # length(participant_levels)
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

  individual_effect <- plots$individual_effect +
    # Make discrete instead of continuous
    ggplot2::scale_fill_continuous(
      breaks = c(0, 1),
      labels = c("No similarity", "Perfect similarity"),
      type = "viridis",
      guide = ggplot2::guide_legend(
        title = "Functional\nConnectome\nSimilarity", reverse = TRUE, order = 1
      )
    ) +
    # These rectangles are purely here to make the legend show up; they aren't
    # actually visible in the plot.
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = Session),
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
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = State),
      data = state_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(case_colours[9], case_colours[3]),
      guide = ggplot2::guide_legend(order = 2)
    )

  plots_legend <- individual_effect |>
    ggpubr::get_legend() |>
    ggpubr::as_ggplot()

  individual_effect <- individual_effect +
    ggplot2::ggtitle("Individual effect") +
    ggplot2::theme(legend.position = "none")

  design <- "ABCGHI
             DEFGHI"

  p1 <- group_effect
  p2 <- session_effect
  p3 <- state_effect
  p4 <- individual_effect
  p5 <- individual_session_effect
  p6 <- individual_state_effect
  p7 <- plots_legend

  patch <- patchwork::wrap_plots(
    A = p1, B = p2, C = p3, D = p4, E = p5, F = p6,
    G = patchwork::plot_spacer(), H = p7, I = patchwork::plot_spacer(),
    design = design,
    guides = "collect",
    widths = c(1, 1, 1, 0.025, 0.2, 0.025)
  )

  patch <- patch &
    ggplot2::theme(
      rect = ggplot2::element_rect(colour = NA_character_, fill = "transparent"),
      plot.background = ggplot2::element_rect(colour = NA_character_, fill = "transparent"),
      panel.background = ggplot2::element_rect(colour = NA_character_, fill = "transparent")
    )

  ggplot2::ggsave(
    filename = filename,
    plot = patch,
    device = ragg::agg_png,
    width = 58.34944,
    height = 25.8233,
    units = "cm",
    dpi = "retina",
    scaling = 1.2,
    bg = "transparent"
  )

}

save_similarity_archetypes_figure(
  filename = here("manuscripts/presentation/figures/archetypes-transparent.png"),
  plots = tar_read(similarity_archetype_plot)
)

# Amplitude Contrasts ----
targets::tar_load(amplitude_similarity_contrasts_delta)
targets::tar_load(amplitude_similarity_contrasts_theta)
targets::tar_load(amplitude_similarity_contrasts_alpha)
targets::tar_load(amplitude_similarity_contrasts_beta)
targets::tar_load(amplitude_similarity_contrasts_gamma)

amplitude_contrasts_all_bands <- list(
  delta = amplitude_similarity_contrasts_delta,
  theta = amplitude_similarity_contrasts_theta,
  alpha = amplitude_similarity_contrasts_alpha,
  beta  = amplitude_similarity_contrasts_beta,
  gamma = amplitude_similarity_contrasts_gamma
)

amplitude_contrasts_all_bands_tidy <- map_dfr(
  amplitude_contrasts_all_bands,
  function(x) {
    x |>
      tidy_contrast_similarity() |>
      dplyr::filter(
        effect_label %in% c(
          "Main effect",
          "Between sessions and states",
          "Within sessions and states")
      )
  },
  .id = "frequency_band"
) |>
  dplyr::mutate(
    frequency_band = stringr::str_to_title(frequency_band),
    frequency_band = factor(
      frequency_band,
      levels = c("Delta", "Theta", "Alpha", "Beta", "Gamma")
    )
  )

plot_similarity_contrasts_all_bands <- function(contrasts) {

  # Plotting
  ggplot2::ggplot(contrasts, ggplot2::aes(x = estimate, y = frequency_band)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linetype = 2, alpha = 0.2) +
    ggdist::stat_interval(
      ggplot2::aes(
        xdist = distributional::dist_student_t(
          df = df, mu = estimate, sigma = std.error
        ),
        colour = effect_direction,
        colour_ramp = stat(level)
      ),
      .width = c(0.66, 0.95)
    ) +
    ggplot2::scale_colour_discrete(
      drop = FALSE,
      type = c(
        "#3aaf85", # Green
        "#1b6ca8", # Blue
        "#cd201f"  # Red
      )
    ) +
    ggdist::scale_colour_ramp_discrete(
      from = "#fafdce",
      range = c(0.33, 1)
    ) +
    # This is a slightly hacky way to get a bar the same height as the interval
    # for the point estimate since the standard ggplot2 geoms aren't well-suited
    # for this.
    ggdist::geom_interval(
      ggplot2::aes(
        # When the estimate is very small the width of the point estimate
        # interval should be smaller so it doesn't overplot the CIs.
        xmin = ifelse(.upper_0.95 - .lower_0.95 < 0.01, estimate - 0.00025, estimate - 0.00075),
        xmax = ifelse(.upper_0.95 - .lower_0.95 < 0.01, estimate + 0.00025, estimate + 0.00075)
      )
    ) +
    # Set the scale of the x-axis so it's easier to compare plots. These values
    # were selected based on the data, so would need to be changed with different
    # data.
    ggplot2::scale_x_continuous(
      breaks = c(-0.04, -0.02, 0, 0.02, 0.04),
      labels = function(x) ifelse(x == 0, 0, x)
    ) +
    # Adding padding makes it hard to centre the secondary x-axis, so don't
    # expand.
    ggplot2::coord_cartesian(xlim = c(-0.05, 0.05), expand = FALSE) +
    ggplot2::labs(
      x = paste0(
        "Difference in functional connectome similarity",
        "\n",
        "(Within participant \U2212 Between participant)"
      ),
      y = "Frequency Band",
      colour = "Effect Direction",
      colour_ramp = "Compatibility\nInterval"
    ) +
    ggplot2::facet_wrap(~ effect_label, ncol = 1) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.line.x.top  = ggplot2::element_blank(),
      axis.text.x.top  = ggplot2::element_text(size = 11, colour = "black")
    )

}

amplitude_contrasts <- plot_similarity_contrasts_all_bands(amplitude_contrasts_all_bands_tidy)

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/amplitude-similarity-contrasts.png"),
  plot = amplitude_contrasts,
  device = ragg::agg_png,
  width = 26.4583,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 2
)

# Amplitude Similarity Matrices ----
targets::tar_load(amplitude_similarity_plot_delta)
targets::tar_load(amplitude_similarity_plot_theta)
targets::tar_load(amplitude_similarity_plot_alpha)
targets::tar_load(amplitude_similarity_plot_beta)
targets::tar_load(amplitude_similarity_plot_gamma)

add_legend_graphics <- function(p) {

  ## The axis-tick colours are being used as annotations for the different
  ## sessions and states. This is a bit hackish, but works for our purposes.
  axis_tick_colours <- colorspace::diverging_hcl(
    12, h = c(299, 135), c = 60, l = c(20, 80), power = c(0.7, 1.3)
  )

  case_colours <- rep(
    c(rev(axis_tick_colours[7:12]), axis_tick_colours[1:6]),
    times = 14
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

  p +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = Session),
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
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = State),
      data = state_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(case_colours[9], case_colours[3]),
      guide = ggplot2::guide_legend(order = 2)
    )

}

amplitude_similarity_plot_delta <- amplitude_similarity_plot_delta |>
  add_legend_graphics() +
  ggplot2::ggtitle("Delta")

amplitude_similarity_plot_theta <- amplitude_similarity_plot_theta |>
  add_legend_graphics() +
  ggplot2::ggtitle("Theta")

amplitude_similarity_plot_alpha <- amplitude_similarity_plot_alpha |>
  add_legend_graphics() +
  ggplot2::ggtitle("Alpha")

amplitude_similarity_plot_beta <- amplitude_similarity_plot_beta |>
  add_legend_graphics() +
  ggplot2::ggtitle("Beta")

amplitude_similarity_plot_gamma <- amplitude_similarity_plot_gamma |>
  add_legend_graphics() +
  ggplot2::ggtitle("Gamma")

amplitude_similarity_plots_patch <- (amplitude_similarity_plot_delta +
  amplitude_similarity_plot_theta +
  amplitude_similarity_plot_alpha
  ) /
  (patchwork::plot_spacer() +
  amplitude_similarity_plot_beta +
  amplitude_similarity_plot_gamma +
  patchwork::plot_spacer() +
  patchwork::plot_layout(
    nrow = 1,
    widths = c(0.5, 1, 1, 0.5)
  )
  ) +
  patchwork::plot_layout(
    nrow = 2,
    guides = "collect"
  )

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/amplitude-similarity-matrices.png"),
  plot = amplitude_similarity_plots_patch,
  device = ragg::agg_png,
  width = 58.34944,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 1.2
)

# Phase similarity matrices ----
targets::tar_load(phase_similarity_plot_delta)
targets::tar_load(phase_similarity_plot_theta)
targets::tar_load(phase_similarity_plot_alpha)
targets::tar_load(phase_similarity_plot_beta)
targets::tar_load(phase_similarity_plot_gamma)

phase_similarity_plot_delta <- phase_similarity_plot_delta |>
  add_legend_graphics() +
  ggplot2::ggtitle("Delta")

phase_similarity_plot_theta <- phase_similarity_plot_theta |>
  add_legend_graphics() +
  ggplot2::ggtitle("Theta")

phase_similarity_plot_alpha <- phase_similarity_plot_alpha |>
  add_legend_graphics() +
  ggplot2::ggtitle("Alpha")

phase_similarity_plot_beta <- phase_similarity_plot_beta |>
  add_legend_graphics() +
  ggplot2::ggtitle("Beta")

phase_similarity_plot_gamma <- phase_similarity_plot_gamma |>
  add_legend_graphics() +
  ggplot2::ggtitle("Gamma")

phase_similarity_plots_patch <- (phase_similarity_plot_delta +
  phase_similarity_plot_theta +
  phase_similarity_plot_alpha
  ) /
  (patchwork::plot_spacer() +
     phase_similarity_plot_beta +
     phase_similarity_plot_gamma +
     patchwork::plot_spacer() +
     patchwork::plot_layout(
       nrow = 1,
       widths = c(0.5, 1, 1, 0.5)
     )
  ) +
  patchwork::plot_layout(
    nrow = 2,
    guides = "collect"
  )

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/phase-similarity-matrices.png"),
  plot = phase_similarity_plots_patch,
  device = ragg::agg_png,
  width = 58.34944,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 1.2
)

# Phase Contrasts ----
targets::tar_load(phase_similarity_contrasts_delta)
targets::tar_load(phase_similarity_contrasts_theta)
targets::tar_load(phase_similarity_contrasts_alpha)
targets::tar_load(phase_similarity_contrasts_beta)
targets::tar_load(phase_similarity_contrasts_gamma)

phase_contrasts_all_bands <- list(
  delta = phase_similarity_contrasts_delta,
  theta = phase_similarity_contrasts_theta,
  alpha = phase_similarity_contrasts_alpha,
  beta  = phase_similarity_contrasts_beta,
  gamma = phase_similarity_contrasts_gamma
)

phase_contrasts_all_bands_tidy <- map_dfr(
  phase_contrasts_all_bands,
  function(x) {
    x |>
      tidy_contrast_similarity() |>
      dplyr::filter(
        effect_label %in% c(
          "Main effect",
          "Between sessions and states",
          "Within sessions and states")
      )
  },
  .id = "frequency_band"
) |>
  dplyr::mutate(
    frequency_band = stringr::str_to_title(frequency_band),
    frequency_band = factor(
      frequency_band,
      levels = c("Delta", "Theta", "Alpha", "Beta", "Gamma")
    )
  )

plot_similarity_contrasts_all_bands <- function(contrasts) {

  # Plotting
  ggplot2::ggplot(contrasts, ggplot2::aes(x = estimate, y = frequency_band)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linetype = 2, alpha = 0.2) +
    ggdist::stat_interval(
      ggplot2::aes(
        xdist = distributional::dist_student_t(
          df = df, mu = estimate, sigma = std.error
        ),
        colour = effect_direction,
        colour_ramp = stat(level)
      ),
      .width = c(0.66, 0.95)
    ) +
    ggplot2::scale_colour_discrete(
      drop = FALSE,
      type = c(
        "#3aaf85", # Green
        "#1b6ca8", # Blue
        "#cd201f"  # Red
      )
    ) +
    ggdist::scale_colour_ramp_discrete(
      from = "#fafdce",
      range = c(0.33, 1)
    ) +
    # This is a slightly hacky way to get a bar the same height as the interval
    # for the point estimate since the standard ggplot2 geoms aren't well-suited
    # for this.
    ggdist::geom_interval(
      ggplot2::aes(
        # When the estimate is very small the width of the point estimate
        # interval should be smaller so it doesn't overplot the CIs.
        xmin = ifelse(.upper_0.95 - .lower_0.95 < 0.01, estimate - 0.00025, estimate - 0.00075),
        xmax = ifelse(.upper_0.95 - .lower_0.95 < 0.01, estimate + 0.00025, estimate + 0.00075)
      )
    ) +
    # Set the scale of the x-axis so it's easier to compare plots. These values
    # were selected based on the data, so would need to be changed with different
    # data.
    ggplot2::scale_x_continuous(
      breaks = c(-0.1, -0.05, 0, 0.05, 0.1),
      labels = function(x) ifelse(x == 0, 0, x)
    ) +
    # Adding padding makes it hard to centre the secondary x-axis, so don't
    # expand.
    ggplot2::coord_cartesian(xlim = c(-0.15, 0.15), expand = FALSE) +
    ggplot2::labs(
      x = paste0(
        "Difference in functional connectome similarity",
        "\n",
        "(Within participant \U2212 Between participant)"
      ),
      y = "Frequency Band",
      colour = "Effect Direction",
      colour_ramp = "Compatibility\nInterval"
    ) +
    ggplot2::facet_wrap(~ effect_label, ncol = 1) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.line.x.top  = ggplot2::element_blank(),
      axis.text.x.top  = ggplot2::element_text(size = 11, colour = "black")
    )

}

phase_contrasts <- plot_similarity_contrasts_all_bands(phase_contrasts_all_bands_tidy)

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/phase-similarity-contrasts.png"),
  plot = phase_contrasts,
  device = ragg::agg_png,
  width = 26.4583,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 2
)

# Phase alpha individual contrasts ----
targets::tar_load(phase_similarity_contrasts_alpha)
targets::tar_load(phase_similarity_subset_contrasts_alpha)

plot_subset_similarity_contrasts <- function(group_object, subset_objects) {

  # Note: The gradient fill in the ridges is currently commented out, as it's a
  # bit busy when also filling the ridges by effect direction, with the colour
  # of a ridge becoming less clear when the effect sizes are very small.

  # Tidy data to prepare for plotting
  contrasts <- tidy_contrast_similarity(group_object) |>
    dplyr::filter(
      effect_label %in% c(
        "Between sessions and states",
        "Within sessions and states")
    ) |>
    dplyr::mutate(effect_label = factor(effect_label))
  subset_contrasts <- tidy_subset_contrast_similarity(subset_objects) |>
    dplyr::filter(
      effect_label %in% c(
        "Between sessions and states",
        "Within sessions and states")
    ) |>
    dplyr::mutate(effect_label = factor(effect_label))

  subset_contrasts <- contrasts |>
    dplyr::select(effect_label, group_estimate = estimate) |>
    dplyr::left_join(subset_contrasts)

  # Plotting
  my_scale <- ggplot2::scale_x_continuous(
    breaks = scales::extended_breaks(n = 3)
  )
  my_scale$train(
    c(min(subset_contrasts$.lower_0.95), max(subset_contrasts$.upper_0.95))
  )
  my_breaks <- my_scale$get_breaks()
  #my_breaks <- ifelse(0 %in% my_breaks, my_breaks, sort(append(my_breaks, 0)))

  p <- ggplot2::ggplot(subset_contrasts) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linetype = 5, alpha = 0.2) +
    # Plot group-level intervals so there's a reference point to compare the
    # sub-models to.
    ggdist::stat_interval(
      ggplot2::aes(
        y = "Group",
        xdist = distributional::dist_student_t(
          df = df, mu = estimate, sigma = std.error
        ),
        colour = effect_direction,
        colour_ramp = stat(level)
      ),
      .width = c(0.66, 0.95),
      data = contrasts
    ) +
    ggplot2::scale_colour_discrete(
      drop = FALSE,
      type = c(
        "#3aaf85", # Green
        "#1b6ca8", # Blue
        "#cd201f"  # Red
      )
    ) +
    ggdist::scale_colour_ramp_discrete(
      from = "#fafdce",
      range = c(0.33, 1)
    ) +
    ggdist::geom_interval(
      ggplot2::aes(
        y = "Group", ymin = "Group", ymax = "Group",
        # When the estimate is very small the width of the point estimate
        # interval should be smaller so it doesn't overplot the CIs.
        xmin = ifelse(.upper_0.95 - .lower_0.95 < 0.01, estimate - 0.00025, estimate - 0.00075),
        xmax = ifelse(.upper_0.95 - .lower_0.95 < 0.01, estimate + 0.00025, estimate + 0.00075)
      ),
      orientation = "horizontal",
      data = contrasts
    ) +
    # Confidence distributions for participant sub-models
    ggdist::stat_slab(
      ggplot2::aes(
        x = estimate,
        y = participant,
        # This throws a warning that the aesthetic is being ignored, but it's
        # needed for the code in `fill_ramp` to work.
        group_estimate = group_estimate,
        xdist = distributional::dist_student_t(df = df, mu = estimate, sigma = std.error),
        fill = effect_direction
        # fill_ramp = stat(abs(group_estimate - x))
      ),
      color = "black",
      size = 0.25,
      height = 2,
      expand = TRUE, # I can't disable this for whatever reason
      trim = FALSE, # I can't disable this for whatever reason
      # fill_type = "gradient"
    ) +
    ggplot2::scale_fill_discrete(
      type = c(
        "#3aaf85", # Green
        "#1b6ca8", # Blue
        "#cd201f"  # Red
      )
    ) +
  # General plot settings
  ggplot2::scale_x_continuous(
    breaks = my_breaks,
    labels = function(x) ifelse(x == 0, 0, x)
  ) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(1.1, 0))) +
    ## These legends already appear in the group-level interval plot, so they
    ## don't need to be repeated here.
    ggplot2::guides(
      fill = "none",
      colour_ramp = "none"
    ) +
    # Keep the empty " " facet to enforce facet alignment. It will be removed
    # manually after building the plot. Make the direction vertical so that
    # the "Main effect" facet is above the " " facet.
    ggplot2::facet_wrap(~ effect_label, ncol = 1) +
    ggplot2::labs(
      title = "Alpha Band Individual Differences",
      x = paste0(
        "Difference in functional connectome similarity",
        "\n",
        "(Within participant y \U2212 Between participant)"
      ),
      y = "Participant",
      colour = "Effect Direction"
    ) +
    ggplot2::theme_classic()

  p

}

phase_similarity_subset_contrasts_plot_alpha <-
  plot_subset_similarity_contrasts(
    phase_similarity_contrasts_alpha,
    phase_similarity_subset_contrasts_alpha
  )

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/phase-similarity-subset-contrasts-alpha-2.png"),
  plot = phase_similarity_subset_contrasts_plot_alpha,
  device = ragg::agg_png,
  width = 26.4583,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 2
)

# Phase alpha patchwork ----
targets::tar_load(phase_connectivity_profile_plot_alpha)
targets::tar_load(phase_similarity_plot_alpha)

phase_alpha_patch <- phase_connectivity_profile_plot_alpha +
  phase_similarity_plot_alpha +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A")

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/phase-alpha.png"),
  plot = phase_alpha_patch,
  device = ragg::agg_png,
  width = 58.34944,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 2
)

# Phase connectivity profiles ----
targets::tar_load(phase_connectivity_profile_plot_delta)
targets::tar_load(phase_connectivity_profile_plot_theta)
targets::tar_load(phase_connectivity_profile_plot_alpha)
targets::tar_load(phase_connectivity_profile_plot_beta)
targets::tar_load(phase_connectivity_profile_plot_gamma)

phase_connectivity_profile_plot_delta <- phase_connectivity_profile_plot_delta +
  ggplot2::ggtitle("Delta")

phase_connectivity_profile_plot_theta <- phase_connectivity_profile_plot_theta +
  ggplot2::ggtitle("Theta")

phase_connectivity_profile_plot_alpha <- phase_connectivity_profile_plot_alpha +
  ggplot2::ggtitle("Alpha")

phase_connectivity_profile_plot_beta <- phase_connectivity_profile_plot_beta +
  ggplot2::ggtitle("Beta")

phase_connectivity_profile_plot_gamma <- phase_connectivity_profile_plot_gamma +
  ggplot2::ggtitle("Gamma")

phase_connectivity_profile_plots_patch <- (
  phase_connectivity_profile_plot_delta +
  phase_connectivity_profile_plot_theta +
  phase_connectivity_profile_plot_alpha
) /
  (patchwork::plot_spacer() +
     phase_connectivity_profile_plot_beta +
     phase_connectivity_profile_plot_gamma +
     patchwork::plot_spacer() +
     patchwork::plot_layout(
       nrow = 1,
       widths = c(0.5, 1, 1, 0.5)
     )
  ) +
  patchwork::plot_layout(
    nrow = 2,
    guides = "collect"
  )

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/phase-connectivity-profiles.png"),
  plot = phase_connectivity_profile_plots_patch,
  device = ragg::agg_png,
  width = 58.34944,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 1.2
)

# Amplitude connectivity profiles ----
targets::tar_load(amplitude_connectivity_profile_plot_delta)
targets::tar_load(amplitude_connectivity_profile_plot_theta)
targets::tar_load(amplitude_connectivity_profile_plot_alpha)
targets::tar_load(amplitude_connectivity_profile_plot_beta)
targets::tar_load(amplitude_connectivity_profile_plot_gamma)

amplitude_connectivity_profile_plot_delta <- amplitude_connectivity_profile_plot_delta +
  ggplot2::ggtitle("Delta")

amplitude_connectivity_profile_plot_theta <- amplitude_connectivity_profile_plot_theta +
  ggplot2::ggtitle("Theta")

amplitude_connectivity_profile_plot_alpha <- amplitude_connectivity_profile_plot_alpha +
  ggplot2::ggtitle("Alpha")

amplitude_connectivity_profile_plot_beta <- amplitude_connectivity_profile_plot_beta +
  ggplot2::ggtitle("Beta")

amplitude_connectivity_profile_plot_gamma <- amplitude_connectivity_profile_plot_gamma +
  ggplot2::ggtitle("Gamma")

amplitude_connectivity_profile_plots_patch <- (
  amplitude_connectivity_profile_plot_delta +
    amplitude_connectivity_profile_plot_theta +
    amplitude_connectivity_profile_plot_alpha
) /
  (patchwork::plot_spacer() +
     amplitude_connectivity_profile_plot_beta +
     amplitude_connectivity_profile_plot_gamma +
     patchwork::plot_spacer() +
     patchwork::plot_layout(
       nrow = 1,
       widths = c(0.5, 1, 1, 0.5)
     )
  ) +
  patchwork::plot_layout(
    nrow = 2,
    guides = "collect"
  )

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/amplitude-connectivity-profiles.png"),
  plot = amplitude_connectivity_profile_plots_patch,
  device = ragg::agg_png,
  width = 58.34944,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 1.2
)

# Amplitude connectivity plots (free fill limits) ----

amplitude_connectivity_profile_plot_delta_free <- amplitude_connectivity_profile_plot_delta +
  ggplot2::ggtitle("Delta") +
  scale_fill_viridis_c(
    option = "cividis",
    limits = NULL,
    guide = guide_colorbar(order = 1)
  )

amplitude_connectivity_profile_plot_theta_free <- amplitude_connectivity_profile_plot_theta +
  ggplot2::ggtitle("Theta") +
  scale_fill_viridis_c(
    option = "cividis",
    limits = NULL,
    guide = guide_colorbar(order = 1)
  )

amplitude_connectivity_profile_plot_alpha_free <- amplitude_connectivity_profile_plot_alpha +
  ggplot2::ggtitle("Alpha") +
  scale_fill_viridis_c(
    option = "cividis",
    limits = NULL,
    guide = guide_colorbar(order = 1)
  )

amplitude_connectivity_profile_plot_beta_free <- amplitude_connectivity_profile_plot_beta +
  ggplot2::ggtitle("Beta") +
  scale_fill_viridis_c(
    option = "cividis",
    limits = NULL,
    guide = guide_colorbar(order = 1)
  )

amplitude_connectivity_profile_plot_gamma_free <- amplitude_connectivity_profile_plot_gamma +
  ggplot2::ggtitle("Gamma") +
  scale_fill_viridis_c(
    option = "cividis",
    limits = NULL,
    guide = guide_colorbar(order = 1)
  )

amplitude_connectivity_profile_plots_free_patch <- (
  amplitude_connectivity_profile_plot_delta_free +
    amplitude_connectivity_profile_plot_theta_free +
    amplitude_connectivity_profile_plot_alpha_free
) /
  (patchwork::plot_spacer() +
     amplitude_connectivity_profile_plot_beta_free +
     amplitude_connectivity_profile_plot_gamma_free +
     patchwork::plot_spacer() +
     patchwork::plot_layout(
       nrow = 1,
       widths = c(0.5, 1, 1, 0.5)
     )
  ) +
  patchwork::plot_layout(
    nrow = 2,
    guides = "auto"
  ) &
  ggplot2::theme(
    legend.key.height = grid::unit(0.38, "cm"),
    legend.key.width  = grid::unit(0.4, "cm")
  )

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/amplitude-connectivity-profiles-free.png"),
  plot = amplitude_connectivity_profile_plots_free_patch,
  device = ragg::agg_png,
  width = 58.34944,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 1.2
)

# Gratton et al. figures ----
gratton_magnitude_plot <- tibble(
  effect = factor(
    c("Other", "Within Individuals", "Between Individuals"),
    levels = c("Other", "Within Individuals", "Between Individuals")
  ),
  value  = c(0.35, 0.52, 0.56)
) |>
  dplyr::mutate(
    value = value / sum(value)
  ) |>
  ggplot(aes(x = 0, y = value, fill = effect)) +
  geom_col(width = 0.1) +
  geom_text(
    aes(label = effect, x = 0.075, hjust = 0),
    position = position_stack(vjust = 0.5)
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  scale_fill_manual(
    values = c("grey", "#3aaf85","#1b6ca8"),
    guide = guide_legend(reverse = TRUE)
  ) +
  coord_cartesian(xlim = c(-.05, .5)) +
  guides(
    y = "axis_truncated"
  ) +
  labs(
    y = "Proportion of Total Network Similarity Effects",
    fill = "Effect",
    caption = "Adapted from Gratton et al., 2018"
  ) +
  theme_classic() +
  theme(
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x  = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_blank(),
    legend.position = "none"
  )

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/gratton-magnitude.png"),
  plot = gratton_magnitude_plot,
  device = ragg::agg_png,
  width = 15.3458,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 2.35
)
