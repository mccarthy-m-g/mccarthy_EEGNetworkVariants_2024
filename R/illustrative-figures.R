library(targets)
library(tidyverse)
library(patchwork)
library(reticulate)
use_miniconda("r-reticulate-mne")
np <- import("numpy")
sp <- import("scipy.signal")

# oscillatory-bands-plot ------------------------------------------------------

# Values come from: @penttonenbuzsaki_NaturalLogarithmicRelationship_2003
oscillatory_bands <- tribble(
  ~ name,       ~ lower,  ~ upper,  ~ centre,
  "slow 4",     1/36.161, 1/13.303, 1/21.933,
  "slow 3",     1/13.303, 1/4.894,  1/8.069,
  "slow 2",     1/4.894,  1/1.800,  1/2.968,
  "slow 1",     1/1.800,  1/0.662,  1/1.092,
  "delta",      1/0.662,  1/0.244,  1/0.402,
  "theta",      1/0.244,  1/0.090,  1/0.148,
  "alpha/beta", 1/0.090,  1/0.033,  1/0.054,
  "gamma",      30,       82.44,    50,
  "fast",       82.44,    224.08,   135.91,
  "ultra fast", 224.08,   609.12,   369.45
) |>
  mutate(
    name = factor(
      name,
      levels = c(
        "slow 4", "slow 3", "slow 2", "slow 1",
        "delta", "theta", "alpha/beta", "gamma", "fast", "ultra fast"
      )
    )
  )

bands_plot <- ggplot(oscillatory_bands, aes(x = centre, y = name)) +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  scale_x_continuous(
    trans = "log",
    n.breaks = 7,
    labels = function(x) round(x, 2),
    sec.axis = sec_axis(
      trans = log,
      name = "log Frequency (Hz)",
      breaks = -4:6
    )
  ) +
  labs(
    x = "Frequency (Hz)",
    y = "Oscillatory Band"
  ) +
  theme_classic()

ggplot2::ggsave(
  filename = "figures/oscillatory-bands-plot.png",
  plot = bands_plot,
  device = ragg::agg_png,
  width = 21.59,
  height = 12,
  units = "cm",
  dpi = "retina",
  scaling = 1
)

# oscillation-plot ------------------------------------------------------------

#' @param time in seconds.
#' @param frequency in Hertz.
#' @param phase in degrees.
#' @param amplitude in arbitrary units.
sine_wave <- function(time, frequency, phase, amplitude = 1) {
  amplitude * sin(2 * pi * frequency * time + (phase / (180 / pi)) )
}

p1 <- ggplot() +
  geom_function(
    aes(linetype = "1 Hz"),
    fun = sine_wave,
    args = list(frequency = 1, phase = 0, amplitude = 1),
    n = 10^2
  ) +
  geom_function(
    aes(linetype = "2 Hz"),
    fun = sine_wave,
    args = list(frequency = 2, phase = 0, amplitude = 1),
    n = 10^2
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s"),
    sec.axis = sec_axis(
      trans = function(.x) .x * 360,
      name = "Phase",
      breaks = seq(0, 360, by = 45),
      labels = function(.x) paste0(.x, "\u00BA")
    )) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    linetype = "Frequency"
  ) +
  theme_classic()

p2 <- ggplot() +
  geom_segment(
    aes(
      x = phase, xend = phase, y = 0, yend = amplitude,
      linetype = "1 Hz", colour = time, alpha = .8
    ),
    data = tibble(
      time = seq(0, 1, by = 0.125/2),
      phase = time * 360,
      amplitude = abs(sine_wave(time, frequency = 1, phase = 0, amplitude = 1))
    )
  ) +
  geom_segment(
    aes(
      x = phase, xend = phase, y = 0, yend = amplitude,
      linetype = "2 Hz", colour = time
    ),
    data = tibble(
      time = seq(0, 1, by = 0.125/2),
      phase = time * 360,
      amplitude = abs(sine_wave(time, frequency = 2, phase = 0, amplitude = 1))
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 360, by = 45),
    expand = c(0, 0),
    limits = c(0, 360),
    labels = function(.x) paste0(.x, "\u00BA")
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .25)) +
  scale_color_viridis_c(labels = function(.x) paste0(.x, "s")) +
  guides(linetype = "none", alpha = "none") +
  coord_polar(start = 270 / (180 / pi), direction = -1, clip = "off") +
  labs(
    x = "Phase",
    y = "Absolute Amplitude (AU)",
    linetype = "Frequency",
    colour = "Time"
  ) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

patch <- p1 + p2 +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

# patch <- p1 +
#   inset_element(
#     p2 & guides(linetype = "none"),
#     left = 0.5, bottom = 0.45, right = 1.2, top = .99
#   ) +
#   plot_layout(
#     guides = "collect"
#   )

ggplot2::ggsave(
  filename = "figures/oscillation-plot.png",
  plot = patch,
  device = ragg::agg_png,
  width = 21.59,
  height = 12,
  units = "cm",
  dpi = "retina",
  scaling = 1
)

# phase-plot ------------------------------------------------------------------

set.seed(666)

# 500 Hz
time <- seq(0, 1, length.out = 500)

# Sine waves + Gaussian white noise
wave_1 <- sine_wave(time, frequency = 5, phase = 145, amplitude = 1) +
  rnorm(time, mean = 0, sd = .25)
wave_2 <- sine_wave(time, frequency = 5, phase = 25, amplitude = 1) +
  rnorm(time, mean = 0, sd = .25)

# Gaussian white noise
wave_3 <- rnorm(time, mean = 0, sd = .49)

# Panel A and B data
panel_A <- tibble(
  time = c(time, time),
  amplitude = c(wave_1, wave_2),
  signal = rep(c("Signal 1", "Signal 2"), each = 500)
)

panel_B <- tibble(
  time = c(time, time),
  amplitude = c(wave_1, wave_3),
  signal = rep(c("Signal 1", "Signal 3"), each = 500)
)

# Time series plots
p3 <- ggplot(panel_A, aes(x = time, y = amplitude)) +
  geom_path() +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s")
  ) +
  facet_wrap(vars(signal), ncol = 1) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    tag = "A"
  ) +
  theme_classic()

p4 <- ggplot(panel_B, aes(x = time, y = amplitude)) +
  geom_path() +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s")
  ) +
  facet_wrap(vars(signal), ncol = 1) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    tag = "B"
  ) +
  theme_classic()

# Phase difference histograms
phase_diff_A <- panel_A |>
  group_by(signal) |>
  mutate(
    hilbert = sp$hilbert(amplitude),
    phase = np$angle(hilbert)
  ) |>
  select(signal, time, phase) |>
  pivot_wider(names_from = signal, values_from = phase) |>
  mutate(
    difference = as.numeric((`Signal 1` - `Signal 2`) * (180 / pi)),
    difference = if_else(difference < 0, difference + 360, difference)
  )

p5 <- ggplot(phase_diff_A, aes(x = difference)) +
  geom_histogram(binwidth = 8, colour = "black", fill = "white") +
  scale_x_continuous(
    breaks = seq(0, 360, by = 45),
    limits = c(0, 360),
    oob = scales::oob_keep,
    labels = function(.x) paste0(.x, "\u00BA")
  ) +
  coord_polar(start = 270 / (180 / pi), direction = -1, clip = "off") +
  labs(x = "Phase Angle Difference", y = ""
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

phase_diff_B <- panel_B |>
  group_by(signal) |>
  mutate(
    hilbert = sp$hilbert(amplitude),
    phase = np$angle(hilbert)
  ) |>
  select(signal, time, phase) |>
  pivot_wider(names_from = signal, values_from = phase) |>
  mutate(
    difference = as.numeric((`Signal 1` - `Signal 3`) * (180 / pi)),
    difference = if_else(difference < 0, difference + 360, difference)
  )

p6 <- ggplot(phase_diff_B, aes(x = difference)) +
  geom_histogram(binwidth = 8, colour = "black", fill = "white") +
  scale_x_continuous(
    breaks = seq(0, 360, by = 45),
    limits = c(0, 360),
    oob = scales::oob_keep,
    labels = function(.x) paste0(.x, "\u00BA")
  ) +
  coord_polar(start = 270 / (180 / pi), direction = -1, clip = "off") +
  labs(x = "Phase Angle Difference", y = ""
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Patchwork
layout <- "
A#B
C#D
"

patch_2 <- p3 + p4 + p5 + p6 + plot_layout(widths = c(1, .1, 1), design = layout)

patch_2

ggplot2::ggsave(
  filename = "figures/phase-plot.png",
  plot = patch_2,
  device = ragg::agg_png,
  width = 21.59,
  height = 18,
  units = "cm",
  dpi = "retina",
  scaling = 1
)

# amplitude-plot --------------------------------------------------------------

time <- seq(0, 1, length.out = 500)
phase <- 90

# Panel A
carrier_wave <- sine_wave(time, frequency = 50, phase = phase, amplitude = 1)
mod_wave <- sine_wave(time, frequency = 3, phase = phase, amplitude = .5)
am_wave <- (1 + mod_wave) * carrier_wave

amplitude_modulation <- tibble(
  time = c(time, time, time),
  amplitude = c(mod_wave, carrier_wave, am_wave),
  signal = rep(
    c("Modulating Signal", "Carrier Signal", "Amplitude Modulated Signal"),
    each = 500
  )
) |>
  mutate(
    signal = factor(signal, levels = c("Modulating Signal", "Carrier Signal", "Amplitude Modulated Signal"))
  )

amplitude_envelope <- amplitude_modulation |>
  filter(signal == "Amplitude Modulated Signal") |>
  mutate(
    hilbert = sp$hilbert(amplitude),
    envelope = abs(hilbert)
  )

amplitude_plot <- ggplot(amplitude_modulation, aes(x = time, y = amplitude)) +
  geom_path(aes(colour = "Real")) +
  geom_path(
    aes(y = envelope, colour = "Envelope"),
    linewidth = 1,
    data = amplitude_envelope
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s")
  ) +
  scale_colour_manual(
    values = c("black", "#cd201f"), breaks = c("Real", "Envelope")
  ) +
  facet_wrap(vars(signal), ncol = 1) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    colour = "Signal"
  ) +
  theme_classic()

# Panel B
aec_sim <- map_df(
  1:1000,
  function(x) {
    # Making waves
    mod_wave <- sine_wave(time, frequency = 3, phase = phase, amplitude = .5) +
      rnorm(length(time), 0, .3)

    carrier_wave_3 <- sine_wave(time, frequency = 50, phase = phase, amplitude = 1) +
      rnorm(length(time), 0, .2)
    carrier_wave_4 <- sine_wave(time, frequency = 30, phase = phase, amplitude = 1) +
      rnorm(length(time), 0, .2)

    am_wave_3 <- (1 + mod_wave) * carrier_wave_3
    am_wave_4 <- (1 + mod_wave) * carrier_wave_4

    aec_tbl_2 <- tibble(
      time = c(time, time),
      amplitude = c(am_wave_3, am_wave_4),
      signal = rep(c("Signal 1", "Signal 2"), each = 500)
    ) |>
      mutate(
        hilbert = sp$hilbert(amplitude),
        envelope = abs(hilbert)
      )

    # Correlations
    amplitude_corr <- cor(
      x = aec_tbl_2 |> filter(signal == "Signal 1") |> pull(amplitude),
      y = aec_tbl_2 |> filter(signal == "Signal 2") |> pull(amplitude)
    )

    envelope_corr <- cor(
      x = aec_tbl_2 |> filter(signal == "Signal 1") |> pull(envelope),
      y = aec_tbl_2 |> filter(signal == "Signal 2") |> pull(envelope)
    )

    # Results
    tibble(
      signal = c("Real", "Envelope"),
      r = c(amplitude_corr, envelope_corr)
    )
  }
)

sim_plot <- ggplot(aec_sim, aes(x = abs(r), fill = signal)) +
  geom_histogram(binwidth = .05) +
  scale_fill_manual(
    values = c("black", "#cd201f"), breaks = c("Real", "Envelope"),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "Absolute Correlation",
    y = "Count",
    fill = "Signal"
  ) +
  theme_classic()

amplitude_patch <- amplitude_plot + sim_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 1, guides = "collect")

ggplot2::ggsave(
  filename = "figures/amplitude-plot.png",
  plot = amplitude_patch,
  device = ragg::agg_png,
  width = 21.59,
  height = 18,
  units = "cm",
  dpi = "retina",
  scaling = 1
)

# contrasts-plot --------------------------------------------------------------

contrast_outcomes <- tibble(
  effect = rep(
    c(
      "Group effect",
      "Session effect",
      "State effect",
      "Individual effect",
      "Individual-session effect",
      "Individual-state effect"
    ),
    each = 9
  ),
  contrast = rep(
    c(
      "Main effect",
      "Between sessions",
      "Within sessions",
      "Between states",
      "Within states",
      "Between sessions and states",
      "Between sessions and within states",
      "Within sessions and between states",
      "Within sessions and states"
    ),
    times = 6
  ),
  value = c(
    # Group effect
    rep(0, times = 9),
    # Session effect
    rep(0, times = 9),
    # State effect
    rep(0, times = 9),
    # Individual effect
    rep(1, times = 9),
    # Individual-session effect
    c(0.5, 0, 1, 0.5, 0.5, 0, 0, 1, 1),
    # Individual-state effect
    c(0.5, 0.5, 0.5, 0, 1, 0, 1, 0, 1)
  )
) |>
  mutate(
    effect = factor(
      effect,
      levels = c(
        "Group effect",
        "Session effect",
        "State effect",
        "Individual effect",
        "Individual-session effect",
        "Individual-state effect"
      )
    ),
    contrast = factor(
      contrast,
      levels = c(
        "Main effect",
        "Between sessions",
        "Within sessions",
        "Between states",
        "Within states",
        "Between sessions and within states",
        "Within sessions and between states",
        "Between sessions and states",
        "Within sessions and states"
      )
    )
  )

contrasts_plot <- ggplot(contrast_outcomes, aes(x = value, y = contrast)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.2) +
  geom_point() +
  coord_cartesian(xlim = c(-1, 1)) +
  facet_wrap(vars(effect), ncol = 2, dir = "v") +
  scale_x_continuous(
    sec.axis = sec_axis(
      ~ .,
      breaks = c(-0.5, 0.5),
      labels = c(
        "Less similar\nwithin participant",
        "More similar\nwithin participant"
      )
    )
  ) +
  labs(
    x = paste0(
      "Difference in functional connectome similarity",
      "\n",
      "(Within participant \U2212 Between participant)"
    ),
    y = "Contrast"
  ) +
  theme_classic() +
  theme(
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank(),
  )

ggplot2::ggsave(
  filename = "figures/contrasts-plot.png",
  plot = contrasts_plot,
  device = ragg::agg_png,
  width = 21.59,
  height = 18,
  units = "cm",
  dpi = "retina",
  scaling = 1
)

# outcome-plot ----------------------------------------------------------------

outcome_plot <- contrast_outcomes |>
  filter(stringr::str_starts(effect, "Individual")) |>
  mutate(
    weight_individual = case_when(
      effect == "Individual effect" ~ 2/3,
      effect == "Individual-session effect" ~ 1/6,
      effect == "Individual-state effect" ~ 1/6
    ),
    weight_session = case_when(
      effect == "Individual effect" ~ 1/3,
      effect == "Individual-session effect" ~ 3/6,
      effect == "Individual-state effect" ~ 1/6
    ),
    weight_state = case_when(
      effect == "Individual effect" ~ 1/3,
      effect == "Individual-session effect" ~ 1/6,
      effect == "Individual-state effect" ~ 3/6
    )
  ) |>
  group_by(contrast) |>
  summarise(
    value_unweighted = mean(value),
    value_individual = weighted.mean(value, weight_individual),
    value_session = weighted.mean(value, weight_session),
    value_state = weighted.mean(value, weight_state)
  ) |>
  pivot_longer(
    cols = c(value_unweighted, value_individual, value_session, value_state)
  ) |>
  mutate(
    name = case_when(
      name == "value_unweighted" ~ "Equally weighted effects",
      name == "value_individual" ~ "Weighted towards Individual effect",
      name == "value_session" ~ "Weighted towards Individual-session effect",
      name == "value_state" ~ "Weighted towards Individual-state effect"
    )
  ) |>
  ggplot(aes(x = value, y = contrast)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.2) +
  geom_point() +
  coord_cartesian(xlim = c(-1, 1)) +
  facet_wrap(vars(name), ncol = 1) +
  labs(
    x = paste0(
      "Difference in functional connectome similarity",
      "\n",
      "(Within participant \U2212 Between participant)"
    ),
    y = "Contrast",
    colour = "Effect Direction",
    colour_ramp = "Compatibility\nInterval"
  ) +
  theme_classic() +
  theme(
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank(),
  )

ggplot2::ggsave(
  filename = "figures/outcome-plot.png",
  plot = outcome_plot,
  device = ragg::agg_png,
  width = 21.59/1.5,
  height = 18,
  units = "cm",
  dpi = "retina",
  scaling = 1
)

# TODO: Delete code after this point after deciding whether or not the above
# figures are good.

# amplitude-plot --------------------------------------------------------------

time <- seq(0, 1, length.out = 500)
phase <- 90

# Panel A
carrier_wave <- sine_wave(time, frequency = 50, phase = phase, amplitude = 1)
mod_wave <- sine_wave(time, frequency = 3, phase = phase, amplitude = .5)
am_wave <- (1 + mod_wave) * carrier_wave

amplitude_modulation <- tibble(
  time = c(time, time, time),
  amplitude = c(mod_wave, carrier_wave, am_wave),
  signal = rep(
    c("Modulating Signal", "Carrier Signal", "Amplitude Modulated Signal"),
    each = 500
  )
) |>
  mutate(
    signal = factor(signal, levels = c("Modulating Signal", "Carrier Signal", "Amplitude Modulated Signal"))
  )

amplitude_envelope <- amplitude_modulation |>
  filter(signal == "Amplitude Modulated Signal") |>
  mutate(
    hilbert = sp$hilbert(amplitude),
    envelope = abs(hilbert)
  )

amplitude_plot <- ggplot(amplitude_modulation, aes(x = time, y = amplitude)) +
  geom_path(aes(colour = "Real")) +
  geom_path(
    aes(y = envelope, colour = "Envelope"),
    linewidth = 1,
    data = amplitude_envelope
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s")
  ) +
  scale_colour_manual(
    values = c("black", "#cd201f"), breaks = c("Real", "Envelope")
  ) +
  facet_wrap(vars(signal), ncol = 1) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    colour = "Signal"
  ) +
  theme_classic()

# Panel B
carrier_wave_2 <- sine_wave(time, frequency = 30, phase = phase, amplitude = 1)
mod_wave_2 <- sine_wave(time, frequency = 3, phase = phase, amplitude = .6)
am_wave_2 <- (1 + mod_wave_2) * carrier_wave_2

aec_tbl <- tibble(
  time = c(time, time),
  amplitude = c(am_wave, am_wave_2),
  signal = rep(c("Signal 1", "Signal 2"), each = 500)
) |>
  mutate(
    hilbert = sp$hilbert(amplitude),
    envelope = abs(hilbert)
  )

aec_plot <- ggplot(aec_tbl, aes(x = time)) +
  geom_path(aes(y = amplitude, colour = "Real")) +
  geom_path(aes(y = envelope, colour = "Envelope"), linewidth = 1) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s")
  ) +
  scale_colour_manual(
    values = c("black", "#cd201f"), breaks = c("Real", "Envelope")
  ) +
  facet_wrap(vars(signal), ncol = 1) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    colour = "Signal"
  ) +
  theme_classic()

env_cor <- cor(
  x = aec_tbl |> filter(signal == "Signal 1") |> pull(envelope),
  y = aec_tbl |> filter(signal == "Signal 2") |> pull(envelope)
)

amp_cor <- cor(
  x = aec_tbl |> filter(signal == "Signal 1") |> pull(amplitude),
  y = aec_tbl |> filter(signal == "Signal 2") |> pull(amplitude)
)

# Patch
amplitude_patch <- amplitude_plot + aec_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

ggplot2::ggsave(
  filename = "figures/amplitude-plot.png",
  plot = amplitude_patch,
  device = ragg::agg_png,
  width = 21.59,
  height = 18,
  units = "cm",
  dpi = "retina",
  scaling = 1
)

# ----

mod_wave <- sine_wave(time, frequency = 3, phase = phase, amplitude = .5) +
  rnorm(length(time), 0, .1)

carrier_wave_3 <- sine_wave(time, frequency = 50, phase = phase, amplitude = 1) +
  rnorm(length(time), 0, .1)
carrier_wave_4 <- sine_wave(time, frequency = 30, phase = phase, amplitude = 1) +
  rnorm(length(time), 0, .1)

am_wave_3 <- (1 + mod_wave) * carrier_wave_3
am_wave_4 <- (1 + mod_wave) * carrier_wave_4

aec_tbl_2 <- tibble(
  time = c(time, time),
  amplitude = c(am_wave_3, am_wave_4),
  signal = rep(c("Signal 1", "Signal 2"), each = 500)
) |>
  mutate(
    hilbert = sp$hilbert(amplitude),
    envelope = abs(hilbert)
  )

cor(
  x = aec_tbl_2 |> filter(signal == "Signal 1") |> pull(envelope),
  y = aec_tbl_2 |> filter(signal == "Signal 2") |> pull(envelope)
)

ggplot(aec_tbl_2, aes(x = time)) +
  geom_path(aes(y = amplitude, colour = "Real")) +
  geom_path(aes(y = envelope, colour = "Envelope"), linewidth = 1) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s")
  ) +
  scale_colour_manual(
    values = c("black", "#cd201f"), breaks = c("Real", "Envelope")
  ) +
  facet_wrap(vars(signal), ncol = 1) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    colour = "Signal"
  ) +
  theme_classic()

aec_sim <- map_df(
  1:1000,
  function(x) {
    # Making waves
    mod_wave <- sine_wave(time, frequency = 3, phase = phase, amplitude = .5) +
      rnorm(length(time), 0, .3)

    carrier_wave_3 <- sine_wave(time, frequency = 50, phase = phase, amplitude = 1) +
      rnorm(length(time), 0, .2)
    carrier_wave_4 <- sine_wave(time, frequency = 30, phase = phase, amplitude = 1) +
      rnorm(length(time), 0, .2)

    am_wave_3 <- (1 + mod_wave) * carrier_wave_3
    am_wave_4 <- (1 + mod_wave) * carrier_wave_4

    aec_tbl_2 <- tibble(
      time = c(time, time),
      amplitude = c(am_wave_3, am_wave_4),
      signal = rep(c("Signal 1", "Signal 2"), each = 500)
    ) |>
      mutate(
        hilbert = sp$hilbert(amplitude),
        envelope = abs(hilbert)
      )

    # Correlations
    amplitude_corr <- cor(
      x = aec_tbl_2 |> filter(signal == "Signal 1") |> pull(amplitude),
      y = aec_tbl_2 |> filter(signal == "Signal 2") |> pull(amplitude)
    )

    envelope_corr <- cor(
      x = aec_tbl_2 |> filter(signal == "Signal 1") |> pull(envelope),
      y = aec_tbl_2 |> filter(signal == "Signal 2") |> pull(envelope)
    )

    # Results
    tibble(
      signal = c("Real", "Envelope"),
      r = c(amplitude_corr, envelope_corr)
    )
  }
)

sim_plot <- ggplot(aec_sim, aes(x = abs(r), fill = signal)) +
  geom_histogram(binwidth = .05) +
  scale_fill_manual(
    values = c("black", "#cd201f"), breaks = c("Real", "Envelope")
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "Absolute Correlation",
    y = "Count",
    fill = "Signal"
  ) +
  theme_classic()


# Garbage? ----

phase_1 <- 145
phase_2 <- 25

p3 <- ggplot() +
  geom_function(
    aes(linetype = "1 Hz"),
    fun = \(.x) sine_wave(.x, frequency = 5, phase = 145, amplitude = 1) + rnorm(.x, mean = 0, sd = .25),
    n = 1000
  ) +
  geom_function(
    aes(linetype = "2 Hz"),
    fun = \(.x) sine_wave(.x, frequency = 5, phase = 25, amplitude = 1) + rnorm(.x, mean = 0, sd = .25),
    n = 1000,
    colour = "red"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s")) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    linetype = "Frequency"
  ) +
  theme_classic()


phase_1 <- 25
phase_2 <- 125
phase_difference <- abs(phase_1 - phase_2)

p3 <- ggplot() +
  geom_point(aes(x = 0.5, y = 1, colour = "white")) +
  geom_function(
    aes(linetype = paste0(phase_1, "\u00BA")),
    fun = sine_wave,
    args = list(frequency = 4, phase = phase_1, amplitude = 1),
    n = 10^3
  ) +
  geom_function(
    aes(linetype = paste0(phase_2, "\u00BA")),
    fun = sine_wave,
    args = list(frequency = 4, phase = phase_2, amplitude = 1),
    n = 10^3
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(.x) paste0(.x, "s")
  ) +
  scale_color_manual(values = "white") +
  guides(
    linetype = guide_legend(order = 1),
    colour = guide_legend(
      order = 2,
      label = FALSE
    )
  ) +
  labs(
    x = "Time",
    y = "Amplitude (AU)",
    linetype = "Phase",
    colour = "Phase Angle\nDifference"
  ) +
  theme_classic()

p3

p4 <- ggplot() +
  geom_segment(
    aes(x = phase_difference, xend = phase_difference, y = 0, yend = 1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 360, by = 90),
    expand = c(0, 0),
    limits = c(0, 360),
    labels = function(.x) paste0(.x, "\u00BA"),
    position = "top"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 1)) +
  coord_polar(start = 270 / (180 / pi), direction = -1, clip = "off") +
  labs(
    x = "",
    y = ""
  ) +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(hjust = 0)
  )

patch_2 <- p3 +
  theme(plot.margin = margin(r = 10)) +
  inset_element(
    p4,
    left = 0.5, bottom = 0.1, right = 1.67, top = .39
  )

ggplot2::ggsave(
  filename = "figures/phase-plot.png",
  plot = patch_2,
  device = ragg::agg_png,
  width = 21.59,
  height = 12,
  units = "cm",
  dpi = "retina",
  scaling = 1
)
