# PLI polar plots ----

phase_angle_differences <- tibble(
  angle = c(
    seq(150, 210, length.out = 10),
    seq(45,  105, length.out = 10)
  ),
  group = factor(
    c(
      rep("PLI = 0",  times = 10),
      rep("PLI = 1",  times = 10)
    ),
    levels = c("PLI = 0", "PLI = 1")
  )
)

phase_angle_difference_plot <- phase_angle_differences |>
  ggplot() +
  geom_vline(aes(xintercept = angle)) +
  facet_wrap(~ group) +
  coord_polar(start = -1.57, direction = -1, clip = "off") +
  scale_x_continuous(
    breaks = seq(0, 360, by = 45),
    expand = c(0, 0),
    limits = c(0, 360),
    labels = function(x) paste0(x, "\u00BA")
  ) +
  labs(
    title = "Example Phase Angle Difference Plots",
    x = "Phase Angle Difference"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    panel.spacing = unit(2, "lines")
  )

# AEC envelope plots ----

ndarray_to_r <- function(x) {

  channel_1 <- x[seq(1, length(x), 3)]
  channel_2 <- x[seq(2, length(x), 3)]

  tibble(
    Channel = c(
      rep("EEG Channel 1", 1001),
      rep("EEG Channel 2", 1001)
    ),
    Sample = rep(1:1001, times = 2),
    Amplitude = c(channel_1, channel_2)
  )

}

# This data comes from the corresponding objects in `scratch.py`
real_data <- ndarray_to_r(py$first_real) |>
  mutate(Signal = "Real", path_size = 0.2)
amplitude_data <- ndarray_to_r(py$first_amplitude) |>
  mutate(Signal = "Envelope", path_size = 0.5)

signal_data <- rbind(real_data, amplitude_data)

signal_cor <- round(
  cor(
    filter(amplitude_data, Channel == "EEG Channel 1")$Amplitude,
    filter(amplitude_data, Channel == "EEG Channel 2")$Amplitude
  ),
  digits = 2
)

amplitude_envelope_plot <- signal_data |>
  ggplot(aes(x = Sample, y = Amplitude, colour = Signal)) +
  geom_path(aes(size = path_size)) +
  scale_size_identity() +
  scale_color_manual(
    values = c("#cd201f", "black")
  ) +
  scale_y_continuous(labels = function(x) x*1000000) +
  facet_wrap(~ Channel) +
  geom_text(
    data = tibble(
      Channel = "EEG Channel 2",
      label = paste0(
        "AEC = 0.07\n", # Estimated in VS Code
        "    r = ", signal_cor
      )
    ),
    mapping = aes(x = Inf, y = Inf, label = label),
    hjust   = 1.1,
    vjust   = 1.25,
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = "Example Amplitude Envelope Plot",
    y = "Amplitude (Microvolts)"
  ) +
  ggplot2::theme_classic()

# patchwork ----
pli_aec_patch <- phase_angle_difference_plot +
  amplitude_envelope_plot +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")

ggplot2::ggsave(
  filename = here("manuscripts/presentation/figures/pli-aec.png"),
  plot = pli_aec_patch,
  device = ragg::agg_png,
  width = 26.4583,
  height = 25.8233,
  units = "cm",
  dpi = "retina",
  scaling = 2
)
