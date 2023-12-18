# Helper functions ------------------------------------------------------------

#' Sum integers from 1 to n
#'
#' The number of observations in the lower triangle of a matrix can be calculated
#' by summing integers from 1 to n.
#'
#' @param n The largest integer.
sum_ascending <- function(n) {
  n * (n + 1) / 2
}

# Prerequisites to determining n_within and n_between -------------------------

# Per participant values.
n_sessions <- 10
n_states <- 5
n_sessions_split_half <- 2

# Figure 3A, S2C: Number of "boxes" on the diagonal of the matrices.
n_participants <- 9

# Figure 3A, S2C: Number of "boxes" on the off-diag, lower-tri of the matrices.
n_participant_pairs <- sum_ascending(n_participants - 1)

# Figure S2C: Per participant number of observations on the diagonal of the
# matrix.
n_recordings <- n_sessions * n_states
# Figure S2C: Total number of observations on the diagonal of the matrix.
n_recordings_total <- n_recordings * n_participants

# MSC09 had less recordings than all other participants. Here we take the difference
# so it can be used to correct n_recordings_total.
n_recordings_MSC09 <- 35
n_recordings_MSC09_diff <- n_recordings - n_recordings_MSC09
n_recordings_total <- n_recordings_total - n_recordings_MSC09_diff

# Figure 3A: Per participant number of observations on the diagonal of the matrix.
n_recordings_split_half <- n_sessions_split_half * n_states
# Figure 3A: Total number of observations on the diagonal of the matrix.
n_recordings_split_half_total <- n_recordings_split_half * n_participants

# MSC09 had less split-half recordings than all other participants. Here we take
# the difference so it can be used to correct n_recordings_split_half_total
n_recordings_MSC09_split_half <- 9
n_recordings_MSC09_diff_split_half <-
  n_recordings_split_half - n_recordings_MSC09_split_half
n_recordings_split_half_total <-
  n_recordings_split_half_total - n_recordings_MSC09_diff_split_half

# Cohen's d for separate sessions ---------------------------------------------

# Means and standard errors were determined from Figure S2D using
# WebPlotDigitizer: https://automeris.io/WebPlotDigitizer/
mean_within <- 0.6270758122743683
se_within <- 0.03249097

mean_between <- 0.3736462093862817
se_between <- 0.01299639

mean_difference <- 0.2534296

# Figure S2C: Number of observations in the lower-tri of the diagonal "boxes" of
# the matrix.
n_within <- sum_ascending(n_recordings - 1) * n_participants
# Correction for MSC09
n_within_MSC09 <- sum_ascending(n_recordings_MSC09 - 1)
n_within_MSC09_diff <- sum_ascending(n_recordings - 1) - n_within_MSC09
n_within <- n_within - n_within_MSC09_diff

# Figure S2C: Number of observations in the off-diagonal "boxes" of the matrix.
# Note that this can also be calculated with:
# n_between <- sum_ascending(n_recordings_total - 1) - n_within
n_between <- n_recordings^2 * n_participant_pairs
# Correction for MSC09
n_between_MSC09 <- n_recordings_MSC09^2 * (n_participants - 1)
n_between_MSC09_diff <- n_recordings^2 * (n_participants - 1) - n_between_MSC09
n_between <- n_between - n_between_MSC09_diff

sd_within <- se_within * sqrt(n_within)
sd_between <- se_between * sqrt(n_between)

sd_pooled <- sqrt(
  sum(c(sd_within^2, sd_between^2))/2
)

d <- mean_difference / sd_pooled

# Cohen's d for split-half sessions -------------------------------------------

# Means and standard errors were determined from Figure S2E using
# WebPlotDigitizer: https://automeris.io/WebPlotDigitizer/
mean_within_split_half <- 1.083941605839416
se_within_split_half <- 0.03941606

mean_between_split_half <- 0.5649635036496349
se_between_split_half <- 0.01313869

mean_difference_split_half <- 0.5189781

# Figure 3A: Number of observations in the lower-tri of the diagonal "boxes" of
# the matrix.
n_within_split_half <-
  sum_ascending(n_recordings_split_half - 1) * n_participants
# Correction for MSC09
n_within_MSC09_split_half <- sum_ascending(n_recordings_MSC09_split_half - 1)
n_within_MSC09_diff_split_half <-
  sum_ascending(n_recordings_split_half - 1) - n_within_MSC09_split_half
n_within_split_half <- n_within_split_half - n_within_MSC09_diff_split_half

# Figure 3A: Number of observations in the off-diagonal "boxes" of the matrix.
n_between_split_half <- n_recordings_split_half^2 * n_participant_pairs
# Correction for MSC09
n_between_MSC09_split_half <-
  n_recordings_MSC09_split_half^2 * (n_participants - 1)
n_between_MSC09_diff_split_half <-
  n_recordings_split_half^2 * (n_participants - 1) - n_between_MSC09_split_half
n_between_split_half <- n_between_split_half - n_between_MSC09_diff_split_half

sd_within_split_half <- se_within_split_half * sqrt(n_within_split_half)
sd_between_split_half <- se_between_split_half * sqrt(n_between_split_half)

sd_pooled_split_half <- sqrt(
  sum(c(sd_within_split_half^2, sd_between_split_half^2))/2
)

d_split_half <- mean_difference_split_half / sd_pooled_split_half
