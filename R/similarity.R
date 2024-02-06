#' Create an ordered vector of participant IDs
#'
#' @param participant_ids A character vector of participant IDs.
#'
#' @return A character vector of ordered participant IDs.
case_order <- function(participant_ids) {

  tidyr::expand_grid(
    participant = participant_ids,
    session_state = factor(
      c(
        "pre_rc1",
        "pre_rc2",
        "post_rc1",
        "post_rc2",
        "fu_rc1",
        "fu_rc2",
        "pre_ro1",
        "pre_ro2",
        "post_ro1",
        "post_ro2",
        "fu_ro1",
        "fu_ro2"
      )
    )
  ) |>
    dplyr::arrange(participant) |>
    dplyr::mutate(label = paste0(participant, "_", session_state)) |>
    purrr::pluck("label")

}

#' Generate possible pairwise comparisons of functional connectomes
#'
#' @param input List of connectivity matrices and their metadata.
#'
#' @return A tibble.
pairwise_comparisons <- function(input) {

  # Get levels for factors
  participant_levels <- input |>
    purrr::map(~ .x$metadata$participant) |>
    unlist() |>
    unique() |>
    sort()
  session_levels <- c("pre", "post", "fu")
  state_levels <- c("rc1", "rc2", "ro1", "ro2")
  label_levels <- case_order(participant_levels)

  # Get the branch and label names that go together so pair names are human
  # readable
  label_names  <- input |>
    purrr::map(~ .x$metadata$case) |>
    unlist()
  label_names <- label_names |>
    factor(levels = label_levels) |>
    sort()
  branch_names <- names(label_names)

  connectivity_matrix_names <- tibble::tibble(
    label_names = label_names,
    branch_names = branch_names
  )

  # Generate all unique pairwise comparisons over branch names. Make sure to
  # allow repeats for self-comparison pairs. This returns the lower triangle and
  # diagonal pairs.
  branch_combinations <- outer(branch_names, branch_names, "paste", sep = "_x_")
  branch_combinations <- stringr::str_split(
    branch_combinations[lower.tri(branch_combinations, diag = TRUE)],
    "_x_"
  )
  branch_combinations <- tibble::tibble(
    x_branch = purrr::map_chr(branch_combinations, 1),
    y_branch = purrr::map_chr(branch_combinations, 2)
  )

  connectivity_matrix_pairs <- branch_combinations |>
    # Using a join ensures the branch and label names will line up
    dplyr::full_join(
      connectivity_matrix_names,
      by = c("x_branch" = "branch_names")
    ) |>
    dplyr::rename(x_label = label_names) |>
    dplyr::full_join(
      connectivity_matrix_names,
      by = c("y_branch" = "branch_names")
    ) |>
    dplyr::rename(y_label = label_names) |>
    dplyr::mutate(
      # Add metadata to assist with analyses and filtering
      x_participant = factor(
        stringr::str_extract(x_label, "^.*?(?=_)"),
        levels = participant_levels
      ),
      x_session = factor(
        stringr::str_extract(x_label, "(?<=_)(.*)(?=_)"),
        levels = session_levels
      ),
      x_state = factor(
        stringr::str_extract(x_label, "([^_]*)$"),
        levels = state_levels
      ),
      y_participant = factor(
        stringr::str_extract(y_label, "^.*?(?=_)"),
        levels = participant_levels
      ),
      y_session = factor(
        stringr::str_extract(y_label, "(?<=_)(.*)(?=_)"),
        levels = session_levels
      ),
      y_state = factor(
        stringr::str_extract(y_label, "([^_]*)$"),
        levels = state_levels
      ),
      within_participant = ifelse(x_participant == y_participant, TRUE, FALSE),
      within_session = ifelse(x_session == y_session, TRUE, FALSE),
      within_state = ifelse(
        stringr::str_detect(x_state, stringr::str_sub(y_state, 1, 2)),
        TRUE,
        FALSE
      ),
      diagonal = ifelse(x_label == y_label, TRUE, FALSE),
      pair_label = paste(x_label, y_label, sep = "_x_"),
      # There are no equivalent participant pairs to collapse here since we're
      # only using the lower triangle of the matrix.
      pair_participant = factor(
        paste(x_participant, y_participant, sep = "_x_")
      ),
      # Relevel the participant pairs so the first 14 levels are within
      # participant and the remaining levels are between participant.
      pair_participant = forcats::fct_relevel(
        pair_participant,
        c(
          "P03_x_P03",
          "P04_x_P04",
          "P06_x_P06",
          "P07_x_P07",
          "P08_x_P08",
          "P09_x_P09",
          "P11_x_P11",
          "P14_x_P14",
          "P16_x_P16",
          "P17_x_P17",
          "P19_x_P19",
          "P20_x_P20",
          "P21_x_P21",
          "P22_x_P22"
        ),
        after = 0
      ),
      pair_session = paste(x_session, y_session, sep = "_x_"),
      # Equivalent session pairs should be collapsed into a single factor so
      # they aren't treated as different factors.
      pair_session = forcats::fct_collapse(
        pair_session,
        pre_x_pre = c("pre_x_pre"),
        post_x_post = c("post_x_post"),
        fu_x_fu = c("fu_x_fu"),
        pre_x_post = c("pre_x_post", "post_x_pre"),
        pre_x_fu = c("pre_x_fu", "fu_x_pre"),
        post_x_fu = c("post_x_fu", "fu_x_post")
      ),
      # Relevel the session pairs so the first 3 levels are within
      # session and the remaining levels are between session
      pair_session = forcats::fct_relevel(
        pair_session,
        c("pre_x_pre", "post_x_post", "fu_x_fu"),
        after = 0
      ),
      pair_state = paste(x_state, y_state, sep = "_x_"),
      # Equivalent state pairs should be collapsed into a single factor so they
      # aren't treated as different factors. It doesn't matter whether the state
      # came before or after the task, so base this on the state name without
      # the number.
      pair_state = forcats::fct_collapse(
        pair_state,
        rc_x_rc = c("rc1_x_rc1", "rc2_x_rc2", "rc1_x_rc2", "rc2_x_rc1"),
        ro_x_ro = c("ro1_x_ro1", "ro2_x_ro2", "ro1_x_ro2", "ro2_x_ro1"),
        rc_x_ro = c(
          "rc1_x_ro1",
          "ro1_x_rc1",
          "rc1_x_ro2",
          "ro2_x_rc1",
          "rc2_x_ro1",
          "ro1_x_rc2",
          "rc2_x_ro2",
          "ro2_x_rc2"
        )
      ),
      # Relevel the state pairs so the first 3 levels are within
      # state and the remaining levels are between session
      pair_state = forcats::fct_relevel(
        pair_state, c("rc_x_rc", "ro_x_ro"), after = 0
      )
    )

  connectivity_matrix_pairs

}

# TODO: Capture warnings. Different approaches shared here:
# https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#' Estimate similarity between list of connectivity matrices
#'
#' @param input The list of connectivity matrices.
#'
#' @return A tibble.
estimate_similarity <- function(input) {

  # Set up pairwise comparisons to prepare for iteration ----
  connectivity_matrix_pairs <- pairwise_comparisons(input)
  x_branch <- connectivity_matrix_pairs$x_branch
  y_branch <- connectivity_matrix_pairs$y_branch
  ## The names here are the key used for joining the connectivity_matrix_pairs
  ## and similarity_results tibbles later in this function.
  names(x_branch) <- connectivity_matrix_pairs$pair_label

  # Estimate similarity for all pairs ----
  similarity_estimates <- purrr::map2(
    x_branch, y_branch, ~{

      x_matrix <- input[[.x]][["connectivity_matrix"]]
      y_matrix <- input[[.y]][["connectivity_matrix"]]

      # Note: The FactoMineR::coeffRV() function centres the matrices in the
      # background.
      # Note 2: This function can give a warning "NaNs produced" for some
      # estimates. When this happens the skewness value will be "0.000000e+00"
      # and the p-value will be "NaN". This is likely from the permutation test
      # failing. This isn't an issue for our study since  we are only concerned
      # with the RV coefficient, not it's p-value
      rv_coeff <- FactoMineR::coeffRV(x_matrix, y_matrix)

      rv_coeff
    }
  )

  # Collect similarity results ----
  similarity_results <- tibble::tibble(
    pair_label = names(similarity_estimates), # Key for joining tibbles
    rv = purrr::map_dbl(similarity_estimates, "rv"),
    rv.std = purrr::map_dbl(similarity_estimates, "rvstd"),
    mean = purrr::map_dbl(similarity_estimates, "mean"),
    variance = purrr::map_dbl(similarity_estimates, "variance"),
    skewness = purrr::map_dbl(similarity_estimates, "skewness"),
    p.value = purrr::map_dbl(similarity_estimates, "p.value")
  )

  similarity_results <- dplyr::left_join(
    connectivity_matrix_pairs,
    similarity_results
  )

  similarity_results

}

#' Make the similarity results symmetric for plotting
#'
#' @param similarity_results A tibble of similarity results.
#'
#' @return A tibble of similarity results, with values for the lower and upper
#'   triangle of a matrix.
make_symmetric <- function(similarity_results) {

  # Only the lower triangle and diagonal are included in the similarity results
  # so we need to make the upper triangle in order to plot a symmetrical matrix.
  similarity_results_upper <- similarity_results |>
    dplyr::filter(diagonal == FALSE) |>
    dplyr::rename_with(~ stringr::str_replace(.x, "x", "z"), starts_with("x")) |>
    dplyr::rename_with(~ stringr::str_replace(.x, "y", "x"), starts_with("y")) |>
    dplyr::rename_with(~ stringr::str_replace(.x, "z", "y"), starts_with("z")) |>
    dplyr::mutate(pair_label = paste(x_label, y_label, sep = "_x_"))

  # Combine the upper, lower, and diagonals to prepare for plotting
  similarity_results_symmetric <- rbind(similarity_results, similarity_results_upper)

  similarity_results_symmetric

}

#' Plot a similarity matrix
#'
#' @param similarity_results A tibble of similarity results.
#' @param estimate The similarity estimate to plot. Options are `rv`.
#'
#' @return A ggplot.
plot_similarity <- function(similarity_results, estimate) {

  # Combine the upper, lower, and diagonals to prepare for plotting
  similarity_results_symmetric <- make_symmetric(similarity_results)

  # Plot

  ## The axis-tick colours are being used as annotations for the different
  ## sessions and states. This is a bit hackish, but works for our purposes.
  axis_tick_colours <- colorspace::diverging_hcl(
    12, h = c(299, 135), c = 60, l = c(20, 80), power = c(0.7, 1.3)
  )

  case_colours <- rep(
    c(rev(axis_tick_colours[7:12]), axis_tick_colours[1:6]),
    times = 14
  )

  similarity_results_symmetric |>
    ggplot2::ggplot(ggplot2::aes(x = x_label, y = y_label, fill = {{estimate}})) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, 1),
      option = "viridis",
      guide  = ggplot2::guide_colourbar(order = 1)
    ) +
    # Expand reduces the spacing between facets
    ggplot2::scale_x_discrete(expand = c(0, 0), position = "top") +
    # Reverse y-axis so the diagonal goes from the upper-left to lower-right
    ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev) +
    # ggplot2 does not natively support secondary axes for discrete variables,
    # but they can be added using ggh4x. The choice to use secondary axes is
    # purely aesthetic, so that the axis titles and facet strips are on opposite
    # sides of the matrix.
    ggplot2::guides(
      x.sec = ggh4x::guide_axis_manual(title = "Functional Connectomes"),
      y.sec = ggh4x::guide_axis_manual(title = "Functional Connectomes")
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "RV"
    ) +
    facet_grid(
      cols = ggplot2::vars(x_participant),
      rows = ggplot2::vars(y_participant),
      scales = "free",
      switch = "y"
    ) +
    ggplot2::theme(
      panel.spacing = grid::unit(0.1, "lines"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      strip.text.y.left = ggplot2::element_text(angle = 0),
      strip.placement = "outside",
      strip.switch.pad.grid = grid::unit(0, "lines"),
      axis.text = ggplot2::element_blank(),
      axis.ticks.x.top = ggplot2::element_line(colour = rev(case_colours), size = 1),
      axis.ticks.x.bottom = ggplot2::element_blank(),
      axis.ticks.y.left = ggplot2::element_line(colour = case_colours, size = 1),
      # axis.title.y.right = ggplot2::element_text(angle = 90),
      axis.ticks.y.right = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(0.5, "lines")
    )

}

#' Plot the within-participant similarity block of a participant
#'
#' @param similarity_results A tibble of similarity results.
#' @param participant A participant ID.
#'
#' @return A ggplot.
plot_similarity_key <- function(similarity_results, participant) {

  similarity_results |>
    dplyr::filter(
      x_participant == participant & y_participant == participant
    ) |>
    plot_similarity(rv) +
    ggplot2::scale_y_discrete(
      expand = c(0, 0),
      limits = rev,
      position = "right"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(),
      axis.ticks  = ggplot2::element_line(),
      legend.position = "none"
    )

}

#' Highlight a within-participant block in the similarity matrix
#'
#' @param similarity_results A tibble of similarity results.
#' @param participant A participant ID.
#'
#' @return A ggplot.
plot_similarity_highlight <- function(similarity_results, participant) {

  # Combine the upper, lower, and diagonals to prepare for plotting
  similarity_results_within <- similarity_results |>
    dplyr::filter(
      x_participant == participant & y_participant == participant
    ) |>
    make_symmetric()

  similarity_results |>
  plot_similarity(estimate = rv) +
    ggplot2::scale_fill_gradient(
      low = "black",
      high = "white",
      limits = c(0, 1),
      guide = "none"
    ) +
    ggplot2::geom_rect(
      data = similarity_results_within,
      fill = NA,
      colour = "#fc9400",
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      size = 2.4
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_raster(data = similarity_results_within, ggplot2::aes(fill = rv)) +
    ggplot2::scale_fill_continuous(limits = c(0, 1), type = "viridis") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(legend.position = "none")

}

#' Plot hypothetical similarity matrices
#'
#' @param similarity_results A tibble of similarity results.
#'
#' @return List of ggplots.
plot_similarity_archetype <- function(similarity_results) {

  group_effect <- similarity_results |>
    dplyr::mutate(rv = 1.0) |>
    plot_similarity(estimate = rv)

  session_effect <- similarity_results |>
    dplyr::mutate(
      rv = dplyr::case_when(
        within_session == TRUE ~ 1.0,
        within_session == FALSE ~ 0.0
      )
    ) |>
    plot_similarity(estimate = rv)

  state_effect <- similarity_results |>
    dplyr::mutate(
      rv = dplyr::case_when(
        within_state == TRUE ~ 1.0,
        within_state == FALSE ~ 0.0
      )
    ) |>
    plot_similarity(estimate = rv)

  individual_effect <- similarity_results |>
    dplyr::mutate(
      rv = dplyr::case_when(
        within_participant == TRUE ~ 1.0,
        within_participant == FALSE ~ 0.0
      )
    ) |>
    plot_similarity(estimate = rv)

  individual_session_effect <- similarity_results |>
    dplyr::mutate(
      rv = dplyr::case_when(
        within_participant == TRUE & within_session == TRUE ~ 1.0,
        within_participant == TRUE & within_session == FALSE ~ 0.0,
        within_participant == FALSE ~ 0.0
      )
    ) |>
    plot_similarity(estimate = rv)

  individual_state_effect <- similarity_results |>
    dplyr::mutate(
      rv = dplyr::case_when(
        within_participant == TRUE & within_state == TRUE ~ 1.0,
        within_participant == TRUE & within_state == FALSE ~ 0.0,
        within_participant == FALSE ~ 0.0
      )
    ) |>
    plot_similarity(estimate = rv)

  list(
    group_effect              = group_effect,
    session_effect            = session_effect,
    state_effect              = state_effect,
    individual_effect         = individual_effect,
    individual_session_effect = individual_session_effect,
    individual_state_effect   = individual_state_effect
  )

}

# TODO: Descriptives
summarize_similarity <- function() {

}

#' Fit a mixed beta regression to similarity results
#'
#' @param formula A regression formula.
#' @param data A tibble of similarity results.
#'
#' @note Since each of the predictors are categorical, it doesn't make sense to
#  include random slopes in this model.
#'
#' @return A glmmTMB model object.
glmmTMB_similarity <- function(formula, data) {

  similarity_results <- data |>
    # Diagonals shouldn't be included in modelling since they're comparing a
    # recording with itself
    dplyr::filter(diagonal == FALSE)

  similarity_fit <- glmmTMB::glmmTMB(
    formula,
    data = similarity_results,
    family = glmmTMB::beta_family(),
    control = glmmTMB::glmmTMBControl(
      # Parallel processing doesn't appear to be supported in a targets
      # pipeline, so make sure to use serial processing to avoid a warning.
      # At least that's what this should do; the warning still appears but
      # it's ignorable so that's okay.
      parallel = 1,
      # The fit in the beta band has a false convergence using the default
      # optimizer. Changing the optimizer resolves the convergence failure
      # and the results are either the same or nearly identical, so this can
      # be brushed off as a false positive failure.
      optimizer = optim,
      optArgs = list(method = "BFGS")
    )
  )

  similarity_fit

}

#' Check uniformity of GL(M)M's residuals
#'
#' `check_uniformity()` checks generalized linear (mixed) models for uniformity
#' of randomized quantile residuals, which can be used to identify typical model
#' misspecification problems, such as over/underdispersion, zero-inflation, and
#' residual spatial and temporal autocorrelation.
#'
#' @param object Fitted model.
#'
#' @details
#'
#' See `vignette("DHARMa")`.
#'
#' @references
#'
#' - Hartig, F., & Lohse, L. (2022). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models (Version 0.4.5). Retrieved from https://CRAN.R-project.org/package=DHARMa
#' - Dunn, P. K., & Smyth, G. K. (1996). Randomized Quantile Residuals. Journal of Computational and Graphical Statistics, 5(3), 236. https://doi.org/10.2307/1390802
#'
#' @return A ggplot.
check_uniformity <- function(object) {

  # Simulated residuals; see vignette("DHARMa")
  simulated_residuals <- DHARMa::simulateResiduals(object)

  dp <- list(min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)
  ggplot2::ggplot(
    tibble::tibble(scaled_residuals = residuals(simulated_residuals)),
    ggplot2::aes(sample = scaled_residuals)
  ) +
    # The confidence band is working here; it's just really tight so it's hard to see.
    qqplotr::stat_qq_band(distribution = "unif", dparams = list(min = 0, max = 1), alpha = .2, conf = .99) +
    qqplotr::stat_qq_line(distribution = "unif", dparams = dp, size = .8, colour = "#3aaf85") +
    qqplotr::stat_qq_point(distribution = "unif", dparams = dp, size = .5, alpha = .05, colour = "#1b6ca8") +
    ggplot2::labs(
      title = "Uniformity of Residuals",
      subtitle = "Dots should fall along the line",
      x = "Uniform Distribution Quantiles",
      y = "Sample Quantiles"
    ) +
    see::theme_lucid()
}

#' Construct a reference grid from the mixed beta regression
#'
#' @param object A glmmTMB model object.
#'
#' @return An emmeans model object.
emmeans_similarity <- function(object) {

  emmeans::emmeans(
    object,
    specs = list(
      main_effect = ~ within_participant,
      # Two-way interactions
      participant_session = ~ within_participant | within_session,
      participant_state = ~ within_participant | within_state,
      # Three-way interactions
      participant_session_state = ~ within_participant | within_state + within_session
    ),
    # Back-transform from the logit scale to the response scale to make
    # interpretation easier. Back-transformation is done using the regrid argument
    # so that the reference grid is reparameterized to the response scale. This
    # ensures that subsequent EMMs and contrasts will be conducted on the response
    # scale. See ?emmeans::regrid for details.
    regrid = "response"
  )

}

#' Tidy emmeans
#'
#' @param object An emmeans model object.
#'
#' @return A tidy tibble.
tidy_emmeans_similarity <- function(object) {

  effect_labels <- c(
    "main_effect",
    "between_sessions",
    "within_sessions",
    "between_states",
    "within_states",
    "between_sessions_states",
    "between_sessions_within_states",
    "within_sessions_between_states",
    "within_sessions_states"
  )

  effect_labels_nice <- c(
    "Main effect",
    "Between sessions",
    "Within sessions",
    "Between states",
    "Within states",
    "Between sessions and states",
    "Between sessions and within states",
    "Within sessions and between states",
    "Within sessions and states"
  )

  object |>
    purrr::map_dfr(broom::tidy, .id = "effect") |>
    dplyr::mutate(
      effect_label = rep(effect_labels, each = 2),
      within_participant = dplyr::case_when(
        within_participant == TRUE  ~ "within_participants",
        within_participant == FALSE ~ "between_participants"
      ),
      effect_label = paste0(within_participant, "_", effect_label),
      effect_label_nice = rep(effect_labels_nice, each = 2)
    )

}

#' Estimate pairwise contrasts using the reference grid
#'
#' @param object An emmeans model object.
#'
#' @return An emmeans model object.
contrast_similarity <- function(object) {

  emmeans::contrast(
    object,
    method = "revpairwise", # Subtracts the between value from the within value
    infer = TRUE
  )

}

#' Tidy emmeans contrasts
#'
#' @param object An emmeans model object.
#'
#' @return A tidy tibble.
tidy_contrast_similarity <- function(object) {

  effect_labels <- c(
    "Main effect",
    "Between sessions",
    "Within sessions",
    "Between states",
    "Within states",
    "Between sessions and states",
    "Between sessions and within states",
    "Within sessions and between states",
    "Within sessions and states"
  )

  contrasts <- object |>
    purrr::map_dfr(broom::tidy, .id = "effect", conf.int = TRUE) |>
    dplyr::mutate(
      effect_label = factor(effect_labels, levels = effect_labels),
      effect_label = forcats::fct_relevel(
        effect_label, "Between sessions and states", after = 7
      ),
      dist = distributional::dist_student_t(
        df = df, mu = estimate, sigma = std.error
      )
    ) |>
    ggdist::point_interval(dist, .width = c(.66, .95)) |>
    tidyr::pivot_wider(names_from = .width, values_from = c(.lower, .upper)) |>
    dplyr::mutate(
      effect_direction = dplyr::case_when(
        .lower_0.95 > 0 & .upper_0.95 > 0 ~ "Positive",
        .lower_0.95 < 0 & .upper_0.95 > 0 ~ "Unresolved",
        .lower_0.95 < 0 & .upper_0.95 < 0 ~ "Negative"
      ),
      effect_direction = factor(
        effect_direction,
        levels = c("Positive", "Unresolved", "Negative")
      )
    )

  contrasts

}

#' Estimate standardized mean differences for contrasts
#'
#' @param model A glmmTMB model object.
#' @param contrasts An emmeans contrasts object.
#'
#' @return An emmeans object.
contrast_similarity_cohens_d <- function(model, contrasts) {

  purrr::map(
    contrasts,
    \(.x) {
      emmeans::eff_size(
        .x,
        # Add variances then take the square root to get total SD of the
        # variance components.
        sigma = sqrt(sum(insight::get_variance_intercept(model))),
        edf = df.residual(model),
        method = "identity"
      )
    }
  )

}

#' Tidy emmeans standardized mean differences
#'
#' @param object An emmeans model object.
#'
#' @return A tidy tibble.
tidy_contrast_similarity_cohens_d <- function(object) {

  effect_labels <- c(
    "Main effect",
    "Between sessions",
    "Within sessions",
    "Between states",
    "Within states",
    "Between sessions and states",
    "Between sessions and within states",
    "Within sessions and between states",
    "Within sessions and states"
  )

  contrasts <- object |>
    purrr::map_dfr(tibble::as_tibble, .id = "effect") |>
    dplyr::mutate(
      effect_label = factor(effect_labels, levels = effect_labels),
      effect_label = forcats::fct_relevel(
        effect_label, "Between sessions and states", after = 7
      )
    ) |>
    dplyr::select(
      effect_label,
      effect,
      smd = effect.size,
      smd_std.error = SE,
      smd_conf.low = lower.CL,
      smd_conf.high = upper.CL
    )

  contrasts

}

#' Plot similarity contrast intervals
#'
#' @param object An emmeans model object.
#'
#' @return A ggplot.
plot_similarity_contrasts <- function(object) {

  # Tidy data to prepare for plotting
  contrasts <- tidy_contrast_similarity(object)

  # Plotting
  ggplot2::ggplot(contrasts, ggplot2::aes(x = estimate, y = effect_label)) +
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
      breaks = c(-0.1, 0, 0.1),
      labels = function(x) ifelse(x == 0, 0, x),
      sec.axis = ggplot2::sec_axis(
        ~ .,
        breaks = c(-0.075, 0.075),
        labels = c(
          "Less similar\nwithin participant",
          "More similar\nwithin participant"
        )
      )
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
      y = "Effect",
      colour = "Effect Direction",
      colour_ramp = "Compatibility\nInterval"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.line.x.top  = ggplot2::element_blank(),
      axis.text.x.top  = ggplot2::element_text(size = 11, colour = "black")
    )

}

#' Subset similarity results by participant
#'
#' @param similarity_results A tibble of similarity results,
#' @param participants A character vector of participant IDs.
#'
#' @return A list of tibbles.
subset_similarity_results <- function(similarity_results, participants) {

  # The similarity results are technically shaped like the lower diagonal of a
  # matrix, which we don't want when subsetting by participant, since it leads
  # factor levels we don't want.
  similarity_results_ltri <- similarity_results
  # As a solution we can make the "matrix" symmetric
  similarity_results_utri <- similarity_results |>
    dplyr::filter(diagonal == FALSE, within_participant == FALSE) |>
    dplyr::rename_with(~ stringr::str_replace(.x, "x", "z"), tidyselect::starts_with("x")) |>
    dplyr::rename_with(~ stringr::str_replace(.x, "y", "x"), tidyselect::starts_with("y")) |>
    dplyr::rename_with(~ stringr::str_replace(.x, "z", "y"), tidyselect::starts_with("z")) |>
    dplyr::mutate(pair_label = paste(x_label, y_label, sep = "_x_"))

  similarity_results_sym <- rbind(
    similarity_results_ltri,
    similarity_results_utri
  )

  # This corresponds to the rows belonging to a given participant in the
  # similarity matrix, including both within and between participant
  # similarities.
  P00_similarities <- purrr::map(
    purrr::set_names(participants),
    ~{
      similarity_results_sym |>
        dplyr::filter(x_participant == .x)
    }
  )

  P00_similarities

}

#' Fit a mixed beta regression to a subset of similarity results
#'
#' @param formula A regression formula.
#' @param data A list of tibbles of similarity results.
#'
#' @note Since each of the predictors are categorical, it doesn't make sense to
#  include random slopes in this model.
#'
#' @return A list of glmmTMB model objects.
subset_glmmTMB_similarity <- function(formula, data) {

  purrr::map(
    data,
    ~{
      glmmTMB_similarity(formula, data = .x)
    }
  )

}

#' Construct a reference grid from the mixed beta regressions
#'
#' @param object A list of glmmTMB model objects.
#'
#' @return A list of emmeans model objects.
subset_emmeans_similarity <- function(objects) {

  purrr::map(
    objects,
    ~{
      emmeans_similarity(.x)
    }
  )

}

#' Estimate pairwise contrasts using the reference grids
#'
#' @param object A list of emmeans model objects.
#'
#' @return A list of emmeans model objects.
subset_contrast_similarity <- function(objects) {

  purrr::map(
    objects,
    ~{
      contrast_similarity(.x)
    }
  )

}

#' Tidy emmeans subset contrasts
#'
#' @param objects A list of emmeans model objects.
#' @param objects A list of glmmTMB model objects. Used for convergence checks.
#'
#' @return A list of tidy tibbles.
tidy_subset_contrast_similarity <- function(
  objects, models = NULL, check_convergence = TRUE
) {

  effect_labels <- c(
    "Main effect",
    "Between sessions",
    "Within sessions",
    "Between states",
    "Within states",
    "Between sessions and states",
    "Between sessions and within states",
    "Within sessions and between states",
    "Within sessions and states"
  )

  purrr::map2_dfr(
    objects, models,
    ~{
      contrasts <- .x |>
        purrr::map_dfr(broom::tidy, .id = "effect") |>
        dplyr::mutate(
          effect_label = factor(effect_labels, levels = effect_labels),
          effect_label = forcats::fct_relevel(
            effect_label, "Between sessions and states", after = 7
          )

        )

      if (check_convergence) {
        contrasts <- dplyr::mutate(
          contrasts,
          converged = performance::check_convergence(.y)
        )
      }

      contrasts
    },
    .id = "participant"
  ) |>
    dplyr::mutate(
      dist = distributional::dist_student_t(
        df = df, mu = estimate, sigma = std.error
      )
    ) |>
    ggdist::point_interval(dist, .width = c(.66, .95)) |>
    tidyr::pivot_wider(names_from = .width, values_from = c(.lower, .upper)) |>
    dplyr::mutate(
      effect_direction = dplyr::case_when(
        .lower_0.95 > 0 & .upper_0.95 > 0 ~ "Positive",
        .lower_0.95 < 0 & .upper_0.95 > 0 ~ "Unresolved",
        .lower_0.95 < 0 & .upper_0.95 < 0 ~ "Negative"
      ),
      effect_direction = factor(
        effect_direction,
        levels = c("Positive", "Unresolved", "Negative")
      )
    )

}

#' Plot subset similarity contrasts
#'
#' @param group_object A tibble of similarity results.
#' @param subset_objects A list of subset similarity results.
#' @param subset_models A list of subset similarity glmmTMB models.
#'
#' @return A grob.
plot_subset_similarity_contrasts <- function(
  group_object, subset_objects, subset_models
) {

  # Note: The gradient fill in the ridges is currently commented out, as it's a
  # bit busy when also filling the ridges by effect direction, with the colour
  # of a ridge becoming less clear when the effect sizes are very small.

  # Tidy data to prepare for plotting
  contrasts <- tidy_contrast_similarity(group_object)
  subset_contrasts <- tidy_subset_contrast_similarity(subset_objects, subset_models)

  subset_contrasts <- contrasts |>
    dplyr::select(effect_label, group_estimate = estimate) |>
    dplyr::left_join(subset_contrasts) |>
    # Since we have 9 contrasts, there's an odd number which makes lining up
    # facets awkward. We'll create a blank factor level, whose facet we will
    # later remove, to work around this.
    dplyr::mutate(
      effect_label = forcats::fct_expand(effect_label, " "),
      effect_label = forcats::fct_relevel(effect_label, "Main effect", " "),
      participant  = factor(participant)
    ) |>
    # Only show results for converged models.
    dplyr::filter(converged)

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
    # ggdist::scale_fill_ramp_continuous(
    #   from = "#f2f5c6", #ramp_colours[close_colour_index],
    #   name = "Distance from\ngroup-level\npoint estimate",
    #   trans = "identity",
    #   # This adjusts how quickly the fill ramps. It should be set to a meaningful
    #   # value such as the maximum distance we consider compatible with the
    #   # group-level estimate.
    #   limits = c(0, 0.05),
    #   # Fills outside the distance limit should be given the max colour
    #   oob = scales::oob_squish,
    #   guide = ggdist::guide_rampbar(
    #     to = "#a0a284", #"gray80",#"white", #ramp_colours[far_colour_index],
    #     title.vjust = 1.5,
    #     order = 2
    #   )
    # ) +
    # General plot settings
    ggplot2::scale_x_continuous(
      breaks = my_breaks,
      labels = function(x) ifelse(x == 0, 0, x)
    ) +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(1.1, 0)),
      drop = FALSE
    ) +
    ## These legends already appear in the group-level interval plot, so they
    ## don't need to be repeated here.
    ggplot2::guides(
      fill = "none",
      colour_ramp = "none",
      # Hacks to keep plot dimensions with an "invisible" legend
      colour = ggplot2::guide_legend(override.aes = list(color = "white"))
    ) +
    # Keep the empty " " facet to enforce facet alignment. It will be removed
    # manually after building the plot. Make the direction vertical so that
    # the "Main effect" facet is above the " " facet.
    ggplot2::facet_wrap(~ effect_label, ncol = 5, drop = FALSE, dir = "v") +
    ggplot2::labs(
      x = paste0(
        "Difference in functional connectome similarity",
        "\n",
        "(Within participant y \U2212 Between participant)"
      ),
      y = "Participant"
    ) +
    ggplot2::theme_classic() +
    # Hacks to keep plot dimensions with an "invisible" legend
    ggplot2::theme(
      legend.title = ggplot2::element_text(color = "white"),
      legend.text = ggplot2::element_text(color = "white"),
      legend.key = ggplot2::element_rect(color = "white", fill = "white")
    )

  # Remove the blank facet. Solution for removing facet modified from:
  # https://stackoverflow.com/a/30372692/16844576
  g <- ggplot2::ggplotGrob(p)

  # Get the grobs that must be removed
  rm_grobs <- g$layout$name %in% c("strip-t-1-2")
  # Remove grobs
  g$grobs[rm_grobs] <- NULL
  g$layout <- g$layout[!rm_grobs, ]
  # Move x-axis from " " facet to "Main effect" facet
  g$layout[g$layout$name == "axis-b-1-2", c("t", "b")] <- c(9, 9)
  # Move y-axis from " " facet to "Within sessions" facet
  g$layout[g$layout$name == "axis-l-2-1", c("l", "r")] <- c(8, 8)

  # Transform back to class ggplot and print
  ggpubr::as_ggplot(g)

}

#' Plot subset similarity contrasts
#'
#' This is the same as the function above, but includes the legend. Yes, the
#' function above should just give an argument to do this (DRY principles and
#' all that), but I need to make targets with and without the legend anyways
#' since the returned object can't be modified after the fact, and this was
#' easier.
#'
#' @param group_object A tibble of similarity results
#' @param subset_objects A list of subset similarity results
#' @param subset_models A list of subset similarity glmmTMB models.
#'
#' @return A grob.
plot_subset_similarity_contrasts_2 <- function(
  group_object, subset_objects, subset_models
) {

  # Tidy data to prepare for plotting
  contrasts <- tidy_contrast_similarity(group_object)
  subset_contrasts <- tidy_subset_contrast_similarity(subset_objects, subset_models)

  subset_contrasts <- contrasts |>
    dplyr::select(effect_label, group_estimate = estimate) |>
    dplyr::left_join(subset_contrasts) |>
    # Since we have 9 contrasts, there's an odd number which makes lining up
    # facets awkward. We'll create a blank factor level, whose facet we will
    # later remove, to work around this.
    dplyr::mutate(
      effect_label = forcats::fct_expand(effect_label, " "),
      effect_label = forcats::fct_relevel(effect_label, "Main effect", " "),
      participant  = factor(participant)
    ) |>
    # Only show results for converged models.
    dplyr::filter(converged)

  # Plotting
  my_scale <- ggplot2::scale_x_continuous(
    breaks = scales::extended_breaks(n = 3)
  )
  my_scale$train(
    c(min(subset_contrasts$.lower_0.95), max(subset_contrasts$.upper_0.95))
  )
  my_breaks <- my_scale$get_breaks()

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
      ),
      color = "black",
      size = 0.25,
      height = 2,
      expand = TRUE, # I can't disable this for whatever reason
      trim = FALSE, # I can't disable this for whatever reason
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
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(1.1, 0)),
      drop = FALSE
    ) +
    ## These legends don't need to appear.
    ggplot2::guides(
      fill = "none",
      colour_ramp = "none"
    ) +
    # Keep the empty " " facet to enforce facet alignment. It will be removed
    # manually after building the plot. Make the direction vertical so that
    # the "Main effect" facet is above the " " facet.
    ggplot2::facet_wrap(~ effect_label, ncol = 5, drop = FALSE, dir = "v") +
    ggplot2::labs(
      x = paste0(
        "Difference in functional connectome similarity",
        "\n",
        "(Within participant y \U2212 Between participant)"
      ),
      y = "Participant",
      colour = "Effect Direction"
    ) +
    ggplot2::theme_classic()

  # Remove the blank facet. Solution for removing facet modified from:
  # https://stackoverflow.com/a/30372692/16844576
  g <- ggplot2::ggplotGrob(p)

  # Get the grobs that must be removed
  rm_grobs <- g$layout$name %in% c("strip-t-1-2")
  # Remove grobs
  g$grobs[rm_grobs] <- NULL
  g$layout <- g$layout[!rm_grobs, ]
  # Move x-axis from " " facet to "Main effect" facet
  g$layout[g$layout$name == "axis-b-1-2", c("t", "b")] <- c(9, 9)
  # Move y-axis from " " facet to "Within sessions" facet
  g$layout[g$layout$name == "axis-l-2-1", c("l", "r")] <- c(8, 8)

  # Transform back to class ggplot and print
  ggpubr::as_ggplot(g)

}

#' Plot subset similarity contrast bar charts
#'
#' @param delta,theta,alpha,beta,gamma A list of subset similarity results for a
#'   given frequency band.
#'
#' @return A ggplot2 object.
plot_subset_similarity_contrast_bars <- function(
  delta = NULL, theta = NULL, alpha = NULL, beta = NULL, gamma = NULL
) {

  subset_contrasts <- purrr::map_dfr(
    list(
      Delta = delta, Theta = theta, Alpha = alpha, Beta = beta, Gamma = gamma
    ),
    function(.x) {

      if (is.null(.x)) {
        .x
      } else {
        .x |>
          tidy_subset_contrast_similarity(
            models = NA, check_convergence = FALSE
          ) |>
          dplyr::select(effect_label, effect_direction) |>
          dplyr::group_by(effect_label) |>
          dplyr::count(effect_direction)
      }

    },
    .id = "frequency_band"
  )

  subset_contrasts |>
    dplyr::mutate(
      frequency_band = factor(
        frequency_band,
        levels = c("Delta", "Theta", "Alpha", "Beta", "Gamma")
      )
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = effect_label, y = n, fill = effect_direction)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(breaks = c(0, 7, 14)) +
    ggplot2::scale_fill_discrete(
      drop = FALSE,
      type = c(
        "#3aaf85", # Green
        "#1b6ca8", # Blue
        "#cd201f"  # Red
      )
    ) +
    ggplot2::facet_wrap(ggplot2::vars(frequency_band), nrow = 1) +
    ggplot2::labs(
      x = "Contrast",
      y = "Count",
      fill = "Effect Direction"
    ) +
    ggplot2::theme_classic()
}

#' Save similarity matrix figure
#'
#' @param filename A character string of the file path.
#' @param plots A similarity matrix ggplot2 object.
#'
#' @return A character string of the file path.
save_similarity_matrix_figure <- function(filename, similarity_plot) {

  # Add legends for session and state.
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

  similarity_plot <- similarity_plot +
    # The axis ticks need to be a tiny bit wider to prevent gaps when scaling
    # the plot up in size.
    ggplot2::theme(
      axis.ticks.x.top = ggplot2::element_line(size = 1.2),
      axis.ticks.y.left = ggplot2::element_line(size = 1.1)
    ) +
    # These rasters/rectangles are purely here to make the new legends show up;
    # they aren't actually visible in the plot.
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = Session),
      data = session_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = session_legend_colours[6:1],
      guide = ggplot2::guide_legend(
        title = "Session-Recording", order = 4
      )
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = State),
      data = state_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(case_colours[9], case_colours[3]),
      guide = ggplot2::guide_legend(order = 3)
    )

  ggplot2::ggsave(
    filename = filename,
    plot = similarity_plot,
    device = ragg::agg_png,
    width = 21.59,
    height = 16,
    units = "cm",
    dpi = "retina",
    scaling = 1
  )

  filename

}

#' Save similarity archetypes figure
#'
#' @param filename A character string of the file path.
#' @param plots A list of similarity archetype plots.
#'
#' @return A character string of the file path.
save_similarity_archetypes_figure <- function(filename, plots) {

  # Hide the legends so we can only show one for the entire patchwork
  group_effect <- plots$group_effect +
    ggplot2::ggtitle("Group effect") +
    ggplot2::theme(legend.position = "none")
  session_effect <- plots$session_effect +
    ggplot2::ggtitle("Group-session effect") +
    ggplot2::theme(legend.position = "none")
  state_effect <- plots$state_effect +
    ggplot2::ggtitle("Group-state effect") +
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
    # These rasters/rectangles are purely here to make the new legends show up;
    # they aren't actually visible in the plot.
    ggnewscale::new_scale_fill() +
    ggplot2::geom_raster(
      ggplot2::aes(x = NA, y = NA, fill = 1),
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = .1),
      type = "viridis",
      guide = ggplot2::guide_legend(
        title = "Similarity Values", reverse = TRUE, order = 2
      )
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = Session),
      data = session_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = session_legend_colours[6:1],
      guide = ggplot2::guide_legend(
        title = "Session-Recording", order = 4
      )
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = State),
      data = state_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(case_colours[9], case_colours[3]),
      guide = ggplot2::guide_legend(order = 3)
    )

  plots_legend <- individual_effect |>
    ggpubr::get_legend() |>
    ggpubr::as_ggplot()

  individual_effect <- individual_effect +
    ggplot2::ggtitle("Individual effect") +
    ggplot2::theme(legend.position = "none")

  # design <- "ABCGHI
  #            DEFGHI"

  design <- "ADG#I
             BEGHI
             CFG#I"

  p1 <- patchwork::wrap_elements(group_effect)
  p2 <- patchwork::wrap_elements(session_effect)
  p3 <- patchwork::wrap_elements(state_effect)
  p4 <- patchwork::wrap_elements(individual_effect)
  p5 <- patchwork::wrap_elements(individual_session_effect)
  p6 <- patchwork::wrap_elements(individual_state_effect)
  p7 <- plots_legend

  patch <- patchwork::wrap_plots(
    A = p1, B = p2, C = p3, D = p4, E = p5, F = p6,
    G = patchwork::plot_spacer(), H = p7, I = patchwork::plot_spacer(),
    design = design,
    guides = "collect",
    widths = c(1, 1, 0.025, 0.2, 0.025)
  )

  ggplot2::ggsave(
    filename = filename,
    plot = patch,
    device = ragg::agg_png,
    width = 21.59,
    height = 26,
    units = "cm",
    dpi = "retina",
    scaling = 0.65
  )

  filename

}

#' Save similarity results figure
#'
#' @param filename A character vector of the file path.
#' @param connectivity_profile The connectivity profile plot
#' @param similarity_matrix The similarity matrix plot
#' @param group_contrasts The group-level contrast plot
#' @param subset_contrasts The individual-level contrast plot
#'
#' @return A character vector of the file path.
save_results_figure <- function(
  filename,
  connectivity_profile,
  similarity_matrix,
  group_contrasts,
  subset_contrasts = NULL
  ) {

  # FIXME: Using wrap_elements() may be unnecessary now that I'm no longer using
  # ggh4x::facet_grid2() for the similarity matrix.
  p1 <- patchwork::wrap_elements(connectivity_profile)
  p2 <- patchwork::wrap_elements(similarity_matrix)
  p3 <- patchwork::wrap_elements(group_contrasts)
  p4 <- patchwork::wrap_elements(subset_contrasts)

  # Include individual-level sub-models
  if (!is.null(subset_contrasts)) {

    design <- "AA
               BC
               DD"

    patch <- patchwork::wrap_plots(A = p1, B = p2, C = p3, D = p4, design = design) +
      patchwork::plot_annotation(tag_levels = "A")

    ggplot2::ggsave(
      filename = filename,
      plot = patch,
      device = ragg::agg_png,
      width = 21.59,
      height = 26,
      units = "cm",
      dpi = "retina",
      scaling = 0.65
    )

  # Do not include individual-level sub-models
  } else {

    design <- "AA
               BC"

    patch <- patchwork::wrap_plots(A = p1, B = p2, C = p3, design = design) +
      patchwork::plot_annotation(tag_levels = "A")

    ggplot2::ggsave(
      filename = filename,
      plot = patch,
      device = ragg::agg_png,
      width = 21.59,
      height = 18,
      units = "cm",
      dpi = "retina",
      scaling = 0.65
    )

  }

  filename

}

save_subset_similarity_contrast_barplot_figure <- function(
  filename,
  phase,
  amplitude
) {

  design <- ("AAAAAAAAAA
              BB########")

  patch <- phase +
    amplitude +
    patchwork::plot_layout(design = design, guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")

  ggplot2::ggsave(
    filename = filename,
    plot = patch,
    device = ragg::agg_png,
    width = 21.59,
    height = 18,
    units = "cm",
    dpi = "retina",
    scaling = 1
  )

  filename

}

#' Save amplitude similarity patchwork
#'
#' @param delta,theta,beta,gamma Similarity matrix for the respective frequency band.
#'
#' @return A character vector of the file path.
save_amplitude_similarity_patchwork_figure <- function(
  delta, theta, beta, gamma
) {

  filename <- "figures/amplitude-similarity-patchwork.png"

  delta <- delta +
    ggplot2::ggtitle("Delta") +
    ggplot2::theme(legend.position = "none")
  theta <- theta +
    ggplot2::ggtitle("Theta") +
    ggplot2::theme(legend.position = "none")
  beta <- beta +
    ggplot2::ggtitle("Beta") +
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

  gamma <- gamma +
    # These rasters/rectangles are purely here to make the new legends show up;
    # they aren't actually visible in the plot.
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = Session),
      data = session_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = session_legend_colours[6:1],
      guide = ggplot2::guide_legend(
        title = "Session-Recording", order = 4
      )
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = State),
      data = state_legend, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(case_colours[9], case_colours[3]),
      guide = ggplot2::guide_legend(order = 3)
    )

  plots_legend <- gamma |>
    ggpubr::get_legend() |>
    ggpubr::as_ggplot()

  gamma <- gamma +
    ggplot2::ggtitle("Gamma") +
    ggplot2::theme(legend.position = "none")

  design <- "ABEFG
             CDEFG"

  p1 <- patchwork::wrap_elements(delta)
  p2 <- patchwork::wrap_elements(theta)
  p3 <- patchwork::wrap_elements(beta)
  p4 <- patchwork::wrap_elements(gamma)

  p7 <- plots_legend

  patch <- patchwork::wrap_plots(
    A = p1, B = p2, C = p3, D = p4,
    E = patchwork::plot_spacer(), F = p7, G = patchwork::plot_spacer(),
    design = design,
    guides = "collect",
    widths = c(1, 1, 0.025, 0.2, 0.025)
  )

  ggplot2::ggsave(
    filename = filename,
    plot = patch,
    device = ragg::agg_png,
    width = 21.59,
    height = 19,
    units = "cm",
    dpi = "retina",
    scaling = 0.65
  )

  filename

}

#' Save similarity matrix patchwork with jet colour scale
#'
#' @param all_arguments Similarity matrix plot for the respective frequency
#'   band.
#'
#' @return A character vector of the file path.
save_jet_similarity_patchwork_figure <- function(
  phase_delta,
  phase_theta,
  phase_alpha,
  phase_beta,
  phase_gamma,
  amplitude_alpha
) {

  filename <- "figures/jet-similarity-patchwork.png"

  jet <- palettes::as_colour(readRDS(here::here("R/jet-colours.rds")))

  phase_delta <- phase_delta +
    palettes::scale_fill_palette_c(
      jet,
      limits = c(0, NA),
      # This is a bit of a hacky way to get the limits of the fill-scale to
      # ignore the 1s on the diagonal of the matrix. They are set to 0 instead
      # of NA because the latter throws a cryptic error that I couldn't be
      # bothered with.
      trans = scales::trans_new(
        name = "remove_ones",
        transform = \(.x) ifelse(.x == 1, 0, .x),
        inverse = \(.x) ifelse(is.na(.x), 1, .x)
      )
    ) +
    ggplot2::ggtitle("Phase, Delta")

  phase_theta <- phase_theta +
    palettes::scale_fill_palette_c(
      jet,
      limits = c(0, NA),
      trans = scales::trans_new(
        name = "remove_ones",
        transform = \(.x) ifelse(.x == 1, 0, .x),
        inverse = \(.x) ifelse(is.na(.x), 1, .x)
      )
    ) +
    ggplot2::ggtitle("Phase, Theta")

  phase_alpha <- phase_alpha +
    palettes::scale_fill_palette_c(
      jet,
      limits = c(0, NA),
      trans = scales::trans_new(
        name = "remove_ones",
        transform = \(.x) ifelse(.x == 1, 0, .x),
        inverse = \(.x) ifelse(is.na(.x), 1, .x)
      )
    ) +
    ggplot2::ggtitle("Phase, Alpha")

  phase_beta <- phase_beta +
    palettes::scale_fill_palette_c(
      jet,
      limits = c(0, NA),
      trans = scales::trans_new(
        name = "remove_ones",
        transform = \(.x) ifelse(.x == 1, 0, .x),
        inverse = \(.x) ifelse(is.na(.x), 1, .x)
      )
    ) +
    ggplot2::ggtitle("Phase, Beta")

  phase_gamma <- phase_gamma +
    palettes::scale_fill_palette_c(
      jet,
      limits = c(0, NA),
      trans = scales::trans_new(
        name = "remove_ones",
        transform = \(.x) ifelse(.x == 1, 0, .x),
        inverse = \(.x) ifelse(is.na(.x), 1, .x)
      )
    ) +
    ggplot2::ggtitle("Phase, Gamma")

  amplitude_alpha <- amplitude_alpha +
    palettes::scale_fill_palette_c(
      jet,
      limits = c(0, NA),
      trans = scales::trans_new(
        name = "remove_ones",
        transform = \(.x) ifelse(.x == 1, 0, .x),
        inverse = \(.x) ifelse(is.na(.x), 1, .x)
      )
    ) +
    ggplot2::ggtitle("Amplitude, Alpha")

  design <- "AB
             CD
             EF"

  patch <- patchwork::wrap_plots(
    A = phase_delta,
    B = phase_theta,
    C = phase_alpha,
    D = phase_beta,
    E = phase_gamma,
    F = amplitude_alpha,
    design = design
  )

  ggplot2::ggsave(
    filename = filename,
    plot = patch,
    device = ragg::agg_png,
    width = 21.59,
    height = 26,
    units = "cm",
    dpi = "retina",
    scaling = 0.65
  )

  filename

}

#' Make contrast results table with significance tests
#'
#' @param emmeans_tidy The tidied emmeans data.
#' @param contrasts_tidy The tidied contrasts data.
#'
#' @return A flextable object.
make_contrast_results_table_nhst <- function(emmeans_tidy, contrasts_tidy) {

  # Prepare data for display
  tidy_data <- emmeans_tidy |>
    dplyr::select(
      effect_label = effect_label_nice, within_participant, response, std.error
    ) |>
    tidyr::pivot_wider(
      names_from = within_participant,
      values_from = c(response, std.error)
    ) |>
    dplyr::relocate(
      effect_label,
      response_within_participants,
      std.error_within_participants,
      response_between_participants,
      std.error_between_participants
    ) |>
    dplyr::left_join(
      dplyr::select(contrasts_tidy, effect_label, estimate:p.value)
    ) |>
    dplyr::relocate(df, .after = statistic) |>
    dplyr::mutate(
      dplyr::across(
        c(where(is.numeric), -df, -p.value),
        function(x) papaja::apa_num(x, digits = 3)
      ),
      dplyr::across(
        dplyr::starts_with("std.error"),
        function(x) paste0("(", x, ")")
      ),
      df = papaja::print_df(df),
      p.value = papaja::print_p(p.value)
    ) |>
    tidyr::unite(
      "mean_se_within_participant",
      dplyr::ends_with("within_participants"),
      sep = " "
    ) |>
    tidyr::unite(
      "mean_se_between_participant",
      dplyr::ends_with("between_participants"),
      sep = " "
    ) |>
    tidyr::unite(
      "mean_difference",
      c(estimate, std.error),
      sep = " "
    ) |>
    tidyr::unite(
      "ci",
      c(conf.low, conf.high),
      sep = ", "
    ) %>%
    dplyr::rename(
      contrast = effect_label,
      "*t*({unique(.$df)})" := statistic,
      "*p*" = p.value
    ) |>
    dplyr::select(-df) |>
    dplyr::mutate(
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
    ) |>
    dplyr::arrange(contrast)

  # Make results table
  tidy_data |>
    flextable::flextable() |>
    # Header
    flextable::set_header_labels(
      contrast = "",
      mean_se_within_participant = "Mean (SE)",
      mean_se_between_participant = "Mean (SE)",
      mean_difference = "ΔRV (SE)",
      ci = "95% CI"
    ) |>
    flextable::add_header_row(
      values = c(
        "Contrast",
        "Within Participant",
        "Between Participant",
        "Mean Difference",
        ""
      ),
      colwidths = c(1, 1, 1, 2, 2)
    ) |>
    ftExtra::colformat_md(part = "header") |>
    # Footer
    flextable::add_footer_lines(
      "SE = Standard Error, CI = Compatibility Interval"
    ) |>
    flextable::footnote(
      i = 1, j = 1,
      value = flextable::as_paragraph(
        "Results are averaged over the levels of session and state."
      ),
      ref_symbols = "1",
      part = "body",
    ) |>
    flextable::footnote(
      i = 2:3, j = 1,
      value = flextable::as_paragraph("Results are averaged over the levels of state."),
      ref_symbols = "2",
      part = "body",
    ) |>
    flextable::footnote(
      i = 4:5, j = 1,
      value = flextable::as_paragraph(
        "Results are averaged over the levels of session."
      ),
      ref_symbols = "3",
      part = "body",
    ) |>
    # Cell alignment
    flextable::align(part = "header", j = 1, align = "left") |>
    flextable::align(part = "body", j = 1, align = "left") |>
    flextable::align(part = "header", j = 2:7, align = "center") |>
    flextable::align(part = "body", j = 2:7, align = "right") |>
    # Borders
    flextable::hline(
      i = 1, j = c(1, 6, 7),
      border = flextable::fp_border_default(width = 0),
      part = "header"
    ) |>
    flextable::hline_bottom(
      border = flextable::fp_border_default(width = 2),
      part = "footer"
    ) |>
    flextable::fix_border_issues(part = "footer") |>
    # Font
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    # Dimensions
    flextable::width(j = 1, width = 3.59, unit = "cm")

}

#' Make contrast results table with standardized mean differences
#'
#' @param model A glmmTMB object.
#' @param emmeans_tidy The tidied emmeans data.
#' @param contrasts_tidy The tidied contrasts data.
#' @param contrasts_cohens_d_tidy The tidied Cohen's d contrasts data.
#'
#' @return A flextable object.
make_contrast_results_table_smd <- function(
  model, emmeans_tidy, contrasts_tidy, contrasts_cohens_d_tidy
) {

  # Prepare data for display
  tidy_data <- emmeans_tidy |>
    dplyr::select(
      effect_label = effect_label_nice, within_participant, response, std.error
    ) |>
    tidyr::pivot_wider(
      names_from = within_participant,
      values_from = c(response, std.error)
    ) |>
    dplyr::relocate(
      effect_label,
      response_within_participants,
      std.error_within_participants,
      response_between_participants,
      std.error_between_participants
    ) |>
    dplyr::left_join(
      dplyr::select(contrasts_tidy, effect_label, estimate:conf.high)
    ) |>
    dplyr::left_join(
      dplyr::select(contrasts_cohens_d_tidy, effect_label, smd:smd_conf.high)
    ) |>
    dplyr::mutate(
      dplyr::across(
        c(where(is.numeric), -df),
        function(x) papaja::apa_num(x, digits = 3)
      ),
      dplyr::across(
        dplyr::contains("std.error"),
        function(x) paste0("(", x, ")")
      ),
      df = papaja::print_df(df)
    ) |>
    tidyr::unite(
      "mean_se_within_participant",
      dplyr::ends_with("within_participants"),
      sep = " "
    ) |>
    tidyr::unite(
      "mean_se_between_participant",
      dplyr::ends_with("between_participants"),
      sep = " "
    ) |>
    tidyr::unite(
      "mean_difference",
      c(estimate, std.error),
      sep = " "
    ) |>
    tidyr::unite(
      "ci",
      c(conf.low, conf.high),
      sep = ", "
    ) |>
    tidyr::unite(
      "standardized_mean_difference",
      c(smd, smd_std.error),
      sep = " "
    ) |>
    tidyr::unite(
      "smd_ci",
      c(smd_conf.low, smd_conf.high),
      sep = ", "
    ) |>
    dplyr::rename(
      contrast = effect_label
    ) |>
    dplyr::select(-df) |>
    dplyr::mutate(
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
    ) |>
    dplyr::arrange(contrast)

  # Make results table
  tidy_data |>
    flextable::flextable() |>
    # Header
    flextable::set_header_labels(
      contrast = "",
      mean_se_within_participant = "Mean (SE)",
      mean_se_between_participant = "Mean (SE)",
      mean_difference = "ΔRV (SE)",
      ci = "95% CI",
      standardized_mean_difference = "d",
      smd_ci = "95% CI"
    ) |>
    flextable::compose(
      part = "header", i = 1, j = 6,
      flextable::as_paragraph(
        flextable::as_i("d"),
        flextable::as_sub(flextable::as_i("T")),
        " (SE)"
      )
    ) |>
    flextable::add_header_row(
      values = c(
        "Contrast",
        "Within Participant",
        "Between Participant",
        "Mean Difference",
        "Standardized Mean Difference"
      ),
      colwidths = c(1, 1, 1, 2, 2)
    ) |>
    #ftExtra::colformat_md(part = "header") |>
    # Footer
    flextable::add_footer_lines(
      "SE = Standard Error, CI = Compatibility Interval"
    ) |>
    flextable::footnote(
      i = 1, j = 1,
      value = flextable::as_paragraph(
        "Results are averaged over the levels of session and state."
      ),
      ref_symbols = "1",
      part = "body",
    ) |>
    flextable::footnote(
      i = 2:3, j = 1,
      value = flextable::as_paragraph("Results are averaged over the levels of state."),
      ref_symbols = "2",
      part = "body",
    ) |>
    flextable::footnote(
      i = 4:5, j = 1,
      value = flextable::as_paragraph(
        "Results are averaged over the levels of session."
      ),
      ref_symbols = "3",
      part = "body",
    ) |>
    flextable::footnote(
      i = 2, j = 6,
      value = flextable::as_paragraph(
        "Mean differences were standardized using the summed standard ",
        "deviation of the variance components from the fitted model (",
        flextable::as_equation("\\sigma_T"),
        " = ",
        papaja::apa_num(
          sqrt(sum(insight::get_variance_intercept(model))), digits = 3
        ),
        ")."
      ),
      ref_symbols = "4",
      part = "header",
    ) |>
    # Cell alignment
    flextable::align(part = "header", j = 1, align = "left") |>
    flextable::align(part = "body", j = 1, align = "left") |>
    flextable::align(part = "header", j = 2:7, align = "center") |>
    flextable::align(part = "body", j = 2:7, align = "right") |>
    # Borders
    flextable::hline(
      i = 1, j = 1,
      border = flextable::fp_border_default(width = 0),
      part = "header"
    ) |>
    flextable::hline_bottom(
      border = flextable::fp_border_default(width = 2),
      part = "footer"
    ) |>
    flextable::fix_border_issues(part = "footer") |>
    # Font
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    # Dimensions
    flextable::width(j = 1, width = 3.59, unit = "cm")

}

#' Make model fit table
#'
#' @param model A glmmTMB object.
#'
#' @return A flextable object.
make_model_fit_table <- function(model, maximal_model = FALSE) {

  model |>
    gtsummary::tbl_regression(
      label = if (maximal_model) {
        list(
          within_participant ~ "Within Participant",
          within_session ~ "Within Session",
          within_state ~ "Within State",
          "pair_participant.sd__(Intercept)" ~ "Participant Pair",
          "x_label.sd__(Intercept)" ~ "X Connectome",
          "y_label.sd__(Intercept)" ~ "Y Connectome"
        )
      } else {
        list(
          within_participant ~ "Within Participant",
          within_session ~ "Within Session",
          within_state ~ "Within State",
          "x_label.sd__(Intercept)" ~ "X Connectome",
          "y_label.sd__(Intercept)" ~ "Y Connectome"
        )
      },
      show_single_row = c(
        within_participant, within_state, within_session,
        "within_participant:within_session",
        "within_participant:within_state",
        "within_session:within_state",
        "within_participant:within_session:within_state"
      ),
      intercept = TRUE,
      estimate_fun = ~ papaja::print_num(.x, digits = 3, na_string = ""),
      pvalue_fun = papaja::print_p,
      tidy_fun = broom.mixed::tidy
    ) |>
    gtsummary::modify_footnote(update = everything() ~ NA, abbreviation = TRUE) |>
    gtsummary::modify_table_body(
      ~ .x |>
        dplyr::mutate(
          Type = c(
            rep("Fixed Effects", times = 8),
            rep("Random Effects", times = ifelse(maximal_model, 3, 2))
          ),
          .before = label
        ) |>
        flextable::as_grouped_data(groups = "Type") |>
        tibble::as_tibble()
    ) |>
    gtsummary::modify_column_unhide(c(Type, std.error)) |>
    gtsummary::add_glance_source_note(
      label = list(
        sigma ~ "Phi coefficient"
      ),
      glance_fun = broom.mixed::glance
    ) |>
    gtsummary::as_flex_table() |>
    # Header and grouped rows
    flextable::set_header_labels(
      label = "Term",
      estimate = "Estimate",
      p.value = "*p*"
    ) |>
    ftExtra::colformat_md(part = "header") |>
    flextable::merge_at(i = 1, part = "body") |>
    flextable::align(i = 1, align = "left") |>
    flextable::merge_at(i = 10, part = "body") |>
    flextable::align(i = 10, align = "left") |>
    # Footnotes
    flextable::footnote(
      j = 3,
      value = flextable::as_paragraph(
        "Estimates for fixed effects and random effects are, respectively, ",
        "coefficients for the mean model (with logit link) and standard ",
        "deviations for the grouping variables (random intercepts)."
      ),
      ref_symbols = "1",
      part = "header",
    ) |>
    # Cell alignment
    flextable::align(part = "header", i = 1, align = "center") |>
    flextable::align(part = "body", j = 3:6, align = "right") |>
    # Borders
    flextable::hline_bottom(
      border = flextable::fp_border_default(width = 1),
      part = "footer"
    ) |>
    flextable::fix_border_issues(part = "footer") |>
    # Font
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    # Dimensions
    flextable::width(j = 1, width = 1.4, unit = "cm") |>
    flextable::width(j = 2, width = 8.5, unit = "cm") |>
    flextable::width(j = 3, width = 1.8, unit = "cm") |>
    flextable::width(j = 4:6, width = 1.6, unit = "cm")

}

#' Make conervgence/singularity table
#'
#' @param all_arguments A glmmTMB object corresponding to the respective
#'   coupling mode, frequency band, and random effects specification.
#'
#' @return A flextable object.
make_convergence_table <- function(
  phase_delta_reduced,
  phase_delta_maximal,
  phase_theta_reduced,
  phase_theta_maximal,
  phase_alpha_reduced,
  phase_alpha_maximal,
  phase_beta_reduced,
  phase_beta_maximal,
  phase_gamma_reduced,
  phase_gamma_maximal,
  amplitude_alpha_reduced,
  amplitude_alpha_maximal,
  phase_delta_hilbert_reduced,
  phase_delta_hilbert_maximal,
  phase_theta_hilbert_reduced,
  phase_theta_hilbert_maximal,
  phase_alpha_hilbert_reduced,
  phase_alpha_hilbert_maximal,
  phase_beta_hilbert_reduced,
  phase_beta_hilbert_maximal,
  phase_gamma_hilbert_reduced,
  phase_gamma_hilbert_maximal
) {

  models <- tibble::tribble(
    ~model,                            ~coupling_mode,              ~freq_band, ~re,
    list(phase_delta_reduced),         "Phase",                     "Delta", "Reduced",
    list(phase_delta_maximal),         "Phase",                     "Delta", "Maximal",
    list(phase_theta_reduced),         "Phase",                     "Theta", "Reduced",
    list(phase_theta_maximal),         "Phase",                     "Theta", "Maximal",
    list(phase_alpha_reduced),         "Phase",                     "Alpha", "Reduced",
    list(phase_alpha_maximal),         "Phase",                     "Alpha", "Maximal",
    list(phase_beta_reduced),          "Phase",                     "Beta",  "Reduced",
    list(phase_beta_maximal),          "Phase",                     "Beta",  "Maximal",
    list(phase_gamma_reduced),         "Phase",                     "Gamma", "Reduced",
    list(phase_gamma_maximal),         "Phase",                     "Gamma", "Maximal",
    list(amplitude_alpha_reduced),     "Amplitude",                 "Alpha", "Reduced",
    list(amplitude_alpha_maximal),     "Amplitude",                 "Alpha", "Maximal",
    list(phase_delta_hilbert_reduced), "Phase (Hilbert transform)", "Delta", "Reduced",
    list(phase_delta_hilbert_maximal), "Phase (Hilbert transform)", "Delta", "Maximal",
    list(phase_theta_hilbert_reduced), "Phase (Hilbert transform)", "Theta", "Reduced",
    list(phase_theta_hilbert_maximal), "Phase (Hilbert transform)", "Theta", "Maximal",
    list(phase_alpha_hilbert_reduced), "Phase (Hilbert transform)", "Alpha", "Reduced",
    list(phase_alpha_hilbert_maximal), "Phase (Hilbert transform)", "Alpha", "Maximal",
    list(phase_beta_hilbert_reduced),  "Phase (Hilbert transform)", "Beta",  "Reduced",
    list(phase_beta_hilbert_maximal),  "Phase (Hilbert transform)", "Beta",  "Maximal",
    list(phase_gamma_hilbert_reduced), "Phase (Hilbert transform)", "Gamma", "Reduced",
    list(phase_gamma_hilbert_maximal), "Phase (Hilbert transform)", "Gamma", "Maximal"
  )

  # Check models
  models |>
    dplyr::rowwise() |>
    dplyr::mutate(
      converged = performance::check_convergence(model[[1]]),
      singular  = performance::check_singularity(model[[1]])
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-model) |>
    # Make table
    tidyr::pivot_wider(
      names_from = re, values_from = c(converged, singular)
    ) |>
    dplyr::relocate(singular_Reduced, .after = converged_Reduced) |>
    flextable::as_grouped_data(groups = "coupling_mode") |>
    flextable::flextable() |>
    # Header and grouped rows
    flextable::set_header_labels(
      coupling_mode = "",
      freq_band = "",
      converged_Reduced = "Converged",
      singular_Reduced = "Singular",
      converged_Maximal = "Converged",
      singular_Maximal = "Singular"
    ) |>
    flextable::add_header_row(
      values = c(
        "Coupling Mode",
        "Frequency Band",
        "Reduced Model",
        "Maximal Model"
      ),
      colwidths = c(1, 1, 2, 2)
    ) |>
    flextable::merge_at(i = 1, part = "body") |>
    flextable::align(i = 1, align = "left") |>
    flextable::merge_at(i = 7, part = "body") |>
    flextable::align(i = 7, align = "left") |>
    flextable::merge_at(i = 9, part = "body") |>
    flextable::align(i = 9, align = "left") |>
    # Cell alignment
    flextable::align(part = "header", i = c(1, 2), align = "center") |>
    flextable::align(part = "body", j = 3:6, align = "center") |>
    # Borders
    flextable::hline(
      i = 1, j = c(1, 2),
      border = flextable::fp_border_default(width = 0),
      part = "header"
    ) |>
    # Font
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = 10, part = "all")

}

#' Save a figure
#'
#' @param filename A character vector of the file path.
#' @param plot A ggplot2 object.
#' @param width,height Numeric. Dimensions in cm.
#' @param scaling Numeric. Plot scaling.
#'
#' @return A character vector of the file path.
save_figure <- function(
  filename,
  plot,
  width = 21.59,
  height,
  scaling = 0.65
) {

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = "cm",
    dpi = "retina",
    scaling = scaling
  )

  filename

}
