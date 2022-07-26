#' Create an ordered vector of participant IDs
#'
#' @param participant_ids
#'
#' @return
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
#' @param similarity_results
#'
#' @return
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
#' @param similarity_results The tibble of similarity results.
#' @param estimate The similarity estimate to plot. Options are `rv`.
#'
#' @return ggplot
plot_similarity <- function(similarity_results, estimate) {

  # Combine the upper, lower, and diagonals to prepare for plotting
  similarity_results_symmetric <- make_symmetric(similarity_results)

  # Plot
  similarity_results_symmetric |>
    ggplot2::ggplot(ggplot2::aes(x = x_label, y = y_label, fill = {{estimate}})) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_continuous(limits = c(0, 1), type = "viridis") +
    # Expand reduces the spacing between facets
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    # Reverse y-axis so the diagonal goes from the upper-left to lower-right
    ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "RV"
    ) +
    ggh4x::facet_grid2(
      cols = ggplot2::vars(x_participant),
      rows = ggplot2::vars(y_participant),
      scales = "free",
      independent = TRUE,
      switch = "y"
    ) +
    ggplot2::theme(
      panel.spacing = grid::unit(0, "lines"),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      strip.text.y.left = ggplot2::element_text(angle = 0),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

}

#' Plot the within-participant similarity block of a participant
#'
#' @param similarity_results
#' @param participant
#'
#' @return
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
#' @param similarity_results
#' @param participant
#'
#' @return
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
#' @param similarity_results
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
#' @param formula
#' @param data
#'
#' @note Since each of the predictors are categorical, it doesn't make sense to
#  include random slopes in this model.
#'
#' @return glmmTMB
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
#' @return ggplot.
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
#' @param object
#'
#' @return
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

#' Title
#'
#' @param object
#'
#' @return
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

  object |>
    purrr::map_dfr(broom::tidy, .id = "effect") |>
    dplyr::mutate(
      effect_label = rep(effect_labels, each = 2),
      within_participant = dplyr::case_when(
        within_participant == TRUE  ~ "within_participants",
        within_participant == FALSE ~ "between_participants"
      ),
      effect_label = paste0(within_participant, "_", effect_label)
    )

}

#' Estimate pairwise contrasts using the reference grid
#'
#' @param object
#'
#' @return
contrast_similarity <- function(object) {

  emmeans::contrast(
    object,
    method = "revpairwise", # Subtracts the between value from the within value
    infer = TRUE
  )

}

#' Title
#'
#' @param object
#'
#' @return A tibble.
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
      dist = distributional::dist_student_t(
        df = df, mu = estimate, sigma = std.error
      )
    ) |>
    ggdist::point_interval(dist, .width = c(.66, .95)) |>
    tidyr::pivot_wider(names_from = .width, values_from = c(.lower, .upper))

  contrasts

}

#' Plot similarity contrast intervals
#'
#' @param object
#'
#' @return
plot_similarity_contrasts <- function(object) {

  # Tidy data to prepare for plotting
  contrasts <- tidy_contrast_similarity(object)

  # Plotting
  interval_pal <- RColorBrewer::brewer.pal(n = 3, "YlGn")[2:3]

  ggplot2::ggplot(contrasts, ggplot2::aes(x = estimate, y = effect_label)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linetype = 2) +
    ggdist::stat_interval(
      ggplot2::aes(
        xdist = distributional::dist_student_t(
          df = df, mu = estimate, sigma = std.error
        )
      ),
      .width = c(0.66, 0.95)
    ) +
    ggplot2::scale_colour_manual(values = interval_pal) +
    # This is a slightly hacky way to get a bar the same height as the interval
    # for the point estimate since the standard ggplot2 geoms aren't well-suited
    # for this.
    ggdist::geom_interval(
      ggplot2::aes(xmin = estimate - 0.001, xmax = estimate + 0.001)
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(-0.1, 0, 0.1),
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
      colour = "Compatibility\nInterval"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.line.x.top  = ggplot2::element_blank(),
      axis.text.x.top  = ggplot2::element_text(size = 11, colour = "black")
    )

}

# TODO: Add documentation
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

subset_glmmTMB_similarity <- function(formula, data) {

  purrr::map(
    data,
    ~{
      glmmTMB_similarity(formula, data = .x)
    }
  )

}

subset_emmeans_similarity <- function(objects) {

  purrr::map(
    objects,
    ~{
      emmeans_similarity(.x)
    }
  )

}

subset_contrast_similarity <- function(objects) {

  purrr::map(
    objects,
    ~{
      contrast_similarity(.x)
    }
  )

}

tidy_subset_contrast_similarity <- function(objects) {

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

  purrr::map_dfr(
    objects,
    ~{
      contrasts <- .x |>
      purrr::map_dfr(broom::tidy, .id = "effect") |>
      dplyr::mutate(
        effect_label = factor(effect_labels, levels = effect_labels)
      )
      contrasts
    },
    .id = "participant"
  )

}

#' Title
#'
#' @param group_object
#' @param subset_objects
#'
#' @return A ggplot
plot_subset_similarity_contrasts <- function(group_object, subset_objects) {

  # Tidy data to prepare for plotting
  contrasts <- tidy_contrast_similarity(group_object)
  subset_contrasts <- tidy_subset_contrast_similarity(subset_objects)

  subset_contrasts <- contrasts |>
    dplyr::select(effect_label, group_estimate = estimate) |>
    dplyr::left_join(subset_contrasts) |>
    # Since we have 9 contrasts, there's an odd number which makes lining up
    # facets awkward. We'll create a blank factor level, whose facet we will
    # later remove, to work around this.
    dplyr::mutate(
      effect_label = forcats::fct_expand(effect_label, " "),
      effect_label = forcats::fct_relevel(effect_label, "Main effect", " ")
    )

  # Plotting
  ramp_colours <- viridis::plasma(9)
  far_colour_index <- 2
  close_colour_index <- length(ramp_colours) - 1

  p <- ggplot2::ggplot(subset_contrasts) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linetype = 5, alpha = 0.2) +
    # Point estimate and 95% interval for the group-level effect
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .lower_0.95,
        xmax = .upper_0.95,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = 0.2,
      data = contrasts
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .lower_0.66,
        xmax = .upper_0.66,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = 0.2,
      data = contrasts
    ) +
    # This is a slightly hacky way to get a consistent width for the point
    # estimate since geom_vline() varies with graphics device size.
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = estimate - 0.001,
        xmax = estimate + 0.001,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = 0.2,
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
        fill_ramp = stat(abs(group_estimate - x))
      ),
      color = "black",
      size = 0.25,
      height = 2,
      expand = TRUE, # I can't disable this for whatever reason
      trim = FALSE, # I can't disable this for whatever reason
      fill = ramp_colours[far_colour_index],
      fill_type = "gradient"
    ) +
    ggdist::scale_fill_ramp_continuous(
      from = ramp_colours[close_colour_index],
      name = "Distance from\ngroup-level\npoint estimate",
      trans = "identity",
      # This adjusts how quickly the fill ramps. It should be set to a meaningful
      # value such as the maximum distance we consider compatible with the
      # group-level estimate.
      limits = c(0, 0.05),
      # Fills outside the distance limit should be given the max colour
      oob = scales::oob_squish,
      guide = ggdist::guide_rampbar(
        to = ramp_colours[far_colour_index],
        title.vjust = 1.5
      )
    ) +
    # Finally, some subtle lines that can be seen when the ridges overlap the
    # interval.
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = .lower_0.95),
      linetype = 5,
      alpha = 0.2,
      data = contrasts
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = .upper_0.95),
      linetype = 5,
      alpha = 0.2,
      data = contrasts
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = estimate),
      linetype = 2,
      alpha = 0.2,
      data = contrasts
    ) +
    # General plot settings
    ggplot2::scale_x_continuous(
      breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)
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
    ggplot2::theme_classic()
  p

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

#' Title
#'
#' @param filename
#' @param similarity_matrix
#' @param group_contrasts
#' @param subset_contrasts
#'
#' @return A character vector
save_results_figure <- function(
  filename,
  connectivity_profile,
  similarity_matrix,
  group_contrasts,
  subset_contrasts
  ) {

  design <- "AA
             BC
             DD"

  p1 <- patchwork::wrap_elements(connectivity_profile)
  p2 <- patchwork::wrap_elements(similarity_matrix)
  p3 <- patchwork::wrap_elements(group_contrasts)
  p4 <- patchwork::wrap_elements(subset_contrasts)

  patch <- patchwork::wrap_plots(A = p1, B = p2, C = p3, D = p4, design = design) +
    patchwork::plot_annotation(tag_levels = "A")

  ggplot2::ggsave(
    filename = filename,
    plot = patch,
    device = ragg::agg_png,
    width = 21.59,
    height = 25,
    units = "cm",
    dpi = "retina",
    scaling = 0.65
  )

  filename

}
