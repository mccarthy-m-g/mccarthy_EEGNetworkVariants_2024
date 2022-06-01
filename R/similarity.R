#' Helper
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

#' Generate pairwise comparisons
#'
#' @param input List of connectivity matrices and their metadata.
#'
#' @return A tibble.
pairwise_comparisons <- function(input) {

  # Get the branch and label names that go together so pair names are human
  # readable
  label_names  <- input |>
    purrr::map(~ .x$metadata$case) |>
    unlist()
  branch_names <- names(label_names)

  connectivity_matrix_names <- tibble::tibble(
    label_names = label_names,
    branch_names = branch_names
  )

  # Generate all unique pairwise comparisons over branch names. Make sure to
  # allow repeats for self-comparison pairs.
  branch_combinations <- gtools::combinations(
    length(branch_names), 2, branch_names, repeats.allowed = TRUE
  )
  colnames(branch_combinations) <- c("x_branch", "y_branch")
  branch_combinations <- tibble::as_tibble(branch_combinations)

  # Get levels for factors
  participant_levels <- label_levels <- input |>
    purrr::map(~ .x$metadata$participant) |>
    unlist() |>
    unique() |>
    sort()
  session_levels <- c("pre", "post", "fu")
  state_levels <- c("rc1", "rc2", "ro1", "ro2")
  label_levels <- case_order(participant_levels)

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
      # Make labels a factor to assist with plotting order
      x_label = factor(x_label, levels = label_levels),
      y_label = factor(y_label, levels = label_levels),
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

#' Helper
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

#' Plot similarity matrix
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
      fill = "RV Coefficient"
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

#' Title
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

#' Title
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

#' Title
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

#' Title
#'
#' @param formula
#' @param data
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

#' Title
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
#' @return
plot_similarity_contrasts <- function(object) {

  # Tidy data to prepare for plotting
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
    purrr::map_dfr(broom::tidy, .id = "effect") |>
      dplyr::mutate(
        effect_label = factor(effect_labels, levels = effect_labels)
      )

  # Plotting
  interval_pal <- RColorBrewer::brewer.pal(n = 3, "YlGn")[2:3]

  ggplot2::ggplot(contrasts, ggplot2::aes(x = estimate, y = effect_label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = estimate, yend = effect_label),
      alpha = 0.05
    ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linetype = 2) +
    ggdist::stat_interval(
      ggplot2::aes(
        xdist = distributional::dist_student_t(
          df = df.residual(similarity_fit), mu = estimate, sigma = std.error
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
        "(Within participant - Between participant)"
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
