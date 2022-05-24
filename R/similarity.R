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

  # Generate all unique pairwise comparisons over branch names and label names
  branch_names <- names(input)
  label_names  <- input |>
    purrr::map(~ .x$metadata$case) |>
    unlist()

  # Make sure to allow repeats for self-comparison pairs
  branch_combinations <- gtools::combinations(
    length(branch_names), 2, branch_names, repeats.allowed = TRUE
  )
  label_combinations <- gtools::combinations(
    length(label_names), 2, label_names, repeats.allowed = TRUE
  )

  # Get levels for labels
  label_levels <- input |>
    purrr::map(~ .x$metadata$participant) |>
    unlist() |>
    unique() |>
    case_order()

  # Add metadata to assist with analyses and filtering
  connectivity_matrix_pairs <- tibble::tibble(
    x_branch = branch_combinations[,1],
    y_branch = branch_combinations[,2],
    x_label = factor(label_combinations[,1], levels = label_levels),
    y_label = factor(label_combinations[,2], levels = label_levels),
    x_participant = stringr::str_extract(x_label, "^.*?(?=_)"),
    x_session = stringr::str_extract(x_label, "(?<=_)(.*)(?=_)"),
    x_state = stringr::str_extract(x_label, "([^_]*)$"),
    y_participant = stringr::str_extract(y_label, "^.*?(?=_)"),
    y_session = stringr::str_extract(y_label, "(?<=_)(.*)(?=_)"),
    y_state = stringr::str_extract(y_label, "([^_]*)$"),
    within_participant = ifelse(x_participant == y_participant, TRUE, FALSE),
    within_session = ifelse(x_session == y_session, TRUE, FALSE),
    within_state = ifelse(
      stringr::str_detect(x_state, stringr::str_sub(y_state, 1, 2)),
      TRUE,
      FALSE
    ),
    diagonal = ifelse(x_label == y_label, TRUE, FALSE),
    pair_label = paste(x_label, y_label, sep = "_x_")
  )

  connectivity_matrix_pairs
}

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
  similarity_results_symmetric <- rbind(similarity_results, similarity_results_upper) |>
    # Ensure levels are in order for plotting
    dplyr::mutate(dplyr::across(where(is.character), forcats::as_factor))

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

# TODO: Add target that maps over participants
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
