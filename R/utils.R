#' Set caption
#'
#' @param x A flextable object.
#' @param caption Character string.
#'
#' @return A flextable object.
set_table_caption_apa <- function(x, caption) {
  flextable::set_caption(
    x,
    caption = flextable::as_paragraph(
      flextable::as_chunk(
        paste0("\n\n", caption, "\n"),
        props = flextable::fp_text_default(
          italic = TRUE, font.family = "Times New Roman", font.size = 12
        )
      )
    ),
    word_stylename = "Table Caption",
    align_with_table = FALSE
  )
}
