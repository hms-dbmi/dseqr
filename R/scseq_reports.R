#' Remove ggplot xaxis title, text, and ticks
#'
#' @return \code{theme}
#'
#' @keywords internal
theme_no_xaxis <- function() {
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank())
}

#' Remove ggplot yaxis title, text, and ticks
#'
#' @return \code{theme}
#'
#' @keywords internal
theme_no_yaxis <- function() {
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())
}

#' Make ggplot axes and text dimgray
#'
#' @param with_nums Include axis ticks/text? Default is TRUE.
#'
#' @return \code{theme}
#'
#' @keywords internal
theme_dimgray <- function(with_nums = TRUE) {

  axis.line <- ggplot2::element_line(size = 0.1, color = 'dimgray')

  if (with_nums) {
    axis.text <- ggplot2::element_text(colour = 'dimgray')
    axis.ticks <- ggplot2::element_line(size = 0.1, color = 'dimgray')

  } else {
    axis.text <- axis.ticks <- ggplot2::element_blank()

  }

  ggplot2::theme(axis.line.y = axis.line,
                 axis.line.x = axis.line,
                 axis.ticks.x = axis.ticks,
                 axis.ticks.y = axis.ticks,
                 axis.text = axis.text,
                 axis.title = ggplot2::element_text(colour = 'dimgray'),
                 text = ggplot2::element_text(colour = 'dimgray'))
}
