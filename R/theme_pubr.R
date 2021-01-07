#'Publication ready theme
#'
#'Create a publication ready theme. From package ggpubr.
#'
#'@param base_size base font size
#'@param base_family base font family
#'@param border logical value. Default is FALSE. If TRUE, add panel border.
#'@param margin logical value. Default is TRUE. If FALSE, reduce plot margin.
#'@param legend character specifying legend position. Allowed values are one of
#'  c("top", "bottom", "left", "right", "none"). Default is "top" side position.
#'  to remove the legend use legend = "none". Legend position can be also
#'  specified using a numeric vector c(x, y).  In this case it is possible to
#'  position the legend inside the plotting area. x and y are the coordinates of
#'  the legend box. Their values should be between 0 and 1. c(0,0) corresponds
#'  to the "bottom left" and c(1,1) corresponds to the "top right" position. For
#'  instance use legend = c(0.8, 0.2).
#'@param x.text.angle Rotation angle of x axis tick labels. Default value is 0.
#'  Use 90 for vertical text.
#'@param flip logical. If TRUE, grid lines are added to y axis instead of x
#'  axis.
#'@importFrom ggplot2 %+replace%
#'
#'@examples
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'    geom_point(aes(color = gear))
#'
#' # Default plot
#' p
#'
#' # Use theme_pubr()
#' p + theme_pubr()
#'
#'
theme_pubr <- function (base_size = 12, base_family = "",
                        border = FALSE, margin = TRUE,
                        legend = c("top", "bottom", "left", "right", "none"),
                        x.text.angle = 0)
{
  half_line <- base_size/2
  if(!is.numeric(legend)) legend <- match.arg(legend)
  if(x.text.angle > 5) xhjust <- 1 else xhjust <- NULL

  if(border){
    panel.border <- ggplot2::element_rect(fill = NA, colour = "black", size = 0.7)
    axis.line <- ggplot2::element_blank()
  }
  else{
    panel.border <- ggplot2::element_blank()
    axis.line = ggplot2::element_line(colour = "black", size = 0.5)
  }


  if(margin)
    plot.margin <- ggplot2::margin(half_line, half_line, half_line,
                          half_line)
  else plot.margin <- unit(c(0.5,0.3,0.3,0.3),"mm")

  .theme <- ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.border = panel.border,
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = axis.line, axis.text = ggplot2::element_text(color = "black"),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = "black", size = 0.7),
          plot.margin = plot.margin,
          legend.position = legend,
          complete = TRUE)

  if(x.text.angle!=0)
    .theme <- .theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x.text.angle, hjust = xhjust))

  .theme
}
