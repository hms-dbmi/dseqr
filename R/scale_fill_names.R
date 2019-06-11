
#' Get a fill scale with named values
#'
#' @param levs Character vector of levels to get a fill scale for.
#'
#' @return ggplot2 \code{scale_fill_manual} with values names from \code{levs}
#' @export
#' @keywords internal
#'
#' @examples
scale_fill_names <- function(levs) {
  nlevs <- length(levs)
  if (nlevs <= 10) {
    values <- head(.get_palette("tableau10medium"), nlevs)
    names(values) <- levs

  } else {
    if (nlevs <= 20) {
      values <- head(.get_palette("tableau20"), nlevs)
      names(values) <- levs

    } else {
      values <- viridisLite::cividis(nlevs)
      names(values) <- levs
    }
  }

  plot_scale <- ggplot2::scale_fill_manual(values = values)
  return(plot_scale)
}

# Function to define color palettes.
.get_palette <- function(palette_name) {
  switch(palette_name,
         tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                       "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                       "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                       "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
         tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                             "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                             "#CDCC5D", "#6DCCDA")
  )
}
