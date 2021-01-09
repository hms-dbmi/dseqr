#' Modified shinyWidgets dropdownButton
#'
#' Button is wrapped in input-group-btn class to allow inline with forms. Also includes title parameter.
#'
#' @inheritParams shinyWidgets::dropdownButton
#' @param title title attribute for button
#'
#' @keywords internal
dropdownButton <- function(..., circle = TRUE, status = "default",
                           size = "default", icon = NULL,
                           label = NULL, tooltip = FALSE,
                           right = FALSE, up = FALSE,
                           circle_class = '',
                           width = NULL, margin = "10px", inline = FALSE, inputId = NULL, title = NULL, style = '') {
  size <- match.arg(arg = size, choices = c("default", "lg", "sm", "xs"))
  if (is.null(inputId)) {
    inputId <- paste0("drop", sample.int(1e9, 1))
  }



  # dropdown content
  html_ul <- list(
    class = paste("dropdown-menu", ifelse(right, "dropdown-menu-right", "")),
    class = "dropdown-shinyWidgets",
    id = paste("dropdown-menu", inputId, sep = "-"),
    style = if (!is.null(width))
      paste0("width: ", htmltools::validateCssUnit(width), ";"),
    `aria-labelledby` = inputId,
    lapply(X = list(...), FUN = htmltools::tags$li,
           style = paste0("margin-left: ", htmltools::validateCssUnit(margin),
                          "; margin-right: ", htmltools::validateCssUnit(margin), ";"))
  )

  # button
  if (circle) {
    html_button <- circleButton(
      inputId = inputId, icon = icon, status = status, size = size,
      class = paste("dropdown-toggle", circle_class),
      `data-toggle` = "dropdown"
    )
  } else {
    html_button <- list(
      class = paste0("btn btn-", status," action-button dropdown-toggle "),
      class = if (size != "default") paste0("btn-", size),
      style = style,
      type = "button",
      id = inputId,
      `data-toggle` = "dropdown",
      `aria-haspopup` = "true",
      `aria-expanded` = "true",
      list(icon, label),
      tags$span(class = "caret")
    )
    html_button <- do.call(htmltools::tags$button, html_button)
  }

  # tooltip
  if (identical(tooltip, TRUE))
    tooltip <- tooltipOptions(title = label)

  if (!is.null(tooltip) && !identical(tooltip, FALSE)) {
    tooltip <- lapply(tooltip, function(x) {
      if (identical(x, TRUE))
        "true"
      else if (identical(x, FALSE))
        "false"
      else x
    })
    tooltipJs <- htmltools::tags$script(
      sprintf(
        "$('#%s').tooltip({ placement: '%s', title: '%s', html: %s, trigger: 'hover' });",
        inputId, tooltip$placement, tooltip$title, tooltip$html
      )
    )
  } else {
    tooltipJs <- ""
  }

  if( inline ) container <- htmltools::tags$span
  else container <- htmltools::tags$div
  dropdownTag <- container(
    class = ifelse(up, "dropup", "dropdown"),
    class = "btn-dropdown-input input-group-btn",
    title = title,
    html_button, id = paste0(inputId, "_state"),
    do.call(htmltools::tags$ul, html_ul),
    tooltipJs
  )
  shinyWidgets:::attachShinyWidgetsDep(dropdownTag, "dropdown")
}
