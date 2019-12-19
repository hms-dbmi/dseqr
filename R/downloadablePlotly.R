#' Plotly with hidden download link
#' @export
#' @keywords internal
downloadablePlotlyUI <- function(id) {
  ns <- NS(id)
  tagList(
    hiddenDownloadUI(ns('dl_plotly')),
    plotly::plotlyOutput(ns('plotly'))
  )
}

#' Logic for plotly with download data button in modebar
#' @param plotly_fun function to generate base plotly
#' @param fname_fun function that returns filename for download
#' @param data_fun function that saves file for download
#' @export
#' @keywords internal
downloadablePlotly <- function(input, output, session, plotly_fun, fname_fun = function(){NULL}, data_fun = function(){NULL}, title = 'Download plot data') {

  # hidden download link
  callModule(hiddenDownload, 'dl_plotly', reactive(input$dl_data), fname_fun, data_fun)


  # wrap to avoid shinytest errors
  output$plotly <- snapshotPreprocessOutput(
    plotly::renderPlotly({

      # reactive counter to trigger plot data download
      # incremented on clock modebar dl button
      ns_val <- isolate(input$dl_data)
      if (is.null(ns_val)) ns_val <- 0

      # id of ns_val to increment on click
      ns_id = session$ns('dl_data')

      # button to add to modebar
      dl_data_button <- list(
        name = title,
        icon = list(
          path = "M216 0h80c13.3 0 24 10.7 24 24v168h87.7c17.8 0 26.7 21.5 14.1 34.1L269.7 378.3c-7.5 7.5-19.8 7.5-27.3 0L90.1 226.1c-12.6-12.6-3.7-34.1 14.1-34.1H192V24c0-13.3 10.7-24 24-24zm296 376v112c0 13.3-10.7 24-24 24H24c-13.3 0-24-10.7-24-24V376c0-13.3 10.7-24 24-24h146.7l49 49c20.1 20.1 52.5 20.1 72.6 0l49-49H488c13.3 0 24 10.7 24 24zm-124 88c0-11-9-20-20-20s-20 9-20 20 9 20 20 20 20-9 20-20zm64 0c0-11-9-20-20-20s-20 9-20 20 9 20 20 20 20-9 20-20z",
          transform = " ",
          width =  564,
          ascent = 564,
          descent = 0
        ),
        click = htmlwidgets::JS(
          paste0("function(gd) {
             Shiny.setInputValue('", ns_id, "', ", ns_val + 1, ");
          }")
        )
      )


      plotly_fun() %>%
        plotly::config(modeBarButtonsToAdd = list(dl_data_button))
    }),
    function(value) { session$ns('downloadable_plotly') }
  )


}



#' Hidden download link
#' @export
#' @keywords internal
hiddenDownloadUI <- function(id) {
  ns <- NS(id)
  downloadLink(ns("downloadData"), '')
}


#' Logic for hidden download link
#' @param check reactive value to trigger download
#' @param fname_fun function that returns filename for download
#' @param data_fun function that saves file for download
#' @export
#' @keywords internal
hiddenDownload <- function(input, output, session, check, fname_fun, data_fun) {

  observeEvent(check(), {
    runjs(paste0("$('#", session$ns('downloadData'), "')[0].click();"))
  })

  output$downloadData <- downloadHandler(
    filename = fname_fun,
    content = data_fun
  )
}
