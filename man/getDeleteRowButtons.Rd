% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules-sc-utils.R
\name{getDeleteRowButtons}
\alias{getDeleteRowButtons}
\title{Get HTML to make delete buttons in datatable rows}
\usage{
getDeleteRowButtons(session, len, title = "Delete file")
}
\arguments{
\item{session}{Shiny session object used for namespacing.}

\item{len}{Numeric number of rows}

\item{title}{Title for buttons. Default is \code{'Delete file'}}
}
\value{
Character vector with length \code{len} with delete buttons. Event
can be observed in a reactive by \code{input$delete_row} and the row clicked
will have have \code{'input$delete_row_i'} where 'i' is the row.
}
\description{
Get HTML to make delete buttons in datatable rows
}
