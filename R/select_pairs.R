

#' Select pairs and replicates for an experiment.
#'
#' Used to select a pairs for paired-end experiments and replicate samples. Please follow prompt, ensuring correct file name matching and
#' end-type experiment identification.
#'
#' @param data_dir Directory with raw fastq.gz RNA-Seq files.
#' @param pdata_path Path to text file with sample annotations. Must be readable by \code{\link[data.table]{fread}}.
#' The first column should contain sample ids that match a single raw rna-seq data file name.
#'
#' @return data.frame of sample annotations with integer columns \code{Pair} (if paired-end) and \code{Replicate} (if added).
#' @export
#'
#' @examples
#'
#' data_dir <- system.file('extdata', 'IBD', package='drugseqr', mustWork = TRUE)
#' pdata_path <- file.path(data_dir, 'Phenotypes.csv')
#'
#' pdata <- select_pairs(data_dir, pdata_path)
#'
#'
select_pairs <- function(data_dir) {

  # setup ----

  group_colors <- c("#C7E9C0", "#C6DBEF", "#FCBBA1", "#FDD0A2", "#BCBDDC", "#D9D9D9", "#F6E8C3", "#DC143C",
                    "#A1D99B", "#9ECAE1", "#FC9272", "#FDAE6B", "#9E9AC8", "#BDBDBD", "#DFC27D", "#FFFFFF",
                    "#C3B091", "#007FFF", "#00FFFF", "#7FFFD4", "#228B22", "#808000", "#7FFF00", "#BFFF00",
                    "#FFD700", "#DAA520", "#FF7F50", "#FA8072","#FC0FC0", "#CC8899", "#E0B0FF", "#B57EDC", "#843179")

  ncolors <- length(group_colors)

  background <- 'url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'


  fastqs <- list.files(data_dir, '.fastq.gz')
  pdata <- tibble::tibble('File Name' = fastqs)

  pdata <- tibble::add_column(pdata, Pair = NA, Replicate = NA, .before = 1)

  # auto-detect if paired
  fastq_id1s <- get_fastq_id1s(file.path(data_dir, fastqs))
  paired <- detect_paired(fastq_id1s)

  # select and mark auto-detected pair type
  end_types <- c('single-ended', 'pair-ended')
  if (paired) end_types <- end_types[c(2, 1)]
  end_types[1] <- paste(end_types[1], '(detected)')

  message("For paired-end experiments:
            - select and add paired rows (usually contain _1 and _2 in filename)\n")

  message("For paired-end and single-ended experiments:
            - confirm that file names (in column 'File Name') are assigned to the correct sample rows
            - confirm that the experiment was correctly identified as single-ended or paired-end
            - select samples to treat as a single library (if any - e.g. same sample sequenced in replicate)\n")

  message("Click 'Done' to continue with kallisto quantification.")


  # things we will update/return to user
  pairs <- reps <- rep(NA, nrow(pdata))


  #  user interface ----

  ui <- miniUI::miniPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(
      shiny::tags$style("#pdata {white-space: nowrap;}"), # table text on 1 line
      shiny::tags$style(".dt-fake-height {height: 1px;}"), # to make 100% height div work
      shiny::tags$style("td.dt-nopad {padding: 0px !important; height: 100%;}"), # td for bg color group column
      shiny::tags$style("td.dt-nopad div {height: 100%; width: 100%; text-align: center;}"), # div inside td for bg color group column
      shiny::tags$style("td.dt-nopad span {display: inline-block; padding: 8px 10px; color: white;}") # span inside div inside td for bg color group column
    ),
    # title bar
    miniUI::gadgetTitleBar(shiny::textOutput('title', inline=TRUE), left = miniUI::miniTitleBarButton("reset", "Reset")),
    miniUI::miniContentPanel(
      shiny::fillCol(flex = c(NA, NA, 1),
                     shiny::selectizeInput('end_type',
                                           'Confirm end-type:',
                                           choices = end_types),
                     shiny::hr(),
                     DT::dataTableOutput("pdata")
      )
    ),
    miniUI::miniButtonBlock(
      shiny::actionButton("pair", "Pair Samples"),
      shiny::actionButton("rep", "Mark Replicates")
    )
  )

  # server ----

  server <- function(input, output, session) {

    # make reactive state value to keep track of ctrl vs test group
    state <- shiny::reactiveValues(pdata = 0)



    # show phenotype data
    output$pdata <- DT::renderDataTable({

      DT::datatable(
        isolate(pdata_r()),
        class = 'cell-border dt-fake-height',
        rownames = FALSE,
        escape = FALSE, # to allow HTML in table
        options = list(
          columnDefs = list(list(className = 'dt-nopad', targets = c(0, 1))),
          scrollY = FALSE,
          paging = FALSE,
          bInfo = 0
        )
      )
    })

    # pdata reactive so that will update Pair/Replicate column
    pdata_r <- shiny::eventReactive(state$pdata, {

      # update pdata Replicate column
      rep_nums <- sort(unique(setdiff(reps, NA)))
      for (rep_num in rep_nums) {
        color <- group_colors[ncolors - rep_num]
        rows <- which(reps == rep_num)
        pdata[rows, 'Replicate'] <<- paste('<div style="background-color:', color, ';"></div>')
      }

      # update pdata Pair column
      if (paired) {
        pair_nums <- sort(unique(setdiff(pairs, NA)))
        for (pair_num in pair_nums) {
          color <- group_colors[pair_num]
          rows <- which(pairs == pair_num)
          pdata[rows, 'Pair'] <<- paste('<div style="background:', color, background, ';"></div>')
        }
      } else {
        pdata[1:nrow(pdata), 'Pair'] <<- NA
      }

      return(pdata)
    })

    # proxy used to replace data
    pdata_proxy <- DT::dataTableProxy("pdata")
    shiny::observe({
      DT::replaceData(pdata_proxy, pdata_r(), rownames = FALSE)
    })


    # click 'Pair Samples' ----

    shiny::observeEvent(input$pair, {

      # get rows
      rows  <- input$pdata_rows_selected

      # check for incomplete/wrong input
      if (validate_pairs(pairs, rows, reps)) {

        # add rows as a pair
        pair_num <- length(unique(setdiff(pairs, NA))) + 1
        pairs[rows] <<- pair_num

        # update states to trigger updates
        state$pdata <- state$pdata + 1
      }
    })

    # click 'Mark Replicates' ----

    shiny::observeEvent(input$rep, {

      # get rows
      rows  <- input$pdata_rows_selected
      if (validate_reps(pairs, rows, reps)) {
        # add rows as replicates
        rep_num <- length(unique(setdiff(reps, NA))) + 1
        reps[rows] <<- rep_num

        # update states to trigger updates
        state$pdata <- state$pdata + 1
      }
    })

    # click 'Done' ----

    shiny::observeEvent(input$done, {
      if (paired & sum(is.na(pdata$Pair)) != 0) {
        message("All samples must be in a pair for paired-end experiments.")
      } else {

        # replace html with integers identifying pairs/replicates
        if (paired) {
          pdata$Pair <<- pairs
        } else {
          pdata$Pair <<- NULL
        }

        if (all(is.na(reps))) {
          pdata$Replicate <<- NULL
        } else {
          pdata$Replicate <<- reps
        }

        shiny::stopApp(pdata)
      }
    })

    # select 'end-type' ----
    shiny::observeEvent(input$end_type, {

      # single-ended experiments are not paired, so disable
      # hide pair column labels if select single-ended

      if (grepl('^single', input$end_type)) {
        shinyjs::disable("pair")
        paired <<- FALSE
        state$pdata <- state$pdata + 1

      } else {
        shinyjs::enable("pair")
        paired <<- TRUE
        state$pdata <- state$pdata + 1
      }
    })


    # click 'Reset' ----

    shiny::observeEvent(input$reset, {
      pairs <<- rep(NA, nrow(pdata))
      reps <<- rep(NA, nrow(pdata))
      pdata$Pair <<- NA
      pdata$Replicate <<- NA

      # remove groups from table and reset control
      state$pdata <- state$pdata + 1
    })
  }

  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::paneViewer())
}


#' Validate sample pairing for pair-ended RNA seq
#'
#' @param pairs Numeric vector of integers and/or \code{NA}. Positions with the same integer
#'   value indicate samples that are paired in a pair-ended experiment.
#' @param rows Numeric vector of integers indicating selected rows.
#' @param reps Numeric vector of integers and/or \code{NA}. Positions with the same integer
#'   value indicate samples that replicates.
#'
#' @return \code{TRUE} if the pairing is valid, otherwise \code{FALSE}.
#' @export
#' @keywords internal
#'
validate_pairs <- function(pairs, rows, reps) {

  rep_rows <- reps[rows]
  all_rep_rows <- all(rep_rows %in% rows, na.rm = TRUE)

  if (length(rows) < 2) {
    msg <- "Select at least two rows to mark as pairs."

  } else if (!anyNA(rep_rows) && length(unique(rep_rows)) < 2) {
    msg <- "Select at least two non-replicate rows to mark as pairs."

  } else if (!anyNA(rep_rows) && !all_rep_rows) {
    msg <- "All replicates must be included in the same pair."

  } else if (length(unique(rep_rows)) > 2) {
    msg <- "A pair must include exactly two replicate groups or one replicate group and one additional sample."

  } else if (!all(is.na(pairs[rows]))) {
    msg <- "Selected row(s) already belong to a pair. Click 'Reset' if you need to start over."

  } else {
    msg <- NULL
  }

  return(msg)
}

#' Validate sample replicates for RNA seq
#'
#' @inheritParams validate_pairs
#'
#' @return \code{TRUE} if the replicate is valid, otherwise \code{FALSE}.
#' @export
#' @keywords internal
validate_reps <- function(pairs, rows, reps) {

  if (length(rows) < 2) {
    msg <- "Select at least two rows to mark as replicates."

  } else if (!all(is.na(reps[rows]))) {
    msg <- "Selected row(s) already belong to a replicate. Click 'Reset' if you need to start over."

  } else if (!all(is.na(pairs[rows]))) {
    msg <- "Replicates must be specified first for paired samples that include replicates. Click 'Reset' if you need to start over."

  } else {
    msg <- NULL
  }

  return(msg)
}


#' Get first sequence identifiers for fastq.gz files
#'
#' Function is used to determine if an experiment is single or paired.
#'
#' @param fastq_paths Character vector of paths to fastq.gz files
#'
#' @return Named character vector of first sequence id lines (start with @) in \code{fastq_paths}.
#'   Names are the \code{fastq_paths}.
#' @keywords internal
#' @export
#'
#' @examples
#' # required fastq.gz files in data-raw/example-data (e.g. IBD example)
#' fastq_paths <- list.files(file.path('data-raw', 'example-data'), '.fastq.gz$', full.names = TRUE)
#' fastq_id1s <- get_fastq_id1s(fastq_paths)
#'
get_fastq_id1s <- function(fastq_paths) {

  # get first line with @ symbol (sequence identifier)
  fastq_id1s <- sapply(fastq_paths, function(f) {
    incon <- gzfile(f)
    while (TRUE) {
      line = readLines(incon, n = 1)
      if (grepl('^@', line)) {
        break
      }
    }
    close(incon)
    return(line)
  })
  return(fastq_id1s)
}


#' Detect if experiment is pair-ended.
#'
#' @param fastq_id1s Character vector of first sequence identifiers from fastq.gz files. Returned from \code{\link{get_fastq_id1s}}.
#'
#' @return boolean indicating if experiement is pair-ended (\code{TRUE}) or single-ended (\code{FALSE}).
#' @export
#' @keywords internal
detect_paired <- function(fastq_id1s) {

  # older illumina sequence identifiers have 1 part
  # newer illumina sequence identifiers have 2 space-seperated parts
  id_parts <- strsplit(fastq_id1s, ' ')
  older <- all(sapply(id_parts, length) == 1)
  newer <- all(sapply(id_parts, length) == 2)

  if (older) {
    # pair is 1 or 2 at end of sequence id after /
    pairs <- gsub('^.+?/([12])$', '\\1', fastq_id1s)

  } else if (newer) {
    # pair is 1 or 2 followed by : followed by N or Y at beginning of second part
    id_parts2 <- sapply(id_parts, `[`, 2)
    pairs <- gsub('^([12]):[YN]:.+$', '\\1', id_parts2)

  } else {
    stop("fastq.gz files don't appear to be from older/newer Illumina software. Please contact package author.")
  }

  # SRA also accepts /1 and /2 at end of read name
  is_sra <- any(grepl('^@SRR\\d+', fastq_id1s))
  if (is_sra) {
    pairs <- gsub('^.+?/([12]$)', '\\1', pairs)
  }

  # paired experiments will have '1' and '2'
  uniq.pairs <- unique(pairs)
  paired <- setequal(c('1', '2'), uniq.pairs)

  # unpaired will have only '1'
  if (!paired) stopifnot(uniq.pairs == '1')

  return(paired)
}

