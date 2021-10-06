send_slack_error <- function() {
    user <- Sys.getenv('SHINYPROXY_USERNAME', 'localhost')
    error <- recover_error()

    url <- readRDS(system.file('extdata/slack.rds', package = 'dseqr'))

    stack <- slackify_stack(error$stack)

    httr::POST(url = url,
               httr::add_headers('Content-Type' = 'application/json'),
               body = sprintf(
                   '{"text": "`%s` \n%s \n\n *user*: %s 🙎"}',
                   error$message,
                   stack,
                   user
               ))

    shinyjs::alert('Sorry about that! \n\n Error has been reported and will be fixed promptly. \n\n (ﾉ´ｰ`)ﾉ')
}

slackify_stack <- function(stack) {

    is.user <- stack$category == 'user'

    res <- stack %>%
        dplyr::select(-category) %>%
        dplyr::mutate(num = paste0('`', num, ':`')) %>%
        tidyr::unite(col = res, sep=' ') %>%
        dplyr::pull(res)

    res[is.user] <- paste0('*', res[is.user], '*')

    paste('>', res, collapse='\n')
}

recover_error <- function ()  {

    # get calls
    calls <- sys.calls()
    from <- 0L

    # get up one frame from last stop() call
    n <- length(calls)
    for (i in rev(seq_len(n))) {

        calli <- calls[[i]]
        fname <- calli[[1L]]
        if ( "stop(e)" %in% deparse(calli)) {
            from <- i - 1
            break
        }
    }

    frame <- sys.frame(from)

    # write to logfile
    getError(frame$e)
}

# adapted from shiny::printError
getError <- function (cond,
                      full = get_devmode_option("shiny.fullstacktrace", FALSE),
                      offset = getOption("shiny.stacktraceoffset", TRUE)) {

    error_msg <- sprintf(
        "Error in %s: %s",
        shiny:::getCallNames(list(conditionCall(cond))),
        conditionMessage(cond)
    )

    should_drop <- !full
    should_strip <- !full
    should_prune <- !full

    stackTraceCalls <- c(
        attr(cond, "deep.stack.trace", exact = TRUE),
        list(attr(cond, "stack.trace", exact = TRUE))
    )

    stackTraceParents <- lapply(stackTraceCalls, attr, which = "parents", exact = TRUE)
    stackTraceCallNames <- lapply(stackTraceCalls, shiny:::getCallNames)
    stackTraceCalls <- lapply(stackTraceCalls, shiny:::offsetSrcrefs, offset = offset)

    # Use dropTrivialFrames logic to remove trailing bits (.handleSimpleError, h)
    if (should_drop) {
        # toKeep is a list of logical vectors, of which elements (stack frames) to keep
        toKeep <- lapply(stackTraceCallNames, shiny:::dropTrivialFrames)
        # We apply the list of logical vector indices to each data structure
        stackTraceCalls <- mapply(stackTraceCalls, FUN = `[`, toKeep, SIMPLIFY = FALSE)
        stackTraceCallNames <- mapply(stackTraceCallNames, FUN = `[`, toKeep, SIMPLIFY = FALSE)
        stackTraceParents <- mapply(stackTraceParents, FUN = `[`, toKeep, SIMPLIFY = FALSE)
    }

    delayedAssign("all_true", {
        # List of logical vectors that are all TRUE, the same shape as
        # stackTraceCallNames. Delay the evaluation so we don't create it unless
        # we need it, but if we need it twice then we don't pay to create it twice.
        lapply(stackTraceCallNames, function(st) {
            rep_len(TRUE, length(st))
        })
    })

    # stripStackTraces and lapply(stackTraceParents, pruneStackTrace) return lists
    # of logical vectors. Use mapply(FUN = `&`) to boolean-and each pair of the
    # logical vectors.
    toShow <- mapply(
        if (should_strip) shiny:::stripStackTraces(stackTraceCallNames) else all_true,
        if (should_prune) lapply(stackTraceParents, shiny:::pruneStackTrace) else all_true,
        FUN = `&`,
        SIMPLIFY = FALSE
    )

    dfs <- mapply(seq_along(stackTraceCalls), rev(stackTraceCalls), rev(stackTraceCallNames), rev(toShow), FUN = function(i, calls, nms, index) {
        data.frame(
            num = rev(which(index)),
            call = rev(nms[index]),
            loc = rev(shiny:::getLocs(calls[index])),
            category = rev(shiny:::getCallCategories(calls[index])),
            stringsAsFactors = FALSE
        )
    }, SIMPLIFY = FALSE)

    res <- list(
        message = error_msg,
        stack = dfs[[1]]
    )

    return(res)
}

toString.data.frame = function (object, ..., digits=NULL, quote=FALSE, right=TRUE, row.names=TRUE) {
    nRows = length(row.names(object));
    if (length(object)==0) {
        return(paste(
            sprintf(ngettext(nRows, "data frame with 0 columns and %d row", "data frame with 0 columns and %d rows")
                    , nRows)
            , "\\n", sep = "")
        );
    } else if (nRows==0) {
        return(gettext("<0 rows> (or 0-length row.names)\\n"));
    } else {
        # get text-formatted version of the data.frame
        m = as.matrix(format.data.frame(object, digits=digits, na.encode=FALSE));
        # define row-names (if required)
        if (isTRUE(row.names)) {
            rowNames = dimnames(object)[[1]];
            if(is.null(rowNames)) {
                # no row header available -> use row numbers
                rowNames = as.character(1:NROW(m));
            }
            # add empty header (used with column headers)
            rowNames = c("", rowNames);
        }
        # add column headers
        m = rbind(dimnames(m)[[2]], m);
        # add row headers
        m = cbind(rowNames, m);
        # max-length per-column
        maxLen = apply(apply(m, c(1,2), stringr::str_length), 2, max, na.rm=TRUE);

        # add right padding
        ##  t is needed because "If each call to FUN returns a vector of length n, then apply returns an array of dimension c(n, dim(X)[MARGIN])"
        m = t(apply(m, 1, stringr::str_pad, width=maxLen, side="right"));
        m = t(apply(m, 1, stringr::str_pad, width=maxLen+3, side="left"));
        # merge columns
        m = apply(m, 1, paste, collapse="");
        # merge rows (and return)
        return(paste(m, collapse="\n"));
    }
}
