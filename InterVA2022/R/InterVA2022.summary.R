#' Summary of the results obtained from InterVA2022 algorithm
#'
#' This function prints the summary message of the fitted results.
#'
#' @param object fitted object from \code{InterVA2022()}
#' @param top number of top CSMF to show
#' @param id the ID of a specific death to show
#' @param InterVA.rule If it is set to "TRUE", only the top 3 causes reported by
#' InterVA2022 is calculated into CSMF as in InterVA2022. The rest of probabilities
#' goes into an extra category "Undetermined". Default set to "TRUE".
#' @param ... not used
#' @keywords InterVA
#' @importFrom grDevices grey.colors 
#' @importFrom graphics barplot par pie 
#' @importFrom utils data write.table
#' @export
#' @examples
#' \dontrun{
#' data(va2022)
#' 
#' ## to get easy-to-read version of causes of death make sure the column
#' ## orders match InterVA2022 standard input this can be monitored by checking
#' ## the warnings of column names
#'
#' sample.output <- InterVA2022(va2022, HIV = "h", Malaria = "l", Covid = "v",
#'     write = FALSE, directory = tempdir(), filename = "VA2022_result",
#'     output = "extended", append = FALSE)
#'
#' summary(sample.output)
#' summary(sample.output, top = 10)
#' summary(sample.output, id = " 1")
#' }
summary.interVA2022 <- function(object, top = 5, id = NULL,
                                InterVA.rule = TRUE, ...){

    if(is.null(object$dev)) {
        data("causetext2022", envir = environment())
    	causetext2022 <- get("causetext2022", envir  = environment())
    	causenames <- causetext2022[4:67, 2]
    	causeindex <- 4:67
    } else {
        InterVA2022 <- FALSE
        causenames <- names(object$VA[[1]]$wholeprob)
        causeindex <- 1:length(causenames)
    }

    out <- NULL
    va <- object$VA
    out$top <- top
    out$N <- length(va)
    out$Malaria <- object$Malaria
    out$HIV <- object$HIV
    out$Covid <- object$Covid

    ## Get population distribution
    ## Initialize the population distribution
    dist <- NULL
    for (i in 1:length(va)) {
        if (!is.null(va[[i]][14])) {
            dist <- rep(0, length(unlist(va[[i]][14])))
            break
        }
    }
    ## determine how many causes from top need to be summarized
    undeter <- 0

    if (is.null(dist)) {
      cat("No va probability found in input")
      return()
      }
    ## Add the probabilities together
    if (!InterVA.rule) {
        for (i in 1:length(va)) {
            if (is.null(va[[i]][14])) {
              undeter = undeter + 1
              next
              }
            this.dist <- unlist(va[[i]][14])
            dist <- dist + this.dist
        }
        ## Normalize the probability for CODs
        if (undeter > 0){
            dist.cod <- c(dist[causeindex], undeter)
            dist.cod <- dist.cod/sum(dist.cod)
            names(dist.cod) <- c(causenames, "Undetermined")
        } else {
            csmf <- dist[causeindex]/sum(dist[causeindex])
            names(csmf) <- causenames
        }
    } else {
        csmf <- CSMF.interVA2022(va)
    }
    csmf <- data.frame(cause = names(csmf), likelihood = csmf)
    rownames(csmf) <- NULL

    if (!is.null(id)) {
        index <- which(object$ID == id)
        if (is.null(index)) {
            stop("Error: provided ID not found")
        } else if (is.null(va[[i]][14])) {
            out$undet <- TRUE
        } else {
            out$undet <- FALSE
            probs.tmp <- object$VA[[index]][14][[1]]
            out$preg <- probs.tmp[1:3]
            out$probs <- probs.tmp[causeindex]
            topcauses <- sort(out$probs, decreasing = TRUE)[1:top]
            out$indiv.top <- data.frame(Cause = names(topcauses))
            out$indiv.top$Likelihood <- topcauses
        }
        out$id.toprint <- id
    } else {
        out$csmf.ordered <- csmf[order(csmf[,2], decreasing = TRUE),]
    }

    out$InterVA.rule <- InterVA.rule
    class(out) <- "interVA2022_summary"
    return(out)
}

#' Print method for summary of the results obtained from InterVA2022 algorithm
#'
#' This function prints the summary message of the fitted results.
#'
#' @param x summary of InterVA2022 results
#' @param ... not used
#' @export
#' @keywords InterVA
print.interVA2022_summary <- function(x, ...) {
    # print single death summary
    if (!is.null(x$id.toprint)) {
        cat(paste0("InterVA2022 fitted top ", x$top, " causes for death ID: ",
                   x$id.toprint, "\n\n"))
        if (x$undet) {
            cat("Cause of death undetermined\n")
        } else {
            x$indiv.top[, 2] <- round(x$indiv.top[, 2], 4)
            print(x$indiv.top, row.names = FALSE, right = FALSE)
        }
    # print population summary
    } else {
        cat(paste("InterVA2022 fitted on", x$N, "deaths\n"))
        if (x$InterVA.rule) {
            cat("CSMF calculated using reported causes by InterVA2022 only\nThe remaining probabilities are assigned to 'Undetermined'\n")
        } else {
            cat("CSMF calculated using distribution over all causes\nwithout 'Undetermined' category\n")
        }
        cat("\n")
        cat(paste("Top", x$top,  "CSMFs:\n"))
        csmf.out.ordered <- x$csmf.ordered[1:x$top, ]
        csmf.out.ordered[, 2] <- round(csmf.out.ordered[, 2], 4)
        print(csmf.out.ordered, right = FALSE, row.names = F)
    }
}
