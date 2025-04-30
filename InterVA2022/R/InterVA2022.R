#' Provide InterVA2022 analysis on the data input.
#'
#' This function implements the InterVA algorithm specifically for data
#' collected using the 2022 WHO verbal autopsy instrument.  It
#' produces individual cause of death (COD) and population cause-specific mortality
#' fractions.  The output is saved in a .csv file specified by user.
#' The calculation is based on the conditional and prior distribution
#' of 64 CODs. The function can also  save the full probability distribution
#' of each individual to file. All information about each individual is
#' saved to a va class object.
#'
#' Be careful if the input file does not match InterVA input format strictly.
#' The function will run normally as long as the number of symptoms are
#' correct. Any inconsistent symptom names will be printed in console as
#' warning. If there is a wrong match of symptom from warning, please change 
#' the input to the correct order.
#'
#' @param Input A matrix input, or data read from csv files in the same format
#' as required by InterVA2022.
#' @param HIV An indicator of the level of prevalence of HIV. The input should
#' be one of the following: "h"(high),"l"(low), or "v"(very low).
#' @param Malaria An indicator of the level of prevalence of Malaria. The input
#' should be one of the following: "h"(high),"l"(low), or "v"(very low).
#' @param Covid An indicator of the level of prevalence of Covid19. The input
#' should be one of the following: "h"(high),"l"(low), or "v"(very low).
#' @param write A logical value indicating whether or not the output (including
#' errors and warnings) will be saved to file.  If the value is set to TRUE, the
#' user must also provide a value for the parameter "directory".
#' @param directory The directory to store the output from InterVA2022. It should
#' either be an existing valid directory, or a new folder to be created. If no
#' path is given and the parameter for "write" is true, then the function stops
#' and and error message is produced.
#' @param filename The file name the user wish to save the output. No extension
#' needed. The output is in .csv format by default.
#' @param output "classic": The same delimited output format as InterVA2022; or
#' "extended": delimited output followed by full distribution of cause of
#' death probability.
#' @param append A logical value indicating whether or not the new output
#' should be appended to the existing file.
#' @param groupcode A logical value indicating whether or not the group code
#' will be included in the output causes.
#' @param sci A data frame that contains the symptom-cause-information (aka
#' Probbase) that InterVA uses to assign a cause of death.
#' @param ... not used
#' @return \item{ID }{ identifier from batch (input) file} \item{MALPREV
#' }{ selected malaria prevalence} \item{HIVPREV }{ selected HIV prevalence}
#' \item{COVIDPREV }{ selected Covid19 prevalence}
#' \item{PREGSTAT }{most likely pregnancy status} \item{PREGLIK }{ likelihood of
#' PREGSTAT} \item{PRMAT }{ likelihood of maternal death} \item{INDET
#' }{ indeterminate outcome} \item{CAUSE1 }{ most likely cause} \item{LIK1 }{
#' likelihood of 1st cause} \item{CAUSE2 }{ second likely cause} \item{LIK2 }{
#' likelihood of 2nd cause} \item{CAUSE3 }{ third likely cause} \item{LIK3 }{
#' likelihood of 3rd cause}
#' \item{wholeprob }{ full distribution of causes of death}
#'
#' @author CDC Foundation
#' @seealso \code{\link{InterVA2022.plot}}
#' @keywords InterVA
#' @export InterVA2022
#' @examples
#'
#' data(va2022)
#' 
#' ## to get easy-to-read version of causes of death make sure the column
#' ## orders match InterVA2022 standard input this can be monitored by checking
#' ## the warnings of column names
#'
#' sample.output1 <- InterVA2022(va2022, HIV = "h", Malaria = "l", 
#'     Covid="v", write = FALSE, directory = tempdir(),
#'     filename = "VA2022_result", output = "extended", append = FALSE)
#'
#' \dontrun{
#' ## to get causes of death with group code for further usage
#' sample.output2 <- InterVA2022(va2022, HIV = "h", Malaria = "l", Covid="v",
#'     write = FALSE, directory = "VA test", filename = "VA2022_result_wt_code",
#'     output = "classic", append = FALSE, groupcode = TRUE)
#'}
#' 
#'
InterVA2022 <- function (Input, HIV, Malaria, Covid,
                         write = TRUE, directory = NULL,
                         filename = "VA2022_result", output = "classic",
                         append = FALSE, groupcode = FALSE, sci = NULL,
                         ...) 
{
  VA2022 <- function(ID, MALPREV, HIVPREV, COVIDPREV, PREGSTAT, PREGLIK,
                     CAUSE1, LIK1, CAUSE2, LIK2, CAUSE3, LIK3, INDET,
                     wholeprob, ...) {
    ID <- ID
    MALPREV <- as.character(MALPREV)
    HIVPREV <- as.character(HIVPREV)
    COVIDPREV <- as.character(COVIDPREV)
    PREGSTAT <- PREGSTAT
    PREGLIK <- PREGLIK
    wholeprob <- wholeprob
    va2022.out <- list(ID = ID, MALPREV = MALPREV, HIVPREV = HIVPREV,
                       COVIDPREV = COVIDPREV, PREGSTAT = PREGSTAT,
                       PREGLIK = PREGLIK, CAUSE1 = CAUSE1, LIK1 = LIK1,
                       CAUSE2 = CAUSE2, LIK2 = LIK2, CAUSE3 = CAUSE3,
                       LIK3 = LIK3, INDET = INDET, wholeprob = wholeprob)
    va2022.out
  }
  
  save.va2022 <- function(x, filename, write) {
    if (!write) {
      return()
    }
    x <- x[-14]
    x <- as.matrix(x)
    filename <- paste(filename, ".csv", sep = "")
    write.table(t(x), file = filename, sep = ",", append = TRUE,
                row.names = FALSE, col.names = FALSE)
  }
  save.va2022.prob <- function(x, filename, write) {
    if (!write) {
      return()
    }
    prob <- unlist(x[14])
    x <- x[-14]
    x <- unlist(c(as.matrix(x), as.matrix(prob)))
    filename <- paste(filename, ".csv", sep = "")
    write.table(t(x), file = filename, sep = ",", append = TRUE,
                row.names = FALSE, col.names = FALSE)
  }
  if (is.null(directory) & write)
    stop("error: please provide a directory (required when write = TRUE)")
  if (is.null(directory)) 
    directory = getwd()
  dir.create(directory, showWarnings = FALSE)
  globle.dir <- getwd()
  setwd(directory)
  
  if (is.null(sci)) {
    data("probbase2022", envir = environment())
    probbase2022 <- get("probbase2022", envir = environment())
    probbase2022 <- as.matrix(probbase2022)
    probbase2022Version <- probbase2022[1, 3]
  }
  if (!is.null(sci)) {
    validSCI <- TRUE
    if (!is.data.frame(sci) & !is.matrix(sci)) validSCI <- FALSE
    if (nrow(sci) != 343) validSCI <- FALSE
    if (ncol(sci) != 71) validSCI <- FALSE
    if (!validSCI) {
      stop("error: invalid sci (must be data frame or matrix with 343 rows and 71 columns).")
    }
    probbase2022 <- as.matrix(sci)
    probbase2022Version <- probbase2022[1, 3]
  }
  message("Using Probbase version:  ", probbase2022Version)
  data("causetext2022", envir = environment())
  causetext2022 <- get("causetext2022", envir = environment())
  if (groupcode) {
    causetext2022 <- causetext2022[, -2]
  } else {
    causetext2022 <- causetext2022[, -3]
  }
  if (write) {
    cat(paste("Error & warning log built for InterVA2022", Sys.time(), "\n"),
        file = "errorlog2022.txt", append = FALSE)
  }
  Input <- as.matrix(Input)
  if (dim(Input)[1] < 1) {
    stop("error: no data input")
  }
  N <- dim(Input)[1]
  S <- dim(Input)[2]
  if (S != dim(probbase2022)[1]) {
    stop("error: invalid data input format. Number of values incorrect")
  }
  if (tolower(colnames(Input)[S]) != "i446o") {
    stop("error: the last variable should be 'i446o'")
  }
  data("va2022", envir = environment())
  va2022 <- get("va2022", envir = environment())
  valabels = colnames(va2022)
  count.changelabel = 0
  for (i in 1:S) {
    if (tolower(colnames(Input)[i]) != tolower(valabels)[i]) {
      warning(paste("Input column '", colnames(Input)[i],
                    "' does not match InterVA2022 standard: '",
                    valabels[i], "'", sep = ""),
              call. = FALSE, immediate. = TRUE)
      count.changelabel = count.changelabel + 1
    }
  }
  if (count.changelabel > 0) {
    warning(paste(count.changelabel,
                  "column names changed in input. \n If the change in ",
                  "undesirable, please change in the input to match ",
                  "standard InterVA2022 input format."),
            call. = FALSE, immediate. = TRUE)
    colnames(Input) <- valabels
  }
  pb_ncol <- ncol(probbase2022)
  pb_nrow <- nrow(probbase2022)
  
  # recode prior
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "I"  ] <- 1
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "A+" ] <- 0.8
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "A"  ] <- 0.5
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "A-" ] <- 0.2
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "B+" ] <- 0.1
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "B"  ] <- 0.05
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "B-" ] <- 0.02
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "C+" ] <- 0.01
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "C"  ] <- 0.005
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "C-" ] <- 0.002
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "D+" ] <- 0.001
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "D"  ] <- 5e-04
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "D-" ] <- 1e-04
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "E"  ] <- 1e-05
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == "N"  ] <- 0
  probbase2022[1, 5:pb_ncol][probbase2022[1, 5:pb_ncol] == ""   ] <- 0
  # recode pregnancy indicators
  probbase2022[, 5:7][probbase2022[, 5:7] == "I"  ] <- 1
  probbase2022[, 5:7][probbase2022[, 5:7] == "A+" ] <- 0.8
  probbase2022[, 5:7][probbase2022[, 5:7] == "A"  ] <- 0.5
  probbase2022[, 5:7][probbase2022[, 5:7] == "A-" ] <- 0.2
  probbase2022[, 5:7][probbase2022[, 5:7] == "B+" ] <- 0.1
  probbase2022[, 5:7][probbase2022[, 5:7] == "B"  ] <- 0.05
  probbase2022[, 5:7][probbase2022[, 5:7] == "B-" ] <- 0.02
  probbase2022[, 5:7][probbase2022[, 5:7] == "C+" ] <- 0.01
  probbase2022[, 5:7][probbase2022[, 5:7] == "C"  ] <- 0.005
  probbase2022[, 5:7][probbase2022[, 5:7] == "C-" ] <- 0.002
  probbase2022[, 5:7][probbase2022[, 5:7] == "D+" ] <- 0.001
  probbase2022[, 5:7][probbase2022[, 5:7] == "D"  ] <- 5e-04
  probbase2022[, 5:7][probbase2022[, 5:7] == "D-" ] <- 1e-04
  probbase2022[, 5:7][probbase2022[, 5:7] == "E"  ] <- 1e-05
  probbase2022[, 5:7][probbase2022[, 5:7] == "N"  ] <- 0
  probbase2022[, 5:7][probbase2022[, 5:7] == ""   ] <- 0
  # recode Pr(S|C)
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "I" ] <- 1
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "A" ] <- 0.8
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "B" ] <- 0.5
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "C" ] <- 0.1
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "D" ] <- 0.01
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "E" ] <- 0.001
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "F" ] <- 1e-04
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "G" ] <- 1e-05
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "H" ] <- 1e-06
  probbase2022[2:pb_nrow, 8:pb_ncol][probbase2022[2:pb_nrow, 8:pb_ncol] == "N" ] <- 0
  
  probbase2022[1, 1:4] <- rep(0, 4)
  Sys_Prior <- as.numeric(probbase2022[1, ])
  D <- length(Sys_Prior)
  HIV <- tolower(HIV)
  Malaria <- tolower(Malaria)
  Covid <- tolower(Covid)
  if (!(HIV %in% c("h", "l", "v")) || 
      !(Malaria %in% c("h","l", "v")) ||
      !(Covid %in% c("h", "l", "v"))) {
    stop("error: the HIV, Malaria, and Covid indicators should be one ",
         "of the three: 'h', 'l', and 'v'")
  }
  if (HIV == "h") Sys_Prior[10] <- 0.05
  if (HIV == "l") Sys_Prior[10] <- 0.005
  if (HIV == "v") Sys_Prior[10] <- 1e-05
  if (Malaria == "h") {
    Sys_Prior[12] <- 0.05
    Sys_Prior[34] <- 0.05
  }
  if (Malaria == "l") {
    Sys_Prior[12] <- 0.005
    Sys_Prior[34] <- 1e-05
  }
  if (Malaria == "v") {
    Sys_Prior[12] <- 1e-05
    Sys_Prior[34] <- 1e-05
  }
  if (Covid == "h") Sys_Prior[20] <- 0.05
  if (Covid == "l") Sys_Prior[20] <- 0.005
  if (Covid == "v") Sys_Prior[20] <- 1e-05
  
  ID.list <- rep(NA, N)
  VAresult <- vector("list", N)
  if (write && append == FALSE) {
    header = c("ID", "MALPREV", "HIVPREV", "COVIDPREV", "PREGSTAT", "PREGLIK", 
               "CAUSE1", "LIK1", "CAUSE2", "LIK2", "CAUSE3", "LIK3",
               "INDET")
    if (output == "extended") 
      header = c(header, as.character(causetext2022[, 2]))
    write.table(t(header), file = paste(filename, ".csv", sep = ""),
                row.names = FALSE, col.names = FALSE, sep = ",")
  }
  nd <- max(1, round(N/100))
  np <- max(1, round(N/10))
  
  if (write) {
    cat(paste("\n\n", "the following records are incomplete and excluded ",
              "from further processing:", "\n\n",
              sep=""), file = "errorlog2022.txt", append = TRUE)
  }
  
  errors <- NULL
  for (i in 1:N) {
    if (i%%nd == 0) {
      cat(".")
    }
    if (i%%np == 0) {
      cat(paste(round(i/N * 100), "% completed\n", sep = ""))
    }
    if (i == N) {
      cat(paste("100% completed\n", sep = ""))
    }
    
    index.current <- as.character(Input[i, 1])
    Input[i, which(toupper(Input[i, ]) == "N")] <- "0"
    Input[i, which(toupper(Input[i, ]) == "Y")] <- "1"
    Input[i, which(Input[i, ] != "1" & Input[i, ] != "0")] <- NA
    input.current <- as.numeric(Input[i, ])
    
    input.current[1] <- 0                      
    if (sum(input.current[5:11], na.rm=TRUE) < 1) {
      if (write) {
        errors <- rbind(errors,
                        paste(index.current,
                              " Error in age indicator: Not Specified "))
      }
      next
    }
    if (sum(is.na(input.current[3:4])) == 2) {
      if (write) {
        errors <- rbind(errors,
                        paste(index.current,
                              " Error in sex indicator: Not Specified "))
      }
      next
    }
    if (sum(input.current[19:343], na.rm=TRUE) < 1) {
      if (write) {
        errors <- rbind(errors,
                        paste(index.current,
                              " Error in indicators: No symptoms specified "))
      }
      next
    }
    
    new.input <- rep(0, S)
    for (y in 2:S) {
      if (!is.na(input.current[y])) {
        if (input.current[y] == 1) {
          new.input[y] <- 1
        }
      }
    }
    
    input.current[input.current==0] <- 1  ## don't think you need this
    input.current[1] <- 0
    input.current[is.na(input.current)] <- 0
    reproductiveAge <- 0
    preg_state      <- " "
    lik.preg        <- " "
    
    if ( (new.input[3] == 1 ||  new.input[4] == 1) &&
        (new.input[16] == 1 || new.input[17] == 1 || new.input[18] == 1)) {
      reproductiveAge <- 1
    }
    
    prob <- Sys_Prior[5:D]
    temp <- which(new.input[2:length(input.current)] == 1)
    for (jj in 1:length(temp)) {
      temp_sub <- temp[jj]
      for (j in 5:D) {
        prob[j - 4] <- prob[j - 4] * as.numeric(probbase2022[temp_sub + 1, j])
      }
      if (sum(prob[1:3]) > 0)
        prob[1:3] <- prob[1:3]/sum(prob[1:3])
      if (sum(prob[4:67]) > 0)
        prob[4:67] <- prob[4:67]/sum(prob[4:67])
    }
    names(prob) <- causetext2022[, 2]
    prob_A <- prob[ 1: 3]
    prob_B <- prob[ 4:67]
    
    ## Determine Preg_State and Likelihood
    if (sum(prob_A) == 0 || reproductiveAge == 0) {
      preg_state <- "n/a"
      lik.preg <- " "
    }
    if (max(prob_A) < 0.1 & reproductiveAge == 1) {
      preg_state <- "indeterminate"
      lik.preg <- " "
    }
    if (which.max(prob_A) == 1 && prob_A[1] >= 0.1 && reproductiveAge == 1) {
      preg_state <- "Not pregnant or recently delivered"
      lik.preg <- as.numeric(round(prob_A[1]/sum(prob_A) * 100))
    }
    if (which.max(prob_A) == 2 && prob_A[2] >= 0.1 && reproductiveAge == 1) {
      preg_state <- "Pregnancy ended within 6 weeks of death"
      lik.preg <- as.numeric(round(prob_A[2]/sum(prob_A) * 100))
    }
    if (which.max(prob_A) == 3 && prob_A[3] >= 0.1 && reproductiveAge == 1) {
      preg_state <- "Pregnant at death"
      lik.preg <- as.numeric(round(prob_A[3]/sum(prob_A) * 100))
    }
    
    ## Determine the output of InterVA
    prob.temp <- prob_B
    if (max(prob.temp) < 0.4) {
      cause1 <- lik1 <- cause2 <- lik2 <- cause3 <- lik3 <- " "
      indet <- 100
    }
    if (max(prob.temp) >= 0.4) {
      lik1 <- round(max(prob.temp) * 100)
      cause1 <- names(prob.temp)[which.max(prob.temp)]
      prob.temp <- prob.temp[-which.max(prob.temp)]
      lik2 <- round(max(prob.temp) * 100)
      cause2 <- names(prob.temp)[which.max(prob.temp)]
      if (max(prob.temp) < 0.5 * max(prob_B))
        lik2 <- cause2 <- " "
      prob.temp <- prob.temp[-which.max(prob.temp)]
      lik3 <- round(max(prob.temp) * 100)
      cause3 <- names(prob.temp)[which.max(prob.temp)]
      if (max(prob.temp) < 0.5 * max(prob_B)) 
        lik3 <- cause3 <- " "
      top3 <- as.numeric(c(lik1, lik2, lik3))
      indet <- round(100 - sum(top3, na.rm=TRUE))
    }
    
    ID.list[i] <- index.current
    VAresult[[i]] <- VA2022(ID = index.current, MALPREV = Malaria,
                            HIVPREV = HIV, COVIDPREV = Covid,
                            PREGSTAT = preg_state, PREGLIK = lik.preg,
                            CAUSE1 = cause1, LIK1 = lik1, CAUSE2 = cause2,
                            LIK2 = lik2, CAUSE3 = cause3, LIK3 = lik3,
                            INDET = indet, wholeprob = c(prob_A, prob_B))
    if (output == "classic") 
      save.va2022(VAresult[[i]], filename = filename, write)
    if (output == "extended") 
      save.va2022.prob(VAresult[[i]], filename = filename, write)
  }
  if (write) {
    cat(errors, file="errorlog2022.txt", sep="\n", append=TRUE)
  }
  
  setwd(globle.dir)
  
  out <- list(ID = ID.list[which(!is.na(ID.list))],
              VA2022 = VAresult[which(!is.na(ID.list))], 
              Malaria = Malaria, HIV = HIV, Covid = Covid)
  class(out) <- "interVA2022"
  return(out)
}
