.onLoad <- function(lib, pkg){
  Rdpack::Rdpack_bibstyles(package = pkg, authors = "LongNames")
  invisible(NULL)
}


#' @title Estimating Multiple Breakpoints for a Sequence of Realizations of
#' Bernoulli Variables
#'
#' @description The iterative procedure estimates structural changes in the
#' success probability of Bernoulli variables. It estimates the number and
#' location of the breakpoints as well as the success probabilities of the
#' sequences between each pair of neighbouring breakpoints.
#' @author Nicolas Froelich
#' @param data A two-column matrix with the location in the first column and
#' the corresponding realizations of the Bernoulli variables in the second
#' column, a vector with the ordered, realizations of the Bernoulli variables or
#' an equivalent data frame. Note that the realizations of the vector
#' respectively the second column of the matrix or the data frame must be zero
#' or one for each element.
#' @param number_bp Number of breakpoints if known a priori. By default, the
#' number of breakpoints is unknown.
#' @param max_bp The maximum number of breakpoints to be estimated (just for the
#' case, where the number of breakpoints is unknown a priori and the chosen
#' information criterion does not stop the procedure before)
#' @param inf_crit Must be one of "BIC" (Bayesian Information Criterion,
#' default), "HQC" (Hannan-Quinn Criterion) or "AIC" (Akaike Information
#' Criterion)
#' @param ext_out If TRUE (default), all function values are stored in the
#' iterative procedure and hidden printed in the output afterwards. This may take
#' additional computing time in large data sets or simulation studies. For the
#' method \link[=plot.mBP]{plot}, the default setting is required.
#' @importFrom Rdpack reprompt
#' @references
#'\insertRef{froelich}{MultipleBreakpoints}
#' @return A list containing the following elements:
#' \item{Breakpoints}{A vector containing the estimated breakpoints in
#' increasing order.}
#' \item{Probabilities}{A vector containing the estimated success
#' probabilities in each class.}
#' \item{Information Criterion}{A vector containing the values of the chosen
#' Information Criterion before the first iteration (thus without a breakpoint)
#' and after each new estimated breakpoint}
#' \item{S}{Only available, if ext_out set to TRUE. A matrix
#' containing the function values, each column representing one iteration}
#' @seealso S3 method \link[=plot.mBP]{plot} for the class "mBP".
#' @examples
#' mBP <- multiple_breakpoints(c(rbinom(1000, 1, 0.5),
#'                             rbinom(1000, 1, 0.1),
#'                             rbinom(1000, 1, 0.2)))
#' plot(mBP)
#' multiple_breakpoints(matrix(c(sort(rnorm(1000)),
#'                             rbinom(500, 1, 0.5),
#'                             rbinom(500, 1, 0.1)), ncol = 2),
#'                             inf_crit = "HQC")
#' multiple_breakpoints(matrix(c(sort(rnorm(1500)),
#'                             rbinom(500, 1, 0.1),
#'                             rbinom(500, 1, 0.3),
#'                             rbinom(500, 1, 0.4)), ncol = 2), number_bp = 2)
#' multiple_breakpoints(matrix(c(1:200, rep(201,5), 202:396,
#'                             rbinom(250,1,0.9), rbinom(150,1,0.75)),
#'                             ncol = 2), number_bp = 1)
#' @import graphics
#' @import stats
#' @export
multiple_breakpoints <- function(data, number_bp = "Unknown", max_bp = 80,
                                inf_crit = "BIC", ext_out = "TRUE") {
  stopifnot("data must be a vector, matrix or data frame" =
              (is.vector(data) | is.data.frame(data) |
                 is.matrix(data)) & !is.list(data))
  stopifnot("inf_crit must be AIC, BIC or HQC" = inf_crit %in% c("AIC", "BIC",
                                                               "HQC"))
  if (is.vector(data)) {
      data <- cbind(seq_len(length(data)), data)
  }
  if (is.data.frame(data)) {
    data <- as.matrix(data)
    if (ncol(data) == 1) data <- cbind(seq_len(length(data)), data)
  }
  if (is.matrix(data)) {
    stopifnot("Vector respectively second column of matrix or data frame is not
              binary" = all(sort(unique(data[, 2])) == c(0, 1)))
      if (dim(data)[2] == 1) {
        data <- cbind(seq_len(length(data)), data)
      }
      data <- data[order(data[, 1]), ]
      number <- nrow(data)
      data2 <- data
      colnames(data) <- c("V1", "V2")
      data <- as.matrix(cbind(aggregate(V2 ~ V1, data = data, sum),
                              aggregate(V2 ~ V1, data = data, length)[, 2]))

        data <- cbind(data, 1 / number *
                        cumsum((data[, 2] - data[, 3] *
                                  mean(data2[, 2])) / sd(data2[, 2])))
        if (ext_out == TRUE) {
          dataexport <- data[, 1]
        }
        l <- limit <- correction <- 0
        stand_dev <-  1
        estimate <- safe_index  <-  ic <- numeric()
          loglik <- sum(data[, 2]) * log(mean(data2[, 2])) +
            (sum(data[, 3]) - sum(data[, 2])) * log(1 - mean(data2[, 2]))
        ic <- - 2 * loglik
        q <- 2
        estimate[1] <- -Inf
         group <- vector(mode = "list")

        if (number_bp == "Unknown") {
          number_bp <- max_bp
          }

        while (limit < number_bp) {
          l <- l + 1
          correction <- correction - 1
          group <- group2 <- vector(mode = "list")
          estimate[length(estimate) + 1] <- data[which.max(abs(data[, 4])), 1]
          estimatesort <- c(sort(estimate), Inf)
          border <- which(data[, 1] <=
                            estimatesort[which(estimate[length(estimate)] ==
                                                 estimatesort) + 1] &
                            data[, 1] >
                            estimatesort[which(estimate[length(estimate)] ==
                                                 estimatesort) - 1])
          if (((quantile(abs(data[border, 4]), type = 1, probs = 0.85) >=
               max(abs(data[, 4])) * 0.95) & length(data[border, 4]) >= 20) &
             max(abs(data[, 4])) > 0) {
            safe_index <- c(safe_index, length(estimate))
            correction <- 3
          } else{
            if (length(safe_index) > 0) {
              estimate <- estimate[-safe_index]
              estimatesort <- c(sort(estimate), Inf)
              safe_index <- numeric()
            }
            loglik <- numeric()
            for (i in 2:length(estimatesort)) {
              group[[i - 1]] <-
                matrix(c(data[data[, 1] <= estimatesort[i] &
                                data[, 1] > estimatesort[i - 1], 2:3]),
                                       ncol = 2)
              if (sum(group[[i - 1]][, 1]) == 0 | sum(group[[i - 1]][, 1]) ==
                  sum(group[[i - 1]][, 2])) {
                loglik[i - 1] <- 0
              } else {
                loglik[i - 1] <- sum(group[[i - 1]][, 1]) *
                  log(weighted.mean(group[[i - 1]][, 1] / group[[i - 1]][, 2],
                                    group[[i - 1]][, 2])) +
                  (sum(group[[i - 1]][, 2]) - sum(group[[i - 1]][, 1])) *
                  log(1 - weighted.mean(group[[i - 1]][, 1] /
                                          group[[i - 1]][, 2],
                                        group[[i - 1]][, 2]))
              }
            }
            if (inf_crit == "BIC") {
              ic <- c(ic, - 2 * sum(loglik) + log(number) *
                        (length(estimate) - 1))
            } else if (inf_crit == "HQC") {
              ic <- c(ic, - 2 * sum(loglik) + 2 * (length(estimate) - 1) *
                        log(log(number)))
            } else if (inf_crit == "AIC") {
              ic <- c(ic, - 2 * sum(loglik) + 2 * (length(estimate) - 1))
            }
            if (number_bp >= max_bp) {
              if (ic[length(ic)] > ic[length(ic) - 1] & (correction < 1)) {
                limit <- Inf
              } else {
                limit <- limit + 1
                if (ext_out == TRUE) {
                  dataexport <- cbind(dataexport, data[, 4])
                }
              }
              correction <- 0
            } else {
              limit <- limit + 1
              if (ext_out == TRUE) {
                dataexport <- cbind(dataexport, data[, 4])
              }
            }
          }
          vektor <- numeric()
          stand_dev <- vector(mode = "list")
          for (j in 2:length(estimatesort)) {
            group[[j - 1]] <-
              matrix(c(data[data[, 1] <= estimatesort[j] &
                            data[, 1] > estimatesort[j - 1], 2:3]),
                                   ncol = 2)
            group2[[j - 1]] <- c(data2[data2[, 1] <= estimatesort[j] &
                                       data2[, 1] > estimatesort[j - 1], 2])
            if (sum(group[[j - 1]][, 2]) > 1) {

              stand_dev[[j - 1]] <- sd(group2[[j - 1]])
              if (stand_dev[[j - 1]] == 0) {
                stand_dev[[j - 1]] <- 1
              }
            } else {
              stand_dev[[j - 1]] <- 1
            }
            vektor <-
              c(vektor,
                1 / number * cumsum((group[[j - 1]][, 1] - group[[j - 1]][, 2] *
                                       weighted.mean(group[[j - 1]][, 1] /
                                                       group[[j - 1]][, 2],
                                                     group[[j - 1]][, 2])) /
                                      stand_dev[[j - 1]]))
          }
          data <- cbind(data[, 1:3], round(vektor, 12))
          q <- q + 1
        }
        if (limit == Inf) {
          estimate <- estimate[-length(estimate)]
          }
        colnames(data) <- c("X", "Y", paste("S_n", 1:(ncol(data) - 2)))
        success <- numeric()
        estimatesort <- c(sort(estimate), Inf)
        for (i in 2:length(estimatesort)) {
          group[[i - 1]] <-
            matrix(c(data[data[, 1] <= estimatesort[i] &
                            data[, 1] > estimatesort[i - 1], 2:3]),
                                 ncol = 2)
          success[i - 1] <- sum(group[[i - 1]][, 1]) / sum(group[[i - 1]][, 2])
        }

        if (ext_out == TRUE) {
          output <- c(list(Breakpoints = sort(estimate[-1]),
                           Probabilities = success,
                           ic = ic, S = dataexport))
          names(output)[3] <- inf_crit
          class(output) <- c("mBP")
          print(output[1:3])
          return(invisible(output))

        } else {
          output <- c(list(Breakpoints = sort(estimate[-1]),
                           Probabilities = success, ic = ic))
          names(output)[3] <- inf_crit
          class(output) <- c("mBP")
          return(output)
        }
  }
}


#' @title Plotting the Results of the \link{multiple_breakpoints} function
#'
#' @description  Plotting the empirical processes, the success probabilities and
#' breakpoints estimated by the \link{multiple_breakpoints} function
#' @author Nicolas Froelich
#' @param x The result of a call to \link{multiple_breakpoints}
#' @param ask logical value. If \code{TRUE} (and the R session is interactive)
#' the user is asked for input, before a new figure is drawn
#' (see \code{\link{devAskNewPage}}).
#' @param ... Further arguments are currently ignored. Only for compatibility
#' with generic functions.
#' @usage \method{plot}{mBP}(x, ask=TRUE, ...)
#' @importFrom Rdpack reprompt
#' @references
#'\insertRef{froelich}{MultipleBreakpoints}
#' @import grDevices
#' @examples
#' mBP <- multiple_breakpoints(matrix(c(sort(rnorm(2000)),
#'                             rbinom(1000, 1, 0.2),
#'                             rbinom(1000, 1, 0.6)), ncol = 2))
#' plot(mBP)
#' @export
plot.mBP <- function(x, ask = TRUE, ...) {
    stopifnot("Plot cannot be drawn. Set ext_out argument to TRUE in
              multiple_breakpoints command" = length(x) != 3)
    devAskNewPage(ask = ask)
    if (length(dim(x[[4]])) != 0) {
      old_par <- par("mar")
      par(mar = c(4.5, 5, 2.5, 0))
       par(omd = c(0, 0.8, 0, 1))

      if (dim(x[[4]])[2] >= 3) {
        plot(x[[4]][, 1], x[[4]][, 2],
             ylim = c(1.1 * min(x[[4]][, -1]), 1.1 * max(x[[4]][, -1])),
             xlab = "x", ylab = expression(paste(S[n]^(m), "(x)", sep = "")),
             type = "l")
        points(c(min(x[[4]][, 1]), max(x[[4]][, 1])), c(0, 0),
               type = "l", lty = "dashed")
        location <- which.max(abs(x[[4]][, 2]))
        segments(x[[4]][location, 1], 0, x1 = x[[4]][location, 1],
                 y1 = x[[4]][location, 2], lty = "dashed")
        for (i in 3:ncol(x[[4]])) {
          points(x[[4]][, 1], x[[4]][, i], col = i - 1, type = "l")
          location <- which.max(abs(x[[4]][, i]))
          segments(x[[4]][location, 1], 0, x1 = x[[4]][location, 1],
                   y1 = x[[4]][location, i], col = i - 1, lty = "dashed")
        }
        legend(par("usr")[2], par("usr")[4], xpd = NA,
               legend = paste("m=", 1:(ncol(x[[4]]) - 1), sep = ""),
               col = 1:(ncol(x[[4]]) - 1), lty = rep(1, (ncol(x[[4]]) - 1)),
               bty = "n", title = "Step")
      } else {
        plot(x[[4]][, 1], x[[4]][, 2], ylim = c(1.1 * min(x[[4]][, -1]),
                                            1.1 * max(x[[4]][, -1])),
             xlab = "x", ylab = expression(paste(S[n]^m, "(x)", sep = "")),
             type = "l")
        points(c(min(x[[4]][, 1]), max(x[[4]][, 1])), c(0, 0),
               type = "l", lty = "dashed")
        location <- which.max(abs(x[[4]][, 2]))
        segments(x[[4]][location, 1], 0, x1 = x[[4]][location, 1],
                 y1 = x[[4]][location, 2], lty = "dashed")
        legend(par("usr")[2], par("usr")[4], xpd = NA,
               legend = paste("m=", 1, sep = ""), col = 1, lty = 1, bty = "n",
               title = "Step")
      }
    }
    par(mar = c(4.5, 5, 2.5, 1))
    if (length(x[[1]]) == 0) {
      plot(0, xlab = "x", ylab = "Succes Probability",
           ylim = c(0, 1.1 * x[[2]]), xlim = c(0, 1), type = "n")
      points(0:1, rep(x[[2]], 2), type = "l")
    } else {
      if (length(x[[1]]) > 1) {
        x_range <- range(x[[1]])
      } else {
        x_range <- c(0.95 * x[[1]], 1.05 * x[[1]])
      }
      add_to_x_range <- 0.05 * diff(x_range)
      xvec <- c(x_range[1] - add_to_x_range, rep(x[[1]], each = 2),
                x_range[2] + add_to_x_range)
      yvec <- rep(x[[2]], each = 2)
      par(las = 1)
      plot(0, xlab = "Breakpoints", ylab = "Succes Probability",
           ylim = c(0, max(1.1 * x[[2]])), xlim = c(min(xvec), max(xvec)),
           type = "n", xaxt = "n", yaxt = "n")
      axis(1, at = x[[1]], labels = round(x[[1]], 3))
      axis(2, at = x[[2]], labels = round(x[[2]], 3))

      for (i in seq(1, (length(xvec) - 1), by = 2)) {
        points(xvec[i:(i + 1)], yvec[i:(i + 1)], type = "l")
      }
      for (i in seq(2, (length(xvec) - 1), by = 2)) {
        points(xvec[i:(i + 1)], c(0, max(yvec[i:(i + 1)])),
               type = "l", lty = "dotted")
      }
    }
    par(mar = old_par, las = 0, ask = FALSE)
}
