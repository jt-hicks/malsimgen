#------------------------------------------------
#' Function that returns values from optional seasonal deterministic model
#'
#' \code{genrandwalk} generate random walk of EIR (recursive fn)
#'
#' @param x position in list
#' @param vol Volatility of the random walk
#' @param randWalk Values in random walk so far
#' @param min_EIR Minimum EIR value allowed
#' @param max_EIR Maximum EIR value allowed
#'
#' @export
#
genrandwalk <- function(x,vol,randWalk,min_EIR,max_EIR) {
  if (x == 0)    return (randWalk)
  else return(genrandwalk(x-1,vol,c(randWalk,max(min_EIR,min(exp(log(randWalk[length(randWalk)])+rnorm(1)*vol),max_EIR))),min_EIR,max_EIR))
}
