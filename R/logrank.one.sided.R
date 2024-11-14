#' Perform one-sided logrank test
#'
#' This functions performs the one-sided logrank test. The standard logrank test in survival package
#' only produces two-sided test. This function can facilitate one-sided logrank test. Larger z statistic indicates
#' better treatment effect.
#'
#' @param time survival time
#' @param event event status (0 = censor, 1 = event)
#' @param group group indicator (0 = control, 1 = experimental arm)
#' @param STRATA strata variable for stratified log-rank test
#'
#' @return
#' \describe{
#' \item{z}{Test statistics z value}
#' \item{p}{one sided p value}
#' }
#'
#' @examples
#'
#' n <- 100
#' time <- c(rexp(n, rate=log(2)/12), rexp(n, rate=log(2)/12*1.2))
#' event <- rep(1, n*2)
#' group <- c(rep(0, n), rep(1, n))
#' STRATA <- rep(c(0,1), n)
#'
#' logrank.one.sided(time, event, group, STRATA)
#'
#' @importFrom survival survdiff Surv strata
#' @export
#'
logrank.one.sided <- function(time, event, group, STRATA=NULL){
  if(is.null(STRATA)){
    lr.test <- survdiff(Surv(time, event) ~ group)
  }else{
    lr.test <- survdiff(Surv(time, event) ~ group + strata(STRATA))
  }
  #convert to z value in correct direction: z>0 means better experimental arm.
  if (is.matrix(lr.test$obs)) {
    otmp <- apply(lr.test$obs, 1, sum)
    etmp <- apply(lr.test$exp, 1, sum)
  } else {
    otmp <- lr.test$obs
    etmp <- lr.test$exp
  }
  better <- as.numeric(otmp[2] < etmp[2])
  sign <- 2*better - 1
  z <- sqrt(lr.test$chisq) * sign
  return(list(z = z, p = 1-pnorm(z), obs = otmp))
}

