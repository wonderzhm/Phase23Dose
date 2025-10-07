#' Simulate a single arm survival data
#'
#' This function will simulate a single arm survival data with flexible enrollment distribution.
#' @param n total number of subjects to be enrolled.
#' @param duration enrollment duration in months; must be an integer.
#' @param hazard hazard rate
#' @param orr objective response rate
#' @param rho the correlation between ORR and survival based on Gaussian copula.
#' @param dropout dropout hazard rate.
#' @param a weight parameter in cumulative enrollment pattern. When \code{a=1}, it is uniform.
#'
#' @return A time-to-event dataset consists of
#' \itemize{
#' \item \code{enterTime}: subject entry time relative to the study starting date
#' \item \code{response}: 0=no response; 1=response
#' \item \code{survTime}: observed event time
#' \item \code{event}: 0=censored; 1=event
#' \item \code{calendarTime}: total time on study, i.e. \code{enterTime} + \code{survTime}.
#' }
#' @importFrom rlang .data
#' @importFrom dplyr mutate arrange
#' @importFrom stats pnorm qnorm rexp
#' @importFrom mvtnorm rmvnorm
#' @export
#'
#' @examples
#' d <- sim_single_arm(n=50)
#' head(d)
sim_single_arm <- function(n = 50, duration = 2, hazard = 0.1, orr = 0.2,
                           rho = 0.7, dropout = 0, a = 1){
  ## Simulate survival times and tumor responses
  z <- rmvnorm(n=n, mean=c(0,0), sigma=matrix(c(1, rho, rho, 1), 2, 2))
  z_surv <- z[,1]
  z_orr <- z[,2]
  surv_time <- -pnorm(q = z_surv, log.p = TRUE)/hazard
  response <- (z_orr <= qnorm(orr))+0
  cen_time <- rexp(n, rate = max(dropout, 1e-10))

  ## Simulate enrollment
  if(duration<=1){
    enroll_time <- runif(n, min=0, max=duration)
  }else{
    nmonth <- ceiling(duration)
    n_per_month <- rep(NA, nmonth) # number of pts per month
    tmp <- 0
    for (i in 1:(nmonth-1)) {
      #ith month: cumulative #pts
      cN0i <- round((i/duration)^a * n)
      n_per_month[i] <- cN0i - tmp
      tmp = cN0i
    }
    n_per_month[nmonth] = n - sum(n_per_month[1:(nmonth-1)])
    enroll_time <- runif(n_per_month[1], min = 0, max = 1)
    if(nmonth>2){
      for(j in 2:(nmonth-1)){
        enroll_time <- c(enroll_time, runif(n_per_month[j], min = (j-1), max = j))
      }
    }
    enroll_time <- c(enroll_time, runif(n_per_month[nmonth], min = (nmonth-1), max = duration))
  }

  ## Combine the data
  d <- data.frame(enterTime = enroll_time, response = response,
                  surv_time = surv_time, cen_time = cen_time) %>%
    mutate(survTime = ifelse(surv_time<=cen_time, surv_time, cen_time)) %>%
    mutate(calendarTime = .data$enterTime + .data$survTime) %>%
    mutate(event = ifelse(surv_time<=cen_time, 1, 0)) %>%
    dplyr::select(-surv_time, -cen_time) %>%
    arrange(enterTime)
  return(d)
}
