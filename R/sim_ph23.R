#' Simulate a seamless phase 2/3 trial
#'
#' This function will simulate a phase 2/3 trial where patients will be randomized to all arms in stage 1
#' and only the selected arm and control arm will go into stage 2. Since we don't know the arm selection in
#' advance, we will simulated patients for all arms in stage 2 but with adjusted accrual rate assuming stage 2
#' only contains the selected arm and the control.
#' @param n1 total number of subjects to be enrolled in stage 1 for control and treatment arms, where
#' the number of treatment arms can be greater than 2.
#' @param n2 total number of subjects to be enrolled in stage 2 for control and the selected treatment arm.
#' @param duration1 enrollment duration in months for stage 1; must be an integer.
#' @param duration2 enrollment duration in months for stage 2; must be an integer.
#' @param enrollment_hold holding period in month after stage 1 enrollment prior to enrollment of stage 2 patients.
#' 0 means seamless enrollment.
#' @param hazard hazard rate for control and treatment arms.
#' @param orr ORR for control and treatment arms.
#' @param rho the correlation between ORR and survival based on Gaussian copula.
#' @param dropout dropout hazard rate for control and treatment arms.
#' @param a1 weight parameter in cumulative enrollment pattern for stage 1. When \code{a=1}, it is uniform.
#' @param a2 weight parameter in cumulative enrollment pattern for stage 2. When \code{a=1}, it is uniform.
#'
#' @return A time-to-event dataset consists of
#' \itemize{
#' \item \code{stage}: 1=phase II; 2=phase III.
#' \item \code{trt}: 0=control; 1=dose 1; 2=dose 2; ...
#' \item \code{enterTime}: subject entry time relative to the study starting date
#' \item \code{response}: 0=no response; 1=response
#' \item \code{survTime}: observed event time
#' \item \code{event}: 0=censored; 1=event
#' \item \code{calendarTime}: total time on study, i.e. \code{enterTime} + \code{survTime}.
#' }
#' @importFrom rlang .data
#' @importFrom dplyr mutate arrange bind_rows
#' @importFrom stats pnorm qnorm rexp
#' @importFrom mvtnorm rmvnorm
#' @export
#'
#' @examples
#' d <- sim_ph23(n1 = c(50, 50, 50, 50), n2 = c(150, 150))
#' head(d)
sim_ph23 <- function(n1 = c(50, 50, 50, 50), n2 = c(150, 150),
                     duration1 = 8, duration2 = 12, enrollment_hold = 5,
                     hazard = c(0.1, 0.09, 0.08, 0.07),
                     orr = c(0.2, 0.25, 0.30, 0.40), rho = 0.7,
                     dropout = c(0, 0, 0, 0), a1 = 1, a2 = 1){
  ## number of treatment arms
  narms <- length(n1)

  ## simulate stage 1 data
  # control arm
  d1 <- sim_single_arm(n = n1[1], duration = duration1, hazard = hazard[1],
                       orr = orr[1], rho = rho, dropout = dropout[1], a = a1) %>%
    mutate(stage = 1, trt = 0, .before = .data$enterTime)
  # treatment arms
  for(i in 2:narms){
    di <- sim_single_arm(n = n1[i], duration = duration1, hazard = hazard[i],
                         orr = orr[i], rho = rho, dropout = dropout[i], a = a1) %>%
      mutate(stage = 1, trt = i-1, .before = .data$enterTime)
    d1 <- bind_rows(d1, di)
  }
  d1 <- d1 %>% arrange(enterTime)

  ## simulate stage 2 data
  # control arm
  d2 <- sim_single_arm(n = n2[1], duration = duration2, hazard = hazard[1],
                       orr = orr[1], rho = rho, dropout = dropout[1], a = a2) %>%
    mutate(enterTime = .data$enterTime + duration1 + enrollment_hold) %>%
    mutate(calendarTime = .data$enterTime + .data$survTime) %>%
    mutate(stage = 2, trt = 0, .before = .data$enterTime)
  # treatment arms
  for(i in 2:narms){
    di <- sim_single_arm(n = n2[2], duration = duration2, hazard = hazard[i],
                         orr = orr[i], rho = rho, dropout = dropout[i], a = a2) %>%
      mutate(enterTime = .data$enterTime + duration1 + enrollment_hold) %>%
      mutate(calendarTime = .data$enterTime + .data$survTime) %>%
      mutate(stage = 2, trt = i-1, .before = .data$enterTime)
    d2 <- bind_rows(d2, di)
  }
  d2 <- d2 %>% arrange(enterTime)

  ## Combine the data
  d <- bind_rows(d1, d2)
  return(d)
}
