#' Get operating characteristics via simulations for the phase 2/3 design
#'
#' @param nsim number of replicates
#' @param test.method method for getting adjusted p-value: \code{"dunnett"} or \code{"simes"}.
#' @param approach approach for combining data:
#' \code{"disjoint"}: disjoint-subjects approach;
#' \code{"ii"}: independent incremental approach.
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
#' @param min_followup minimum follow-up in month for stage 1 before dose selection. It should be smaller than
#' \code{enrollment_hold}. For example, a reasonable choice is \code{min_followup = 4} with \code{enrollment_hold=5}
#' which gives us at least 4 months follow-up for stage 1 patients and 1 month for dose selection.
#' @param targetEventsIA_all target number of events at IA for all subjects from both stage 1 and 2.
#' @param targetEventsFA_all target number of events at FA for all subjects from both stage 1 and 2.
#' @param w weights parameter for stage 1 data: for independent incremental approach, it will be length of 1;
#' for the disjoint-subjects approach, it will be length of 2 for IA and FA.
#' @param bound_z z-scale rejection boundaries and IA and FA.
#' @param alpha Type I error, always one-sided.
#' @param update_bound whether to re-calculate the FA rejection boundary using the updated
#' correlation between IA and FA test statistics.
#' @param nonselected_max_followup The final data cutoff time for non-selected arms.
#' The default \code{NULL} means the non-selected arms will not be followed up after dose selection.
#' @param dose_selection_endpoint either "ORR" or "Survival"
#'
#' @return It returns a matrix with each row corresponding to the analysis results for each trial,
#' where the analysis results include: "Selected dose", "IA time", "IA reject", "FA reject", "FA time",
#' "Study duration", "Total sample size", "observed correlation", and "actual FA boundary".
#' @importFrom dplyr filter group_by summarise
#' @importFrom rlang .data
#' @importFrom gsDesign gsDesign sfLDOF
#' @export
#'
#' @examples
#' res <- getOC_ph23(nsim=10)
#' apply(res, 2, mean, na.rm=TRUE)
getOC_ph23 <- function(nsim = 1000, test.method = "dunnett",
                       approach = "disjoint", n1 = c(50, 50, 50, 50), n2 = c(150, 150),
                       duration1 = 8, duration2 = 12, enrollment_hold = 5,
                       hazard = c(0.1, 0.09, 0.08, 0.07),
                       orr = c(0.2, 0.25, 0.30, 0.40), rho = 0.7,
                       dropout = c(0, 0, 0, 0), a1 = 1, a2 = 1,
                       min_followup = 4,
                       targetEventsIA_all = 240, targetEventsFA_all = 320,
                       w = NULL, bound_z = c(2.44, 2),
                       alpha = 0.025, update_bound = TRUE,
                       nonselected_max_followup = NULL,
                       dose_selection_endpoint = "ORR"){
  ## rename boundaries
  bound_z_IA <- bound_z[1]
  bound_z_FA <- bound_z[2]

  ## Start simulation
  results.names <- c("Selected dose", "IA time", "IA reject", "FA reject", "FA time",
                     "Study duration", "Total sample size", "observed correlation",
                     "actual FA boundary")
  results <- matrix(NA, nrow=nsim, ncol=length(results.names))
  colnames(results) <- results.names

  for(sim in 1:nsim){
    ## Simulate survival times and tumor responses
    num_trt <- length(hazard)-1
    dat <- sim_ph23(n1 = n1, n2 = n2, duration1 = duration1, duration2 = duration2,
                    enrollment_hold = enrollment_hold, hazard = hazard,
                    orr = orr, rho = rho, dropout = dropout, a1 = a1, a2 = a2)

    ## dose selection
    d1 <- dat %>% filter(.data$stage==1)
    IAd_time <- duration1 + min_followup
    if(dose_selection_endpoint == "ORR"){
      orr_est <- d1 %>% filter(.data$trt!=0) %>%
        group_by(.data$trt) %>% summarise(orr = mean(.data$response))
      selected <- orr_est$trt[which.max(orr_est$orr)]
    }else{
      d1IAd <- cut_by_date(d1, cut_time = IAd_time)
      pvalues <- rep(NA, num_trt)
      for(i in 1:num_trt){
        d1IAi <- d1IAd %>% filter(.data$trt%in%c(0, i))
        res <- logrank.one.sided(time = d1IAi$survTimeCut, event = d1IAi$eventCut,
                                 group = as.factor(d1IAi$trt))
        pvalues[i] <- res$p
      }
      selected <- which.min(pvalues)
    }

    ## non-selected arms will be cut at IAd_time or nonselected_max_followup (if specified)
    if(is.null(nonselected_max_followup)){
      dat2 <- cut_nonselected_arms(dat = dat, selected = selected,
                                   nonselected_max_followup = IAd_time)
    }else{
      dat2 <- cut_nonselected_arms(dat = dat, selected = selected,
                                   nonselected_max_followup = nonselected_max_followup)
    }

    ## IA
    if(approach=="ii"){
      res_IA <- getZ1(dat = dat2, IAd_time = IAd_time, w = w[1], selected = selected, targetEvents = targetEventsIA_all,
                      test.method = test.method)
    }else{
      res_IA <- getZstat(dat = dat2, w = w[1], selected = selected, targetEvents = targetEventsIA_all,
                         test.method = test.method)
      obsEventsIA_1 <- res_IA$obsEvents_stage1
    }
    z_IA_tilde <- res_IA$Z_tilde
    z_IA <- res_IA$Z
    obsEventsIA <- res_IA$obsEvents_all
    IA_time <- res_IA$cut_time
    IA_sample_size <- res_IA$sample_size
    w_IA <- res_IA$w
    # decision
    IA_reject <- (z_IA_tilde >= bound_z_IA)

    ## FA
    if(IA_reject){
      FA_reject <- IA_reject
      FA_time <- NA
      Study_duration <- IA_time
      total_sample_size <-  IA_sample_size
      obs_cor <- NA
      if(update_bound) bound_z_FA_obs <- NA
    }else{
      if(approach=="disjoint"){
        res_FA <- getZstat(dat = dat2, w = w[2], selected = selected, targetEvents = targetEventsFA_all,
                           test.method = test.method)
        FA_time <- res_FA$cut_time
        obsEventsFA_1 <- res_FA$obsEvents_stage1
        obsEventsFA <- res_FA$obsEvents_all
        FA_sample_size <- res_FA$sample_size
        z_FA_tilde <- res_FA$Z_tilde
        w_FA <- res_FA$w
        obs_cor <- w_IA*w_FA*sqrt(obsEventsIA_1/obsEventsFA_1) +
          sqrt(1-w_IA^2)*sqrt(1-w_FA^2)*sqrt((obsEventsIA-obsEventsIA_1)/(obsEventsFA-obsEventsFA_1))
      }else{
        # calculate incremental statistic
        d <- dat2 %>% filter(.data$trt%in%c(0, selected))
        dFA <- cut_by_event(d, targetEvents = targetEventsFA_all)
        FA_time <- dFA$calendarCutoff[1]
        obsEventsFA <- sum(dFA$eventCut)
        res <- logrank.one.sided(time = dFA$survTimeCut, event = dFA$eventCut,
                                 group = as.factor(dFA$trt))
        z_FA <- res$z
        z_FA_tilde <- z_FA + sqrt(obsEventsIA/obsEventsFA)*(z_IA_tilde-z_IA)
        FA_sample_size <- nrow(dFA) + sum((n1[-1])[-selected])
        obs_cor <- sqrt(obsEventsIA/obsEventsFA)
      }
      if(update_bound){
        obs_info_fraction <- obs_cor^2
        du <- gsDesign(k = 2, test.type = 1, alpha = alpha, sfu = sfLDOF,
                       n.I = c(targetEventsIA_all, targetEventsIA_all/obs_info_fraction),
                       maxn.IPlan = targetEventsFA_all)
        bound_z_FA_obs <- du$upper$bound[2]
        FA_reject <- (z_FA_tilde >= bound_z_FA_obs)
      }else{
        FA_reject <- (z_FA_tilde >= bound_z_FA)
      }
      Study_duration <- FA_time
      total_sample_size <-  FA_sample_size
    }

    ## Save results
    results[sim, 1] <- selected
    results[sim, 2] <- IA_time
    results[sim, 3] <- IA_reject
    results[sim, 4] <- FA_reject
    results[sim, 5] <- FA_time
    results[sim, 6] <- Study_duration
    results[sim, 7] <- total_sample_size
    results[sim, 8] <- obs_cor
    if(update_bound){
      results[sim, 9] <- bound_z_FA_obs
    }else{
      results[sim, 9] <- bound_z_FA
    }
  }
  return(results)
}
