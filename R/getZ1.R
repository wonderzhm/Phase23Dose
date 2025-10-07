#' Get test statistic value using independent incremental approach
#'
#' @param dat A time-to-event dataset returned from \link{sim_ph23}.
#' @param IAd_time the data cutoff time for dose selection.
#' @param w weight parameter for stage 1 data.
#' @param selected The selected dose level.
#' @param targetEvents The target number of events from both stages at which analysis is to be made.
#' @param test.method The method for getting adjusted p-value: \code{"dunnett"} or \code{"simes"}.
#'
#' @return It returns a list including \code{z} the test statistic, \code{obsEvents_stage1} the
#' observed number of events from stage 1 and \code{cut_time} the data cut time.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' d <- sim_ph23(n1 = c(20, 20, 20, 20), n2 = c(30, 30), duration1 = 4, enrollment_hold = 5)
#' getZ1(d, IAd_time = 13, selected = 2, targetEvents = 30, test.method="dunnett")
getZ1 <- function(dat, IAd_time, w = NULL, selected = 1, targetEvents, test.method = "dunnett"){
  num_trt <- length(unique(dat$trt)) - 1
  # stage 1:
  d1 <- dat %>% filter(.data$stage==1);
  d1IAd <- cut_by_date(d1, cut_time = IAd_time)
  pvalues <- rep(NA, num_trt)
  for(i in 1:num_trt){
    d1IAi <- d1IAd %>% filter(.data$trt%in%c(0, i))
    res <- logrank.one.sided(time = d1IAi$survTimeCut, event = d1IAi$eventCut,
                             group = as.factor(d1IAi$trt))
    pvalues[i] <- res$p
  }
  names(pvalues) <- paste("H", 1:num_trt, sep = "")
  g <- matrix(1/(num_trt-1), nrow = num_trt, ncol = num_trt)
  diag(g) <- 0
  rownames(g) <- colnames(g) <- names(pvalues)
  weights <- rep(1, num_trt)/num_trt
  cr <- matrix(0.5, num_trt, num_trt)
  diag(cr) <- 1
  if(test.method == "dunnett"){
    pall <- getPvals.dunnett(g=g, w=weights, cr=cr, p=pvalues, adjusted = FALSE)
    pselected <- pall[!is.na(pall[, selected]),]
    padjusted <- apply(pselected, 1, min, na.rm = TRUE)
  }else{
    padjusted <- getPvals.simes(g=g, w=weights, p=pvalues, selected = selected)
  }
  d1IAd_selected <- d1IAd %>% filter(.data$trt%in%c(0, selected))
  obsEventsIAd <- sum(d1IAd_selected$eventCut)
  Zd_tilde <- qnorm(1-max(padjusted))
  Zd <- qnorm(1-pvalues[selected])
  # stage 2: include selected arm and control from both stages
  d <- dat %>% filter(.data$trt%in%c(0, selected))
  dIA1 <- cut_by_event(d, targetEvents = targetEvents)
  IA1_time <- dIA1$calendarCutoff[1]
  obsEventsIA1 <- sum(dIA1$eventCut)
  res <- logrank.one.sided(time = dIA1$survTimeCut, event = dIA1$eventCut,
                           group = as.factor(dIA1$trt))
  Z1 <- res$z
  # combine using independent incremental
  if(is.null(w)) w <- sqrt(obsEventsIAd/obsEventsIA1)
  Z1_tilde <- w*Zd_tilde + sqrt(1-w^2)*(sqrt(obsEventsIA1)*Z1-sqrt(obsEventsIAd)*Zd)/
    sqrt(obsEventsIA1-obsEventsIAd)
  #Z1_tilde <- Z1 + sqrt(obsEventsIAd/obsEventsIA1)*(Zd_tilde-Zd)
  # sample size
  d2IA1 <- dIA1 %>% filter(.data$stage==2)
  sample_size <- nrow(d1IAd) + nrow(d2IA1)
  # resulted test statistic
  return(list(Z_tilde = Z1_tilde, Z = Z1, obsEvents_all = obsEventsIA1,
              w = w, cut_time = IA1_time, sample_size = sample_size))
}
