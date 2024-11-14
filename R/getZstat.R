#' Get combination test statistic value using disjoint-subjects approach
#'
#' @param dat A time-to-event dataset returned from \link{sim_ph23}.
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
#' d <- sim_ph23(n1_per_arm = 20, n2_per_arm = 30)
#' getZstat(d, w=0.8, selected = 2, targetEvents = 25, test.method="dunnett")
getZstat <- function(dat, w = NULL, selected, targetEvents, test.method = "dunnett"){
  num_trt <- length(unique(dat$trt)) - 1
  # cut data
  d <- dat %>% filter(.data$trt%in%c(0, selected))
  dIA <- cut_by_event(d, targetEvents = targetEvents)
  res <- logrank.one.sided(time = dIA$survTimeCut, event = dIA$eventCut,
                           group = as.factor(dIA$trt))
  ZIA <- res$z # non-adjusted Z statistic
  IA_time <- dIA$calendarCutoff[1]
  obsEventsIA <- sum(dIA$eventCut)
  # stage 1:
  d1 <- dat %>% filter(.data$stage==1)
  d1IA <- cut_by_date(d1, cut_time = IA_time)
  pvalues <- rep(NA, num_trt)
  for(i in 1:num_trt){
    d1IAi <- d1IA %>% filter(.data$trt%in%c(0, i))
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
  d1IAselected <- d1IA %>% filter(.data$trt%in%c(0, selected))
  obsEventsIA_1 <- sum(d1IAselected$eventCut)
  ZIA_1 <- qnorm(1-max(padjusted)) # adjusted Z statistic for stage 1 subjects
  # stage 2:
  d2IA <- dIA %>% filter(.data$stage==2)
  res <- logrank.one.sided(time = d2IA$survTimeCut, event = d2IA$eventCut,
                           group = as.factor(d2IA$trt))
  ZIA_2 <- res$z # Z statistics for stage 2 subjects, no adjustment needed
  # p-value combination
  if(is.null(w)) w <- sqrt(obsEventsIA_1/obsEventsIA)
  ZIA_tilde <- w*ZIA_1 + sqrt(1-w^2)*ZIA_2
  # sample size
  sample_size <- nrow(d1IA) + nrow(d2IA)
  # resulted test statistic
  return(list(Z_tilde = ZIA_tilde, Z = ZIA, obsEvents_stage1 = obsEventsIA_1,
              obsEvents_all = obsEventsIA, w = w,
              cut_time = IA_time, sample_size = sample_size))
}
