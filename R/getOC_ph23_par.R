#' Get operating characteristics via simulations for the phase 2/3 design
#'
#' @param ncore number of cores for parallel computing.
#' @param seed random seed for reproducibility.
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
#' @importFrom parallel makeCluster detectCores clusterCall clusterExport
#' @importFrom parallel parLapply stopCluster clusterSetRNGStream
#' @importFrom doParallel registerDoParallel
#' @export
#'
#' @examples
#' res <- getOC_ph23_par(ncore = 5, nsim=10)
#' apply(res, 2, mean, na.rm=TRUE)
getOC_ph23_par <- function(ncore = 2, seed = 2024, nsim = 10, test.method = "dunnett",
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
  ## Start simulation
  if(ncore > detectCores()) stop("ncore is too big!")
  cl <- makeCluster(ncore, type = "PSOCK")
  registerDoParallel(cl)
  ## Export Functions to the Cluster
  tmp1 <- clusterCall(cl, function() {library(Phase23Dose)})
  ## Function input
  # Set one master seed; each worker gets an independent L'Ecuyer-CMRG substream
  clusterSetRNGStream(cl, iseed = seed)
  results <- parLapply(
    cl, rep(nsim, ncore), fun = getOC_ph23, test.method = test.method, approach = approach,
    n1 = n1, n2 = n2, duration1 = duration1, duration2 = duration2,
    enrollment_hold = enrollment_hold, hazard = hazard,
    orr = orr, rho = rho, dropout = dropout, a1 = a1, a2 = a2,
    min_followup = min_followup, targetEventsIA_all = targetEventsIA_all,
    targetEventsFA_all = targetEventsFA_all,
    w = w, bound_z = bound_z, alpha = alpha, update_bound = update_bound,
    nonselected_max_followup = nonselected_max_followup,
    dose_selection_endpoint = dose_selection_endpoint)
  ## Close Cluster
  stopCluster(cl)
  res <- do.call(rbind, results) # combine results
  return(res)
}
