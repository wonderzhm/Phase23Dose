#' @importFrom gMCP generateWeights simes.test
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats qnorm

## non-selected arms will be cut at dose selection time
cut_nonselected_arms <- function(dat, selected, nonselected_max_followup){
  d1 <- dat %>% filter(.data$stage==1)
  d1nonselected <- d1 %>% filter(!(.data$trt%in%c(0, selected)))
  dremaining <- dat %>% filter((.data$trt%in%c(0, selected)))
  d1nonselected_cut <- cut_by_date(d1nonselected, nonselected_max_followup)
  d1nonselected$survTime <- d1nonselected_cut$survTimeCut
  d1nonselected$event <- d1nonselected_cut$eventCut
  d1nonselected$calendarTime <- d1nonselected$enterTime + d1nonselected$survTime
  dat_cut <- dplyr::bind_rows(d1nonselected, dremaining)
  return(dat_cut)
}


# Get p-values for simes' method
getPvals.simes <- function(g, w, p, selected){
  n <- length(p)
  weis <- generateWeights(g, w)
  res <- rep(NA, nrow(weis))
  for(i in 1:nrow(weis)){
    res[i] <- simes.test(pvalues = p, weights = weis[i, (1:n)+n])
  }
  return(res[weis[,selected]==1])
}

# Get p-values for dunnett method
getPvals.dunnett <- function(g,w,cr,p,adjusted=TRUE,upscale=FALSE){#, alternatives="less"){
  hint <- generateWeights(g,w)
  res <- t(apply(hint,1,Pvals.dunnett,p=p,cr=cr,upscale=upscale))#, alternatives=alternatives))
  if(adjusted){
    return(ad.p(res))
  } else {
    return(res)
  }
}

## At the moment hypotheses that are not tested at all get an adj. p-value of 1
ad.p <- function(P){
  p.ad <- rep(NA,ncol(P))
  for(i in 1:ncol(P)){
    out <- apply(P[!is.na(P[,i]),],1,min,na.rm=T)
    p.ad[i] <- ifelse(length(out)>0,max(out),1)
  }
  return(p.ad)
}

## pvalues for dunnett method
Pvals.dunnett <- function(h,cr,p,upscale, alternatives="less") {
  #  if(a > .5){
  #    stop("alpha levels above .5 are not supported")
  #  }
  n <- length(h)
  I <- h[1:(n/2)]
  w <- h[((n/2)+1):n]
  hw <- sapply(w,function(x) !isTRUE(all.equal(x,0)))
  e <- which(I>0 & hw)
  zb <- rep(NA,n/2)
  if(length(e) == 0){
    return(zb)
  }
  zb[e] <- p.dunnet.ph23(p[e],cr[e,e],w[e],upscale, alternatives=alternatives)
  zb[which(I>0 & !hw)] <- 1
  return(zb)
}

p.dunnet.ph23 <- function(p,cr,w,upscale, alternatives="less"){
  if(length(cr)>1){
    conn <- conn.comp.ph23(cr)
  } else {
    conn <- 1
  }
  twosided <- alternatives==rep("two.sided", length(w))
  lconn <- sapply(conn,length)
  conn <- lapply(conn,as.numeric)
  e <- sapply(1:length(p),function(i){
    sum(sapply(conn,function(edx){
      if(length(edx)>1){
        if (upscale=="o3") {
          return((1-pmvnorm(
            lower=ifelse(twosided[edx],qnorm(pmin(1,(w[edx]*p[i]/(w[i]*sum(w))))/2),-Inf),
            upper=ifelse(twosided[edx],qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]*sum(w))))/2),qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]*sum(w)))))),
            corr=cr[edx,edx], seed = 20240716, abseps=10^-5)))
        } else {
          return((1-pmvnorm(
            lower=ifelse(twosided[edx],qnorm(pmin(1,(w[edx]*p[i]/(w[i])))/2),-Inf),
            upper=ifelse(twosided[edx],qnorm(1-pmin(1,(w[edx]*p[i]/(w[i])))/2),qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]))))),
            corr=cr[edx,edx], seed = 20240716,abseps=10^-5))/ifelse(upscale,1,sum(w)))
        }
      } else {
        if(upscale=="o3" || !upscale){
          return((w[edx]*p[i]/(w[i]*sum(w))))
        } else {
          return((w[edx]*p[i]/(w[i])))
        }
      }
    }))})

  e <- pmin(e,1)
  e
}

conn.comp.ph23 <- function(m){
  N <- 1:ncol(m)
  M <- numeric(0)
  out <- list()
  while(length(N)>0){
    Q <- setdiff(N,M)[1]
    while(length(Q)>0){
      w <- Q[1]
      M <- c(M,w)
      Q <- setdiff(unique(c(Q,which(!is.na(m[w,])))),M)
    }
    out <- c(out,list(M))
    N <- setdiff(N,M)
    M <- numeric(0)
  }
  return(out)
}
