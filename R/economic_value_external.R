#' @title
#' General function to compute an economic value for carcass fat, carcass conformation or carcass weight.
#' @description
#' This function creates an economic value for either carcass fat and carcass conformation or carcass weight. It uses the "internal finction using thresholds". First it distinguishes between carcass conformation/carcass fat and carcass weight. Carcass weight has thresholds as input parameter, but carcass conformation and carcass fat have population frequencies per class as input parameters. There the thresholds are computed first from the frequencies using the characteristics of normal distributions.
#' @param pvec_class_freq  discrete empirical distribution of population for given trait over different classes
#' @param pvec_threshold   thresholds under continuous distribution given by pricing system
#' @param pvec_price       vector of prices for given classes
#' @param pn_mean          current population phenotypic mean
#' @param pn_sd            current population phenotypic standard deviation
#' @param pn_delta_mean    small change of current phenotypic population mean
#' @return economic value
#' @export compute_economic_value
compute_economic_value <- function( pvec_class_freq = NULL,
                                    pvec_threshold  = NULL,
                                    pvec_price,
                                    pn_mean,
                                    pn_sd,
                                    pn_delta_mean ){
  # general format requirements for input variables
  if(length(pn_mean)>1){
    cat("pn_mean must be one number")
    return()
  }
  if(length(pn_sd)>1){
    cat("pn_sd must be one number")
    return()
  }
  if(length(pn_delta_mean)>1){
    cat("pn_delta_mean must be one number")
    return()
  }
  if(pn_sd<0){
    cat("pn_sd must be positive")
    return()
  }
  if(is.null(pvec_class_freq)) {
    # case of carcass weight
    # computation of economic value for carcass weight
    # special format requirements for input variables
    if(length(pvec_threshold)!=length(pvec_price)-1){
      cat("pvec_threshold must be one vector element shorter than pvec_price ")
      return()
    }

    # fit threshold vector by adding infinity to its tail
    pvec_threshold <- c(pvec_threshold,Inf)
    # compute economic value
    ev_result <- compute_internal_economic_value(pvec_threshold = pvec_threshold,
                                                 pn_sd = pn_sd, pn_mean = pn_mean, pvec_price = pvec_price,
                                                 pn_delta_mean = pn_delta_mean)
    return(ev_result)



  } else if(is.null(pvec_threshold)){
    # case of cf* and cc*
    # special format requirements for input variables
    if(length(pvec_class_freq)!=length(pvec_price)){
      cat("pvec_class_freq must be the same length as pvec_price. ")
      return()
    }
    # define thresholds within normal distribution
    vec_cumsum <- cumsum(pvec_class_freq)
    if (vec_cumsum[length(vec_cumsum)] > 1) {
      vec_cumsum[length(vec_cumsum)] <- 1
    }
    pvec_threshold<- qnorm(vec_cumsum,sd=pn_sd,mean=pn_mean, lower.tail = T)
    # compute economic value
    ev_result <- compute_internal_economic_value(pvec_threshold = pvec_threshold,
                                                 pn_sd = pn_sd, pn_mean = pn_mean, pvec_price = pvec_price,
                                                 pn_delta_mean = pn_delta_mean)

    return(ev_result)
  }

}
